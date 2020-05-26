#include <iostream>
#include <fstream>
#include <queue>
#include <unordered_set>
#include <numeric>
#include <random>
#include <htslib/sam.h>
#include <htslib/hts.h>


#include "config.h"
#include "cluster.h"
#include "libs/cptl_stl.h"
#include "ks-test.h"

std::string workdir;
std::mutex mtx;
std::mt19937 rng {std::random_device{}()};

config_t config;

int MAX_READ_IS;

const int MAX_BUFFER_SIZE = 100;

std::unordered_map<std::string, int> contig_name2tid;
std::unordered_map<std::string, size_t> contig_name2len;
std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;

std::queue<open_samFile_t*> bam_pool;

//alglib::real_1d_array population_del, population_ins;
std::vector<double> population_del_v, population_ins_v;
double population_del_mean, population_ins_mean;

open_samFile_t* get_bam_reader(std::string bam_fname) {
    if (!bam_pool.empty()) {
        open_samFile_t* o = bam_pool.front();
        bam_pool.pop();
        return o;
    }

    open_samFile_t* o = open_samFile(bam_fname.c_str());
    return o;
}

void release_bam_reader(open_samFile_t* reader) {
    bam_pool.push(reader);
}


int POP_SIZE = 10000;

int max_size_to_test() {
    return config.max_is - config.min_is;
}

bool is_stat_testable(prediction_t* pred) {
    if (pred->sv_type == SV_TYPES.DEL || pred->sv_type == SV_TYPES.DUP || pred->sv_type == SV_TYPES.INS) {
        return pred->len() <= config.max_is + max_size_to_test();
    } else {
        return false;
    }
}

double mean(std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double variance(std::vector<double>& v, double m = NAN) {
    if (m == NAN) m = mean(v);
    double acc = 0;
    for (double d : v) {
        acc += (d-m)*(d-m);
    }
    return acc/v.size();
}

bam1_t* add_to_queue(std::deque<bam1_t*>& q, bam1_t* o, int size_limit) {
    bam1_t* t = NULL;
    while (q.size() >= size_limit) {
        t = q.front();
        q.pop_front();
    }
    q.push_back(o);
    return t;
}

std::vector<double> gen_population_for_dels(std::string bam_fname) {
    open_samFile_t* open_sam = open_samFile(bam_fname.c_str());
    samFile* bam_file = open_sam->file;
    hts_idx_t *idx = open_sam->idx;
    bam_hdr_t* header = open_sam->header;

    std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
    std::vector<double> population;
    bam1_t* read = bam_init1();
    while (population.size() < POP_SIZE*100) {
        char contig[1000]; int pos;
        rnd_pos_fin >> contig >> pos;
        if (pos < 2*config.max_is || pos >= contig_name2len[std::string(contig)]-2*config.max_is) continue;

        char region[1000];
        sprintf(region, "%s:%d-%d", contig, pos-2*config.max_is, pos);
        hts_itr_t* iter = sam_itr_querys(idx, header, region);
        while (sam_itr_next(bam_file, iter, read) >= 0) {
            if (!is_valid(read, false) || bam_is_rev(read)) continue;

            if (is_samechr(read) && !is_samestr(read) && !is_outward(read, config.min_is) &&
                    read->core.isize > 0 && read->core.isize < MAX_READ_IS) {
                int start = read->core.pos+read->core.l_qseq/2, end = get_mate_endpos(read)-read->core.l_qseq/2;
                if (start >= end) continue;

                if (start <= pos && pos <= end) {
                    population.push_back((double) read->core.isize);
                }
            }
        }
        sam_itr_destroy(iter);
    }
    close_samFile(open_sam);
    bam_destroy1(read);

    std::shuffle(std::begin(population), std::end(population), rng);
    population.resize(POP_SIZE);

    return population;
}

std::vector<double> gen_population_for_ins(std::string bam_fname) {
    open_samFile_t* open_sam = open_samFile(bam_fname.c_str());
    samFile* bam_file = open_sam->file;
    hts_idx_t *idx = open_sam->idx;
    bam_hdr_t* header = open_sam->header;

    std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
    std::vector<double> population;
    bam1_t* read = bam_init1();
    while (population.size() < POP_SIZE) {
        char contig[1000]; int pos;
        rnd_pos_fin >> contig >> pos;

        char region[1000];
        sprintf(region, "%s:%d-%d", contig, pos-200, pos);
        hts_itr_t* iter = sam_itr_querys(idx, header, region);
        while (sam_itr_next(bam_file, iter, read) >= 0) {
            if (!is_valid(read, false) || bam_is_rev(read)) continue;

            if (is_samechr(read) && !is_samestr(read) && !is_outward(read, config.min_is) &&
                abs(read->core.isize) < MAX_READ_IS) {
                if (read->core.isize > 0 && read->core.pos > pos-200 && read->core.pos <= pos) {
                    population.push_back((double) read->core.isize);
                } else if (read->core.isize < 0 && bam_endpos(read) > pos-200 && bam_endpos(read) <= pos) {
                    population.push_back((double) -read->core.isize);
                }
            }
        }
        sam_itr_destroy(iter);
    }
    close_samFile(open_sam);
    bam_destroy1(read);

    std::shuffle(std::begin(population), std::end(population), rng);
    population.resize(POP_SIZE);

    return population;
}

void make_contig_name2tid(std::string bam_fname) {
    open_samFile_t* open_sam = open_samFile(bam_fname.c_str());
    samFile* bam_file = open_sam->file;
    bam_hdr_t* header = open_sam->header;

    for (int i = 0; i < header->n_targets; i++) {
        contig_name2tid[std::string(header->target_name[i])] = i;
        contig_name2len[std::string(header->target_name[i])] = header->target_len[i];
    }
    close_samFile(open_sam);
}

void stat_testing(int id, std::string bam_fname, prediction_t* pred) {

    mtx.lock();
    open_samFile_t* open_sam = get_bam_reader(bam_fname);
    mtx.unlock();

    std::vector<double> sample;
    std::string contig = contig_id2name[pred->bp1.contig_id];
    if (pred->sv_type == SV_TYPES.DEL) {
        char region[1000];
        int midpoint = (pred->bp1.pos()+pred->bp2.pos())/2;
        sprintf(region, "%s:%d-%d", contig.c_str(), midpoint-MAX_READ_IS, midpoint+MAX_READ_IS);

        std::deque<bam1_t*> two_way_buffer, forward_buffer;
        std::unordered_set<std::string> r1_ok;

        hts_itr_t* iter = sam_itr_querys(open_sam->idx, open_sam->header, region);
        bam1_t* read = bam_init1();

        int i = 0;
        while (sam_itr_next(open_sam->file, iter, read) >= 0 && i < MAX_BUFFER_SIZE-1) {
            if (is_valid(read, false)) {
                bam1_t *read2 = bam_dup1(read);
                add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
                add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
                i++;
            }
        }

        while (!forward_buffer.empty()) {
            while (sam_itr_next(open_sam->file, iter, read) >= 0) {
                if (is_valid(read, false)) {
                    bam1_t *read2 = bam_dup1(read);
                    bam1_t* to_destroy = add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
                    bam_destroy1(to_destroy);
                    add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
                    break;
                }
            }

            bam1_t* read = forward_buffer.front();
            forward_buffer.pop_front();

            std::string qname(bam_get_qname(read));
            if (bam_is_rev(read) && r1_ok.count(qname) > 0
                && check_SNP(read, two_way_buffer, config.avg_depth)) {
                sample.push_back((double) abs(read->core.isize));
            }

            if (is_samechr(read) && !is_samestr(read) && !is_outward(read, config.min_is) &&
                read->core.isize > 0 && read->core.isize < MAX_READ_IS) {
                int start = read->core.pos + read->core.l_qseq / 2, end =
                        get_mate_endpos(read) - read->core.l_qseq / 2;
                if (start >= end) continue;

                if (start <= midpoint && midpoint <= end && check_SNP(read, two_way_buffer, config.avg_depth)) {
                    r1_ok.insert(qname);
                }
            }

        }

        for (bam1_t* r : two_way_buffer) {
            bam_destroy1(r);
        }
        bam_destroy1(read);
        sam_itr_destroy(iter);
    } else {
        char region[1000];
        int left_pos = pred->bp1.pos(), right_pos = pred->bp2.pos();
        sprintf(region, "%s:%d-%d", contig.c_str(), left_pos-200, right_pos+200);

        hts_itr_t* iter = sam_itr_querys(open_sam->idx, open_sam->header, region);
        bam1_t* read = bam_init1();
        while (sam_itr_next(open_sam->file, iter, read) >= 0) {
            if (!is_valid(read, false)) continue;

            if (is_samechr(read) && !is_samestr(read) && !is_outward(read, config.min_is) &&
                abs(read->core.isize) < MAX_READ_IS) {
                if (read->core.isize > 0 && read->core.pos > left_pos-200 && read->core.pos <= left_pos) {
                    sample.push_back((double) read->core.isize);
                } else if (read->core.isize < 0 && bam_endpos(read) > right_pos && bam_endpos(read) <= right_pos+200) {
                    sample.push_back((double) -read->core.isize);
                }
            }
        }
        bam_destroy1(read);
        sam_itr_destroy(iter);
    }

    mtx.lock();
    release_bam_reader(open_sam);
    mtx.unlock();

    if (sample.size() < 5) {
        pred->pval = 1.0;
        return;
    }

    double pval;
    if (pred->sv_type == SV_TYPES.DEL) {
        pval = ks_test(population_del_v, sample);
    } else if (pred->sv_type == SV_TYPES.INS || pred->sv_type == SV_TYPES.DUP) {
        pval = ks_test(population_ins_v, sample);
    }

    pred->pval = pval;

    double sample_mean = mean(sample);
    double pop_mean = pred->sv_type == SV_TYPES.DEL ? population_del_mean : population_ins_mean;
    pred->size = sample_mean - pop_mean;

    for (int i = 0; i < sample.size(); i++) sample[i] -= pred->size;
    if (pred->sv_type == SV_TYPES.DEL) {
        pval = ks_test(population_del_v, sample);
    } else if (pred->sv_type == SV_TYPES.INS || pred->sv_type == SV_TYPES.DUP) {
        pred->size = -pred->size;
        pval = ks_test(population_ins_v, sample);
    }
    pred->shift_pval = pval;

    double var = variance(sample, sample_mean);
    pred->conf_ival = 1.96*sqrt(var/sample.size());
}

void find_spanning(breakpoint_t& bp, std::string& bam_fname) {
    std::string contig = contig_id2name[bp.contig_id];

    mtx.lock();
    open_samFile_t* open_sam = get_bam_reader(bam_fname);
    mtx.unlock();

    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), bp.pos()-2*config.max_is, bp.pos());

    hts_itr_t* iter = sam_itr_querys(open_sam->idx, open_sam->header, region);
    bam1_t* read = bam_init1();
    while (sam_itr_next(open_sam->file, iter, read) >= 0) {
        if (!is_valid(read, false)) continue;

        if (!bam_is_rev(read) && is_samechr(read) && !is_samestr(read) && !is_outward(read, config.min_is) &&
            read->core.isize >= config.min_is && read->core.isize <= config.max_is) {
            if ((bp.dir == 'R' && bp.pos() >= read->core.pos && bp.pos() < read->core.mpos) ||
                (bp.dir == 'L' && bp.pos() > bam_endpos(read) && bp.pos() <= get_mate_endpos(read))) {
                bp.spanning_pairs++;
            }
        }
        if (read->core.pos+10 < bp.pos() && bp.pos() < bam_endpos(read)-10) {
            bp.spanning_reads++;
        }
    }

    mtx.lock();
    release_bam_reader(open_sam);
    mtx.unlock();
};

void odds_ratio(int id, std::string bam_fname, prediction_t* pred) {
    find_spanning(pred->bp1, bam_fname);
    find_spanning(pred->bp2, bam_fname);
}

int main(int argc, char* argv[]) {
    std::string bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";

    config = parse_config(workdir + "/config.txt");
    MAX_READ_IS = config.max_is + max_size_to_test();

    std::vector<prediction_t*> preds;

    std::ifstream predictions_fin(workspace + "/predictions.raw");
    std::string line;
    while (predictions_fin >> line) {
        prediction_t* pred = new prediction_t(line);
        if ((pred->disc_pairs+pred->bp1.sc_reads > 1 && pred->disc_pairs+pred->bp2.sc_reads > 1) ||
                is_stat_testable(pred)) {
            preds.push_back(pred);
        } else {
            delete pred;
        }
    }

    make_contig_name2tid(bam_fname);
    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
        contig_name2id[contig_name] = contig_id;
    }

    std::cout << "Loading simple repeats from " << config.simple_rep_fname << std::endl;
    std::ifstream repeat_fin(config.simple_rep_fname);
    simple_repeat_t last_rep;
    while (getline(repeat_fin, line)) {
        if (line[0] != '#') {
            simple_repeat_t rep(line);
            if (last_rep.chr == rep.chr && last_rep.start <= rep.start && last_rep.end >= rep.end) // skip if contained in previous
                continue;
            last_rep = rep;

            if (rep.end-rep.start < config.max_is/* && rep.end-rep.start >= config.min_sv_len*/) {
//                cluster_t* c = new cluster_t(anchor_t('R', contig_name2id[rep.chr], rep.start-50, rep.start, 0),
    //                                             anchor_t('L', contig_name2id[rep.chr], rep.end, rep.end+50, 0),
    //                                            DISC_TYPES.SI, 0);
//                prediction_t* pred = new prediction_t(c, DISC_TYPES.SI);
//                preds.push_back(pred);

//                c = new cluster_t(anchor_t('R', contig_name2id[rep.chr], rep.start-50, rep.start, 0),
//                                             anchor_t('L', contig_name2id[rep.chr], rep.end, rep.end+50, 0),
//                                             DISC_TYPES.LI, 0);
//                pred = new prediction_t(c, DISC_TYPES.LI);
//                preds.push_back(pred);
            }
        }
    }
    std::cout << "Repeats loaded" << std::endl;

    population_del_v = gen_population_for_dels(bam_fname);
    population_del_mean = mean(population_del_v);

    population_ins_v = gen_population_for_ins(bam_fname);
    population_ins_mean = mean(population_ins_v);

    std::cout << "Population built." << std::endl;

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (prediction_t* pred : preds) {
        if (is_stat_testable(pred)) {
            std::future<void> future1 = thread_pool.push(stat_testing, bam_fname, pred);
            futures.push_back(std::move(future1));
            std::future<void> future2 = thread_pool.push(odds_ratio, bam_fname, pred);
            futures.push_back(std::move(future2));
        } else {
            std::future<void> future = thread_pool.push(odds_ratio, bam_fname, pred);
            futures.push_back(std::move(future));
        }
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    std::ofstream predictions_filter_out(workdir + "/predictions.info");
    for (prediction_t* pred : preds) {
        predictions_filter_out << pred->to_str() << "\n";
        delete pred;
    }
    predictions_filter_out.close();

    while (!bam_pool.empty()) {
        close_samFile(bam_pool.front());
        bam_pool.pop();
    }
}
