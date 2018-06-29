#include <iostream>
#include <fstream>
#include <unordered_set>

#include "config.h"
#include "cluster.h"

config_t config;

std::unordered_map<std::string, IntervalTree<simple_repeat_t> *> repeat_trees;
std::vector<std::string> contig_id2name;

int support(prediction_t& pred) {
    return pred.disc_pairs + std::min(pred.bp1.sc_reads, pred.bp2.sc_reads);
}

double ptn_score(prediction_t& pred) {
//    return std::max(double(pred.disc_pairs+pred.bp1.sc_reads)/std::max(pred.bp1.spanning_reads,1),
//                    double(pred.disc_pairs+pred.bp2.sc_reads)/std::max(pred.bp2.spanning_reads,1));
    return double(pred.disc_pairs+pred.bp1.sc_reads)/std::max(pred.bp1.spanning_reads,1);
}

std::string formatted_print(prediction_t& pred) {
    char str[10000];
    if (pred.bp1.sc_reads > 0 && pred.sv_type != SV_TYPES.INS) {
        sprintf(str, "ID=%d BP1=%c:%s:%d BP2=%c:%s:%d %s DISC=%d SC=%d,%d SPANNING=%d,%d PVAL=%lf EST_SIZE=%d SHIFT-PVAL=%lf",
                pred.id,
                pred.bp1.dir, contig_id2name[pred.bp1.contig_id].c_str(), pred.bp1.pos(),
                pred.bp2.dir, contig_id2name[pred.bp2.contig_id].c_str(), pred.bp2.pos(),
                svt_to_str(pred.sv_type).c_str(), pred.disc_pairs, pred.bp1.sc_reads, pred.bp2.sc_reads,
                pred.bp1.spanning_reads, pred.bp2.spanning_reads, pred.pval, pred.get_size(), pred.shift_pval);
    } else {
        sprintf(str, "ID=%d BP1=%c:%s:%d BP2=%c:%s:%d %s DISC=%d SC=%d,%d SPANNING=%d,%d PVAL=%lf EST_SIZE=%d:%d SHIFT-PVAL=%lf",
                pred.id,
                pred.bp1.dir, contig_id2name[pred.bp1.contig_id].c_str(), pred.bp1.pos(),
                pred.bp2.dir, contig_id2name[pred.bp2.contig_id].c_str(), pred.bp2.pos(),
                svt_to_str(pred.sv_type).c_str(), pred.disc_pairs, pred.bp1.sc_reads, pred.bp2.sc_reads,
                pred.bp1.spanning_reads, pred.bp2.spanning_reads, pred.pval, pred.size - pred.conf_ival,
                pred.size + pred.conf_ival, pred.shift_pval);
    }
    return str;
}

int overlap(int s1, int e1, int s2, int e2) {
    return std::max(0, std::min(e1, e2)-std::max(s1, s2));
}

std::unordered_set<int> in_reps_ids;

bool operator < (const prediction_t& p1, const prediction_t& p2) {
    bool in_rep1 = in_reps_ids.count(p1.id);
    bool in_rep2 = in_reps_ids.count(p2.id);
    if (in_rep1 != in_rep2) return in_rep1 < in_rep2;
    if (p1.bp1.sc_reads*p1.bp2.sc_reads != p2.bp1.sc_reads*p2.bp2.sc_reads) {
        return p1.bp1.sc_reads*p1.bp2.sc_reads > p2.bp1.sc_reads*p2.bp2.sc_reads;
    }
    if (p1.disc_pairs != p2.disc_pairs) {
        return p1.disc_pairs > p2.disc_pairs;
    }
    if (p1.pval != p2.pval) {
        return p1.pval < p2.pval;
    }
    return p1.id < p2.id;
}

int main(int argc, char* argv[]) {
    std::string workdir = argv[1];
    bool no_filter = std::string(argv[2]) == "no-filter";
    double alpha, ptn_ratio;
    int minsize;
    if (!no_filter) {
        alpha = argc > 2 ? std::stod(argv[2]) : 0.01;
        ptn_ratio = argc > 3 ? std::stod(argv[3]) : 0.33;
        minsize = argc > 4 ? std::stoi(argv[4]) : 50;

        config = parse_config(workdir + "/config.txt");
        if (argc > 5) {
            repeat_trees = parse_simple_repeats(argv[5], config.min_sv_len);
        }
    }


    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
    }

    std::vector<prediction_t> retained;

    std::ifstream predictions_fin(workdir + "/predictions.info");
    std::string line;
    while (predictions_fin >> line) {
        prediction_t pred(line);
        if (pred.sv_type == SV_TYPES.INS && pred.bp1.pos() > pred.bp2.pos()) std::swap(pred.bp1, pred.bp2); //TODO: temporary fix

        if (no_filter) {
            retained.push_back(pred);
            continue;
        }

        if (pred.pval < 0) { // odds ratio
            if (support(pred) >= config.avg_depth/6 && ptn_score(pred) > 0.33) {
                retained.push_back(pred);
            }
        } else { // stat testing
            // if no rep file was given, then we assume all SV are inside repetitive regions
            bool in_rep = true;
            std::string contig_name = contig_id2name[pred.bp1.contig_id];
            if (repeat_trees[contig_name] != NULL) {
                in_rep = false;
                std::vector<Interval<simple_repeat_t> > reps = repeat_trees[contig_name]->findOverlapping(pred.bp1.pos(), pred.bp2.pos());
                for (Interval<simple_repeat_t> rep : reps) {
                    if (pred.sv_type != SV_TYPES.DEL ||
                            overlap(rep.start, rep.stop, pred.bp1.pos(), pred.bp2.pos()) >= config.min_sv_len) {
                        in_reps_ids.insert(pred.id);
                        in_rep = true;
                    }
                }
            }

            if (!in_rep && pred.disc_pairs < std::max(3, config.avg_depth/8) && pred.bp1.sc_reads+pred.bp2.sc_reads == 0) continue;

            // if it's a deletion and it has enough SC support, accept automatically
            // it seems for insertions SC support is not so reliable, so only accept if in unique region
            int min_disc = (pred.sv_type == SV_TYPES.DEL || pred.sv_type == SV_TYPES.INS || pred.sv_type == SV_TYPES.NOV)
                           && pred.len() < 100 ? 0 : 1;
            if ((pred.sv_type == SV_TYPES.DEL || !in_rep) && pred.disc_pairs >= min_disc
                && pred.bp1.sc_reads >= 5 && pred.bp2.sc_reads >= 5) {
                retained.push_back(pred);
                continue;
            }

            if (pred.shift_pval <= alpha) pred.size *= 2;
            if (pred.pval <= alpha && pred.get_size() >= minsize) {
                if (pred.sv_type == SV_TYPES.DEL && pred.bp1.sc_reads+pred.bp2.sc_reads > 0 && pred.disc_pairs == 0 &&
                        abs(pred.size-(pred.bp2.pos()-pred.bp1.pos())) > pred.conf_ival) continue;
                retained.push_back(pred);
            }
        }
    }

    std::sort(retained.begin(), retained.end());
    for (prediction_t pred : retained) {
        std::cout << formatted_print(pred) << " " << in_reps_ids.count(pred.id) << std::endl;
    }
}
