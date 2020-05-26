#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "sam_utils.h"
#include "config.h"
#include "libs/cptl_stl.h"

config_t config;
std::unordered_map<std::string, int> contig_name2id;
std::vector<int> tid_to_contig_id;

std::mutex mtx;

const int MAX_BUFFER_SIZE = 100;

std::string workdir;

bam1_t* add_to_queue(std::deque<bam1_t*>& q, bam1_t* o, int size_limit) {
    bam1_t* t = NULL;
    while (q.size() >= size_limit) {
        t = q.front();
        q.pop_front();
    }
    q.push_back(o);
    return t;
}

samFile* get_writer(std::string name, bam_hdr_t* header) {
    samFile* writer = sam_open((workdir + "/workspace/" + name).c_str(), "wb");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + (workdir + name);
    }
    return writer;
}

void categorize(int id, std::string contig, std::string bam_fname, int target_len) {

    samFile* bam_file = sam_open(bam_fname.c_str(), "r");
    if (bam_file == NULL) {
        throw "Unable to open BAM file.";
    }

    hts_idx_t* idx = sam_index_load(bam_file, bam_fname.c_str());
    if (idx == NULL) {
        throw "Unable to open BAM index.";
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == NULL) {
        throw "Unable to open BAM header.";
    }

    int contig_id = contig_name2id[contig];
    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), 1, target_len);

    mtx.lock();
    std::cout << "Categorizing " << region << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(idx, header, region);
    bam1_t* read = bam_init1();

    std::unordered_map<std::string, std::pair<bam1_t*, disc_type_t> > pairs;
    std::unordered_map<disc_type_t, samFile*> writers;

    samFile* clip_writer = get_writer(std::to_string(contig_id) + "-CLIP.bam", header);

    samFile* ow_writer = get_writer(std::to_string(contig_id) + "-OW.bam", header);
    writers[DISC_TYPES.OW] = ow_writer;
    samFile* li_writer = get_writer(std::to_string(contig_id) + "-LI.bam", header);
    writers[DISC_TYPES.LI] = li_writer;
    samFile* si_writer = get_writer(std::to_string(contig_id) + "-SI.bam", header);
    writers[DISC_TYPES.SI] = si_writer;

    std::deque<bam1_t*> two_way_buffer, forward_buffer;

    int i = 0;
    while (sam_itr_next(bam_file, iter, read) >= 0 && i < MAX_BUFFER_SIZE-1) {
        if (is_valid(read, false)) {
            bam1_t *read2 = bam_dup1(read);
            add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
            add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
            i++;
        }
    }

    while (!forward_buffer.empty()) {
        while (sam_itr_next(bam_file, iter, read) >= 0) {
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

        int64_t mq = get_mq(read);
        if (read->core.qual < MIN_MAPQ && mq < MIN_MAPQ) continue;

        // clipped read
        if (read->core.qual >= MIN_MAPQ && is_clipped(read) && check_SNP(read, two_way_buffer, config.avg_depth)) {
            int ok = sam_write1(clip_writer, header, read);
            if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
        }

        // we accept one of the mates having mapq 0 only if they are on different chromosomes
        if (read->core.qual < MIN_MAPQ || mq < MIN_MAPQ || is_dc_pair(read) || is_mate_unmapped(read)) {
            continue;
        }

        std::string qname = std::string(bam_get_qname(read));
        if (pairs.count(qname) && check_SNP(read, two_way_buffer, config.avg_depth)) {
            auto p = pairs[qname];
            int ok = sam_write1(writers[p.second], header, is_first_in_pair(read, p.second) ? read : p.first);
            if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
            bam_destroy1(p.first);
            pairs.erase(qname);
        }

        if (read->core.pos < read->core.mpos && is_outward(read, config.min_is)) {
            // only check ss for leftmost in pair - if second, then it should already be stored in the map "pairs"
            if (check_SNP(read, two_way_buffer, config.avg_depth)) {
                pairs[qname] = std::make_pair(bam_dup1(read), DISC_TYPES.OW);
            }
        } else if (read->core.pos < read->core.mpos && read->core.isize > config.max_is
                   && read->core.mpos-bam_endpos(read) > config.min_is-2*config.read_len && is_inward(read, config.min_is)) {
            if (check_SNP(read, two_way_buffer, config.avg_depth)) {
                pairs[qname] = std::make_pair(bam_dup1(read), DISC_TYPES.LI);
            }
        } else if (read->core.pos < read->core.mpos && read->core.isize < config.min_is && is_inward(read, config.min_is)) {
            if (check_SNP(read, two_way_buffer, config.avg_depth)) {
                pairs[qname] = std::make_pair(bam_dup1(read), DISC_TYPES.SI);
            }
        }
    }

    sam_close(clip_writer);

    sam_close(ow_writer);
    sam_close(si_writer);
    sam_close(li_writer);

    sam_close(bam_file);
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];
    workdir = std::string(argv[2]);
    std::string workspace = workdir + "/workspace";

    config = parse_config(workdir + "/config.txt");

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_name2id[contig_name] = contig_id;
    }

    ctpl::thread_pool thread_pool(config.threads);

    samFile* bam_file = sam_open(bam_fname.c_str(), "r");
    if (bam_file == NULL) {
        std::cerr << "Unable to open BAM file." << std::endl;
        return -1;
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    tid_to_contig_id.resize(header->n_targets);
    std::vector<std::future<void> > futures;
    for (int i = 0; i < header->n_targets; i++) {
        std::future<void> future = thread_pool.push(categorize, header->target_name[i], bam_fname, header->target_len[i]);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cout << s << std::endl;
        }
    }
}
