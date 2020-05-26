#include <iostream>
#include <fstream>
#include <vector>
#include <htslib/sam.h>

#include "libs/cptl_stl.h"
#include "config.h"
#include "sam_utils.h"
#include "cluster.h"

config_t config;
std::mutex mtx;
std::string workdir;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;

bam_hdr_t* header;

void categorize(int id, int contig_id, std::string& clip_fname) {
    mtx.lock();
    std::cout << "Categorizing SC for " << contig_id << " (" << contig_id2name[contig_id] << ")" << std::endl;
    mtx.unlock();

    samFile* bam_file = sam_open(clip_fname.c_str(), "r");
    if (bam_file == NULL) {
        throw "Unable to open BAM file.";
    }

    hts_idx_t* idx = sam_index_load(bam_file, clip_fname.c_str());
    if (idx == NULL) {
        throw "Unable to open BAM index.";
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == NULL) {
        throw "Unable to open BAM header.";
    }

    std::ofstream li_writer(workdir + "/workspace/" + std::to_string(contig_id) + "-LI.txt");
    std::ofstream ow_writer(workdir + "/workspace/" + std::to_string(contig_id) + "-OW.txt");

    std::string contig = contig_id2name[contig_id];
    hts_itr_t* iter = sam_itr_querys(idx, header, contig.c_str());
    bam1_t* read = bam_init1();

    while (sam_itr_next(bam_file, iter, read) >= 0) {
        // clip must nearly all (>80%) aligned
        if (!is_valid(read, false) || read->core.qual < MIN_MAPQ || bam_endpos(read)-read->core.pos < read->core.l_qseq*0.8) continue;

        int anchor_contig_id, anchor_start, anchor_end; char anchor_dir; int anchor_sc_reads;
        char* seq_name = bam_get_qname(read);
        sscanf(seq_name, "%d_%d_%d_%c_%d", &anchor_contig_id, &anchor_start, &anchor_end, &anchor_dir, &anchor_sc_reads);

        auto opp_dir = [] (char dir) { return dir == 'L' ? 'R' : 'L'; };
        anchor_t a_anchor(anchor_dir, anchor_contig_id, anchor_start, anchor_end, anchor_sc_reads);
        anchor_t a_clip(bam_is_rev(read) ? anchor_dir : opp_dir(anchor_dir),
                        contig_id, read->core.pos, bam_endpos(read), 1);

        if (bam_is_rev(read)) continue; // since we are not predicting inversions and inverted transpositions

        std::ofstream* writer;
        bool anchor_first;
        disc_type_t dt;
        if ((anchor_dir == 'L' && anchor_start-bam_endpos(read) >= config.min_sv_len) ||
                   (anchor_dir == 'R' && read->core.pos-anchor_end >= config.min_sv_len)) {
            anchor_first = (anchor_dir == 'R');
            dt = DISC_TYPES.LI;
            writer = &li_writer;
        } else if ((anchor_dir == 'L' && bam_endpos(read)-anchor_start >= config.min_sv_len) ||
                   (anchor_dir == 'R' && anchor_end-read->core.pos >= config.min_sv_len)) {
            anchor_first = (anchor_dir == 'L');
            dt = DISC_TYPES.OW;
            writer = &ow_writer;
        } else continue;

        cluster_t cluster(anchor_first ? a_anchor : a_clip, anchor_first ? a_clip : a_anchor, dt, 0);
        *writer << cluster.to_str() << "\n";
    }

    li_writer.close();
    ow_writer.close();

    bam_destroy1(read);

    hts_itr_destroy(iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);

    sam_close(bam_file);
}

int main(int argc, char* argv[]) {

    workdir = std::string(argv[1]);
    std::string workspace = workdir + "/workspace";
    std::string clip_fname = workspace + "/CLIPS.sorted.bam";
    samFile* clip_file = sam_open(clip_fname.c_str(), "r");

    int code = sam_index_build(clip_fname.c_str(), 0);
    if (code != 0) {
        throw "Cannot read " + clip_fname;
    }

    header = sam_hdr_read(clip_file);
    if (header == NULL) {
        throw "Unable to open BAM header.";
    }

    config = parse_config(workdir + "/config.txt");

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
        contig_name2id[contig_name] = contig_id;
    }

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (int contig_id = 1; contig_id <= header->n_targets; contig_id++) {
        std::future<void> future = thread_pool.push(categorize, contig_id, clip_fname);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    sam_close(clip_file);
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