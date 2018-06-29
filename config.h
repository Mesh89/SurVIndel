#ifndef SURVEYOR_CONFIG_H
#define SURVEYOR_CONFIG_H

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "libs/IntervalTree.h"

int MIN_MAPQ = 1;
int MIN_CLIP_LEN = 5;
int MIN_CLIP_CONSENSUS_LEN = 15;
int MIN_DC_MAPQ = 20;
double MAX_SEQ_ERROR = 0.04;

int MAX_READ_SUPPORTED = 10000;

struct config_t {
    int threads;
    std::string rmsk_fname = "", simple_rep_fname = "";
    int read_len;
    int avg_depth;
    int max_is, min_is;
    int min_sv_len, max_sc_dist;
};

struct repeat_t {
    std::string chr;
    int start, end;
    std::string type;

    repeat_t() {}
    repeat_t(std::string& line) {
        char temp[100];
        std::stringstream ss(line);

        for (int i = 0; i < 5; i++) {
            ss >> temp;
        }
        ss >> chr >> start >> end;

        ss >> temp >> temp;
        ss >> type >> temp;
        type += "-" + std::string(temp);
    }
};

struct simple_repeat_t {
    std::string chr;
    int start, end;
    int period;

    simple_repeat_t() {}
    simple_repeat_t(std::string& line) {
        char temp[100];
        std::stringstream ss(line);

        ss >> temp;
        ss >> chr >> start >> end;
        ss >> temp;
        ss >> period;
    }
};


// return a map: key is the contig name, value is an interval tree of the repeats for that contig
std::unordered_map<std::string, IntervalTree<repeat_t>* > parse_repeats(std::string fname) {
    std::unordered_map<std::string, IntervalTree<repeat_t>* > repeat_trees;
    std::ifstream rep_fin(fname);
    if (!rep_fin.is_open()) {
        std::cerr << "WARN: Could not open rmsk file" << std::endl;
    } else {
        std::unordered_map<std::string, std::vector<Interval<repeat_t> > > v;
        std::string line;
        while (getline(rep_fin, line)) {
            if (line[0] != '#') {
                repeat_t rep(line);
                v[rep.chr].push_back(Interval<repeat_t>(rep.start, rep.end, rep));
            }
        }

        for (auto& k : v) {
            repeat_trees[k.first] = new IntervalTree<repeat_t>(k.second);
        }
    }
    return repeat_trees;
}

// return a map: key is the contig name, value is an interval tree of the repeats for that contig
std::unordered_map<std::string, IntervalTree<simple_repeat_t>* > parse_simple_repeats(std::string fname, int min_sv_len) {
    std::unordered_map<std::string, IntervalTree<simple_repeat_t>* > repeat_trees;
    std::ifstream rep_fin(fname);
    if (!rep_fin.is_open()) {
        std::cerr << "WARN: Could not open simpleRepeats file" << std::endl;
    } else {
        std::unordered_map<std::string, std::vector<Interval<simple_repeat_t> > > v;
        std::string line;
        while (getline(rep_fin, line)) {
            if (line[0] != '#') {
                simple_repeat_t rep(line);
                if (rep.end-rep.start+1 < min_sv_len) continue;
                v[rep.chr].push_back(Interval<simple_repeat_t>(rep.start, rep.end, rep));
            }
        }

        for (auto& k : v) {
            repeat_trees[k.first] = new IntervalTree<simple_repeat_t>(k.second);
        }
    }
    return repeat_trees;
}

config_t parse_config(std::string file) {
    std::unordered_map<std::string, std::string> config_params;
    std::ifstream fin(file);
    std::string name, value;
    while (fin >> name >> value) {
        config_params[name] = value;
    }
    fin.close();

    config_t config;
    config.threads = stoi(config_params["threads"]);
    if (config_params.count("rmsk")) {
        config.rmsk_fname = config_params["rmsk"];
    }
    if (config_params.count("simple_rep")) {
        config.simple_rep_fname = config_params["simple_rep"];
    }
    config.read_len = stoi(config_params["read_len"]);
    config.avg_depth = stoi(config_params["avg_depth"]);
    config.min_is = stoi(config_params["min_is"]);
    config.max_is = stoi(config_params["max_is"]);
    config.min_sv_len = stoi(config_params["min_sv_len"]);
    config.max_sc_dist = stoi(config_params["max_sc_dist"]);
    return config;
};

#endif //SURVEYOR_CONFIG_H
