//
// Created by knirschl on 06.07.23.
//

#ifndef BA_META_H
#define BA_META_H

#include <vector>
#include <unordered_map>
#include <string>

std::vector<std::string> idx2leafname{}; // map node array index to locus name ("12_0_1")
std::vector<int> leaf_indices{}; // leaf array indices
std::unordered_map<std::string, std::string> leafname2groupname{}; // map locus name to species name
std::unordered_map<std::string, int> groupname2id{}; // map species name ("12") to id
std::unordered_map<std::string, int> leafname2matidx{};
std::unordered_map<std::string, int> groupname2matidx{};

std::string get_name(int node) {
    if (node < 0 || node >= idx2leafname.size()) {
        return "";
    }
    return idx2leafname[node];
}

std::string get_name_or_idx(int node) {
    if (node < 0 || node >= idx2leafname.size()) {
        return std::to_string(node);
    }
    return idx2leafname[node];

}

#endif //BA_META_H
