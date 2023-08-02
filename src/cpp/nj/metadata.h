//
// Created by knirschl on 02.08.23.
//

#ifndef BA_METADATA_H
#define BA_METADATA_H

#include <unordered_map>
#include <string>
#include <vector>

/**
 * Metadata (mappings) that are valid globally.
 */
struct global_metadata {
    // map species name ("12") to covered_groups id
    std::unordered_map<std::string, int> groupname2id{};
    // map species name to index in species tree matrix
    std::unordered_map<std::string, int> groupname2matidx{};
} glob_mdata;

/**
 * Metadata that are only valid for one tree.
 */
struct metadata {
    // the root of the tree
    int root{};
    // map node array index to locus name ("12_0_1")
    std::vector<std::string> idx2nodename{};
    // leaf array indices
    std::vector<int> leaf_indices{};
    // map locus name to species name
    std::unordered_map<std::string, std::string> leafname2groupname{};
    // map locus name to index in alignment/corrected matrix
    std::unordered_map<std::string, int> leafname2matidx{};

    /**
     * Returns the name of a node in the tree.
     *
     * @param node  the index of the node in the node array of the tree
     * @return the name of the node or an empty string
     */
    std::string get_name(int node) const {
        if (node < 0 || node >= idx2nodename.size()) {
            return "";
        }
        return idx2nodename[node];
    }

    /**
     * Returns the name of a node in the tree or its index as a string if the
     * node has no name.
     *
     * @param node  the index of the node in the node array of the tree
     * @return the name of the node or {@param node} as string
     */
    std::string get_name_or_idx(int node) const {
        if (node < 0 || node >= idx2nodename.size()) {
            return std::to_string(node);
        }
        return idx2nodename[node];

    }
};

#endif //BA_METADATA_H
