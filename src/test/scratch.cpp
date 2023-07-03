using dist_t = double;
/**
 * (TODO APro rerooting: tag_and_root, make_edge_list, reroot, reverse_branch)
 * TODO mad rerooting
 * TODO get_pairs_to_modify
 *
 * TODO NJ/parse rewrite to support this tree layout
 * TODO parse mapping
 *
 * TODO encapsulate in struct Tree/TaggedTree/AnnotatedTree?
 */


/*
 * https://github.com/SchadePaul/MasterThesis/blob/master/njst/src/node.h
 * https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp
 */

#include <vector>
#include <unordered_map>
#include <string>
#include <format>

std::vector<std::string> idx2leafname{}; // map node array index to locus name ("12_0_1")
std::unordered_map<std::string, std::string> leafname2groupname{}; // map locus name to species name
std::unordered_map<std::string, int> groupname2id{}; // map species name ("12") to id

// https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp#L55
struct DynamicBitset {
private:
    int size{};
    std::vector<u_int64_t> vec;

public:
    DynamicBitset() {}

    DynamicBitset(int sz) : size(sz), vec((sz + 63) / 64) {}

    void set(int i) {
        if (i >= size) {
            size = i + 1;
            if ((size + 63) / 64 > vec.size()) {
                vec.resize((size + 63) / 64);
            }
        }
        vec[i / 64] |= (1LL << (i % 64));
    }

    DynamicBitset operator|(const DynamicBitset &b) const {
        if (size < b.size) return b | *this;
        DynamicBitset res(size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] | b.vec[i];
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            res.vec[i] = vec[i];
        }
        return res;
    }

    DynamicBitset operator&(const DynamicBitset &b) const {
        if (size < b.size) return b & *this;
        DynamicBitset res(b.size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] & b.vec[i];
        }
        return res;
    }

    DynamicBitset operator^(const DynamicBitset &b) const {
        if (size < b.size) return b ^ *this;
        DynamicBitset res(size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] ^ b.vec[i];
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            res.vec[i] = vec[i];
        }
        return res;
    }

    DynamicBitset operator-(const DynamicBitset &b) const {
        DynamicBitset res(size);
        for (int i = 0; i < vec.size(); i++) {
            if (i < b.vec.size()) res.vec[i] = vec[i] & ~b.vec[i];
            else res.vec[i] = vec[i];
        }
        return res;
    }

    DynamicBitset operator~() const {
        DynamicBitset res(size);
        for (int i = 0; i < vec.size(); i++) {
            res.vec[i] = ~vec[i];
        }
        return res;
    }

    bool operator==(const DynamicBitset &b) const {
        if (size < b.size) return b == *this;
        for (int i = 0; i < b.vec.size(); i++) {
            if (vec[i] != b.vec[i]) return false;
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            if (vec[i] != 0) return false;
        }
        return true;
    }

    bool operator!=(const DynamicBitset &b) const {
        return !(*this == b);
    }

    bool isDisjointTo(const DynamicBitset &b) const {
        if (size < b.size) return b.isDisjointTo(*this);
        for (int i = 0; i < b.vec.size(); i++) {
            if ((vec[i] & b.vec[i]) != 0) return false;
        }
        return true;
    }
};

struct Node {
    int idx{-1};
    int left_child_idx{-1};
    int right_child_idx{-1};
    int parent_idx{-1};
    int score{-1};
    bool is_dup{};
    bool is_leaf{}; // can be incorrect
    DynamicBitset covered_groups{};
    dist_t branch_length{-1};
};

std::vector<Node> tree;

std::string to_string(int cur) {
    if (tree[cur].is_leaf) {
        return idx2leafname.at(cur);
    }
    return std::format("{}({}, {})", idx2leafname.at(cur), to_string(tree[cur].left_child_idx),
                       to_string(tree[cur].right_child_idx));
}

std::string print_info(int cur) {
    if (tree[cur].is_leaf) {
        return std::format("LEAF {}[par= {}, dup= {}]", idx2leafname.at(cur),
                           tree[cur].parent_idx != -1 ? idx2leafname.at(tree[cur].parent_idx) : "root",
                           tree[cur].is_dup);
    }
    return std::format("NODE {}[par= {}, l= {}, r= {}, dup= {}];    {};    {}", idx2leafname.at(cur),
                       tree[cur].parent_idx != -1 ? idx2leafname.at(tree[cur].parent_idx) : "root",
                       idx2leafname.at(tree[cur].left_child_idx), idx2leafname.at(tree[cur].right_child_idx),
                       tree[cur].is_dup,
                       print_info(tree[cur].left_child_idx), print_info(tree[cur].right_child_idx));
}

auto make_leafs(std::vector<std::string> leafnames) {
    tree = std::vector<Node>{leafnames.size()};
    for (auto &leafname: leafnames) {
        int leaf_idx = tree.size();
        tree.emplace_back();
        tree[leaf_idx].idx = leaf_idx;
        tree[leaf_idx].is_leaf = true;
        tree[leaf_idx].score = 0;
        idx2leafname[leaf_idx] = leafname;
        // give associated groupname an id or use existing one to set bitset
        std::string groupname = leafname2groupname.at(leafname);
        if (!(groupname2id.contains(groupname))) {
            groupname2id.emplace(groupname, groupname2id.size());
        }
        tree[leaf_idx].covered_groups.set(groupname2id.at(groupname));
    }

    return tree;
}

int tag_APro(int cur) {
    /* TODO
     * if (tree[cur].score != -1) {
        return tree[cur].score;
    }
    */
    if (tree[cur].is_leaf) {
        return tree[cur].score = 0;
    }
    tree[cur].score = tag_APro(tree[cur].left_child_idx) + tag_APro(tree[cur].right_child_idx);
    DynamicBitset &covered_left = tree[tree[cur].left_child_idx].covered_groups;
    DynamicBitset &covered_right = tree[tree[cur].right_child_idx].covered_groups;
    tree[cur].covered_groups = covered_left | covered_right;

    // from astral-pro project
    if (!covered_left.isDisjointTo(covered_right)) {
        tree[cur].is_dup = true;
        tree[cur].score++;
        if (tree[cur].covered_groups != covered_left) {
            tree[cur].score++;
        }
        if (tree[cur].covered_groups != covered_right) {
            tree[cur].score++;
        }
    } else {
        tree[cur].is_dup = false; // TODO remove?
    }
    //return tree[cur].score;

    // my implementation of pauls thesis
    /*
    if (covered_left.isDisjointTo(covered_right)) {
        // speciation
        tree[cur].is_dup = false;
        return tree[cur].score;
    }
    tree[cur].is_dup = true;
    if (covered_left == tree[cur].covered_groups || covered_right == tree[cur].covered_groups) {
        if (covered_left == covered_right) {
            return tree[cur].score++;
        }
        return tree[cur].score += 2;
    }
    return tree[cur].score += 3;
    */
}