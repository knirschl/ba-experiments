using dist_t = double;
/**
 * TODO mad rerooting
 * TODO get_pairs_to_modify -> lca (efficient for many queries: https://www.baeldung.com/cs/tree-lowest-common-ancestor)
 *
 * TODO NJ/parse rewrite to support this tree layout
 * TODO parse mapping
 *
 * TODO encapsulate in struct Tree/TaggedTree/AnnotatedTree?
 */


/*
 * https://github.com/SchadePaul/MasterThesis/blob/master/njst/src/tree.c
 * https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp
 */

#include <vector>
#include <unordered_map>
#include <string>
#include <format>
#include <limits>

std::vector<std::string> idx2leafname{}; // map node array index to locus name ("12_0_1")
std::vector<int> leaf_indices{}; // leaf array indices
std::unordered_map<std::string, std::string> leafname2groupname{}; // map locus name to species name
std::unordered_map<std::string, int> groupname2id{}; // map species name ("12") to id

/**
 * Implementation of a dynamic bitset by ASTER.
 * @link https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp#L55
 */
struct DynamicBitset {
private:
    int size{};
    std::vector<u_int64_t> vec;

public:
    DynamicBitset() {}

    DynamicBitset(int sz) : size(sz), vec((sz + 63) / 64) {}

    /**
     * Set the i-th bit.
     *
     * @param i position of the bit to set
     */
    void set(int i) {
        if (i >= size) {
            size = i + 1;
            if ((size + 63) / 64 > vec.size()) {
                vec.resize((size + 63) / 64);
            }
        }
        vec[i / 64] |= (1LL << (i % 64));
    }

    /**
     * ORs two bitsets without changing them. The smaller bitset gets virtually padded with zeros at
     * the higher bits.
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the bigger one
     */
    DynamicBitset operator|(const DynamicBitset &b) const {
        if (size < b.size) {
            return b | *this;
        }
        DynamicBitset res(size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] | b.vec[i];
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            res.vec[i] = vec[i];
        }
        return res;
    }

    /**
     * ANDs two bitsets without changing them. The bigger bitset gets virtually truncated at the
     * higher bits.
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the smaller one
     */
    DynamicBitset operator&(const DynamicBitset &b) const {
        if (size < b.size) {
            return b & *this;
        }
        DynamicBitset res(b.size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] & b.vec[i];
        }
        return res;
    }

    /**
     * XORs two bitsets without changing them. All bits at greater positions than the size of the
     * smaller bitset are copied from the bigger bitset.
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the bigger one
     */
    DynamicBitset operator^(const DynamicBitset &b) const {
        if (size < b.size) {
            return b ^ *this;
        }
        DynamicBitset res(size);
        for (int i = 0; i < b.vec.size(); i++) {
            res.vec[i] = vec[i] ^ b.vec[i];
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            res.vec[i] = vec[i];
        }
        return res;
    }

    /**
     * Subtracts two bitsets without changing them. Subtracting two bits from one another is defined
     * as ANDing with the NEGated subtrahend.
     * @code a - b := a & ~b
     *
     * @param b the other bitset
     * @return a new DynamicBitset of the size of the minuend (this)
     */
    DynamicBitset operator-(const DynamicBitset &b) const {
        DynamicBitset res(size);
        for (int i = 0; i < vec.size(); i++) {
            if (i < b.vec.size()) {
                res.vec[i] = vec[i] & ~b.vec[i];
            } else {
                res.vec[i] = vec[i];
            }
        }
        return res;
    }

    /**
     * NEGates a bitset without changing it.
     *
     * @return a new DynamicBitset of the size of this bitset
     */
    DynamicBitset operator~() const {
        DynamicBitset res(size);
        for (int i = 0; i < vec.size(); i++) {
            res.vec[i] = ~vec[i];
        }
        return res;
    }

    /**
     * Compares two bitsets for equality. Two bitsets are equal if they have the same bit at all
     * positions and the bigger bitset consists of only zeros at every position greater than the
     * size of the smaller bitset.
     *
     * @param b the other bitset
     * @return true if the two DynamicBitsets are equal under the above definition, false otherwise
     */
    bool operator==(const DynamicBitset &b) const {
        if (size < b.size) {
            return b == *this;
        }
        for (int i = 0; i < b.vec.size(); i++) {
            if (vec[i] != b.vec[i]) {
                return false;
            }
        }
        for (int i = b.vec.size(); i < vec.size(); i++) {
            if (vec[i] != 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Compares two bitsets for inequality. Two bitsets are inequal if they have a different bit at
     * all positions and the bigger bitset consists of only ones at every position greater than the
     * size of the smaller bitset.
     *
     * @param b the other bitset
     * @return true if the two DynamicBitsets are not equal under the above definition, false
     *          otherwise
     */
    bool operator!=(const DynamicBitset &b) const {
        return !(*this == b);
    }

    /**
     * Checks if two bitsets are disjoint to one another. Two bitsets are disjoint if at every
     * position only one or none of the bitsets has a bit set.
     *
     * @param b the other bitset
     * @return true if the two DynamicBitsets are disjoint to one another under the above
     *          definition, false otherwise
     */
    bool isDisjointTo(const DynamicBitset &b) const {
        if (size < b.size) {
            return b.isDisjointTo(*this);
        }
        for (int i = 0; i < b.vec.size(); i++) {
            if ((vec[i] & b.vec[i]) != 0) {
                return false;
            }
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

std::string get_name_or_idx(int node) {
    if (node < 0 || node >= idx2leafname.size()) {
        return std::to_string(node);
    }
    return idx2leafname[node];

}

std::string to_string(int cur) {
    if (tree[cur].is_leaf) {
        return idx2leafname.at(cur);
    }
    return std::format("{}({}, {})", get_name_or_idx(cur), to_string(tree[cur].left_child_idx),
                       to_string(tree[cur].right_child_idx));
}

std::string print_info(int cur) {
    if (tree[cur].is_leaf) {
        return std::format("LEAF {}[par= {}, dup= {}]", idx2leafname.at(cur), get_name_or_idx(tree[cur].parent_idx),
                           tree[cur].is_dup);
    }
    return std::format("NODE {}[par= {}, l= {}, r= {}, dup= {}];    {};    {}",
                       get_name_or_idx(cur), get_name_or_idx(tree[cur].parent_idx),
                       get_name_or_idx(tree[cur].left_child_idx),
                       get_name_or_idx(tree[cur].right_child_idx), tree[cur].is_dup,
                       print_info(tree[cur].left_child_idx), print_info(tree[cur].right_child_idx));
}

int make_node() {
    int node_idx = tree.size();
    tree.emplace_back();
    return node_idx;
}

int make_node(int node_idx, int left, int right) {
    if (node_idx >= tree.size() || node_idx < 0) {
        node_idx = make_node();
    }
    tree[node_idx].idx = node_idx;
    tree[node_idx].left_child_idx = left;
    tree[node_idx].right_child_idx = right;
    //tree[node_idx].covered_groups = tree[left].covered_groups | tree[right].covered_groups;

    tree[left].parent_idx = node_idx;
    tree[right].parent_idx = node_idx;

    return node_idx;
}

int make_node(int left, int right) {
    return make_node(-1, left, right);
}

auto make_leafs(const std::vector<std::string> &leafnames) {
    tree = std::vector<Node>{leafnames.size()};
    for (int leaf_idx{}; leaf_idx < leafnames.size(); leaf_idx++) {
        tree[leaf_idx].idx = leaf_idx;
        tree[leaf_idx].is_leaf = true;
        tree[leaf_idx].score = 0;
        if (idx2leafname.size() <= leaf_idx) {
            idx2leafname.resize(leaf_idx + 1);
        }
        idx2leafname[leaf_idx] = leafnames.at(leaf_idx);
        leaf_indices.push_back(leaf_idx);
        // give associated groupname an id or use existing one to set bitset
        std::string groupname = leafname2groupname.at(leafnames.at(leaf_idx));
        if (!(groupname2id.contains(groupname))) {
            groupname2id.emplace(groupname, groupname2id.size());
        }
        tree[leaf_idx].covered_groups.set(groupname2id.at(groupname));
    }

    return tree;
}

auto make_tree(const std::vector<std::string>& leafnames) {
    tree = make_leafs(leafnames);
    int i{};
    while (i + 1 < tree.size()) {
        make_node(i, i + 1);
        i += 2;
    }

    return tree;
}

void reverse_branch(int new_parent, bool overwrite_left) {
    Node &new_parent_node{tree[new_parent]};
    int old_parent{new_parent_node.parent_idx};
    Node &old_parent_node{tree[old_parent]};
    if (old_parent_node.parent_idx == -1) {
        int old_sibling{
                old_parent_node.left_child_idx == new_parent ? old_parent_node.right_child_idx
                                                             : old_parent_node.left_child_idx};
        tree[old_sibling].parent_idx = new_parent;
        if (overwrite_left) {
            new_parent_node.left_child_idx = old_sibling;
        } else {
            new_parent_node.right_child_idx = old_sibling;
        }
        return;
    }
    if (overwrite_left) {
        new_parent_node.left_child_idx = old_parent;
    } else {
        new_parent_node.right_child_idx = old_parent;
    }
    reverse_branch(old_parent, old_parent_node.left_child_idx == new_parent);
    old_parent_node.parent_idx = new_parent;
}

int reroot(int root_child) {
    int old_parent{tree[root_child].parent_idx};
    if (old_parent == -1) {
        // already root
        return root_child;
    }
    int root{make_node()};
    bool overwrite_left{tree[old_parent].left_child_idx == root_child};
    reverse_branch(old_parent, overwrite_left);
    if (overwrite_left) {
        make_node(root, root_child, old_parent);
    } else {
        make_node(root, old_parent, root_child);
    }

    return root;
}

auto make_quadratic_double_matrix(const size_t size) {
    std::vector<std::vector<double>> mat{size};
    for (auto &row: mat) {
        row.resize(size);
    }
    return mat;
}

void node_to_node_dist(int root, std::vector<std::vector<double>> &dist_mat) {

}

void leaf_to_leaf_relative_deviation(int cur, std::vector<std::vector<double>> const &dist_mat,
                                     std::vector<std::vector<double>> &rel_dev_mat) {
    // TODO
}

int rms(int root, int cur, std::vector<std::vector<double>> const &dist_mat,
        std::vector<std::vector<double>> const &rel_dev_mat, double &rho) {
    // TODO
    return 0;
}

void mad_root(int root, int top_id, double top_rho) {
    // TODO
}

/**
 * @link https://github.com/SchadePaul/MasterThesis/blob/master/njst/src/tree.c#L857
 *
 * @param root
 * @return
 */
int mad(int root) {
    std::vector<std::vector<double>> nn_dist_matrix{make_quadratic_double_matrix(tree.size())};
    std::vector<std::vector<double>> ll_rel_dev_matrix{
            make_quadratic_double_matrix(leaf_indices.size())};

    // calculate all node to node distance
    node_to_node_dist(root, nn_dist_matrix);
    // calculate relative deviation
    leaf_to_leaf_relative_deviation(root, nn_dist_matrix, ll_rel_dev_matrix);

    // get scores for all possible roots

    double score{};
    double best_score{std::numeric_limits<double>::max()};
    double rho{};
    double best_rho{2};
    int best_idx{-1};

    int cur{root};
    while (1) {
        while (tree[cur].left_child_idx >= 0) {
            cur = tree[cur].left_child_idx;
            score = rms(root, cur, nn_dist_matrix, ll_rel_dev_matrix, rho);
            if (score < best_score) {
                best_score = score;
                best_rho = rho;
                best_idx = cur;
            }
        }
        while (tree[cur].parent_idx >= 0 && tree[tree[cur].parent_idx].right_child_idx == cur) {
            cur = tree[cur].parent_idx;
        }
        if (cur == root) {
            break;
        }
        cur = tree[tree[cur].parent_idx].right_child_idx;
        score = rms(root, cur, nn_dist_matrix, ll_rel_dev_matrix, rho);
        if (score < best_score) {
            best_score = score;
            best_rho = rho;
            best_idx = cur;
        }
    }

    mad_root(root, best_idx, best_rho);
    return root;
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
    return tree[cur].score;

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

int root_APro() {
    int best_score{std::numeric_limits<int>::max()};
    int best_root{-1};
    int best_counter{};

    for (int root{}; root < tree.size(); root++) {
        int score{tag_APro(root)};
        if (score == best_score) {
            best_counter++;
            if (rand() % best_counter == 0) {
                best_root = root;
            }
        } else if (score < best_score) {
            best_score = score;
            best_root = root;
            best_counter = 1;
        }
    }

    //return best_root;
    return reroot(best_root);
}

/**
 * Computes the depth of the node in the tree.
 * @link https://www.geeksforgeeks.org/lowest-common-ancestor-in-a-binary-tree-using-parent-pointer/
 * @param node
 * @return
 */
int depth(int node) {
    int d{-1};
    while (node >= 0) {
        d++;
        node = tree[node].parent_idx;
    }

    return d;
}

/**
 * Easy but not efficient for many queries Lowest Common Ancestor algorithm.
 * Time: O(h)
 * Space: O(1)
 * @link https://www.geeksforgeeks.org/lowest-common-ancestor-in-a-binary-tree-using-parent-pointer/
 *
 * @param low_node
 * @param high_node
 * @return
 */
int lca(int low_node, int high_node) {
    // find difference of depths of the two nodes
    int depth_difference{depth(low_node) - depth(high_node)};
    // if high_node is deeper, swap low_node and high_node
    if (depth_difference < 0) {
        int tmp{low_node};
        low_node = high_node;
        high_node = tmp;
        depth_difference *= -1;
    }

    // move low_node up until it reaches the same level as high_node
    while (depth_difference--) {
        low_node = tree[low_node].parent_idx;
    }

    while (low_node >= 0 && high_node >= 0) {
        if (low_node == high_node) {
            return low_node;
        }
        low_node = tree[low_node].parent_idx;
        high_node = tree[high_node].parent_idx;
    }

    // no common ancestor
    return -1;
}

/**
 * Computes all leaf pairs whose lowest common ancestor was tagged as a speciation event.
 *
 * @return
 */
std::vector<std::pair<int, int>> get_speciation_pairs() {
    std::vector<std::pair<int, int>> pairs{};
    for (int i{}; i < leaf_indices.size(); i++) {
        int leaf1_id{tree[leaf_indices[i]].idx};
        for (int j{i + 1}; j < leaf_indices.size(); j++) {
            int leaf2_id{tree[leaf_indices[j]].idx};
            if (!tree[lca(leaf1_id, leaf2_id)].is_dup) {
                pairs.emplace_back(leaf1_id, leaf2_id);
            }
        }
    }

    return pairs;
}

#include <iostream>
#include <sstream>

std::string to_string(std::vector<std::pair<int, int>> ps) {
    std::ostringstream oss;
    for (auto& p : ps) {
        oss << std::format("[{},{}], ", p.first, p.second);
    }

    return oss.str();
}

/** ######
 *  TESTS
 * #######
 * make                             : Correct
 * mad                              : not tested
 *  node_to_node_dist               : not tested
 *  leaf_to_leaf_relative_deviation : not tested
 *  rms                             : not tested
 *  mad_root                        : not tested
 * tag_APro                         : not tested
 * root_APro                        : not tested
 * reroot                           : Correct
 *  reverse_branch                  : Correct
 * lca                              : Correct
 *  depth                           : Correct
 * get_speciation_pairs             : Correct
 *
 */
int main() {
    leafname2groupname = {{"a", "A"},
                          {"b", "B"},
                          {"c", "C"},
                          {"d", "A"},
                          {"e", "E"},};
    std::vector<std::string> leafnames = {"a", "b", "c", "d", "e"};
    tree = make_tree(leafnames);
    std::cout << to_string(tree.size() - 1);
    std::cout << "\n" << print_info(tree.size() - 1);
    std::cout << "\n" << to_string(get_speciation_pairs());
    tree[8].is_dup = true;
    std::cout << "\n" << to_string(get_speciation_pairs());
    tree[7].is_dup = true;
    std::cout << "\n" << to_string(get_speciation_pairs());
}