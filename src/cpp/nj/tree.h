//
// Created by knirschl on 02.05.23.
//

#ifndef BA_TREE_H
#define BA_TREE_H

#include <omp.h>
#include <vector>
#include <string>
#include <format>
#include "../misc/DynamicBitset.h"
#include "../misc/meta.h"
#include "../misc/common_types.h"

/**
 * Representation of a tree as a vector.
 */
struct Tree {
    /*
     * Public and private definitions
     */
private:
    /**
     * Internal representation of a node.
     */
    struct Node {
        /**
         * Indices in `tree` (node array)
         */
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
public:
    std::vector<Node> tree;
    metadata mdata;

private:

    /**
     * Adds an "empty" node to the tree.
     *
     * @return the index of the new node in the tree vector
     */
    int make_node() {
        int node_idx = tree.size();
        tree.emplace_back();
        tree[node_idx].idx = node_idx;
        return node_idx;
    }

    /**
     * Adds a node to the tree with `left` and `right` as its children.
     *
     * @param left  the left child
     * @param right the right child
     * @return the index of the new node in the tree vector
     */
    int make_node(int left, int right) {
        return make_node(-1, left, -1, right, -1);
    }

    /**
     * Resets the node at index `node_idx` or adds a new one at the next free.
     *
     * @param node_idx  the preferred index of the new node
     * @param left      the left child
     * @param bl_left   the branch length from the left child to the new node
     * @param right     the right child
     * @param bl_right  the branch length from the right child to the new node
     * @return the index of the new node in the tree vector
     */
    int make_node(int node_idx, int left, dist_t bl_left, int right, dist_t bl_right) {
        if (node_idx >= tree.size() || node_idx < 0) {
            node_idx = make_node();
        }
        tree[node_idx].idx = node_idx;
        tree[node_idx].left_child_idx = left;
        tree[node_idx].right_child_idx = right;
        // TODO needed?
        //tree[node_idx].covered_groups = tree[left].covered_groups | tree[right].covered_groups;

        tree[left].parent_idx = node_idx;
        tree[left].branch_length = bl_left;
        tree[right].parent_idx = node_idx;
        tree[right].branch_length = bl_right;

        return node_idx;
    }

    /**
     * Builds a tree from a set of leafs. Overwrites `Tree::tree`.
     *
     * @param leafnames the names of the leafs
     */
    void make_tree(const std::vector<std::string> &leafnames) {
        make_leafs(leafnames);
        int i{};
        while (i + 1 < tree.size()) {
            mdata.root = make_node(i, i + 1);
            i += 2;
        }
    }

    static vector_t<std::string> nwk_split(const std::string &str) {
        vector_t<std::string> vec;
        size_t from{};
        for (size_t i{}; i < str.length(); i++) {
            if (str[i] == '\t' || str[i] == ' ' || str[i] == ';') {
                continue;
            } else if (str[i] == '(' || str[i] == ')' || str[i] == ':' || str[i] == ',') {
                if (size_t count = i - from) {
                    vec.emplace_back(str.substr(from, count));
                }
                vec.push_back(std::string() + str[i]);
                from = i + 1;
            }
        }
        if (size_t count = str.length() - 1 - from) { // -1 to skip ';'
            vec.emplace_back(str.substr(from, count));
        }

        return vec;
    }

    void nwk_name_bl(const vector_t<std::string> &vec, int idx, int &pos) {
        if (vec[pos] != ":") {
            if (mdata.idx2nodename.size() <= idx) {
                mdata.idx2nodename.resize(pos + 1);
            }
            mdata.idx2nodename[idx] = vec[pos++];
        }
        if (vec[pos] == ":") {
            tree[idx].branch_length = stod(vec[pos + 1]);
            pos += 2;
        }
    }

    /**
     * @link https://github.com/ANGSD/newick
     * @param vec
     * @param pos
     * @return
     */
    int nwk_parse(const vector_t<std::string> &vec, int &pos) {
        int ndi{make_node()};
        if(pos > vec.size()) {
            return ndi;
        }
        // leaf case
        if (vec[pos] != "(") {
            tree[ndi].is_leaf = true;
            tree[ndi].score = 0;
            nwk_name_bl(vec, ndi, pos);
            mdata.leaf_indices.push_back(ndi);
            // give associated groupname an id or use existing one to set bitset
            std::string groupname = mdata.leafname2groupname[mdata.idx2nodename[ndi]];
            if (!(glob_mdata.groupname2id.contains(groupname))) {
                glob_mdata.groupname2id.emplace(groupname, glob_mdata.groupname2id.size());
            }
            tree[ndi].covered_groups.set(glob_mdata.groupname2id[groupname]);
            return ndi;
        }
        // recursive cases
        if (vec[pos] == "(") {
            while (1) {
                int nndi = nwk_parse(vec, ++pos);
                tree[nndi].parent_idx = ndi;
                if (tree[ndi].left_child_idx == -1) {
                    tree[ndi].left_child_idx = nndi;
                } else {
                    tree[ndi].right_child_idx = nndi;
                }
                // as only two children are supported this could be moved to if (continue) and else (break)
                // this way, future expansion is a bit easier
                if (vec[pos] == ",") {
                    continue;
                }
                else {
                    break;
                }
            }
        }
        // no name or branch length
        if (vec[pos] == ")") {
            if (vec[++pos] == ")")
                return ndi;
        }
        nwk_name_bl(vec, ndi, pos);
        return ndi;
    }

    void make_tree_nwk(const std::string &newick) {
        vector_t<std::string> vec{nwk_split(newick)};
        std::cout << vec_to_string(vec) << "\n";
        mdata.idx2nodename.resize(vec.size() / 4); // probably enough
        int z{0};
        nwk_parse(vec, z);
    }

public:
    /**
     * Constructs an empty tree.
     */
    Tree() : tree{} {}

    /**
     * Constructs a tree from a set of leafs.
     *
     * @param leafnames the names of the leafs
     */
    explicit Tree(const std::vector<std::string> &leafnames, metadata mdata) : mdata{
            std::move(mdata)} {
        make_tree(leafnames);
    }

    explicit Tree(const std::string &newick, metadata mdata) : mdata{std::move(mdata)} {
        make_tree_nwk(newick);
    }

    explicit Tree(metadata mdata) : mdata{std::move(mdata)} {}

    /**
     * Creates n = leafnames.size() many leaf nodes and set mappings accordingly.
     * Overwrites `Tree::tree`.
     *
     * @param leafnames names of the leafs
     */
    void make_leafs(const std::vector<std::string> &leafnames) {
        const size_t size_ln{leafnames.size()};
        // reset structs
        mdata.leaf_indices.clear();
        mdata.idx2nodename.clear();
        mdata.idx2nodename.resize(size_ln);
        tree = std::vector<Node>{size_ln};
        mdata.root = -1;

        for (int leaf_idx{}; leaf_idx < size_ln; leaf_idx++) {
            tree[leaf_idx].idx = leaf_idx;
            tree[leaf_idx].is_leaf = true;
            tree[leaf_idx].score = 0;
            //if (idx2nodename.size() <= leaf_idx) {
            //    idx2nodename.resize(leaf_idx + 1);
            //}
            mdata.idx2nodename[leaf_idx] = leafnames[leaf_idx];
            mdata.leaf_indices.push_back(leaf_idx);
            // give associated groupname an id or use existing one to set bitset
            std::string groupname = mdata.leafname2groupname[leafnames[leaf_idx]];
//#pragma omp critical
            {
                if (!(glob_mdata.groupname2id.contains(groupname))) {
                    glob_mdata.groupname2id.emplace(groupname, glob_mdata.groupname2id.size());
                }
            };
            tree[leaf_idx].covered_groups.set(glob_mdata.groupname2id[groupname]);
        }
    }

    /**
     * Adds a node to the tree with `left` and `right` as its children.
     *
     * @param left      the left child
     * @param bl_left   the branch length from the left child to the new node
     * @param right     the right child
     * @param bl_right  the branch length from the right child to the new node
     * @return the index of the new node in the tree vector
     */
    int make_node(int left, dist_t bl_left, int right, dist_t bl_right) {
        return make_node(-1, left, bl_left, right, bl_right);
    }

    /**
     * Reverses all branches from `new_parent` upwards until the root. The parent of `new_parent`
     * will be set as the left or right child of `new_child` depending on the `overwrite_left` flag.
     * The parent of the parent of `new_parent` will be set to `new_parent` and the root is skipped.
     *
     * @param new_parent        the new parent (old child)
     * @param overwrite_left    set the old parent as left or right child of the new parent
     */
    void reverse_branch(int new_parent, bool overwrite_left) {
        Node &new_parent_node{tree[new_parent]};
        int old_parent{new_parent_node.parent_idx};
        Node &old_parent_node{tree[old_parent]};
        if (old_parent_node.parent_idx == -1) {
            int old_sibling{
                    old_parent_node.left_child_idx == new_parent ? old_parent_node.right_child_idx
                                                                 : old_parent_node.left_child_idx};
            tree[old_sibling].parent_idx = new_parent;
            tree[old_sibling].branch_length += new_parent_node.branch_length;
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
        old_parent_node.branch_length = new_parent_node.branch_length;
    }

    /**
     * Inserts a root along the branch (`root_child`, parent(`root_child`)).
     * If `root_child` or its parent is already the root, nothing changes. Else, a new root node is
     * inserted.
     *
     * @param root_child one child of the new root
     * @return the index of the new root in the tree vector
     */
    int reroot(int root_child) {
        int old_parent{tree[root_child].parent_idx};
        if (old_parent == -1) {
            // already root
            return root_child;
        }
        int old_parent_parent{tree[old_parent].parent_idx};
        if (old_parent_parent == -1) {
            // old_parent is root
            return old_parent;
        }
        int new_root{make_node()};
        bool overwrite_left{tree[old_parent].left_child_idx == root_child};
        reverse_branch(old_parent, overwrite_left);
        dist_t new_bl{tree[root_child].branch_length / 2};
        if (overwrite_left) {
            make_node(new_root, root_child, new_bl, old_parent, new_bl);
        } else {
            make_node(new_root, old_parent, new_bl, root_child, new_bl);
        }

        mdata.root = new_root;
        return mdata.root;
    }

    /**
     *
     * https://github.com/SchadePaul/MasterThesis/blob/master/njst/src/tree.c
     * https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp
     *
     * @param cur
     * @return
     */
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
            tree[cur].is_dup = false;
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

    /**
     *
     * https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp
     *
     * @return
     */
    int reroot_APro() {
        int best_score{std::numeric_limits<int>::max()};
        int best_root{-1};
        int best_counter{};

        std::vector<Node> backup_tree{tree};

        for (int test_root{}; test_root < tree.size(); test_root++) {
            // TODO Make it work without rerooting the tree all the time
            reroot(test_root);
            int score{tag_APro(test_root)};
            /*if (score == best_score) {
                best_counter++;
                if (rand() % best_counter == 0) {
                    best_root = test_root;
                }
            } else */if (score < best_score) {
                best_score = score;
                best_root = test_root;
                best_counter = 1;
            }
            tree = backup_tree;
        }

        best_root = reroot(best_root);
        tag_APro(best_root);
        return best_root;
    }

    // TODO impl
    int reroot_MAD() {
        return -1;
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
    * Computes all leaf pairs whose lowest common ancestor was not tagged as a duplication event.
    *
    * @return
    */
    std::vector<std::pair<int, int>> get_speciation_pairs() {
        std::vector<std::pair<int, int>> pairs{};
        for (int i{}; i < mdata.leaf_indices.size(); i++) {
            int leaf1_id{tree[mdata.leaf_indices[i]].idx};
            for (int j{i + 1}; j < mdata.leaf_indices.size(); j++) {
                int leaf2_id{tree[mdata.leaf_indices[j]].idx};
                if (!tree[lca(leaf1_id, leaf2_id)].is_dup) {
                    pairs.emplace_back(leaf1_id, leaf2_id);
                }
            }
        }
        return pairs;
    }

    [[nodiscard]] std::string to_string(int cur) const {
        if (tree[cur].is_leaf) {
            return mdata.idx2nodename[cur];
        }
        int lc{tree[cur].left_child_idx};
        int rc{tree[cur].right_child_idx};
        return std::format("({}:{}{}, {}:{}{})",
                           to_string(lc), tree[lc].is_leaf ? "" : mdata.get_name(lc),
                           tree[lc].branch_length,
                           to_string(rc), tree[rc].is_leaf ? "" : mdata.get_name(rc),
                           tree[rc].branch_length);
    }

    [[nodiscard]] std::string to_string() const {
        std::string out{"{}"};
        for (int i{1}; i < size(tree); i++) {
            if (tree[i].parent_idx < 0) {
            out += ", {" + to_string(i) + "}";
            }
        }
        return out;
    }

    [[nodiscard]] std::string to_newick() const {
        return std::format("{};", to_string(mdata.root));
    }

    [[nodiscard]] std::string node_info() const {
        std::string out{"{}"};
        for (int i{1}; i < size(tree); i++) {
            if (tree[i].parent_idx < 0) {
                out += ", {" + node_info(i) + "}";
            }
        }
        return out;
    }

    [[nodiscard]] std::string node_info(int cur) const {
        if (tree[cur].is_leaf) {
            return std::format("LEAF {}[par= {}, bl= {}, dup= {}]", mdata.idx2nodename[cur],
                               mdata.get_name_or_idx(tree[cur].parent_idx), tree[cur].branch_length,
                               tree[cur].is_dup);
        }
        return std::format("NODE {}[par= {}, l= {}, r= {}, bl= {}, dup= {}];    {};    {}",
                           mdata.get_name_or_idx(cur),
                           mdata.get_name_or_idx(tree[cur].parent_idx),
                           mdata.get_name_or_idx(tree[cur].left_child_idx),
                           mdata.get_name_or_idx(tree[cur].right_child_idx),
                           tree[cur].branch_length, tree[cur].is_dup,
                           node_info(tree[cur].left_child_idx),
                           node_info(tree[cur].right_child_idx));
    }
};

#endif //BA_TREE_H
