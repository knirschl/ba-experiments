using dist_t = double;
/**
 * TODO mad rerooting
 * TODO get_pairs_to_modify -> lca (efficient for many queries: https://www.baeldung.com/cs/tree-lowest-common-ancestor)
 *
 * TODO NJ/parse rewrite to support this tree layout
 * TODO parse mapping
 *
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