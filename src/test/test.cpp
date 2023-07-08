//
// Created by knirschl on 08.07.23.
//

#include <iostream>
#include <sstream>
#include "../cpp/nj/tree.h"

Tree make_tree() {
    leafname2groupname = {{"a", "A"},
                          {"b", "B"},
                          {"c", "C"},
                          {"d", "A"},
                          {"e", "E"},};
    Tree tree{{"a", "b", "c", "d", "e"}};
    for (auto &n: tree.tree) {
        n.branch_length = 1;
    }

    return tree;
}

template<typename U, typename V>
std::string to_string(std::pair<U, V> p) {
    return std::format("[{},{}]", p.first, p.second);
}

template<typename U>
std::string to_string(std::vector<U> vs) {
    std::ostringstream oss;
    for (auto &v: vs) {
        oss << std::format("{}, ", to_string(v));
    }

    return oss.str();
}

template<typename U, typename V>
std::string to_string(std::unordered_map<U, V> m) {
    std::ostringstream oss;
    for (auto &e: m) {
        oss << std::format("{}, ", to_string(e));
    }

    return oss.str();
}

int main() {
    Tree test_tree{make_tree()};
    std::cout << test_tree.to_newick() << "\n" << test_tree.node_info(test_tree.root) << "\n";
    test_tree.reroot_APro();
    std::cout << test_tree.to_newick() << "\n" << test_tree.node_info(test_tree.root) << "\n";
    std::cout << to_string(test_tree.get_speciation_pairs()) << "\n";

    return 0;
}