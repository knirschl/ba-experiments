//
// Created by knirschl on 02.05.23.
//

#ifndef BA_TREE_H
#define BA_TREE_H

#include <string>
#include <format>
#include <memory>
#include <iostream>

using dist_t = double;

struct NTree {
    virtual std::string to_string() = 0;
};

struct NLeaf : NTree {
    std::string id;

    std::string to_string() override {
        return id;
    }
};

struct NNode : NTree {
    std::shared_ptr<NTree> left{};
    std::shared_ptr<NTree> right{};
    dist_t branch_length_left;
    dist_t branch_length_right;

    std::string to_string() override {
        /* C++20
        return std::format("({}:{},{}:{})",
                           left->to_string(),
                           branch_length_left,
                           right->to_string(),
                           branch_length_right);
        */
        return "("
            + left->to_string() + ":" + std::to_string(branch_length_left)
            + "," + right->to_string() + ":" + std::to_string(branch_length_right)
            + ")";
    }
};

std::string tree_to_string(NTree& tree) {
     // C++20
    return std::format("({});", tree.to_string());


    return "(" + tree.to_string() + ")";
}

void print_tree(NTree& tree) {
    std::cout << tree_to_string(tree) << "\n";
}

#endif //BA_TREE_H
