//
// Created by knirschl on 02.05.23.
//

#ifndef BA_TREE_H
#define BA_TREE_H

#include <string>
#include <format>
#include <memory>
#include <utility>
#include "../common_types.h"

struct NTree {
    virtual std::string to_string() const = 0;
};

struct NLeaf : NTree {
    std::string id;

    explicit NLeaf(std::string id) {
        this->id = std::move(id);
    }

    std::string to_string() const override {
        return id;
    }
};

template<typename T>
struct NNode : NTree {
    std::shared_ptr<NTree> left{};
    std::shared_ptr<NTree> right{};
    T branch_length_left;
    T branch_length_right;

    NNode(std::shared_ptr<NTree> n_left, T bl_left, std::shared_ptr<NTree> n_right, T bl_right) {
        left = std::move(n_left);
        branch_length_left = bl_left;
        right = std::move(n_right);
        branch_length_right = bl_right;
    }

    std::string to_string() const override {
        return std::format("({}:{},{}:{})",
                           left->to_string(),
                           branch_length_left,
                           right->to_string(),
                           branch_length_right);
    }
};

template<typename T>
std::shared_ptr<NTree> join(std::shared_ptr<NTree> const& n_left, T const& bl_left,
                            std::shared_ptr<NTree> const& n_right, T const& bl_right) {
    NNode<T> parent{n_left, bl_left, n_right, bl_right};
    return std::make_shared<NNode<T>>(parent);
}

std::string to_newick(NTree const& tree) {
    return std::format("{};", tree.to_string());
}

#endif //BA_TREE_H
