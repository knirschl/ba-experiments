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

struct TNode {
    std::string id;
    std::shared_ptr<TNode> par{};
    std::shared_ptr<TNode> left{};
    std::shared_ptr<TNode> right{};

    explicit TNode(std::string id) {
        this->id = std::move(id);
    }

    std::string to_string() {
        if (left == nullptr && right == nullptr) {
            // leaf
            return id;
        }
        return std::format("{}({}, {})",
                           id,
                           left->to_string(),
                           right->to_string());
    }

    std::string print_info() {
        if (left == nullptr && right == nullptr) {
            return std::format("LEAF {}[par= {}]", id, par->id);
        }
        return std::format("NODE {}[par= {}, l= {}, r= {}];    {};    {}",
                           id, par != nullptr ? par->id : "root", left->id, right->id,
                           left->print_info(), right->print_info());
    }
};

struct Tree {
    std::string id;
    std::shared_ptr<Tree> parent{};
    virtual std::string to_string() const = 0;
    virtual std::string print_info() = 0;
};

struct Leaf : Tree {

    explicit Leaf(std::string id) {
        this->id = std::move(id);
    }

    std::string to_string() const override {
        return id;
    }

    std::string print_info() override {
        return std::format("LEAF {}[par= {}]", id, parent->id);
    }
};

template<typename T>
struct Node : Tree {
    std::shared_ptr<Tree> left{};
    std::shared_ptr<Tree> right{};
    T branch_length_left;
    T branch_length_right;

    explicit Node(std::string id) {
        this->id = std::move(id);
    }

    Node(std::shared_ptr<Tree> n_left, T bl_left, std::shared_ptr<Tree> n_right, T bl_right) {
        left = std::move(n_left);
        branch_length_left = bl_left;
        right = std::move(n_right);
        branch_length_right = bl_right;
    }

    std::string to_string() const override {
        return std::format("{}({}, {})",
                           id,
                           left->to_string(),
                           right->to_string());
    }

    std::string print_info() override {
        return std::format("NODE {}[par= {}, l= {}, r= {}];    {};    {}",
                    id, parent != NULL ? parent->id : "root", left->id, right->id,
                    left->print_info(), right->print_info());
    }
};

template<typename T>
std::shared_ptr<Tree> join(std::shared_ptr<Tree> const& n_left, T const& bl_left,
                            std::shared_ptr<Tree> const& n_right, T const& bl_right) {
    Node<T> parent{n_left, bl_left, n_right, bl_right};
    return std::make_shared<Node<T>>(parent);
}

std::string to_newick(Tree const& tree) {
    return std::format("{};", tree.to_string());
}

#endif //BA_TREE_H
