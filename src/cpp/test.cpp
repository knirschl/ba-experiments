//
// Created by knirschl on 29.06.23.
//

#include <iostream>
#include <algorithm>
#include "nj/tree.h"

void reverse_branch(std::shared_ptr<TNode> &new_parent,
                    std::shared_ptr<TNode> &old_parent,
                    bool overwrite_left_child) {
    if (old_parent->par == NULL) {
        // if next parent is root skip root
        std::shared_ptr<TNode> old_par_child{};
        if (old_parent->left == new_parent) {
            old_par_child = old_parent->right;
        } else {
            old_par_child = old_parent->left;
        }
        old_par_child->par = new_parent;
        if (overwrite_left_child) {
            new_parent->left = old_par_child;
        } else {
            new_parent->right = old_par_child;
        }
        return;
    }

    if (overwrite_left_child) {
        new_parent->left = old_parent;
    } else {
        new_parent->right = old_parent;
    }
    reverse_branch(old_parent, old_parent->par, old_parent->left == new_parent);
    old_parent->par = new_parent;
}

std::shared_ptr<TNode> reroot(std::shared_ptr<TNode> one,
                              std::shared_ptr<TNode> two) {
    if (one->par == NULL) {
        // one is old root
        return one;
    } else if (two->par == NULL) {
        // two is old root
        return two;
    } else if (one->par == two) {
        // two is parent of one
        bool is_one_left_child = two->left == one;
        reverse_branch(two, two->par, is_one_left_child);
    } else {
        // one is parent of two
        bool is_two_left_child = one->left == two;
        reverse_branch(one, one->par, is_two_left_child);
    }

    std::shared_ptr<TNode> new_root = std::make_shared<TNode>("00");
    new_root->left = one;
    one->par = new_root;
    new_root->right = two;
    two->par = new_root;

    return new_root;
}

std::vector<std::pair<std::shared_ptr<TNode>, std::shared_ptr<TNode>>> make_edge_list(std::shared_ptr<TNode> root) {
    std::vector<std::pair<std::shared_ptr<TNode>, std::shared_ptr<TNode>>> edge_list{};
    if (root->is_leaf()) {
        return edge_list;
    }
    edge_list.emplace_back(root, root->left);
    edge_list.emplace_back(root, root->right);
    auto left_edge_list = make_edge_list(root->left);
    std::for_each(left_edge_list.begin(), left_edge_list.end(),
                  [&edge_list](auto &edge) { edge_list.push_back(edge); });
    auto right_edge_list = make_edge_list(root->right);
    std::for_each(right_edge_list.begin(), right_edge_list.end(),
                  [&edge_list](auto &edge) { edge_list.push_back(edge); });

    return edge_list;
}

int tag(std::shared_ptr<TNode> &u) {
    if (u->is_leaf()) {
        return u->score = 0;
    }

    auto ul = u->left;
    auto ur = u->right;
    u->species_map = ul->species_map | ur->species_map;
    u->score = tag(ul) + tag(ur);
/*
    if (ul->species_map == ur->species_map) {
        // N(l) = N(r)
        u->dup = true;
        return u->score += 1;
    } else if ((ul->species_map & ur->species_map) == std::min(ul->species_map, ur->species_map)) {
        // N(l) subset N(r) or N(r) subset N(l)
        u->dup = true;
        return u->score += 2;
    } else if (ul->species_map & ur->species_map) {
        // N(l) cut N(r) != 0
        u->dup = true;
        return u->score += 3;
    } else {
        u->dup = false;
        return u->score;
    }*/


    if ((ul->species_map & ur->species_map) == 0) {
        // speciation
        u->dup = false;
        return u->score;
    }

    u->dup = true;
    if (ul->species_map == u->species_map || ur->species_map  == u->species_map) {
        if (ul->species_map == ur->species_map) {
            return u->score += 1;
        }
        return u->score += 2;
    }
    return u->score += 3;
}

/**
 * ASTAL-Pro implementation.
 * @param root
 * @return
 */
std::shared_ptr<TNode> tag_and_root(std::shared_ptr<TNode> &root) {
    int s{INT_MAX};
    std::shared_ptr<TNode> r{root};
    for (auto &edge: make_edge_list(r)) {
        std::shared_ptr<TNode> re = reroot(edge.first, edge.second);
        int se = tag(re);
        if (se < s) {
            r = re;
            s = se;
        }
    }
    reroot(r->left, r->right);
    tag(r);
    return r;
}

int main() {
    std::shared_ptr<TNode> a = std::make_shared<TNode>("a");
    std::shared_ptr<TNode> b = std::make_shared<TNode>("b");
    std::shared_ptr<TNode> d = std::make_shared<TNode>("d");
    std::shared_ptr<TNode> c = std::make_shared<TNode>("c");
    c->species_map = 1 << 1;
    std::shared_ptr<TNode> g1 = std::make_shared<TNode>("g");
    g1->species_map = 1 << 2;
    std::shared_ptr<TNode> f = std::make_shared<TNode>("f");
    f->species_map = 1 << 3;
    std::shared_ptr<TNode> g2 = std::make_shared<TNode>("g");
    g2->species_map = 1 << 2;

    a->left = b;
    b->par = a;
    a->right = c;
    c->par = a;
    b->left = d;
    d->par = b;
    b->right = g1;
    g1->par = b;
    d->left = f;
    f->par = d;
    d->right = g2;
    g2->par = d;

    std::cout << "##############\n  Start tree\n##############\n" << a->to_string() << "\n";
    std::cout << a->print_info() << "\n\n";

    tag(a);
    std::cout << "#################\n  tagged tree\n#################\n" << a->to_string() << "\n"
              << a->print_info() << "\n\n";

    std::shared_ptr<TNode> rerooted = tag_and_root(a);
    std::cout << "#################\n  tag-rooted tree\n#################\n" << rerooted->to_string() << "\n"
              << rerooted->print_info() << "\n\n";

    return 0;
}