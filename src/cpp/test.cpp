//
// Created by knirschl on 29.06.23.
//

#include <iostream>
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

int main() {
    std::shared_ptr<TNode> a = std::make_shared<TNode>("a");
    std::shared_ptr<TNode> b = std::make_shared<TNode>("b");
    std::shared_ptr<TNode> d = std::make_shared<TNode>("d");
    std::shared_ptr<TNode> c = std::make_shared<TNode>("c");
    std::shared_ptr<TNode> e = std::make_shared<TNode>("e");
    std::shared_ptr<TNode> f = std::make_shared<TNode>("f");
    std::shared_ptr<TNode> g = std::make_shared<TNode>("g");

    a->left = b;
    b->par = a;
    a->right = c;
    c->par = a;
    b->left = d;
    d->par = b;
    b->right = e;
    e->par = b;
    d->left = f;
    f->par = d;
    d->right = g;
    g->par = d;

    std::cout << "##############\n  Start tree\n##############\n" << a->to_string() << "\n";
    std::cout << a->print_info() << "\n\n";

    std::shared_ptr<TNode> rerooted = reroot(d, b);

    std::cout << "#################\n  Rerooted tree\n#################\n" << rerooted->to_string() << "\n"
              << rerooted->print_info() << "\n\n";

    return 0;
}