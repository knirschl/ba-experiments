#include <utility>
#include <array>

//
// Created by knirschl on 18.04.23.
//

#ifndef BA_NEIGHBORJOINING_H
#define BA_NEIGHBORJOINING_H

class Node {
public:
    virtual std::string toString() {
        return "Abstract Node";
    }
};

class Leaf : public Node {
private:
    std::string m_name;

public:
    explicit Leaf(std::string name) {
        m_name = std::move(name);
    }
    std::string toString() override {
        return m_name;
    }
};

class Tree : public Node {
private:
    Node* m_n1 {};
    double m_branchLength1 {};
    Node* m_n2 {};
    double m_branchLength2 {};

public:
    Tree(Node* n1, double bl1, Node* n2, double bl2) {
        m_n1 = n1;
        m_branchLength1 = bl1;
        m_n2 = n2;
        m_branchLength2 = bl2;
    }

    std::string toString() override {
        return "[(" + m_n1->toString() + "," + std::to_string(m_branchLength1) + "),("
            + m_n2->toString() + "," + std::to_string(m_branchLength2) + ")]";
    }
};

template<std::size_t N>
Node* neighborJoining(std::array<std::array<double, N>, N> distMatrix, std::array<Node*, N> trees);

#endif //BA_NEIGHBORJOINING_H
