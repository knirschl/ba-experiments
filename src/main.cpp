//
// Created by knirschl on 18.04.23.
//

#include <iostream>
#include <array>
#include "NeighborJoining.h"

int main() {
    // example from de.wiki
    std::array<std::array<double, 4>, 4> distMatrixDE {
            std::array<double, 4> {0, 3, 14, 12},
            std::array<double, 4> {3, 0, 13, 11},
            std::array<double, 4> {14, 13, 0, 4},
            std::array<double, 4> {12, 11, 4, 0}};

    Leaf mensch{"mensch"};
    Leaf maus{"maus"};
    Leaf rose{"rose"};
    Leaf tulpe{"tulpe"};
    std::array<Node*, 4> treesDE{&mensch, &maus, &rose, &tulpe};

    Node* treeDE = neighborJoining<4>(distMatrixDE, treesDE);

/*
    // example from en.wiki
    std::array<std::array<double, 5>, 5> distMatrixEN {
            std::array<double, 5> {0, 5, 9, 9, 8}, // why is std::array<> only needed in first line?
            {5, 0, 10, 10, 9},
            {9, 10, 0, 8, 7},
            {9, 10, 8, 0, 3},
            {8, 9, 7, 3, 0}};

    Leaf a{"a"};
    Leaf b{"e"};
    Leaf c{"c"};
    Leaf d{"d"};
    Leaf e{"e"};
    std::array<Node*, 5> treesEN{&a, &b, &c, &d, &e};

    Node* treeEN = neighborJoining<5>(distMatrixEN, treesEN);
    */
    // output
    std::cout << "Neighbor-joined treeDE: " << treeDE->toString() << std::endl;
    //std::cout << "Neighbor-joined treeEN: " << treeEN->toString() << std::endl;
    return 0;
}
