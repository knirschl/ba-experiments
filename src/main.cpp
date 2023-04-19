//
// Created by knirschl on 18.04.23.
//

#include <iostream>
#include "NeighborJoining.h"

int main() {
    // example from de.wiki
    double distMatrixDE[4 * 4] = {
            0, 3, 14, 12,
            3, 0, 13, 11,
            14, 13, 0, 4,
            12, 11, 4, 0};

    Node* treesDE[4] = {new Leaf {"Mensch"}, new Leaf {"Maus"},
                     new Leaf {"Rose"}, new Leaf {"Tulpe"}};

    Node* treeDE = neighborJoining(distMatrixDE, treesDE, 4);


    // example from en.wiki
    double distMatrixEN[5 * 5] = {
            0, 5, 9, 9, 8,
            5, 0, 10, 10, 9,
            9, 10, 0, 8, 7,
            9, 10, 8, 0, 3,
            8, 9, 7, 3, 0};

    Node* treesEN[5] = {new Leaf {"a"}, new Leaf {"b"},
                      new Leaf {"c"}, new Leaf {"d"},
                      new Leaf {"e"}};

    Node* treeEN = neighborJoining(distMatrixEN, treesEN, 5);


    // output
    std::cout << "Neighbor-joined treeDE: " << treeDE->toString() << std::endl;
    std::cout << "Neighbor-joined treeEN: " << treeEN->toString() << std::endl;
    return 0;
}
