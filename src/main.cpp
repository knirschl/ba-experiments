//
// Created by knirschl on 18.04.23.
//

#include <iostream>
#include <array>
#include <memory>
#include "NeighborJoining.h"

int main() {
    // example from de.wiki
    std::array<std::array<double, 4>, 4> distMatrixDE {
            std::array<double, 4> {0, 3, 14, 12},
            std::array<double, 4> {3, 0, 13, 11},
            std::array<double, 4> {14, 13, 0, 4},
            std::array<double, 4> {12, 11, 4, 0}};
    std::array<std::shared_ptr<Node>, 4> treesDE{std::make_shared<Leaf>("mensch"),
                                                 std::make_shared<Leaf>("maus"),
                                                 std::make_shared<Leaf>("rose"),
                                                 std::make_shared<Leaf>("tulpe")};

    std::shared_ptr<Node> treeDE = neighborJoining<4>(distMatrixDE, treesDE);

    // example from en.wiki
    std::array<std::array<double, 5>, 5> distMatrixEN {
            std::array<double, 5> {0, 5, 9, 9, 8},
            std::array<double, 5> {5, 0, 10, 10, 9},
            std::array<double, 5> {9, 10, 0, 8, 7},
            std::array<double, 5> {9, 10, 8, 0, 3},
            std::array<double, 5> {8, 9, 7, 3, 0}};
    std::array<std::shared_ptr<Node>, 5> treesEN{std::make_shared<Leaf>("a"),
                                                 std::make_shared<Leaf>("b"),
                                                 std::make_shared<Leaf>("c"),
                                                 std::make_shared<Leaf>("d"),
                                                 std::make_shared<Leaf>("e")};

    std::shared_ptr<Node> treeEN = neighborJoining<5>(distMatrixEN, treesEN);

    // output
    std::cout << "Neighbor-joined treeDE: " << treeDE->toString() << std::endl;
    std::cout << "Neighbor-joined treeEN: " << treeEN->toString() << std::endl;
    return 0;
}
