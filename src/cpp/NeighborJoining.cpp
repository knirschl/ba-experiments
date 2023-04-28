//
// Created by knirschl on 18.04.23.
//

#include <memory>
#include <utility>
#include "NeighborJoining.h"

/**
 * Calculate the branch lengths in the inferred tree.
 * v_i,u = (d_i,j + r_i - r_j) / 2
 * v_j,u = d_i,j - v_i,u
 *
 * @param dist
 * @param r1
 * @param r2
 * @return pair containing the two branch lengths, first corresponds to branch of tree that belongs to r1
 */
std::pair<double, double> calculateBranchLengths(double dist, double r1, double r2) {
    return std::pair<double, double> {
            (dist + r1 - r2) / 2.0,
            dist - ((dist + r1 - r2) / 2.0)
    };
}

/**
 * Joins two trees by adding a common parent.
 *
 * @param n1
 * @param bl1
 * @param n2
 * @param bl2
 * @return parent Tree node with branches of length bl1 and bl2 to the to Nodes
 */
std::shared_ptr<Tree> join(const std::shared_ptr<Node>& n1, double bl1, const std::shared_ptr<Node>& n2, double bl2) {
    return std::make_shared<Tree>(n1, bl1, n2, bl2);
}
