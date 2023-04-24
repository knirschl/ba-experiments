//
// Created by knirschl on 18.04.23.
//

#include <iostream>
#include <limits>
#include "NeighborJoining.h"

// debug func
void printMatrix(const double* distMatrix, int size) {
    std::cout << "Printing matrix:" << std::endl;
    for (int i{0}; i < size; ++i) {
        for (int j{0}; j < size; ++j) {
            std::cout << *(distMatrix + i * size + j) << " ";
        }
        std::cout << std::endl;
    }
}

// debug func
void printTrees(Node** trees, int size) {
    std::cout << "Printing trees:" << std::endl;
    for (int i{0}; i < size; ++i) {
        std::cout << trees[i]->toString() << std::endl;
    }
}

/**
 * Calculate mean distance from every taxon to every taxon.
 * r_i = 1/(N-2) * sum{k=1,N}(d_i,k)
 *
 * @tparam N
 * @param distMatrix
 * @return array containing r_i at index i
 */
template<std::size_t N>
std::array<double, N> meanDistance(std::array<std::array<double, N>, N> distMatrix) {
    std::array<double, N> r;
    for (std::size_t i{0}; i < N; ++i) {
        double sumK{0};
        for (std::size_t k{0}; k < N; ++k) {
            sumK += distMatrix[i][k];
        }
        r[i] = (1.0 / (N - 2)) * sumK; // 1.0 for double division
    }

    return r;
}

/**
 * Calculate intermediate matrix.
 * M_i,j = d_i,j - (r_i + r_j)
 *
 * @tparam N
 * @param distMatrix
 * @param r 
 * @return intermediate matrix
 */
template<std::size_t N>
std::array<std::array<double, N>, N> interMatrix(std::array<std::array<double, N>, N> distMatrix,
                                                 std::array<double, N> r) {
    std::array<std::array<double, N>, N> iM;
    for (std::size_t i{0}; i < N; ++i) {
        for (std::size_t j{0}; j < N; ++j) {
            iM[i][j] = distMatrix[i][j] - (r[i] + r[j]);
        }
    }

    return iM;
}

/**
 * Find the minimum in the matrix.
 * 
 * @param matrix
 * @param size 
 * @return pair containing the indices describing the position of the first minimum
 */
template<std::size_t N>
std::pair<std::size_t, std::size_t> findMin(std::array<std::array<double, N>, N> matrix) {
    double min = std::numeric_limits<double>::infinity();
    std::pair<std::size_t, std::size_t> minDistPos;
    for (std::size_t i{0}; i < N; ++i) {
        for (std::size_t j{0}; j < N; ++j) {
            if (i == j) {
                continue;
            }
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                minDistPos.first = i;
                minDistPos.second = j;
            }
        }
    }

    return minDistPos;
}

/**
 * Calculate new distance matrix.
 * Removes rows/columns defined in minDistPos from distMatrix and adds new row/column at the end storing distances of the
 * newly merged neighbors/trees to everything else.
 *
 * @tparam N
 * @param distMatrix
 * @param minDistPos 
 * @return matrix of size size-1
 */
template<std::size_t N>
std::array<std::array<double, N - 1>, N - 1> calculateNewMatrix(std::array<std::array<double, N>, N> distMatrix,
                                                        std::pair<std::size_t, std::size_t> minDistPos) {
    std::array<std::array<double, N - 1>, N - 1> newMatrix;
    // copy old values of rows/columns that won't be changed
    std::size_t xOff{0};
    for (std::size_t i{ 0 }; i < N; ++i) {
        // skip row if one of the merged ones
        if (i == minDistPos.first || i == minDistPos.second) {
            xOff++;
        } else {
            std::size_t yOff{0};
            for (std::size_t j{0}; j < N; ++j) {
                // skip column if one of the merged ones
                if (j == minDistPos.first || j == minDistPos.second) {
                    yOff++;
                } else {
                    // index at newMatrix depends on how many rows/columns were already skipped
                    newMatrix[i - xOff][j - yOff] = distMatrix[i][j];
                }
            }
        }
    }

    // calculate new distances from merged clusters to everything else, put in last row/column
    std::size_t off{0};
    for (std::size_t k{0}; k < N - 2; ++k) {
        // calculate offset (how much needs to be skipped) via black magic
        while (k + off == minDistPos.first || k + off == minDistPos.second) {
            off++;
        }
        /*
         * d_u,k = (d_i,k + d_j,k - d_i,j) / 2
         * u = (i,j) the merged clusters, k is the current cluster
         */
        newMatrix[k][N - 2] = (distMatrix[minDistPos.first][k + off]
                                + distMatrix[minDistPos.second][k + off]
                                - distMatrix[minDistPos.first][minDistPos.second])
                                        / 2.0;
        newMatrix[N - 1][k] = (distMatrix[minDistPos.first][k + off]
                               + distMatrix[minDistPos.second][k + off]
                               - distMatrix[minDistPos.first][minDistPos.second])
                                       / 2.0;
    }
    newMatrix[N - 1][N - 1] = 0; // set new diagonal entry to 0

    return newMatrix;
}

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
Tree* join(Node* n1, double bl1, Node* n2, double bl2) {
    return new Tree {n1, bl1, n2, bl2};
}

/**
 * Joins two trees by adding a common parent.
 *
 * @param n1
 * @param n2
 * @param bl
 * @return parent Tree node with branches of length bl.first and bl.second to the to Nodes
 */
Tree* join(Node* n1, Node* n2, std::pair<double, double> bl) {
    return join(n1, bl.first, n2, bl.second);
}

/**
 * Joins two trees by adding a common parent.
 *
 * @param n1
 * @param n2
 * @param bl
 * @return parent Tree node with branches of length bl to the to Nodes
 */
Tree* join(Node* n1, Node* n2, double bl) {
    return join(n1, bl, n2, bl);
}

/**
 * Copies unjoined trees as they do not get changed.
 * 
 * @param trees 
 * @param minDistPos 
 * @param size 
 * @return array of size size-1 where the last entry is not set
 */
template<std::size_t N>
std::array<Node*, N - 1> copyUnjoinedTrees(std::array<Node*, N> trees, std::pair<std::size_t, std::size_t> minDistPos) {
    std::array<Node*, N - 1> copy;
    std::size_t off{0};
    for (std::size_t i{0}; i < N; ++i) {
        // skip and add offset if this tree is one of the two that will be merged
        if (i == minDistPos.first || i == minDistPos.second) {
            ++off;
        } else {
            copy[i - off] = trees[i];
        }
    }
    
    return copy;
}

/**
 * Neighbor-Joining after {@link https://de.wikipedia.org/wiki/Neighbor-Joining-Algorithmus}.
 *
 * @param distMatrix
 * @param trees
 * @param size
 * @return pointer to a node, in most cases this is the root of the inferred tree
 */
template<std::size_t N>
Node* neighborJoining(std::array<std::array<double, N>, N> distMatrix, std::array<Node*, N> trees) {
    //printTrees(trees, size);
    //printMatrix(distMatrix, size);
    if (N > 2) {
        // calculations
        std::array<double, N> r = meanDistance<N>(distMatrix);
        //std::cout << "r = [" << r[0] << ", " << r[1] << ", " << r[2] << ", " << r[3] << "]" << std::endl;
        //printMatrix(iM, size);
        std::pair<std::size_t, std::size_t> minDistPos = findMin<N>(interMatrix<N>(distMatrix, r));
        //std::cout << "Joining " << trees[minDistPos[0]]->toString() << "(idx=" << minDistPos[0] << ",bl=" << blMin[0] << ") and "
        //    << trees[minDistPos[1]]->toString() << "(idx=" << minDistPos[1] << ",bl=" << blMin[1]<< ")" << std::endl;
        std::array<std::array<double, N - 1>, N - 1> newDistMatrix = calculateNewMatrix<N>(distMatrix, minDistPos);

        // tree building
        std::pair<double, double> branchLengths = calculateBranchLengths(
                distMatrix[minDistPos.first][minDistPos.second],
                r[minDistPos.first],
                r[minDistPos.second]);
        std::array<Node*, N - 1> newTrees = copyUnjoinedTrees<N>(trees, minDistPos);
        // add joined trees (neighbors)
        newTrees[N - 2] = join(trees[minDistPos.first], trees[minDistPos.second], branchLengths);

        return neighborJoining<N - 1>(newDistMatrix, newTrees);
    } else if (N == 2) {
        // base case
        //std::cout << "Joining last two subtrees with distance " << distMatrix[1] << std::endl;
        // connecting branch defined by two branches of same length to parent -> "rooted" tree...
        return join(trees[0], trees[1], distMatrix[1]);
    } else {
        // edge case, only one or none clusters
        return trees[0];
    }
}
