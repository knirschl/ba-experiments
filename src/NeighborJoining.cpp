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
 * Caluclate mean distance from every taxon to every taxon.
 * r_i = 1/(N-2) * sum{k=1,N}(d_i,k)
 *
 * @param distMatrix 
 * @param size 
 * @return array containing r_i at index i
 */
double* meanDistance(const double* distMatrix, int size) {
    double* r = new double[size];
    for (int i{0}; i < size; ++i) {
        double sumK{0};
        for (int k{0}; k < size; ++k) {
            sumK += *(distMatrix + i * size + k);
        }
        r[i] = (1.0 / (size - 2)) * sumK; // 1.0 for double division
    }

    return r;
}

/**
 * Calculate intermediate matrix.
 * M_i,j = d_i,j - (r_i + r_j)
 * 
 * @param distMatrix 
 * @param size 
 * @param r 
 * @return intermediate matrix as 1D array
 */
double* interMatrix(const double* distMatrix, int size, const double* r) {
    double* iM = new double[size * size];
    for (int i{0}; i < size; ++i) {
        for (int j{0}; j < size; ++j) {
            iM[i * size + j] = distMatrix[i * size + j] - (r[i] + r[j]);
        }
    }

    return iM;
}

/**
 * Find the minimum in the matrix.
 * 
 * @param matrix
 * @param size 
 * @return array containing the indices describing the position of the first minimum
 */
int* findMin(const double* matrix, int size) {
    double min = std::numeric_limits<double>::infinity();
    int* minDistPos = new int[2];
    for (int i{0}; i < size; ++i) {
        for (int j{0}; j < size; ++j) {
            if (i == j) {
                continue;
            }
            if (matrix[i * size + j] < min) {
                min = matrix[i * size + j];
                minDistPos[0] = i;
                minDistPos[1] = j;
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
 * @param distMatrix 
 * @param size 
 * @param minDistPos 
 * @return matrix as 1D array of size size-1
 */
double* calculateNewMatrix(const double* distMatrix, int size, const int minDistPos[2]) {
    double* newMatrix = new double[(size - 1) * (size - 1)];
    // copy old values of rows/columns that won't be changed
    int xOff{0};
    for (int i{ 0 }; i < size; ++i) {
        // skip row if one of the merged ones
        if (i == minDistPos[0] || i == minDistPos[1]) {
            xOff++;
        } else {
            int yOff{0};
            for (int j{0}; j < size; ++j) {
                // skip column if one of the merged ones
                if (j == minDistPos[0] || j == minDistPos[1]) {
                    yOff++;
                } else {
                    // index at newMatrix depends on how many rows/columns were already skipped
                    newMatrix[(i - xOff) * (size - 1) + (j - yOff)] = distMatrix[i * size + j];
                }
            }
        }
    }

    // calculate new distances from merged clusters to everything else, put in last row/column
    int off{0};
    for (int k{0}; k < size - 2; ++k) {
        // calculate offset (how much needs to be skipped) via black magic
        while (k + off == minDistPos[0] || k + off == minDistPos[1]) {
            off++;
        }
        /*
         * d_u,k = (d_i,k + d_j,k - d_i,j) / 2
         * u = (i,j) the merged clusters, k is the current cluster
         */
        newMatrix[k * (size - 1) + size - 2] = (distMatrix[minDistPos[0] * size + k + off]
                                                + distMatrix[minDistPos[1] * size + k + off]
                                                - distMatrix[minDistPos[0] * size + minDistPos[1]])
                                               / 2.0;
        newMatrix[(size - 1) * (size - 2) + k] = (distMatrix[minDistPos[0] * size + k + off]
                                                  + distMatrix[minDistPos[1] * size + k + off]
                                                  - distMatrix[minDistPos[0] * size + minDistPos[1]])
                                                 / 2.0;
    }
    newMatrix[(size - 1) * (size - 1) - 1] = 0; // set new diagonal entry to 0

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
 * @return array containging the two branch lengths, index 0 corresponds to branch of tree that belongs to r1
 */
double* calculateBranchLengths(double dist, double r1, double r2) {
    double* bl = new double[2];
    bl[0] = (dist + r1 - r2) / 2.0;
    bl[1] = dist - bl[0];

    return bl;
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
 * Copies unjoined trees as they do not get changed.
 * 
 * @param trees 
 * @param minDistPos 
 * @param size 
 * @return array of size size-1 where the last entrie is not set
 */
Node** copyUnjoinedTrees(Node** trees, int* minDistPos, int size) {
    Node** copy = new Node*[size - 1];
    int off{0};
    for (int i{0}; i < size; ++i) {
        // skip and add offset if this tree is one of the two that will be merged
        if (i == minDistPos[0] || i == minDistPos[1]) {
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
Node* neighborJoining(double* distMatrix, Node** trees, int size) {
    //printTrees(trees, size);
    //printMatrix(distMatrix, size);
    if (size > 2) {
        // calculations
        double* r = meanDistance(distMatrix, size);
        //std::cout << "r = [" << r[0] << ", " << r[1] << ", " << r[2] << ", " << r[3] << "]" << std::endl;
        //printMatrix(iM, size);
        int* minDistPos = findMin(interMatrix(distMatrix, size, r), size);
        //std::cout << "Joining " << trees[minDistPos[0]]->toString() << "(idx=" << minDistPos[0] << ",bl=" << blMin[0] << ") and "
        //    << trees[minDistPos[1]]->toString() << "(idx=" << minDistPos[1] << ",bl=" << blMin[1]<< ")" << std::endl;
        double* newDistMatrix = calculateNewMatrix(distMatrix, size, minDistPos);

        // tree building
        double* branchLengths = calculateBranchLengths(distMatrix[minDistPos[0] * size + minDistPos[1]],
                                                       r[minDistPos[0]], r[minDistPos[1]]);
        Node** newTrees = copyUnjoinedTrees(trees, minDistPos, size);
        // add joined trees (neighbors)
        newTrees[size - 2] = join(trees[minDistPos[0]], branchLengths[0],
                                  trees[minDistPos[1]], branchLengths[1]);

        return neighborJoining(newDistMatrix, newTrees, size - 1);
    } else if (size == 2) {
        // base case
        //std::cout << "Joining last two subtrees with distance " << distMatrix[1] << std::endl;
        // connecting branch defined by two branches of same length to parent -> "rooted" tree...
        return join(trees[0], distMatrix[1], trees[1], distMatrix[1]);
    } else {
        // edge case, only one or none clusters
        return *trees;
    }
}
