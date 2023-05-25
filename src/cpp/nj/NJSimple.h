//
// Created by knirschl on 22.05.23.
//

#ifndef BA_NJSIMPLE_H
#define BA_NJSIMPLE_H

#include "tree.h"
#include "../common_types.h"

template<typename T>
vector_t<T> meanDistance(matrix_t<T> const& distMatrix) {
    const size_t n = distMatrix.size();
    vector_t<T> r;
    for (size_t i{0}; i < n; ++i) {
        T sumK{0};
        for (size_t k{0}; k < n; ++k) {
            sumK += distMatrix[i][k];
        }
        r.push_back((1.0 / (n - 2)) * sumK); // 1.0 for double division
    }

    return r;
}

template<typename T>
matrix_t<T> interMatrix(matrix_t<T> const& distMatrix, vector_t<T> const& r) {
    matrix_t<T> iM;
    const size_t n = distMatrix.size();
    for (size_t i{0}; i < n; ++i) {
        vector_t<T> row;
        for (size_t j{0}; j < n; ++j) {
            row.push_back(distMatrix[i][j] - (r[i] + r[j]));
        }
        iM.push_back(row);
    }

    return iM;
}

template<typename T>
std::pair<size_t, size_t> findMin(matrix_t<T> const& matrix) {
    T min = std::numeric_limits<T>::infinity();
    const size_t n = matrix.size();
    std::pair<size_t, size_t> minDistPos;
    for (size_t i{0}; i < n; ++i) {
        for (size_t j{0}; j < n; ++j) {
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

template<typename T>
matrix_t<T> calculateNewMatrix(matrix_t<T> const& distMatrix, std::pair<size_t, size_t> const& minDistPos) {
    const size_t n = distMatrix.size();
    matrix_t<T> newMatrix{};
    newMatrix.resize(n - 1);
    // pre-initialize to prevent creating one too many
    for (size_t i{ 0 }; i < n - 1; ++i) {
        newMatrix[i].resize(n - 1);
    }
    // copy old values of rows/columns that won't be changed
    size_t xOff{ 0 };
    for (size_t i{ 0 }; i < n; ++i) {
        // skip row if one of the merged ones
        if (i == minDistPos.first || i == minDistPos.second) {
            xOff++;
        } else {
            size_t yOff{ 0 };
            for (size_t j{ 0 }; j < n; ++j) {
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
    size_t off{0};
    for (size_t k{0}; k < n - 2; ++k) {
        // calculate offset (how much needs to be skipped) via black magic
        while (k + off == minDistPos.first || k + off == minDistPos.second) {
            off++;
        }
        /*
         * d_u,k = (d_i,k + d_j,k - d_i,j) / 2
         * u = (i,j) the merged clusters, k is the current cluster
         */
        newMatrix[k][n - 2] = (distMatrix[minDistPos.first][k + off]
                               + distMatrix[minDistPos.second][k + off]
                               - distMatrix[minDistPos.first][minDistPos.second])
                              / 2.0;
        newMatrix[n - 2][k] = (distMatrix[minDistPos.first][k + off]
                               + distMatrix[minDistPos.second][k + off]
                               - distMatrix[minDistPos.first][minDistPos.second])
                              / 2.0;
    }
    newMatrix[n - 2][n - 2] = 0; // set new diagonal entry to 0

    return newMatrix;
}

// const& to be sdure if bigger datatypes than simple numbers are used
template<typename T>
std::pair<T, T> calculateBranchLengths(T const& dist, T const& r1, T const& r2) {
    return std::pair<T, T> {
            (dist + r1 - r2) / 2.0,
            dist - ((dist + r1 - r2) / 2.0)
    };
}

vector_t<std::shared_ptr<NTree>> copyUnjoinedTrees(vector_t<std::shared_ptr<NTree>> const& trees,
                                                   std::pair<size_t, size_t> const& minDistPos) {
    vector_t<std::shared_ptr<NTree>> copy;
    const size_t n = trees.size();
    for (size_t i{0}; i < n; ++i) {
        // skip and add offset if this tree is one of the two that will be merged
        if (i == minDistPos.first || i == minDistPos.second) {
            continue;
        }
        copy.push_back(trees[i]);
    }

    return copy;
}

template<typename T>
std::shared_ptr<NTree> neighborJoining(matrix_t<T> const& distMatrix, vector_t<std::shared_ptr<NTree>> const& trees) {
    const size_t n = distMatrix.size();
    /*static_assert([n, distMatrix]() {
        bool quadratic = true;
        for (auto row : distMatrix) {
            quadratic &= row.size() == n;
        }
        return quadratic;
    }, "distMatrix needs to be quadratic");
    static_assert(n == trees.size(), "Same amount of leafs as rows in matrix needed");*/
    if (n > 2) {
        // calculations
        vector_t<T> r = meanDistance(distMatrix);
        std::pair<size_t, size_t> minDistPos = findMin(interMatrix(distMatrix, r));
        matrix_t<T> newDistMatrix = calculateNewMatrix(distMatrix, minDistPos);

        //tree building
        std::pair<T, T> branchLengths = calculateBranchLengths(
                distMatrix[minDistPos.first][minDistPos.second],
                r[minDistPos.first],
                r[minDistPos.second]);
        vector_t<std::shared_ptr<NTree>> newTrees = copyUnjoinedTrees(trees, minDistPos);
        // add joined trees (neighbors)
        newTrees.push_back(join(trees[minDistPos.first], branchLengths.first,
                                             trees[minDistPos.second], branchLengths.second));

        return neighborJoining(newDistMatrix, newTrees);
    } else if (n == 2) {
        // base case
        // connecting branch defined by two branches of same length to parent -> "rooted" tree...
        return join(trees[0], distMatrix[0][1], trees[1], distMatrix[0][1]);
    } else {
        // edge case, only one or none clusters
        return trees[0];
    }
}

#endif //BA_NJSIMPLE_H
