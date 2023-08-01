//
// Created by knirschl on 20.05.23.
//

#ifndef BA_COMMON_TYPES_H
#define BA_COMMON_TYPES_H

#include <vector>

template<typename T>
using vector_t = std::vector<T>;

template<typename T>
using matrix_t = std::vector<vector_t<T>>;

using dist_t = double;
using dist_vector_t = vector_t<dist_t>;
using dist_matrix_t = std::vector<dist_vector_t>;

/**
  * Scales two matrices with independent scales and adds them together. Writes the result to a new
  * matrix and returns that one. The type parameter of the new matrix is the resulting type of
  * (T + S).
  *
 * @tparam T
 * @tparam S
 * @param left
 * @param right
 * @param scale_l
 * @param scale_r
 * @return the sum of the two scaled matrices
 */
template<typename T, typename S>
matrix_t<decltype(std::declval<T>() + std::declval<S>())>
mat_scaled_add(matrix_t<T> const &left, matrix_t<S> const &right,
               T const &scale_l = 1.0, S const &scale_r = 1.0) {
    matrix_t<decltype(std::declval<T>() + std::declval<S>())> result;
    result.resize(left.size());
    for (int i{}; i < left.size(); i++) {
        result[i].resize(left[i].size());
        for (int j{}; j < left[i].size(); j++) {
            result[i][j] = scale_l * left[i][j] + scale_r * right[i][j];
        }
    }
    return result;
}

std::string vec_to_string(const vector_t<double> &vec) {
    return "[" + accumulate(vec.begin() + 1, vec.end(),
                            std::to_string(vec[0]), [](const std::string &acc, double b) {
                return acc + ", " + std::to_string(b);
            }) + "]";
}

std::string vec_to_string(const vector_t<std::string> &vec) {
    return "[" + accumulate(vec.begin() + 1, vec.end(),
                            vec[0], [](const std::string &acc, std::string b) {
                return acc + ", " + b;
            }) + "]";
}

std::string mat_to_string(const matrix_t<double> &mat) {
    return "[" + accumulate(mat.begin() + 1, mat.end(),
                            vec_to_string(mat[0]), [](const std::string &acc, vector_t<double> b) {
                return acc + ",\n" + vec_to_string(b);
            }) + "]";
}

#endif //BA_COMMON_TYPES_H
