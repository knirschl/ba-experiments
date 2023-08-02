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
