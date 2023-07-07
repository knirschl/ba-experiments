//
// Created by knirschl on 20.05.23.
//

#ifndef BA_COMMON_TYPES_H
#define BA_COMMON_TYPES_H

#include <vector>

template <typename T>
using vector_t = std::vector<T>;

template <typename T>
using matrix_t = std::vector<vector_t<T>>;

using dist_t = double;
using dist_vector_t = vector_t<dist_t>;
using dist_matrix_t = std::vector<dist_vector_t>;

#endif //BA_COMMON_TYPES_H
