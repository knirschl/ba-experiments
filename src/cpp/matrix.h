//
// Created by knirschl on 02.05.23.
//

#ifndef BA_MATRIX_H
#define BA_MATRIX_H

#include <vector>
#include <cstddef>
#include <array>

template<typename T>
std::vector<std::vector<T>> make_zero_matrix(size_t rows, size_t columns) {
    std::vector<std::vector<T>> result{rows};
    for (auto c : result) {
        c.resize(columns);
    }
    return result;
}

template<typename T, size_t N, size_t M>
std::vector<std::vector<T>> make_matrix_from_array(std::array<std::array<T, M>, N> const& arr) {
    std::vector<std::vector<T>> result{N};
    for (int row{}; row < N; ++row) {
        result[row].resize(M);
        for (int col{}; col < M; ++col) {
            result[row][col] = arr[row][col];
        }
    }
    return result;
}

// structs ?

#endif //BA_MATRIX_H
