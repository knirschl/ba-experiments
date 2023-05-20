//
// Created by knirschl on 02.05.23.
//

#ifndef BA_MATRIX_H
#define BA_MATRIX_H

#include <vector>
#include <functional>
#include <string>
#include <numeric>
#include "../common_types.h"

template<typename T>
matrix_t<T> make_zero_matrix(size_t rows, size_t columns) {
    matrix_t<T> result{rows};
    for (auto c : result) {
        c.resize(columns);
    }
    return result;
}

/**
 * Adds two vectors. Writes the result to the given vector (param 3).
 * @tparam T 
 * @param left 
 * @param right 
 * @param result 
 */
template<typename T>
void vecadd(vector_t<T>& left, vector_t<T>& right, vector_t<T>& result) {
    std::transform(left.begin(), left.end(), right.begin(),
                   std::back_inserter(result),
                   [](T a, T b) { return a + b; });
}

/**
 * Adds two vectors. Writes the result to a new vector and returns that one.
 * @tparam T
 * @param left 
 * @param right 
 * @return 
 */
template<typename T>
vector_t<T> vecadd(vector_t<T>& left, vector_t<T>& right) {
    vector_t<T> result;
    std::transform(left.begin(), left.end(), right.begin(),
                   std::back_inserter(result),
                   [](T a, T b) { return a + b; });
    return result;
}

template<typename R, typename S>
void matadd(matrix_t<R>& left, matrix_t<S> const& right) {
    //std::transform(left.begin(), left.end(), right.begin(), std::back_inserter(left), std::bit_and<>());
}

/**
 * Scales a vector by a scalar. Writes the result to the given vector (param 3).
 * Type is the same as the given vector.
 * @tparam R
 * @tparam S
 * @param vec
 * @param scalar
 * @param result
 */
template<typename R, typename S>
void vecscale(vector_t<R>& vec, const S scalar, vector_t<R>& result) {
    std::transform(vec.begin(), vec.end(), std::back_inserter(result), [scalar](R e) { return e * scalar; });
}

/**
 * Scales a vector by a scalar. Writes the result to a new vector and returns that one.
 * Type is the same as the given vector.
 * @tparam R
 * @tparam S
 * @param vec
 * @param scalar
 * @return
 */
template<typename R, typename S>
vector_t<R> vecscale(vector_t<R>& vec, const S scalar) {
    vector_t<R> result;
    std::transform(vec.begin(), vec.end(), std::back_inserter(result), [scalar](R e) { return e * scalar; });
    return result;
}

/**
 * Scales a matrix by a scalar. Overwrites the matrix.
 * @tparam R
 * @tparam S
 * @param mat
 * @param scalar
 */
template<typename R, typename S>
void matscale(matrix_t<R>& mat, const S scalar) {
    //std::transform(matrix.begin(), matrix.end(), std::back_inserter(matrix), std::multiplies<>());
    std::for_each(mat.begin(), mat.end(),
                   [scalar](vector_t<R>& e) {
        e = vecscale(e, scalar);
    });
}

/**
 * Returns a string representation of a vector.
 * @tparam T
 * @param vec
 * @return
 */
template<typename T>
std::string vecstr(vector_t<T> vec) {
    std::string str = std::accumulate(vec.begin() + 1, vec.end(), "[" + std::to_string(vec[0]),
                                      [](const std::string& acc, T e) {
        return acc + ", " + std::to_string(e);
    });
    return str + "]";
}

/**
 * Returns a string representation of a matrix
 * @tparam T
 * @param mat
 * @return
 */
template<typename T>
std::string matstr(matrix_t<T> mat) {
    std::string str = std::accumulate(mat.begin() + 1, mat.end(), "[" + vecstr(mat[0]),
                                      [](const std::string& acc, vector_t<T> e) {
        return acc + ",\n " + vecstr<T>(e);
    });
    return str + "]";
}

#endif //BA_MATRIX_H
