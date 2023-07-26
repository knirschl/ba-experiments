//
// Created by knirschl on 25.05.23.
//

#ifndef BA_WRITE_FILE_H
#define BA_WRITE_FILE_H

#include <fstream>
#include <iostream>
#include "../nj/tree.h"
#include "../misc/common_types.h"

bool write_newick(Tree const &tree, std::ofstream &writer) {
    if (writer.is_open() && writer.good()) {
        writer << tree.to_newick() << '\n';
    }
    return !writer.fail();
}

bool write_newick(Tree const &tree, std::string const &out_file) {
    std::ofstream streamed_file{out_file};
    bool is_written = write_newick(tree, streamed_file);
    streamed_file.close();
    return is_written;
}

bool write_phylip(const dist_matrix_t &matrix, const vector_t<std::string> &names, std::ofstream &writer) {
    if (names.size() != matrix.size()) {
        std::cout << "Matrix and names don't match in size.\n";
    }
    const size_t ms{matrix.size()};
    if (writer.is_open() && writer.good()) {
        writer << ms << '\n';
        for (int row{}; row < ms; row++) {
            if (names[row].length() > 10) {
                // not really possible because of parse_file.cpp
                return false;
            }
            // one-time padding
            writer << std::setw(10) << names[row];
            for (auto &col: matrix[row]) {
                writer << col << ' ';
            }
            writer << '\n';
        }
    }
    return !writer.fail();
}

bool write_phylip(const dist_matrix_t &matrix, const vector_t<std::string> &names, const std::string &out_file) {
    std::ofstream streamed_file{out_file};
    // number precision after comma << keep word on left, fill chars to right
    streamed_file << std::setprecision(17) << std::left;
    bool is_written = write_phylip(matrix, names, streamed_file);
    streamed_file.close();
    return is_written;
}

#endif //BA_WRITE_FILE_H
