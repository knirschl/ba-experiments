//
// Created by balin on 20.05.23.
//

#ifndef BA_PARSE_PHYLIP_H
#define BA_PARSE_PHYLIP_H

#include <vector>
#include <fstream>
#include <string>
#include <regex>
#include <iostream>
#include <string>
#include "../common_types.h"

// non-catching outer group, then one catching group for X.XXXXXX and one for X.XXXXe-XX
static const std::string DOUBLE_REGEX= "(?:(0.[0-9]+)|([1-9].[0-9]+e-[0-9]{2}))";
static const std::string WHITESPACE_DOUBLE_REGEX = "\\s*" + DOUBLE_REGEX;
static const std::string IDENT_REGEX = "([0-9]+)";

template<typename T>
matrix_t<T> parse(std::ifstream& stream) {
    matrix_t<T> matrix {};
    vector_t<std::string> ids {};
    std::string line {};
    int seq_count {};
    // first line := size
    if (stream.is_open() && stream.good()) {
        std::getline(stream, line);
        seq_count = std::stoi(line);
        matrix.resize(seq_count);
        for (auto& c : matrix) { // ref needed to update actually column and not a copy
            c.resize(seq_count);
        }
        ids.resize(seq_count);
    }
    // build regex depending on seq_count
    std::string dist_regex_str {};
    for (int i{}; i < seq_count; i++) {
        dist_regex_str += WHITESPACE_DOUBLE_REGEX;
    }
    std::regex line_regex{IDENT_REGEX + dist_regex_str};
    std::smatch results;

    size_t outer_counter{};
    while (stream) {
        size_t inner_counter{};
        std::getline(stream, line);
        // match regex
        if (std::regex_match(line, results, line_regex)) {
            ids[outer_counter] = results[1]; // first group := id
            for (int i{2}; i < results.size(); i++) {
                if (results[i].matched) {
                    // check to convert only matched groups (either X.XXX or X.XXe-XX)
                    matrix[outer_counter][inner_counter++] = std::stold(results[i]);
                }
            }
        }
        outer_counter++;
    }

    return matrix;
}

template<typename T>
matrix_t<T> parse_from_file(std::string const& file) {
    std::ifstream streamed_file {file};
    return parse<T>(streamed_file);
}

#endif //BA_PARSE_PHYLIP_H
