//
// Created by knirschl on 20.05.23.
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
#include <sstream>

using namespace std::literals;

// anything but whitespace, capturing
static const std::string IDENT_REGEX_STR{R"((\S+))"};
// 0, -0.00000, 0.12312
static const std::string ZERO_POINT_DBL_REGEX_STR{R"(([+\-]?0(?:\.\d+)?))"};
// e00, e12, e+01, e-99, non-capturing
static const std::string POWER10_REGEX_STR{R"((?:e[+\-]?\d{2}))"};
// 1, 213123, -12, capturing
static const std::string NON_ZERO_INT_REGEX_STR{R"([+\-]?[1-9]\d*)"};
// 12321, 12321.3213, 3123.213e-23
static const std::string INT_FLP_SCI_REGEX_STR{"(" + NON_ZERO_INT_REGEX_STR + "(?:\\.\\d+" + POWER10_REGEX_STR + "?)?)"};
// non-capturing outer group, then one catching group for X.XXXXXX and one for X.XXXXe+-XX
static const std::string MATRIX_ELEMENT_REGEX_STR{"(?:" + ZERO_POINT_DBL_REGEX_STR + "|" + INT_FLP_SCI_REGEX_STR + ")"};
static const std::string WHITESPACE_MAT_ELEMENT_REGEX_STR{R"(\s*)" + MATRIX_ELEMENT_REGEX_STR};

template<typename T>
std::pair<matrix_t<T>, vector_t<std::string>> parse(std::ifstream& reader) {
    matrix_t<T> matrix {};
    vector_t<std::string> ids {};
    std::string line {};
    int seq_count {};
    // first line := size
    if (reader.is_open() && reader.good()) {
        std::getline(reader, line);
        seq_count = std::stoi(line);
        matrix.resize(seq_count);
        for (auto& c : matrix) { // ref needed to update actually column and not a copy
            c.resize(seq_count);
        }
        ids.resize(seq_count);
    }
    // build regex depending on seq_count
    std::stringstream line_regex_stream;
    line_regex_stream << IDENT_REGEX_STR << "(?:"s << WHITESPACE_MAT_ELEMENT_REGEX_STR << ")"s << "{"s << seq_count << "}"s;
    std::regex line_regex{line_regex_stream.str()};
    std::regex whitespace_mat_element_regex{WHITESPACE_MAT_ELEMENT_REGEX_STR};
    std::smatch results;

    size_t outer_counter{};
    while (reader) {
        size_t inner_counter{};
        std::getline(reader, line);
        // match regex
        if (std::regex_match(line, results, line_regex)) {
            ids[outer_counter] = results[1]; // first group := id
            line.erase(line.begin(), line.begin() + 10); // remove id
            while (std::regex_search(line, results, whitespace_mat_element_regex)) {
                for (int i{1}; i < results.size(); i++) {
                    if (results[i].matched) {
                        // check to convert only matched group (either X.XXX or X.XXe-XX)
                        matrix[outer_counter][inner_counter++] = std::stold(results[i]);
                        break;
                    }
                }
                line = results.suffix();
            }
        }
        outer_counter++;
    }

    return std::pair<matrix_t<T>, vector_t<std::string>>{std::move(matrix), std::move(ids)};
}

template<typename T>
std::pair<matrix_t<T>, vector_t<std::string>> parse_from_file(std::string const& file) {
    std::ifstream streamed_file{file};
    auto read = parse<T>(streamed_file);
    streamed_file.close();
    return read;
}

#endif //BA_PARSE_PHYLIP_H
