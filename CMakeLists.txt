cmake_minimum_required(VERSION 3.24)
project(BA)

set(CMAKE_CXX_STANDARD 20)

add_executable(BA src/cpp/main.cpp src/cpp/nj/tree.h src/cpp/nj/matrix.h src/cpp/io/parse_phylip.h src/cpp/common_types.h src/cpp/nj/NJSimple.h src/cpp/io/write_tree.h src/cpp/io/parse_cli.h)
