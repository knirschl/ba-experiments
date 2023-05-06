cmake_minimum_required(VERSION 3.24)
project(BA)

set(CMAKE_CXX_STANDARD 23)

add_executable(BA src/cpp/main.cpp src/cpp/tree.h src/cpp/matrix.h)
