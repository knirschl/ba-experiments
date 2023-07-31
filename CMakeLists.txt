cmake_minimum_required(VERSION 3.24)
project(BA)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(BA src/cpp/main.cpp)
find_package(OpenMP)
target_link_libraries(BA PRIVATE OpenMP::OpenMP_CXX)