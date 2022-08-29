cmake_minimum_required(VERSION 3.20)
project(MiQsDD)

set(CMAKE_CXX_STANDARD 17)

add_executable(MiQsDD main.cpp /include/ test/test_pkg.cpp)
