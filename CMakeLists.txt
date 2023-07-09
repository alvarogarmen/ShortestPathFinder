#Put everything, including the graph files under the same directory
cmake_minimum_required(VERSION 3.25)
project(praktikum)

set(CMAKE_CXX_STANDARD 17)

add_executable(praktikum main.cpp argtable3.c DijkstraBidirectional.hh Experiments.hh)
