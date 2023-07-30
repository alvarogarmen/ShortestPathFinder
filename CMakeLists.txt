cmake_minimum_required(VERSION 3.16)
project(praktikum)

set(CMAKE_BUILD_TYPE Release)

# Set optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

add_executable(praktikum main.cpp argtable3.c)
