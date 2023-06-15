#Put everything, including the graph files under the same directory
cmake_minimum_required(VERSION 3.25)
project(praktikum)

set(CMAKE_CXX_STANDARD 17)

add_executable(praktikum main.cpp argtable3.h argtable3.c Dijkstra.hh Graph.hh Node.hh readgraph.hh AStar.hh ALT.hh ALT-Landmarks.hh Dijkstra_Saving.hh APQ_Saving.hh AStar_Saving.hh Helper.hh Dumpster.h ALT_Saving.h ProcessInput.hh ALT_Bidirectional.hh HighwayHierarchiesStar.hh AStarBidirectional.hh)
