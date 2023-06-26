//
// Created by alvar on 05/06/2023.
//
/*
#ifndef PRAKTIKUM_DUMPSTER_H
#define PRAKTIKUM_DUMPSTER_H

#include <cmath>
#include <unordered_set>
#include "Graph.hh"
#include "APQ.hh"

/*
double AStarBidirectionalSaving(Graph myGraph, double sourceNode, double targetNode, const std::string& exploredFileName, const std::string& pathFileName) {
    APQ forwardAPQ = APQ();
    APQ backwardAPQ = APQ();
    std::set<double> forwardVisited;
    std::set<double> backwardVisited;
    std::vector<double> forwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> backwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> forwardPath(myGraph.nodes.size(), -1);
    std::vector<double> backwardPath(myGraph.nodes.size(), -1);

    forwardAPQ.insertNode(sourceNode - 1, 0);
    backwardAPQ.insertNode(targetNode - 1, 0);
    forwardDist[sourceNode - 1] = 0;
    backwardDist[targetNode - 1] = 0;

    double bestPath = INT_MAX;

    std::ofstream exploredNodeFile(exploredFileName);
    std::ofstream pathFile(pathFileName);

    while (!forwardAPQ.isEmpty() && !backwardAPQ.isEmpty()) {
        double forwardCurrentNode = forwardAPQ.popMin();
        double backwardCurrentNode = backwardAPQ.popMin();

        forwardVisited.insert(forwardCurrentNode);
        backwardVisited.insert(backwardCurrentNode);



        double forwardStartEdge = (forwardCurrentNode > 0) ? myGraph.edgeStarts[forwardCurrentNode - 1] + 1 : 0;
        double forwardEndEdge = myGraph.edgeStarts[forwardCurrentNode];

        double backwardStartEdge = (backwardCurrentNode > 0) ? myGraph.edgeStarts[backwardCurrentNode - 1] + 1 : 0;
        double backwardEndEdge = myGraph.edgeStarts[backwardCurrentNode];

        // Forward search
        for (double forwardEdgeIndex = forwardStartEdge; forwardEdgeIndex <= forwardEndEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardCurrentNode], myGraph.nodes[forwardEdge]);

            if (forwardDist[forwardCurrentNode] + forwardWeight < forwardDist[forwardEdge]) {
                forwardDist[forwardEdge] = forwardDist[forwardCurrentNode] + forwardWeight;
                forwardPath[forwardEdge] = forwardCurrentNode;

                if (backwardVisited.count(forwardEdge) > 0) {
                    double pathCost = forwardDist[forwardEdge] + backwardDist[forwardEdge];
                    if (pathCost < bestPath) {
                        bestPath = pathCost;
                    }
                }

                double forwardH = distance(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode - 1]);
                double forwardF = forwardDist[forwardEdge] + forwardH;

                if (forwardAPQ.contains(forwardEdge)) {
                    forwardAPQ.decreaseKey(forwardEdge, forwardF);
                } else {
                    forwardAPQ.insertNode(forwardEdge, forwardF);
                }
            }
        }

        // Backward search
        for (double backwardEdgeIndex = backwardStartEdge; backwardEdgeIndex <= backwardEndEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardCurrentNode], myGraph.nodes[backwardEdge]);

            if (backwardDist[backwardCurrentNode] + backwardWeight < backwardDist[backwardEdge]) {
                backwardDist[backwardEdge] = backwardDist[backwardCurrentNode] + backwardWeight;
                backwardPath[backwardEdge] = backwardCurrentNode;

                if (forwardVisited.count(backwardEdge) > 0) {
                    double pathCost = backwardDist[backwardEdge] + forwardDist[backwardEdge];
                    if (pathCost < bestPath) {
                        bestPath = pathCost;
                    }
                }

                double backwardH = distance(myGraph.nodes[backwardEdge], myGraph.nodes[sourceNode - 1]);
                double backwardF = backwardDist[backwardEdge] + backwardH;

                if (backwardAPQ.contains(backwardEdge)) {
                    backwardAPQ.decreaseKey(backwardEdge, backwardF);
                } else {
                    backwardAPQ.insertNode(backwardEdge, backwardF);
                }
            }
        }
        //Check if shortest path already found
        if (bestPath != INT_MAX) {
            break;
        }
    }

    std::cout << "All available edges relaxed" << std::endl;
    if (bestPath == INT_MAX) {
        return -1;
    }

    double currentNode = targetNode - 1;
    while (currentNode != sourceNode - 1) {
        pathFile << myGraph.nodes[currentNode].coordinateX << " " << myGraph.nodes[currentNode].coordinateY << "\n";
        currentNode = forwardPath[currentNode];
    }
    pathFile << myGraph.nodes[sourceNode - 1].coordinateX << " " << myGraph.nodes[sourceNode - 1].coordinateY << "\n";

    exploredNodeFile.close();
    pathFile.close();

    return bestPath;
}*/