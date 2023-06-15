#ifndef PRAKTIKUM_ASTARBIDIRECTIONAL_HH
#define PRAKTIKUM_ASTARBIDIRECTIONAL_HH
#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>
double ManhattanDistance(Node source, Node target) {                             //Calculate Manhattan distance
    return std::abs(static_cast<double>(source.coordinateX - target.coordinateX))
                + std::abs(static_cast<double>(source.coordinateY - target.coordinateY));
}
double AStarBidirectional(Graph myGraph, double sourceNode, double targetNode) {
    APQ forwardAPQ = APQ();
    APQ backwardAPQ = APQ();
    std::set<double> forwardVisited;
    std::set<double> backwardVisited;
    std::vector<double> forwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> backwardDist(myGraph.nodes.size(), INT_MAX);

    forwardAPQ.insertNode(sourceNode - 1, 0);
    backwardAPQ.insertNode(targetNode - 1, 0);
    forwardDist[sourceNode - 1] = 0;
    backwardDist[targetNode - 1] = 0;

    while (!forwardAPQ.isEmpty() && !backwardAPQ.isEmpty()) {
        double forwardNode = forwardAPQ.getMin().first;
        forwardVisited.insert(forwardNode);
        forwardAPQ.popMin();

        double backwardNode = backwardAPQ.getMin().first;
        backwardVisited.insert(backwardNode);
        backwardAPQ.popMin();

        double forwardStartEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double forwardEndEdge = myGraph.edgeStarts[forwardNode];

        for (double edgeIndex = forwardStartEdge; edgeIndex <= forwardEndEdge; edgeIndex++) {
            double forwardEdge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (forwardDist[forwardNode] + weight < forwardDist[forwardEdge]) {
                forwardDist[forwardEdge] = forwardDist[forwardNode] + weight;

                if (forwardVisited.find(forwardEdge) == forwardVisited.end()) {
                    double h = distance(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode - 1]);
                    double f = forwardDist[forwardEdge] + h;

                    if (forwardAPQ.contains(forwardEdge)) {
                        forwardAPQ.decreaseKey(forwardEdge, f);
                    } else {
                        forwardAPQ.insertNode(forwardEdge, f);
                    }
                }
            }
            if (backwardVisited.find(forwardEdge)!= backwardVisited.end()){
                std::cout<<"forwardEdge Weight:"<<weight<<std::endl;
                double forwardDistToNode = forwardDist[forwardEdge];
                double backwardDistToNode = backwardDist[forwardEdge];
                return forwardDistToNode + backwardDistToNode;
            }
        }

        double backwardStartEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double backwardEndEdge = myGraph.edgeStarts[backwardNode];

        for (double edgeIndex = backwardStartEdge; edgeIndex <= backwardEndEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[backwardNode], myGraph.nodes[edge]);

            if (backwardDist[backwardNode] + weight < backwardDist[edge]) {
                backwardDist[edge] = backwardDist[backwardNode] + weight;

                if (backwardVisited.find(edge) == backwardVisited.end()) {
                    double h = ManhattanDistance(myGraph.nodes[edge], myGraph.nodes[sourceNode - 1]);
                    double f = backwardDist[edge] + h;

                    if (backwardAPQ.contains(edge)) {
                        backwardAPQ.decreaseKey(edge, f);
                    } else {
                        backwardAPQ.insertNode(edge, f);
                    }
                }
                if(forwardVisited.find(edge)!=forwardVisited.end()){
                    std::cout<<"backwardEdge"<<weight<<std::endl;
                    double forwardDistToNode = forwardDist[edge];
                    double backwardDistToNode = backwardDist[edge];
                    return forwardDistToNode + backwardDistToNode;
                }
            }
        }

        // Check for meeting point
        for (double node : forwardVisited) {
            if (backwardVisited.find(node) != backwardVisited.end()) {
                double forwardDistToNode = forwardDist[node];
                double backwardDistToNode = backwardDist[node];
                return forwardDistToNode + backwardDistToNode;
            }
        }
    }

    std::cout << "Paths did not meet" << std::endl;
    return -1;
}
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
        double forwardCurrentNode = forwardAPQ.getMin().first;
        double backwardCurrentNode = backwardAPQ.getMin().first;

        forwardVisited.insert(forwardCurrentNode);
        backwardVisited.insert(backwardCurrentNode);

        forwardAPQ.popMin();
        backwardAPQ.popMin();

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
}
#endif