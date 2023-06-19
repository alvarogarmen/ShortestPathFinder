#ifndef PRAKTIKUM_ASTARBIDIRECTIONAL_HH
#define PRAKTIKUM_ASTARBIDIRECTIONAL_HH
#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>

double AStarBidirectional(Graph& myGraph, double& sourceNode, double& targetNode) {
    APQ apqForward = APQ();
    APQ apqBackward = APQ();
    std::set<double> visitedForward;
    std::set<double> visitedBackward;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    double bestPath = INT_MAX;
    double meetingNode = -1;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.getMin().first;
        double backwardNode = apqBackward.getMin().first;

        visitedForward.insert(forwardNode);
        visitedBackward.insert(backwardNode);

        apqForward.popMin();
        apqBackward.popMin();

        if (forwardNode == backwardNode) {
            meetingNode = forwardNode;
            break;
        }

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                double forwardH = distance(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode - 1]);
                double forwardF = distForward[forwardEdge] + forwardH;

                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight < distBackward[backwardEdge]) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                double backwardH = distance(myGraph.nodes[backwardEdge], myGraph.nodes[sourceNode - 1]);
                double backwardF = distBackward[backwardEdge] + backwardH;

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
    }

    if (meetingNode != -1) {
        double shortestPath = distForward[meetingNode] + distBackward[meetingNode];
        if (shortestPath < bestPath)
            bestPath = shortestPath;
    }

    if (bestPath == INT_MAX) {
        return -1;
    }

    return bestPath;
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