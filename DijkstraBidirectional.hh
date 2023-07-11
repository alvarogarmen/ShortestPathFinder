//
// Created by alvar on 09/07/2023.
//

#ifndef PRAKTIKUM_DIJKSTRABIDIRECTIONAL_HH
#define PRAKTIKUM_DIJKSTRABIDIRECTIONAL_HH

#include <set>
#include <fstream>
#include "Graph.hh"
#include "APQ.hh"
#include "Dijkstra.hh"

double BidirectionalDijkstra(Graph& myGraph, double& sourceNode, double& targetNode) {
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<bool> visitedForward(myGraph.nodes.size());
    std::vector<bool> visitedBackward(myGraph.nodes.size());
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    double bestPath = INT_MAX;
    double meetingNode = -1;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        visitedForward[forwardNode]=true;
        visitedBackward[backwardNode]=true;



        if (visitedForward[backwardNode] && visitedBackward[forwardNode]) {
            return bestPath;
        }

        if (distForward[forwardNode] + distBackward[forwardNode] < bestPath) {
            bestPath = distForward[forwardNode] + distBackward[forwardNode];
            meetingNode = forwardNode;
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

                double forwardF = distForward[forwardEdge];

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

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;


                double backwardF = distBackward[backwardEdge];

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

double BidirectionalDijkstraSaving(Graph& myGraph, double& sourceNode, double& targetNode, std::string exploredFileName) {
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::set<double> visitedForward;
    std::set<double> visitedBackward;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> forwardPath(myGraph.nodes.size(), -1);
    std::vector<double> backwardPath(myGraph.nodes.size(), -1);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    double bestPath = INT_MAX;
    double meetingNode = -1;
    std::ofstream exploredNodeFile("experiments/"+exploredFileName);
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        visitedForward.insert(forwardNode);
        exploredNodeFile << myGraph.getNode(forwardNode).coordinateX << " " << myGraph.getNode(forwardNode).coordinateY << std::endl;  // Write explored node to the file

        visitedBackward.insert(backwardNode);
        exploredNodeFile << myGraph.getNode(backwardNode).coordinateX << " " << myGraph.getNode(backwardNode).coordinateY << std::endl;  // Write explored node to the file




        if (visitedForward.find(backwardNode)!=visitedForward.end() && visitedBackward.find(forwardNode)!=visitedBackward.end()) {
            exploredNodeFile.close();
            return bestPath;
        }

        if (distForward[forwardNode] + distBackward[forwardNode] < bestPath) {
            bestPath = distForward[forwardNode] + distBackward[forwardNode];
            meetingNode = forwardNode;
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
                forwardPath[forwardEdge] = forwardNode;

                double forwardH = 0;
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

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                double backwardH = 0;
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

    exploredNodeFile.close();

    return bestPath;
}

double BidirectionalDijkstraSearchSpace(Graph& myGraph, double& sourceNode, double& targetNode) {
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::set<double> visitedForward;
    std::set<double> visitedBackward;
    std::vector<double> visited;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> forwardPath(myGraph.nodes.size(), -1);
    std::vector<double> backwardPath(myGraph.nodes.size(), -1);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;
    visited.push_back(sourceNode);
    visited.push_back(targetNode);
    double bestPath = INT_MAX;
    double meetingNode = -1;
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        visitedForward.insert(forwardNode);

        visitedBackward.insert(backwardNode);




        if (visitedForward.find(backwardNode)!=visitedForward.end() && visitedBackward.find(forwardNode)!=visitedBackward.end()) {
            return visited.size();
        }

        if (distForward[forwardNode] + distBackward[forwardNode] < bestPath) {
            bestPath = distForward[forwardNode] + distBackward[forwardNode];
            meetingNode = forwardNode;
        }

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            visited.push_back(forwardEdge);
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                forwardPath[forwardEdge] = forwardNode;

                double forwardH = 0;
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
            visited.push_back(backwardEdge);
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                double backwardH = 0;
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
        return visited.size();
    }


    return visited.size();
}

#endif //PRAKTIKUM_DIJKSTRABIDIRECTIONAL_HH
