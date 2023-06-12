//
// Created by alvar on 12/06/2023.
//

#ifndef PRAKTIKUM_ALT_BIDIRECTIONAL_HH
#define PRAKTIKUM_ALT_BIDIRECTIONAL_HH

#include <unordered_map>
#include <set>
#include "Graph.hh"
#include "APQ.hh"
#include "Dijkstra.hh"
#include "ALT.hh"

double ALTBidirectional(Graph myGraph, double sourceNode, double targetNode, std::unordered_map<int, double> landmarkDistances) {
    APQ apqForward = APQ();  // Priority queue for forward search
    APQ apqBackward = APQ(); // Priority queue for backward search
    std::set<double> visitedForward;  // Set of visited nodes for forward search
    std::set<double> visitedBackward; // Set of visited nodes for backward search
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);  // Distances from source for forward search
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX); // Distances from target for backward search

    // Initialize the search
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    double meetingNode = -1; // Node where forward and backward searches meet

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        // Perform one step of forward search
        double currentForwardNode = apqForward.getMin().first;
        visitedForward.insert(currentForwardNode);
        apqForward.popMin();

        double startForwardEdge = (currentForwardNode > 0) ? myGraph.edgeStarts[currentForwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[currentForwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[currentForwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[currentForwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[currentForwardNode] + forwardWeight;

                if (visitedBackward.find(forwardEdge) != visitedBackward.end()) {
                    // Node visited in both searches, we have a meeting node
                    meetingNode = forwardEdge;
                    break;
                }

                double forwardH = estimate(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode - 1],
                                           landmarkDistances[forwardEdge], landmarkDistances[targetNode - 1]);
                double forwardF = distForward[forwardEdge] + forwardH;

                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        if (meetingNode != -1)
            break;

        // Perform one step of backward search
        double currentBackwardNode = apqBackward.getMin().first;
        visitedBackward.insert(currentBackwardNode);
        apqBackward.popMin();

        double startBackwardEdge = (currentBackwardNode > 0) ? myGraph.edgeStarts[currentBackwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[currentBackwardNode];

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[currentBackwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[currentBackwardNode] + backwardWeight < distBackward[backwardEdge]) {
                distBackward[backwardEdge] = distBackward[currentBackwardNode] + backwardWeight;

                if (visitedForward.find(backwardEdge) != visitedForward.end()) {
                    // Node visited in both searches, we have a meeting node
                    meetingNode = backwardEdge;
                    break;
                }

                double backwardH = estimate(myGraph.nodes[backwardEdge], myGraph.nodes[sourceNode - 1],
                                            landmarkDistances[backwardEdge], landmarkDistances[sourceNode - 1]);
                double backwardF = distBackward[backwardEdge] + backwardH;

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }

        if (meetingNode != -1)
            break;
    }

    if (meetingNode == -1) {
        // No path found
        std::cout << "No path found" << std::endl;
        return -1;
    }

    // Calculate the shortest path distance
    double shortestPathDistance = INT_MAX;
    for (double node : visitedForward) {
        if (visitedBackward.find(node) != visitedBackward.end()) {
            double pathDistance = distForward[node] + distBackward[node];
            if (pathDistance < shortestPathDistance)
                shortestPathDistance = pathDistance;
        }
    }

    return shortestPathDistance;
}

#endif //PRAKTIKUM_ALT_BIDIRECTIONAL_HH
