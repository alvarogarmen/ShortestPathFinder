#ifndef UNTITLED_DIJKSTRA_SAVING_HH
#define UNTITLED_DIJKSTRA_SAVING_HH

#include "Graph.hh"
#include <climits>
#include "APQ_Saving.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>
#include <fstream>  // Added file handling libraries

double DijkstraSaving(Graph myGraph, double sourceNode, double targetNode) {
    APQSaving apq;
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    std::ofstream exploredFile;  // File stream for explored nodes
    std::ofstream pathFile;      // File stream for path
    exploredFile.open("explored_nodes");  // Open explored nodes file
    pathFile.open("path");                // Open path file

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();
        visited.insert(currentNode);
        exploredFile << myGraph.getNode(currentNode).coordinateX << " " << myGraph.getNode(currentNode).coordinateY << std::endl;  // Write explored node to the file


        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            // Calculate the weight as the Euclidean distance between the nodes
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);
            // Relax the edge if a shorter path is found
            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;
                // Break if we have reached our destination
                if (edge == targetNode - 1) {
                    std::vector<double> path;
                    path.push_back(edge);
                    double prevNode = currentNode;
                    while (prevNode != -1 && prevNode != sourceNode-1) {
                        path.push_back(prevNode);
                        prevNode = apq.getPrev(prevNode);
                    }
                    path.push_back(prevNode);

                    // Write the path to the file in reverse order
                    for (int i = 0; i < path.size()-1; i++) {
                        pathFile << myGraph.getNode(path[i]).coordinateX << " " << myGraph.getNode(path[i]).coordinateY
                        << " "<<path[i]<< " " << std::endl;
                    }

                    exploredFile.close();  // Close the files before returning
                    pathFile.close();
                    return dist[edge];
                }

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge], currentNode);
                } else {
                    apq.insertNode(edge, dist[edge], currentNode);
                }
            }
        }
    }

    exploredFile.close();  // Close the files
    pathFile.close();

    std::cout << "All available edges relaxed" << std::endl;

    if (dist[targetNode - 1] == INT_MAX) {
        return -1;
    }
    return dist[targetNode - 1];
}



double BidirectionalDijkstraSaving(Graph& myGraph, double& sourceNode, double& targetNode, std::string exploredFileName) {
    APQ apqForward = APQ();
    APQ apqBackward = APQ();
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
    std::ofstream exploredNodeFile(exploredFileName);
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


#endif
