// Created by alvar on 29/05/2023.
//

#ifndef UNTITLED_ASTAR_HH
#define UNTITLED_ASTAR_HH

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include "AStarBidirectional.hh"
#include <cmath>
#include <set>

double AStar(Graph& myGraph, double& sourceNode, double& targetNode) {
    APQ apq = APQ(myGraph.nodeCount);
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0);
    dist[sourceNode - 1] = 0;

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();

        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;

                if (edge == targetNode - 1) {
                    return dist[edge];
                }

                double h = distance(myGraph.nodes[edge], myGraph.nodes[targetNode - 1]);
                double f = dist[edge] + h;

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, f);
                } else {
                    apq.insertNode(edge, f);
                }
            }
        }
    }

    std::cout << "All available edges relaxed" << std::endl;
    if (dist[targetNode - 1] == INT_MAX) {
        return -1;
    }
    return dist[targetNode - 1];
}


double AStarSaving(Graph myGraph, double sourceNode, double targetNode) {
    APQSaving apq = APQSaving();
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    std::ofstream exploredFile("experiments/AStar_explored_nodes");  // File stream for explored nodes
    std::ofstream pathFile("experiments/AStar_path");      // File stream for path

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();
        visited.insert(currentNode);
        exploredFile << myGraph.getNode(currentNode).coordinateX << " " << myGraph.getNode(currentNode).coordinateY << std::endl;  // Write explored node to the file


        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;

                if (edge == targetNode - 1) {
                    if (edge == targetNode - 1) {
                        std::vector<double> path;
                        path.push_back(edge);
                        double prevNode = currentNode;
                        while (prevNode != -1 && prevNode != sourceNode - 1) {
                            path.push_back(prevNode);
                            prevNode = apq.getPrev(prevNode);
                        }
                        path.push_back(prevNode);

                        // Write the path to the file in reverse order
                        for (int i = 0; i < path.size() - 1; i++) {
                            pathFile << myGraph.getNode(path[i]).coordinateX << " "
                                     << myGraph.getNode(path[i]).coordinateY
                                     << " " << path[i] << " " << std::endl;
                        }

                        exploredFile.close();  // Close the files before returning
                        pathFile.close();
                        return dist[edge];
                    }
                }

                double h = distance(myGraph.nodes[edge], myGraph.nodes[targetNode - 1]);
                double f = dist[edge] + h;

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, f, currentNode);
                }
                else {
                    apq.insertNode(edge, f, currentNode);
                }
            }
        }
    }

    std::cout << "All available edges relaxed" << std::endl;
    if (dist[targetNode - 1] == INT_MAX) {
        return -1;
    }
    return dist[targetNode - 1];
}


double AStarSearchSpace(Graph myGraph, double sourceNode, double targetNode) {
    APQSaving apq = APQSaving();
    std::vector<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;
    visited.push_back(sourceNode);

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();
        visited.push_back(currentNode);


        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;

                if (edge == targetNode - 1) {
                    if (edge == targetNode - 1) {
                        return visited.size();
                    }
                }

                double h = distance(myGraph.nodes[edge], myGraph.nodes[targetNode - 1]);
                double f = dist[edge] + h;

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, f, currentNode);
                }
                else {
                    apq.insertNode(edge, f, currentNode);
                }
            }
        }
    }

    std::cout << "All available edges relaxed" << std::endl;
    if (dist[targetNode - 1] == INT_MAX) {
        return -1;
    }
    return dist[targetNode - 1];
}

#endif
