// Created by alvar on 01/06/2023.
//

#ifndef UNTITLED_DIJKSTRA_SAVING_HH
#define UNTITLED_DIJKSTRA_SAVING_HH

#include "Graph.hh"
#include <climits>
#include "APQ_Saving.hh"
#include <cmath>
#include <set>
#include <fstream>  // Added file handling libraries

double distance(Node source, Node target) {
    return sqrt(std::pow(static_cast<double>(source.coordinateX - target.coordinateX), 2.0)
                + std::pow(static_cast<double>(source.coordinateY - target.coordinateY), 2.0));
}

double DijkstraSaving(Graph myGraph, double sourceNode, double targetNode) {
    APQ apq = APQ();
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0);
    dist[sourceNode - 1] = 0;

    std::ofstream exploredFile;  // File stream for explored nodes
    std::ofstream pathFile;      // File stream for path
    exploredFile.open("explored_nodes.txt");  // Open explored nodes file
    pathFile.open("path.txt");                // Open path file

    while (!apq.isEmpty()) {
        double currentNode = apq.getMin().first;
        visited.insert(currentNode);
        exploredFile << currentNode << " ";  // Write explored node to the file
        apq.popMin();

        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;

                if (edge == targetNode - 1) {
                    std::vector<double> path;
                    path.push_back(edge);
                    double prevNode = currentNode;
                    while (prevNode != sourceNode - 1) {
                        path.push_back(prevNode);
                        prevNode = apq.getPrev(prevNode);
                    }
                    path.push_back(prevNode);

                    for (int i = path.size() - 1; i >= 0; i--) {
                        pathFile << path[i] << " ";
                    }

                    exploredFile.close();  // Close the files before returning
                    pathFile.close();
                    return dist[edge];
                }

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge]);
                } else {
                    apq.insertNode(edge, dist[edge]);
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

#endif
