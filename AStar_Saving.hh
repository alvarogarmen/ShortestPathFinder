// Created by alvar on 04/06/2023.
//

#ifndef UNTITLED_ASTAR_SAVING_HH
#define UNTITLED_ASTAR_SAVING_HH

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>
#include <fstream>

double AStarSaving(Graph myGraph, double sourceNode, double targetNode) {
    APQSaving apq = APQSaving();
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    std::ofstream exploredFile;  // File stream for explored nodes
    std::ofstream pathFile;      // File stream for path
    exploredFile.open("AStar_explored_nodes.txt");  // Open explored nodes file
    pathFile.open("AStar_path.txt");                // Open path file

    while (!apq.isEmpty()) {
        double currentNode = apq.getMin().first;
        visited.insert(currentNode);
        exploredFile << myGraph.getNode(currentNode).coordinateX << " " << myGraph.getNode(currentNode).coordinateY << std::endl;  // Write explored node to the file
        apq.popMin();

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

#endif
