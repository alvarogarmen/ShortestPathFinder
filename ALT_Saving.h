//
// Created by alvar on 05/06/2023.
//

#ifndef PRAKTIKUM_ALT_SAVING_H
#define PRAKTIKUM_ALT_SAVING_H

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include "Dijkstra_Saving.hh"
#include "ALT.hh"
#include <cmath>
#include <set>
#include <fstream>

double ALTSaving(Graph myGraph, double sourceNode, double targetNode, const std::vector<std::vector<double>>& potentials, std::string exploredNodes, std::string path) {
    APQSaving apq = APQSaving();
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    std::ofstream exploredFile;  // File stream for explored nodes
    std::ofstream pathFile;      // File stream for path
    exploredFile.open(exploredNodes);  // Open explored nodes file
    pathFile.open(path);                // Open path file

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
                    std::vector<double> path;
                    path.push_back(edge);
                    double prevNode = currentNode;
                    while (prevNode != -1 && prevNode != sourceNode - 1) {
                        path.push_back(prevNode);
                        prevNode = apq.getPrev(prevNode);
                    }
                    path.push_back(prevNode);

                    // Write the path to the file in reverse order
                    for (int i = path.size() - 1; i >= 0; i--) {
                        pathFile << myGraph.getNode(path[i]).coordinateX << " "
                                 << myGraph.getNode(path[i]).coordinateY
                                 << " " << path[i] << " " << std::endl;
                    }

                    exploredFile.close();  // Close the files before returning
                    pathFile.close();
                    return dist[edge];
                }

                double f = dist[edge] + potentials[edge][targetNode - 1];

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, f, currentNode);
                } else {
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

#endif //PRAKTIKUM_ALT_SAVING_H
