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

#endif
