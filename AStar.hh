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
    APQ apq = APQ();
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

#endif
