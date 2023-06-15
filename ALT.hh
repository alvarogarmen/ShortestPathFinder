// Created by alvaro on 29/05/2023.
//

#ifndef UNTITLED_ALT_HH
#define UNTITLED_ALT_HH

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>

double estimate(double source, double target, const std::vector<std::vector<double>>& potentials) {
    double potential;

    for (int i = 0; i < potentials[0].size(); i++) {
        double potentialPlus = potentials[target][i] - potentials[source][i];
        double potentialMinus = potentials[source][i] - potentials[target][i];
        potential = (potentialPlus > potentialMinus) ? potentialPlus : potentialMinus;
    }

    return potential;
}


double ALT(Graph myGraph, double sourceNode, double targetNode, const std::vector<std::vector<double>>& potentials) {
    APQ apq = APQ();
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0);
    dist[sourceNode - 1] = 0;

    while (!apq.isEmpty()) {
        double currentNode = apq.getMin().first;
        visited.insert(currentNode);
        apq.popMin();

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

                double h = estimate(edge + 1, targetNode, potentials);
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
