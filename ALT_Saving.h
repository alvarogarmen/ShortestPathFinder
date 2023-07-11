//
// Created by alvar on 05/06/2023.
//

#ifndef PRAKTIKUM_ALT_SAVING_H
#define PRAKTIKUM_ALT_SAVING_H

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include "ALT.hh"
#include <cmath>
#include <set>
#include <fstream>

double ALTSaving(Graph myGraph, double sourceNode, double targetNode, const std::vector<std::vector<double>>& potentials,
                 std::string exploredNodes, std::vector<double>& landmarks) {
    std::vector<double> usefulLandmarks = findUsefulLandmarks(myGraph, sourceNode, targetNode, landmarks);
    APQSaving apq = APQSaving();
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityDist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    std::ofstream exploredFile("experiments/"+exploredNodes);  // File stream for explored nodes
    std::ofstream usefulFile("experiments/"+exploredNodes+"_UsefulLandmarks");
    for (auto i : usefulLandmarks){
        usefulFile << myGraph.getNode(landmarks[i]).coordinateX << " " << myGraph.getNode(landmarks[i]).coordinateY <<std::endl;
    }
    usefulFile.close();

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();
        visited.insert(currentNode);

        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;

                if (edge == targetNode -1){
                    break;
                }
                priorityDist[edge] = priorityDist[currentNode] - estimate(currentNode, targetNode - 1, potentials, usefulLandmarks) + weight + estimate(edge, targetNode - 1, potentials, usefulLandmarks);

                double f = priorityDist[edge];

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, f, currentNode);
                } else {
                    apq.insertNode(edge, f, currentNode);
                }
            }
        }
    }
    for (auto currentNode : visited){
        exploredFile << myGraph.getNode(currentNode).coordinateX << " " << myGraph.getNode(currentNode).coordinateY << std::endl;  // Write explored node to the file
    }
    return dist[targetNode-1];
}


double ALTSearchSpace(Graph myGraph, double sourceNode, double targetNode, const std::vector<std::vector<double>>& potentials,
                      std::vector<double>& landmarks) {
    std::vector<double> usefulLandmarks = findUsefulLandmarks(myGraph, sourceNode, targetNode, landmarks);
    APQSaving apq = APQSaving();
    std::vector<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;
    visited.push_back(sourceNode);

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();


        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            visited.push_back(edge);
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;

                if (edge == targetNode - 1) {
                    return visited.size();
                }

                double h = estimate(edge, targetNode-1, potentials, usefulLandmarks);
                double f = dist[edge] + h;

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
        return visited.size();
    }
    return visited.size();
}

#endif //PRAKTIKUM_ALT_SAVING_H
