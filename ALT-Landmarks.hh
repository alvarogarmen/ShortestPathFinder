//
// Created by alvar on 29/05/2023.
//

#ifndef PRAKTIKUM_ALT_LANDMARKS_HH
#define PRAKTIKUM_ALT_LANDMARKS_HH

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <unordered_set>

std::unordered_map<int, double> computeLandmarkDistancesRandom(Graph myGraph, int numLandmarks) {
    //Initialize the landmarkDistances
    std::unordered_map<int, double> landmarkDistances;
    // Select random landmark nodes
    std::vector<Node> nodes = myGraph.getNodes();
    std::random_shuffle(nodes.begin(), nodes.end());
    std::vector<Node> landmarks(nodes.begin(), nodes.begin() + numLandmarks);

    // Compute distances from landmarks to all other nodes using Dijkstra's algorithm
    for (Node landmark : landmarks) {
        std::vector<double> dist(myGraph.nodes.size(), INT_MAX);
        std::set<double> visited;
        APQ apq;

        apq.insertNode(landmark.nodeId - 1, 0);
        dist[landmark.nodeId - 1] = 0;

        while (!apq.isEmpty()) {
            double currentNode = apq.getMin().first;
            visited.insert(currentNode);
            apq.popMin();

            double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
            double endEdge = myGraph.edgeStarts[currentNode];
f
            for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
                double edge = myGraph.edges[edgeIndex] - 1;
                double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

                if (dist[currentNode] + weight < dist[edge]) {
                    dist[edge] = dist[currentNode] + weight;

                    if (visited.find(edge) == visited.end()) {
                        apq.insertNode(edge, dist[edge]);
                    }
                }
            }
        }

        // Store the landmark distances
        for (int i = 0; i < dist.size(); i++) {
            landmarkDistances[i + 1] += dist[i];
        }
    }

    // Normalize the landmark distances
    for (auto& entry : landmarkDistances) {
        entry.second /= numLandmarks;
    }
    return landmarkDistances;
}

std::vector<int> computeFurthestLandmarks(Graph myGraph, int numLandmarks) {
    std::vector<int> landmarks;
    std::vector<Node> nodes = myGraph.getNodes();
    std::unordered_set<int> selectedNodes;  // To keep track of selected landmarks

    // Randomly select the first landmark node
    int firstLandmark = rand() % nodes.size();
    landmarks.push_back(firstLandmark);
    selectedNodes.insert(firstLandmark);

    // Select the remaining landmarks
    for (int i = 1; i < numLandmarks; i++) {
        double maxDistance = 0;
        int furthestNode = -1;

        // Find the node that is furthest from the already selected landmarks
        for (int j = 0; j < nodes.size(); j++) {
            if (selectedNodes.find(j) == selectedNodes.end()) {
                double minDistance = INT_MAX;

                // Calculate the minimum distance to already selected landmarks
                for (int landmark : landmarks) {
                    double landmarkDistance = distance(myGraph.nodes[landmark], myGraph.nodes[j]);
                    if (landmarkDistance < minDistance) {
                        minDistance = landmarkDistance;
                    }
                }

                // Update the furthest node if it has the maximum minimum distance
                if (minDistance > maxDistance) {
                    maxDistance = minDistance;
                    furthestNode = j;
                }
            }
        }

        // Add the furthest node as a landmark
        landmarks.push_back(furthestNode);
        selectedNodes.insert(furthestNode);
    }
    return landmarks;
}
#endif //PRAKTIKUM_ALT_LANDMARKS_HH
