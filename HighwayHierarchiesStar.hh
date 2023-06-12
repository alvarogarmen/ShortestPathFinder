//
// Created by alvar on 12/06/2023.
//

#ifndef PRAKTIKUM_HIGHWAYHIERARCHIESSTAR_HH
#define PRAKTIKUM_HIGHWAYHIERARCHIESSTAR_HH

#include <unordered_set>
#include <queue>
#include <algorithm>
#include "Graph.hh"
#include "Dijkstra.hh"

void preprocessHighwayHierarchiesStar(Graph& myGraph) {
    // Step 1: Compute node levels
    std::vector<int> nodeLevels(myGraph.nodes.size(), 0);
    std::vector<std::unordered_set<int>> lowerLevelNeighbors(myGraph.nodes.size());

    for (int node = 0; node < myGraph.nodes.size(); node++) {
        for (int neighbor : myGraph.getNeighbors(node)) {
            if (myGraph.nodes[neighbor].coordinateY < myGraph.nodes[node].coordinateY)
                lowerLevelNeighbors[node].insert(neighbor);
        }
    }

    std::queue<int> nodeQueue;
    for (int node = 0; node < myGraph.nodes.size(); node++) {
        if (lowerLevelNeighbors[node].empty()) {
            nodeQueue.push(node);
        }
    }

    while (!nodeQueue.empty()) {
        int node = nodeQueue.front();
        nodeQueue.pop();

        for (int neighbor : myGraph.getNeighbors(node)) {
            lowerLevelNeighbors[neighbor].erase(node);
            if (lowerLevelNeighbors[neighbor].empty()) {
                nodeLevels[neighbor] = nodeLevels[node] + 1;
                nodeQueue.push(neighbor);
            }
        }
    }

    // Step 2: Create hierarchy and compute shortcuts
    for (int level = 1; level < *std::max_element(nodeLevels.begin(), nodeLevels.end()); level++) {
        for (int node = 0; node < myGraph.nodes.size(); node++) {
            if (nodeLevels[node] == level) {
                for (int neighbor : myGraph.getNeighbors(node)) {
                    if (nodeLevels[neighbor] > level) {
                        double shortcutDistance = distance(myGraph.nodes[node], myGraph.nodes[neighbor]);
                        myGraph.addEdge(node, neighbor, shortcutDistance);
                    }
                }
            }
        }
    }
}

double highwayHierarchiesStar(Graph myGraph, double sourceNode, double targetNode) {
    preprocessHighwayHierarchiesStar(myGraph);

    std::vector<double> distances(myGraph.nodes.size(), std::numeric_limits<double>::infinity());
    std::vector<bool> visited(myGraph.nodes.size(), false);
    std::vector<double> forwardDistances(myGraph.nodes.size(), std::numeric_limits<double>::infinity());
    std::vector<double> backwardDistances(myGraph.nodes.size(), std::numeric_limits<double>::infinity());

    forwardDistances[sourceNode - 1] = 0.0;
    backwardDistances[targetNode - 1] = 0.0;

    while (true) {
        // Forward search
        double forwardNode = -1;
        double minForwardDistance = std::numeric_limits<double>::infinity();

        for (int node = 0; node < myGraph.nodes.size(); node++) {
            if (!visited[node] && forwardDistances[node] < minForwardDistance) {
                forwardNode = node;
                minForwardDistance = forwardDistances[node];
            }
        }

        if (forwardNode == -1)
            break;

        visited[forwardNode] = true;

        std::vector<double> forwardNeighbors = myGraph.getNeighbors(forwardNode);

        for (double neighbor : forwardNeighbors) {
            double weight = distance(myGraph.nodes[forwardNode], myGraph.nodes[neighbor]);

            if (forwardDistances[forwardNode] + weight < forwardDistances[neighbor]) {
                forwardDistances[neighbor] = forwardDistances[forwardNode] + weight;
            }
        }

        // Backward search
        double backwardNode = -1;
        double minBackwardDistance = std::numeric_limits<double>::infinity();

        for (int node = 0; node < myGraph.nodes.size(); node++) {
            if (!visited[node] && backwardDistances[node] < minBackwardDistance) {
                backwardNode = node;
                minBackwardDistance = backwardDistances[node];
            }
        }

        if (backwardNode == -1)
            break;

        visited[backwardNode] = true;

        std::vector<double> backwardNeighbors = myGraph.getNeighbors(backwardNode);

        for (double neighbor : backwardNeighbors) {
            double weight = distance(myGraph.nodes[backwardNode], myGraph.nodes[neighbor]);

            if (backwardDistances[backwardNode] + weight < backwardDistances[neighbor]) {
                backwardDistances[neighbor] = backwardDistances[backwardNode] + weight;
            }
        }
    }

    double shortestDistance = forwardDistances[sourceNode - 1] + backwardDistances[targetNode - 1];

    return shortestDistance;
}

std::vector<double> highwayHierarchiesStarPath(Graph myGraph, double sourceNode, double targetNode) {
    preprocessHighwayHierarchiesStar(myGraph);

    std::vector<double> distances(myGraph.nodes.size(), std::numeric_limits<double>::infinity());
    std::vector<bool> visited(myGraph.nodes.size(), false);
    std::vector<double> forwardDistances(myGraph.nodes.size(), std::numeric_limits<double>::infinity());
    std::vector<double> backwardDistances(myGraph.nodes.size(), std::numeric_limits<double>::infinity());

    forwardDistances[sourceNode - 1] = 0.0;
    backwardDistances[targetNode - 1] = 0.0;

    while (true) {
        // Forward search
        double forwardNode = -1;
        double minForwardDistance = std::numeric_limits<double>::infinity();

        for (int node = 0; node < myGraph.nodes.size(); node++) {
            if (!visited[node] && forwardDistances[node] < minForwardDistance) {
                forwardNode = node;
                minForwardDistance = forwardDistances[node];
            }
        }

        if (forwardNode == -1)
            break;

        visited[forwardNode] = true;

        std::vector<double> forwardNeighbors = myGraph.getNeighbors(forwardNode);

        for (double neighbor : forwardNeighbors) {
            double weight = distance(myGraph.nodes[forwardNode], myGraph.nodes[neighbor]);

            if (forwardDistances[forwardNode] + weight < forwardDistances[neighbor]) {
                forwardDistances[neighbor] = forwardDistances[forwardNode] + weight;
            }
        }

        // Backward search
        double backwardNode = -1;
        double minBackwardDistance = std::numeric_limits<double>::infinity();

        for (int node = 0; node < myGraph.nodes.size(); node++) {
            if (!visited[node] && backwardDistances[node] < minBackwardDistance) {
                backwardNode = node;
                minBackwardDistance = backwardDistances[node];
            }
        }

        if (backwardNode == -1)
            break;

        visited[backwardNode] = true;

        std::vector<double> backwardNeighbors = myGraph.getNeighbors(backwardNode);

        for (double neighbor : backwardNeighbors) {
            double weight = distance(myGraph.nodes[backwardNode], myGraph.nodes[neighbor]);

            if (backwardDistances[backwardNode] + weight < backwardDistances[neighbor]) {
                backwardDistances[neighbor] = backwardDistances[backwardNode] + weight;
            }
        }
    }

    double shortestDistance = forwardDistances[sourceNode - 1] + backwardDistances[targetNode - 1];
    std::vector<double> path;

    if (shortestDistance < std::numeric_limits<double>::infinity()) {
        double currentNode = sourceNode - 1;
        path.push_back(currentNode);

        while (currentNode != targetNode - 1) {
            std::vector<double> neighbors = myGraph.getNeighbors(currentNode);

            for (double neighbor : neighbors) {
                double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[neighbor]);

                if (backwardDistances[neighbor] + weight == forwardDistances[currentNode]) {
                    currentNode = neighbor;
                    path.push_back(currentNode);
                    break;
                }
            }
        }
    }

    return path;
}


#endif //PRAKTIKUM_HIGHWAYHIERARCHIESSTAR_HH
