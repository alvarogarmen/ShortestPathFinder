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
#include <random>

std::vector<std::vector<double>> precomputePotentialsEuclidian(Graph myGraph, const std::vector<double>& landmarks){
    std::vector<std::vector<double>> potentials(myGraph.getNodes().size(), std::vector<double>(landmarks.size()));
    std::vector<std::vector<double>> potentialsInverted((landmarks.size()), std::vector<double>(myGraph.getNodes().size()));

    for (double i = 0; i<landmarks.size(); i++){
        for(int j = 0; j<potentialsInverted.size(); j++){
            potentials[i]= DijkstraToALL(myGraph, landmarks[i]);
            //potentials[i][j]=distance(myGraph.getNode(i), myGraph.getNode(landmarks[j]));

        }
    }

    return potentialsInverted;
}
std::vector<std::vector<double>> precomputeLandmarkToNodes(Graph myGraph, const std::vector<double>& landmarks){
    std::vector<std::vector<double>> landmarkToNodes(landmarks.size(), std::vector<double>(myGraph.getNodes().size()));

    int i = 0;
    //First compute all landmarks to all nodes and vice versa
    for(double landmark : landmarks){
        std::cout<<"Iteration LandmarkToNodes: "<<i<<std::endl;
        landmarkToNodes[i] = DijkstraToALL(myGraph, landmark);
        i++;
    }
    return landmarkToNodes;
}
std::vector<std::vector<double>> precomputeNodeToLandmarks(Graph myGraph, const std::vector<std::vector<double>> landmarkToNodes){
    std::vector<std::vector<double>> nodesToLandmarks(landmarkToNodes[0].size(), std::vector<double>(landmarkToNodes.size()));
    //First compute all landmarks to all nodes and vice versa
    for(int i = 0; i<landmarkToNodes.size(); i++){
        for (int j = 0; j<landmarkToNodes[0].size();j++){
            nodesToLandmarks[j][i]=landmarkToNodes[i][j];
        }

    }
    return nodesToLandmarks;
}
std::vector<std::vector<double>> precomputePotentials(Graph myGraph, const std::vector<double>& landmarks) {
    std::vector<std::vector<double>> potentials;

    // Precompute landmark-to-nodes potentials
    std::vector<std::vector<double>> landmarkToNodes = precomputeLandmarkToNodes(myGraph, landmarks);

    // Precompute node-to-landmarks potentials
    std::vector<std::vector<double>> nodeToLandmarks = precomputeNodeToLandmarks(myGraph, landmarkToNodes);

    // Precompute combined potentials
    for (Node node : myGraph.getNodes()) {
        std::vector<double> nodePotentials;
        for (double landmark : landmarks) {
            double potentialPlus = landmarkToNodes[landmark][node.nodeId];
            double potentialMinus = nodeToLandmarks[node.nodeId][landmark];
            double combinedPotential = std::max(potentialPlus, potentialMinus);
            nodePotentials.push_back(combinedPotential);
        }
        potentials.push_back(nodePotentials);
    }

    return potentials;
}



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
            double currentNode = apq.popMin();
            visited.insert(currentNode);

            double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
            double endEdge = myGraph.edgeStarts[currentNode];

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

std::vector<double> computeFurthestLandmarks(Graph myGraph, int numLandmarks, std::string filename) {
    std::vector<double> landmarks;
    std::vector<Node> nodes = myGraph.getNodes();
    std::unordered_set<double> selectedNodes;  // To keep track of selected landmarks

    // Select the first landmark node as the most southwestern one
    int firstLandmark = 0;
    for (Node node : myGraph.nodes){
        if (node.coordinateX+node.coordinateY < myGraph.nodes[firstLandmark].coordinateX+myGraph.nodes[firstLandmark].coordinateY){
            firstLandmark = node.nodeId;
        }
    }
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
    //save the landmarks to a file
    std::ofstream file(filename);

    for (double landmark : landmarks) {
        file << landmark << std::endl;
    }
    file.close();
    std::cout<<"Writing landmarks to "<<filename+"coord"<<std::endl;
    std::ofstream file2(filename+"coord");

    for (double landmark : landmarks) {
        file2 << myGraph.getNode(landmark).coordinateX << " "<<myGraph.getNode(landmark).coordinateY << std::endl;
    }
    file2.close();
    return landmarks;
}

std::vector<double> selectLandmarks(Graph myGraph, double numLandmarks) {
    std::vector<double> landmarks;

    for (double i = 0; i < numLandmarks; i++) {
        double maxDegreeNode = myGraph.findMaxDegreeNode();
        landmarks.push_back(maxDegreeNode);

        // Update the outDegree of the selected landmark and its neighbors to avoid them in the next iteration
        myGraph.getNode(maxDegreeNode).outDegree = 0;
        for (double neighborId : myGraph.getNode(maxDegreeNode).neighbors) {
            myGraph.getNode(neighborId).outDegree--;
        }
    }

    return landmarks;
}

std::vector<double> avoidLandmarkSelection(Graph myGraph, int numLandmarks) {
    std::vector<double> landmarks;
    std::vector<Node> nodes = myGraph.getNodes();
    std::unordered_set<double> selectedNodes;  // To keep track of selected landmarks

    // Randomly select the first landmark node
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> dist(0, nodes.size() - 1);
    int firstLandmark = dist(rng);
    landmarks.push_back(firstLandmark);
    selectedNodes.insert(firstLandmark);

    // Select the remaining landmarks
    for (int i = 1; i < numLandmarks; i++) {
        // Calculate the weights and sizes of nodes
        std::vector<double> weights(nodes.size(), 0.0);
        std::vector<double> sizes(nodes.size(), 0.0);

        // Grow a shortest-path tree Tr from a random node r
        int r = dist(rng);
        std::queue<int> nodeQueue;
        nodeQueue.push(r);
        std::vector<bool> visited(nodes.size(), false);
        visited[r] = true;

        while (!nodeQueue.empty()) {
            int currentNode = nodeQueue.front();
            nodeQueue.pop();

            // Calculate the lower bound d(v, r) for each node v
            double lowerBound = std::numeric_limits<double>::infinity();
            for (double landmark : landmarks) {
                double landmarkDistance = distance(nodes[landmark], nodes[currentNode]);
                if (landmarkDistance < lowerBound) {
                    lowerBound = landmarkDistance;
                }
            }

            // Calculate the weight of each node v
            double nodeWeight = distance(nodes[currentNode], nodes[r]) - lowerBound;
            weights[currentNode] = nodeWeight;

            // Calculate the size of each node v
            sizes[currentNode] = nodeWeight;

            // Traverse the neighbors of the current node
            for (int neighbor : myGraph.getNeighbors(currentNode)) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    nodeQueue.push(neighbor);
                }
            }
        }

        // Set the size of nodes containing a landmark to zero
        for (double landmark : landmarks) {
            sizes[landmark] = 0.0;
        }

        // Pick the node with maximum size as the next landmark
        double maxNodeSize = -1.0;
        int nextLandmark = -1;
        for (int node = 0; node < nodes.size(); node++) {
            double nodeSize = sizes[node];
            if (nodeSize > maxNodeSize) {
                maxNodeSize = nodeSize;
                nextLandmark = node;
            }
        }

        // Add the next landmark to the set
        landmarks.push_back(nextLandmark);
        selectedNodes.insert(nextLandmark);
    }

    return landmarks;
}
#endif //PRAKTIKUM_ALT_LANDMARKS_HH
