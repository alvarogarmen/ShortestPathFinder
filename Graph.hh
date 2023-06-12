//
// Created by alvar on 26/04/2023.
//

#ifndef UNTITLED_GRAPH_HH
#define UNTITLED_GRAPH_HH

#include <vector>
#include "Node.hh"

class Graph{
public:
    Graph();
    void insertNode(Node node);
    void insertEdge(double node);
    std::vector<Node> nodes;            //Vector with all the Node objects
    std::vector<double> edges;          //Vector with all the Edge destinations
    std::vector<double> edgeStarts;     //Vector with all the outdegrees of the Nodes
    double nodeCount;                   //I don't actually use this
    double edgeCount;                   //Nor this. nodes.size() does the trick
    Node& getNode(double nodeId);       //Getter by reference. Useful when reading the files
    std::vector<Node> getNodes();
    double findMaxDegreeNode();
    std::vector<double> getNeighbors(double nodeId);
    void addEdge(double node, double target, double shortcutDistance);

};
Graph::Graph() {}

void Graph::insertNode(Node node) {
    this->nodes.push_back(node);
}

void Graph::insertEdge(double edge) {
    this->edges.push_back(edge);
}

Node& Graph::getNode(double nodeId) {
    return this->nodes.at(nodeId);
}
std::vector<Node> Graph::getNodes() {
    return this->nodes;
}

std::vector<double> Graph::getNeighbors(double nodeId) {
    std::vector<double> neighbors;
    double startEdge = (nodeId > 0) ? edgeStarts[nodeId - 1] + 1 : 0;
    double endEdge = edgeStarts[nodeId];

    for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
        double neighborId = edges[edgeIndex] - 1;
        neighbors.push_back(neighborId);
    }

    return neighbors;
}

double Graph::findMaxDegreeNode() {
    double maxDegreeNode = 0;
    double maxDegree = std::numeric_limits<double>::min();

    for (auto& node : this->nodes) {
        if (node.outDegree > maxDegree) {
            maxDegreeNode = node.nodeId;
            maxDegree = node.outDegree;
        }
    }

    return maxDegreeNode;
}
void Graph::addEdge(double node, double target, double shortcutDistance) {
    // Update the edgeStarts vector to reflect the new edge
    double nodeIndex = node - 1;
    double outDegree = (nodeIndex < edgeStarts.size()) ? edgeStarts[nodeIndex + 1] - edgeStarts[nodeIndex] : edges.size();
    edgeStarts.insert(edgeStarts.begin() + nodeIndex + 1, edgeStarts[nodeIndex] + outDegree + 1);

    // Insert the target node and shortcutDistance into the edges vector
    edges.insert(edges.begin() + edgeStarts[nodeIndex] + outDegree, target);

    // Update the outDegree of the node
    if (nodeIndex < nodes.size())
        nodes[nodeIndex].outDegree++;

    // Update the nodeCount and edgeCount
    nodeCount = nodes.size();
    edgeCount = edges.size();
}
#endif