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
#endif