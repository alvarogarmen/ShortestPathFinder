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
    void buildEdgeStarts();
    std::vector<Node> nodes;
    std::vector<double> edges;
    std::vector<double> edgeStarts;
    double nodeCount;
    double edgeCount;
    Node& getNode(double nodeId);
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
void Graph::buildEdgeStarts() {
    for (Node node:this->nodes){
        this->edgeStarts.push_back(node.start);
    }
}

#endif