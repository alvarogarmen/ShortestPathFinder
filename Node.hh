//
// Created by alvar on 03/05/2023.
//

#ifndef UNTITLED_NODE_HH
#define UNTITLED_NODE_HH

#include <vector>

class Node{
public:
    Node(){                             //A bunch of different constructors, needed at different points of reading the graph
        this->coordinateX=-1;
        this->coordinateY=-1;
        this->outDegree = -1;
        this->nodeId = -1;
    };
    Node(double degree, double coordinateX, double coordinateY, double nodeId) {
        this->outDegree = degree;
        this->coordinateX = coordinateX;
        this->coordinateY = coordinateY;
        this->nodeId = nodeId;
    }
    Node(double CoordinateX, double CoordinateY){
        this->coordinateX=CoordinateX;
        this->coordinateY=CoordinateY;
        this->outDegree = -1;
        this->nodeId = -1;
    }
    Node(double CoordinateX, double CoordinateY, double NodeId){
        this->coordinateX = CoordinateX;
        this->coordinateY = CoordinateY;
        this->nodeId = NodeId;
        this->outDegree = 0;
    }
    void setStart(double newStart) {
        this->outDegree=newStart;
    }
    void addNeighbor(double neighbor){
        this->neighbors.push_back(neighbor);
        this->outDegree++;
    }
    double coordinateX, coordinateY, nodeId;
    double outDegree;
    std::vector<double> neighbors;

};

#endif






