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
        this->start = -1;
        this->nodeId = -1;
    };
    Node(double degree, double coordinateX, double coordinateY, double nodeId) {
        this->start = degree;
        this->coordinateX = coordinateX;
        this->coordinateY = coordinateY;
        this->nodeId = nodeId;
    }
    Node(double CoordinateX, double CoordinateY){
        this->coordinateX=CoordinateX;
        this->coordinateY=CoordinateY;
        this->start = -1;
        this->nodeId = -1;
    }
    Node(double CoordinateX, double CoordinateY, double NodeId){
        this->coordinateX = CoordinateX;
        this->coordinateY = CoordinateY;
        this->nodeId = NodeId;
        this->start = 0;
    }
    void setStart(double newStart) {
        this->start=newStart;
    }
    double coordinateX, coordinateY, nodeId;
    double start;                                   //Start is the outdegree of the node

};

#endif






