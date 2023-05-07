//
// Created by alvar on 03/05/2023.
//

#ifndef UNTITLED_NODE_HH
#define UNTITLED_NODE_HH

#include <vector>

class Node{
public:
    Node(){
        this->coordinateX=-1;
        this->coordinateY=-1;
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
    }
    Node(double CoordinateX, double CoordinateY, double NodeId){
        this->coordinateX = CoordinateX;
        this->coordinateY = CoordinateY;
        this->nodeId = NodeId;
        this->start = 0;
    }
    double getStart() {
        return this->start;
    }
    void setStart(double newStart) {
        this->start=newStart;
    }
    std::vector<double> getCoordinates() {
        return std::vector<double> {this->coordinateX, this->coordinateY};
    }
    double coordinateX, coordinateY, nodeId;
    double start;

};

#endif






