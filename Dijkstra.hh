//
// Created by alvar on 26/04/2023.
//

#ifndef UNTITLED_DIJKSTRA_HH
#define UNTITLED_DIJKSTRA_HH

#include "Graph.hh"
#include <limits.h>
#include "APQ.hh"
#include <cmath>
#include <set>

double distance(Node source, Node target) {                             //Calculate Euclidian distance
    return sqrt(std::pow(source.coordinateX - target.coordinateX, 2)
          + pow(source.coordinateY - target.coordinateY, 2));
}
double Dijkstra(Graph myGraph, double sourceNode, double targetNode){
    //Get a priority queue
    APQ apq = APQ();

    //Keep this to know if we have visited a node or not
    std::set<int> visited;

    //Initialize the distances to infinity
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    //Push the source node with distance 0 into the APQ
    apq.insertNode(sourceNode, 0);
    dist.at(sourceNode-1)=0;

    //Main loop
    while(!apq.isEmpty()){
        std::cout<<"Main Loop"<<std::endl;
        //Active node. Note that if its already in the set, inserting doesn't do anything
        int currentNode = apq.getMin().first;       //The getter doesn't pop the minimum element
        visited.insert(currentNode);
        apq.popMin();                               //It is popped here

        int startEdge = (currentNode-1 > 0) ? myGraph.edgeStarts[currentNode-2] : 0;    //Ternary operation: if we start with the node 0, we'll start with
        int endEdge = myGraph.edgeStarts[currentNode];                                  //the first edge, otherwise, it is taken from the adjacency vector
        //Look for edges to relax
        for (int edgeIndex = startEdge; edgeIndex < endEdge; edgeIndex++) {             //This selects only the relevant edges
            int edge = myGraph.edges[edgeIndex];
            // Calculate the weight as the Euclidean distance between the nodes
            std::cout<<"With current Node: "<< currentNode-1<<" and destination Node: "<< edge-1<< ".Weight: ";
            double weight = distance(myGraph.nodes[currentNode-1], myGraph.nodes[edge-1]);
            std::cout<<weight<<std::endl;
            // Relax the edge if a shorter path is found
            if (dist[currentNode-1] + weight < dist[edge-1]) {
                dist[edge-1] = dist[currentNode-1] + weight;
                //Decrease key operation if the node is already in APQ
                if(apq.contains(edge-1)) {
                    apq.decreaseKey(edge-1, dist[edge-1]);
                    std::cout << "Key decreased" << std::endl;
                }
                //Otherwise, insert it into the APQ
                else{apq.insertNode(edge-1, dist[edge-1]);}

                //Break if we have reached our destination
                if (edge == targetNode){
                    std::cout<<"edgeIndex==targetNode arrived"<<std::endl;
                    return dist[edge-1];
                }
            }
            //In case we don't relax and the element is not in the APQ already, put it in
            if (visited.find(edge-1) == visited.end()){
                std::cout<<"Insertion in the apq"<<std::endl;
                apq.insertNode(edge-1, dist[edge-1]);
            }
        }
    }
    std::cout<<"All edges relaxed"<<std::endl;
    if (dist[targetNode-1]==INT_MAX){
        return -1;                      //To symbolize that it is not reached
    }
    return dist[targetNode-1];          //I think it works, but I haven't tested it yet
}

#endif