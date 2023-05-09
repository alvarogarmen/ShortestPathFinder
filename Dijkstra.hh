//
// Created by alvar on 26/04/2023.
//

#ifndef UNTITLED_DIJKSTRA_HH
#define UNTITLED_DIJKSTRA_HH

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include <cmath>
#include <set>

double distance(Node source, Node target) {                             //Calculate Euclidian distance
    return sqrt(std::pow(static_cast<double>(source.coordinateX - target.coordinateX), 2.0)
                + std::pow(static_cast<double>(source.coordinateY - target.coordinateY), 2.0));
}
double Dijkstra(Graph myGraph, double sourceNode, double targetNode){
    //Get a priority queue
    APQ apq = APQ();


    //Keep this to know if we have visited a node or not
    std::set<double> visited;

    //Initialize the distances to infinity
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    //Push the source node with distance 0 into the APQ
    apq.insertNode(sourceNode-1, 0);
    dist[sourceNode-1]=0;
    //Main loop
    while(!apq.isEmpty()){
        //Active node. Note that if its already in the set, inserting doesn't do anything
        double currentNode = apq.getMin().first;       //The getter doesn't pop the minimum element
        visited.insert(currentNode);
        apq.popMin();                               //It is popped here

        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode-1]+1 : 0;    //Ternary operation: if we start with the node 0, we'll start with
        double endEdge = myGraph.edgeStarts[currentNode];                                  //the first edge, otherwise, it is taken from the adjacency vector
        //Look for edges to relax
        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {             //This selects only the relevant edges
            double edge = myGraph.edges[edgeIndex]-1;
            // Calculate the weight as the Euclidean distance between the nodes
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);
            // Relax the edge if a shorter path is found
            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;
                //Break if we have reached our destination
                if (edge == targetNode-1){
                    return dist[edge];
                }
                //Decrease key operation if the node is already in APQ
                if(apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge]);
                }
                //Otherwise, insert it into the APQ
                else{apq.insertNode(edge, dist[edge]);}


            }
            //In case we don't relax and the element is not in the APQ already, put it in
            /*if (visited.find(edge-1) == visited.end()){
                apq.insertNode(edge-1, dist[edge-1]);
            }*/
        }
    }
    std::cout<<"All available edges relaxed"<<std::endl;
    if (dist[targetNode-1]==INT_MAX){
        return -1;                      //To symbolize that it is not reached
    }
    return dist[targetNode-1];          //I think it works, but I haven't tested it yet
}

#endif