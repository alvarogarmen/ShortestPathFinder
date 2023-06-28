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

double distance(Node& source, Node& target) {                             //Calculate Euclidian distance
    return sqrt(std::pow(static_cast<double>(source.coordinateX - target.coordinateX), 2.0)
                + std::pow(static_cast<double>(source.coordinateY - target.coordinateY), 2.0));
}
double Dijkstra(Graph& myGraph, double& sourceNode, double& targetNode){
    if (sourceNode==targetNode){
        return 0;
    }
    //Get a priority queue
    APQ apq = APQ();

    //Initialize the distances to infinity
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    //Push the source node with distance 0 into the APQ
    apq.insertNode(sourceNode-1, 0);
    dist[sourceNode-1]=0;
    //Main loop
    while(!apq.isEmpty()){
        //Active node. Note that if its already in the set, inserting doesn't do anything
        double currentNode = apq.popMin();

        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode-1]+1 : 0;    //Ternary operation: if we outDegree with the node 0, we'll outDegree with
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
                    std::cout<<"Break early"<<std::endl;
                    return dist[edge];
                }
                //Decrease key operation if the node is already in APQ
                if(apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge]);
                }
                //Otherwise, insert it into the APQ
                else{apq.insertNode(edge, dist[edge]);}


            }

        }
    }
    std::cout<<"All available edges relaxed, with source: "<< sourceNode<<", target: "<<targetNode<<" and distance: "<< dist[targetNode-1]<<std::endl;
    if (dist[targetNode-1]==INT_MAX){
        std::cout<<"WHAT"<<std::endl;
        return -1;                      //To symbolize that it is not reached
    }
    return dist[targetNode-1];          //I think it works, but I haven't tested it yet
}

std::vector<double> DijkstraToALL(Graph myGraph, double sourceNode) {
    //Get a priority queue
    APQ apq = APQ();

    //Initialize the distances to infinity
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    //Push the source node with distance 0 into the APQ
    apq.insertNode(sourceNode - 1, 0);
    dist[sourceNode - 1] = 0;
    //Main loop
    while (!apq.isEmpty()) {
        //Active node. Note that if its already in the set, inserting doesn't do anything
        double currentNode = apq.popMin();

        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1
                                             : 0;    //Ternary operation: if we outDegree with the node 0, we'll outDegree with
        double endEdge = myGraph.edgeStarts[currentNode];                                  //the first edge, otherwise, it is taken from the adjacency vector
        //Look for edges to relax
        for (double edgeIndex = startEdge;
             edgeIndex <= endEdge; edgeIndex++) {             //This selects only the relevant edges
            double edge = myGraph.edges[edgeIndex] - 1;
            // Calculate the weight as the Euclidean distance between the nodes
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);
            // Relax the edge if a shorter path is found
            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;

                //Decrease key operation if the node is already in APQ
                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge]);
                }
                    //Otherwise, insert it into the APQ
                else { apq.insertNode(edge, dist[edge]); }


            }

        }
    }
    std::cout << "All available edges relaxed, with source: " << sourceNode << std::endl;


    return dist;
}

double BidirectionalDijkstra(Graph& myGraph, double& sourceNode, double& targetNode) {
    APQ apqForward = APQ();
    APQ apqBackward = APQ();
    std::set<double> visitedForward;
    std::set<double> visitedBackward;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    double bestPath = INT_MAX;
    double meetingNode = -1;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        visitedForward.insert(forwardNode);
        visitedBackward.insert(backwardNode);



        if (visitedForward.find(backwardNode)!=visitedForward.end() && visitedBackward.find(forwardNode)!=visitedBackward.end()) {
            return bestPath;
        }

        if (distForward[forwardNode] + distBackward[forwardNode] < bestPath) {
            bestPath = distForward[forwardNode] + distBackward[forwardNode];
            meetingNode = forwardNode;
        }

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                double forwardF = distForward[forwardEdge];

                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;


                double backwardF = distBackward[backwardEdge];

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
    }

    if (meetingNode != -1) {
        double shortestPath = distForward[meetingNode] + distBackward[meetingNode];
        if (shortestPath < bestPath)
            bestPath = shortestPath;
    }

    if (bestPath == INT_MAX) {
        return -1;
    }

    return bestPath;
}



#endif