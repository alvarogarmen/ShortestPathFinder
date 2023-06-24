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
            //In case we don't relax and the element is not in the APQ already, put it in
            /*if (visited.find(edge-1) == visited.end()){
                apq.insertNode(edge-1, dist[edge-1]);
            }*/
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

double BidirectionalDijkstra(Graph myGraph, double sourceNode, double targetNode) {
    // Get priority queues for forward and backward searches
    APQ forwardAPQ = APQ();
    APQ backwardAPQ = APQ();

    // Keep track of visited nodes in forward and backward searches
    std::set<double> forwardVisited;
    std::set<double> backwardVisited;

    // Initialize distances to infinity for both searches
    std::vector<double> forwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> backwardDist(myGraph.nodes.size(), INT_MAX);

    // Push the source node with distance 0 into the forward APQ
    forwardAPQ.insertNode(sourceNode - 1, 0);
    forwardDist[sourceNode - 1] = 0;

    // Push the target node with distance 0 into the backward APQ
    backwardAPQ.insertNode(targetNode - 1, 0);
    backwardDist[targetNode - 1] = 0;

    double shortestPath = INT_MAX;

    // Main loop
    while (!forwardAPQ.isEmpty() && !backwardAPQ.isEmpty()) {
        // Active node for forward search
        double forwardNode = forwardAPQ.popMin();
        forwardVisited.insert(forwardNode);

        // Active node for backward search
        double backwardNode = backwardAPQ.popMin();
        backwardVisited.insert(backwardNode);

        // Start and end edges for forward search
        double forwardStartEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double forwardEndEdge = myGraph.edgeStarts[forwardNode];

        // Start and end edges for backward search
        double backwardStartEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double backwardEndEdge = myGraph.edgeStarts[backwardNode];

        // Relax edges for forward search
        for (double forwardEdgeIndex = forwardStartEdge; forwardEdgeIndex <= forwardEndEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (forwardDist[forwardNode] + forwardWeight < forwardDist[forwardEdge]) {
                forwardDist[forwardEdge] = forwardDist[forwardNode] + forwardWeight;

                // Check if the node is already visited in the backward search
                if (backwardVisited.find(forwardEdge) != backwardVisited.end()) {
                    double pathLength = forwardDist[forwardEdge] + backwardDist[forwardEdge];
                    if (pathLength < shortestPath) {
                        shortestPath = pathLength;
                        // Return immediately if shortest path is found
                        return shortestPath;
                    }
                }

                if (forwardAPQ.contains(forwardEdge)) {
                    forwardAPQ.decreaseKey(forwardEdge, forwardDist[forwardEdge]);
                } else {
                    forwardAPQ.insertNode(forwardEdge, forwardDist[forwardEdge]);
                }
            }
        }

        // Relax edges for backward search
        for (double backwardEdgeIndex = backwardStartEdge; backwardEdgeIndex <= backwardEndEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (backwardDist[backwardNode] + backwardWeight < backwardDist[backwardEdge]) {
                backwardDist[backwardEdge] = backwardDist[backwardNode] + backwardWeight;

                // Check if the node is already visited in the forward search
                if (forwardVisited.find(backwardEdge) != forwardVisited.end()) {
                    double pathLength = forwardDist[backwardEdge] + backwardDist[backwardEdge];
                    if (pathLength < shortestPath) {
                        shortestPath = pathLength;
                        // Return immediately if shortest path is found
                        return shortestPath;
                    }
                }

                if (backwardAPQ.contains(backwardEdge)) {
                    backwardAPQ.decreaseKey(backwardEdge, backwardDist[backwardEdge]);
                } else {
                    backwardAPQ.insertNode(backwardEdge, backwardDist[backwardEdge]);
                }
            }
        }
    }

    std::cout << "No path found from source: " << sourceNode << " to target: " << targetNode << std::endl;
    return -1; // To symbolize that there is no path
}


#endif