//
// Created by alvar on 26/04/2023.
//

#ifndef UNTITLED_DIJKSTRA_HH
#define UNTITLED_DIJKSTRA_HH

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "APQ_Saving.hh"
#include <cmath>
#include <set>

double distance(const Node& source, const Node& target) {                             //Calculate Euclidian distance
    return sqrt(std::pow(static_cast<double>(source.coordinateX - target.coordinateX), 2.0)
                + std::pow(static_cast<double>(source.coordinateY - target.coordinateY), 2.0));
}
double Dijkstra(const Graph& myGraph, double& sourceNode, double& targetNode){
    if (sourceNode==targetNode){
        return 0;
    }
    //Get a priority queue
    APQ apq = APQ(myGraph.nodeCount);

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
    if (dist[targetNode-1]==INT_MAX){
        std::cout<<"WHAT"<<std::endl;
        return -1;                      //To symbolize that it is not reached
    }
    return dist[targetNode-1];          //I think it works, but I haven't tested it yet
}

double DijkstraSaving(Graph myGraph, double sourceNode, double targetNode) {
    APQSaving apq;
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    std::ofstream exploredFile("experiments/explored_nodes");  // File stream for explored nodes
    std::ofstream pathFile("experiments/path");      // File stream for path


    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();
        visited.insert(currentNode);
        exploredFile << myGraph.getNode(currentNode).coordinateX << " " << myGraph.getNode(currentNode).coordinateY << std::endl;  // Write explored node to the file


        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            // Calculate the weight as the Euclidean distance between the nodes
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);
            // Relax the edge if a shorter path is found
            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;
                // Break if we have reached our destination
                if (edge == targetNode - 1) {
                    std::vector<double> path;
                    path.push_back(edge);
                    double prevNode = currentNode;
                    while (prevNode != -1 && prevNode != sourceNode-1) {
                        path.push_back(prevNode);
                        prevNode = apq.getPrev(prevNode);
                    }
                    path.push_back(prevNode);

                    // Write the path to the file in reverse order
                    for (int i = 0; i < path.size()-1; i++) {
                        pathFile << myGraph.getNode(path[i]).coordinateX << " " << myGraph.getNode(path[i]).coordinateY
                                 << " "<<path[i]<< " " << std::endl;
                    }

                    exploredFile.close();  // Close the files before returning
                    pathFile.close();
                    return dist[edge];
                }

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge], currentNode);
                } else {
                    apq.insertNode(edge, dist[edge], currentNode);
                }
            }
        }
    }

    exploredFile.close();  // Close the files
    pathFile.close();

    std::cout << "All available edges relaxed" << std::endl;

    if (dist[targetNode - 1] == INT_MAX) {
        return -1;
    }
    return dist[targetNode - 1];
}

double DijkstraRank(Graph myGraph, double sourceNode, double targetNode) {
    APQSaving apq;
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();
        visited.insert(currentNode);

        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            // Calculate the weight as the Euclidean distance between the nodes
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);
            // Relax the edge if a shorter path is found
            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;
                // Break if we have reached our destination
                if (edge == targetNode - 1) {
                    std::vector<double> path;
                    path.push_back(edge);
                    double prevNode = currentNode;
                    while (prevNode != -1 && prevNode != sourceNode-1) {
                        path.push_back(prevNode);
                        prevNode = apq.getPrev(prevNode);
                    }
                    path.push_back(prevNode);

                    return path.size()-1;
                }

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge], currentNode);
                } else {
                    apq.insertNode(edge, dist[edge], currentNode);
                }
            }
        }
    }
    if (dist[targetNode - 1] == INT_MAX) {
        return -1;
    }
    return -1;
}

std::vector<double> DijkstraToALL(Graph myGraph, double sourceNode) {
    //Get a priority queue
    APQ apq = APQ(myGraph.nodeCount);

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


    return dist;
}

double DijkstraSearchSpace(Graph myGraph, double sourceNode, double targetNode) {
    APQSaving apq;
    std::vector<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    apq.insertNode(sourceNode - 1, 0, -1);
    dist[sourceNode - 1] = 0;

    visited.push_back(sourceNode);
    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();


        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            visited.push_back(edge);
            // Calculate the weight as the Euclidean distance between the nodes
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);
            // Relax the edge if a shorter path is found
            if (dist[currentNode] + weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;
                // Break if we have reached our destination
                if (edge == targetNode - 1) {



                    return visited.size();
                }

                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, dist[edge], currentNode);
                } else {
                    apq.insertNode(edge, dist[edge], currentNode);
                }
            }
        }
    }

    std::cout << "All available edges relaxed" << std::endl;

    if (dist[targetNode - 1] == INT_MAX) {
        return visited.size();
    }
    return visited.size();
}




#endif