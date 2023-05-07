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

//TODO: implement distance measurement
double distance(Node source, Node target) {
    return sqrt(std::pow(source.coordinateX - target.coordinateX, 2)
          + pow(source.coordinateY - target.coordinateY, 2));
}
double Dijkstra(Graph myGraph, double sourceNode, double targetNode){

    for (int i :myGraph.edgeStarts){
        std::cout<<i<<std::endl;
    }

    APQ apq = APQ();
    std::set<int> visited;
    //Initialize the distances to infinity

    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    //Push the source node with distance 0 into the APQ
    apq.insertNode(sourceNode, 0);
    std::cout<<sourceNode<<std::endl;
    dist.at(sourceNode-1)=0;
    //Main loop
    while(!apq.isEmpty()){
        std::cout<<"Loop"<<std::endl;
        int currentNode = apq.getMin().first;
        visited.insert(currentNode);
        apq.popMin();

        int startEdge = (currentNode-1 > 0) ? myGraph.edgeStarts[currentNode-2] : 0;
        int endEdge = myGraph.edgeStarts[currentNode];
        for (int edgeIndex = startEdge; edgeIndex < endEdge; edgeIndex++) {
            std::cout<<"EdgeIndex: "<<edgeIndex<<std::endl;
            int edge = myGraph.edges[edgeIndex];
            std::cout<<"Current edge "<< edge<<std::endl;
            // Calculate the weight as the Euclidean distance between the nodes
            std::cout<<"With current Node: "<< currentNode-1<<" and destination Node: "<< edge-1<< ".Weight: ";
            double weight = distance(myGraph.nodes[currentNode-1], myGraph.nodes[edge-1]);
            std::cout<<weight<<std::endl;
            // Relax the edge if a shorter path is found
            if (dist[currentNode-1] + weight < dist[edge-1]) {
                dist[edge-1] = dist[currentNode-1] + weight;
                for(auto g:dist){
                    std::cout<<g<<", ";
                }
                //std::cout<<"Relaxation"<<std::endl;
                if(apq.contains(edge-1)) {
                    apq.decreaseKey(edge-1, dist[edge-1]);
                    std::cout << "Key decreased" << std::endl;
                }
                else{apq.insertNode(edge-1, dist[edge-1]);}
                std::cout<<edge-1<<"son iguales?"<<targetNode-1<<std::endl;
                if (edge == targetNode){
                    std::cout<<"edgeIndex==targetNode arrived"<<std::endl;
                    return dist[edge-1];
                }
            }
            if (visited.find(edge-1) == visited.end()){
                std::cout<<"Insertion in the apq"<<std::endl;
                apq.insertNode(edge-1, dist[edge-1]);
            }
        }
    }
    std::cout<<"All edges relaxed"<<std::endl;
    for (auto i:dist){
        std::cout<<i<<std::endl;
    }
    return dist[targetNode-1];
}

#endif

/*for (int i = 0; i < size; i++){
            for (int j = myGraph.nodes.at(i).start; j<myGraph.edges.size(); j++){
                if(dist.at(u)+ distance(myGraph.nodes.at(u), myGraph.nodes.at(i))<dist.at(i)) {
                    dist.at(i) = dist.at(u) + distance(myGraph.nodes.at(u), myGraph.nodes.at(i));
                    apq.insertNode(i, dist.at(i));
                    if (myGraph.nodes.at(i).nodeId==targetNode-1){


                        //break;
                    }std::cout<<"Line 51"<<std::endl;
                }std::cout<<"Line 52"<<std::endl;

            }std::cout<<"Line 54"<<std::endl;
        }


    }
    std::cout<<targetNode<< "line 58" <<std::endl;
    //std::cout<<"Node "<<targetNode<< " Distance from Source:"<<dist.at(targetNode-1)<<std::endl;
    //std::cout<<"Distance check"<<distance(myGraph.getNode(1), myGraph.getNode(2))<<std::endl;

    return dist.at(targetNode-1);
*/