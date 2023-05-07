//
// Created by alvar on 26/04/2023.
//

#ifndef UNTITLED_DIJKSTRA_HH
#define UNTITLED_DIJKSTRA_HH

#include "Graph.hh"
#include <limits.h>
#include "APQ.hh"
#include <cmath>

//TODO: implement distance measurement
double distance(Node target, Node source) {
    return sqrt(std::pow(source.coordinateX - target.coordinateX, 2) - pow(source.coordinateY - target.coordinateY, 2));
}
double Dijkstra(Graph myGraph, double sourceNode, double targetNode){


    APQ apq = APQ();

    //Initialize the distances to infinity
    double size = myGraph.nodes.size();
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);

    //Push the source node with distance 0 into the APQ
    apq.insertNode(sourceNode, 0);
    dist.at(sourceNode-1)=0;
    std::cout<<std::endl;
    for(int z=0; z<9;z++){

        //std::cout<<dist[z]<<"-"<<std::endl;
        std::cout<<"Node"<<z<<std::endl;
        std::cout<<myGraph.getNode(z).coordinateX<<myGraph.getNode(z).coordinateY;
    }
    std::cout<<"LOOK"<<distance(myGraph.getNode(1), myGraph.getNode(0));
    //Main loop
    while(!apq.isEmpty()){
        double u = apq.getMin().second;
        apq.popMin();
        for (int i = 0; i < size; i++){
            for (int j = myGraph.nodes.at(i).start; j<myGraph.edges.size(); j++){
                if(dist.at(u)+ distance(myGraph.nodes.at(u), myGraph.nodes.at(i))<dist.at(i)) {
                    dist.at(i) = dist.at(u) + distance(myGraph.nodes.at(u), myGraph.nodes.at(i));
                    apq.insertNode(i, dist.at(i));
                    if (myGraph.nodes.at(i).nodeId==targetNode-1){
                        std::cout<<"break"<<std::endl;
                        std::cout<<dist.at(targetNode-1);
                        break;
                    }
                }

            }
        }


    }
    //std::cout<<"Node "<<targetNode<< " Distance from Source:"<<dist.at(targetNode-1)<<std::endl;
    //std::cout<<"Distance check"<<distance(myGraph.getNode(1), myGraph.getNode(2))<<std::endl;

    return dist.at(targetNode-1);

}

#endif