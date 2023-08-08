//
// Created by alvar on 29/07/2023.
//

#ifndef PRAKTIKUM_ALTBIDIRECTIONALSTOPPING_HH
#define PRAKTIKUM_ALTBIDIRECTIONALSTOPPING_HH

#include "Graph.hh"
#include "APQ.hh"
#include "Dijkstra.hh"
#include "ALT.hh"

double ALTBidirectionalStopping(const Graph& myGraph, double sourceNode, double targetNode, const std::vector<std::vector<double>>& potentials) {
    APQ forwardAPQ = APQ(myGraph.nodeCount);
    APQ backwardAPQ = APQ(myGraph.nodeCount);
    std::vector<double> forwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> backwardDist(myGraph.nodes.size(), INT_MAX);

    forwardAPQ.insertNode(sourceNode - 1, 0);
    backwardAPQ.insertNode(targetNode - 1, 0);
    forwardDist[sourceNode - 1] = 0;
    backwardDist[targetNode - 1] = 0;

    double bestPath = INT_MAX;

    while (!forwardAPQ.isEmpty() && !backwardAPQ.isEmpty()) {
        // Perform a forward search step
        double forwardNode = forwardAPQ.popMin();
        if(backwardDist[forwardNode] < INT_MAX){
            return forwardDist[forwardNode]+backwardDist[forwardNode];}

        double forwardStartEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double forwardEndEdge = myGraph.edgeStarts[forwardNode];

        for (double edgeIndex = forwardStartEdge; edgeIndex <= forwardEndEdge; edgeIndex++) {
            double forwardEdge = myGraph.edges[edgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (forwardDist[forwardNode] + forwardWeight < forwardDist[forwardEdge]) {
                forwardDist[forwardEdge] = forwardDist[forwardNode] + forwardWeight;
                double forwardF = forwardDist[forwardEdge] + estimate(forwardEdge, targetNode-1, potentials);
                if (forwardAPQ.contains(forwardEdge)) {
                    forwardAPQ.decreaseKey(forwardEdge, forwardF);
                } else {
                    forwardAPQ.insertNode(forwardEdge, forwardF);
                }



            }
        }

        // Perform a backward search step
        double backwardNode = backwardAPQ.popMin();
        if(forwardDist[backwardNode] < INT_MAX){
            return forwardDist[backwardNode]+backwardDist[backwardNode];}


        double backwardStartEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double backwardEndEdge = myGraph.edgeStarts[backwardNode];

        for (double edgeIndex = backwardStartEdge; edgeIndex <= backwardEndEdge; edgeIndex++) {
            double backwardEdge = myGraph.edges[edgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);
            if (backwardDist[backwardNode] + backwardWeight < backwardDist[backwardEdge]) {
                backwardDist[backwardEdge] = backwardDist[backwardNode] + backwardWeight;
                double backwardF = backwardDist[backwardEdge] + estimate(backwardEdge, sourceNode-1, potentials);
                if (backwardAPQ.contains(backwardEdge)) {
                    backwardAPQ.decreaseKey(backwardEdge, backwardF);
                } else {
                    backwardAPQ.insertNode(backwardEdge, backwardF);
                }


            }
        }




    }

    return bestPath;
}

double ALTBidirectionalStoppingSearchSpace(Graph& myGraph, double& sourceNode, double& targetNode,
                                   std::vector<std::vector<double>>& potentials) {

    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> visited;
    double minPath = INT_MAX;
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    visited.push_back(sourceNode);
    visited.push_back(targetNode);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        visited.push_back(forwardNode);

        double backwardNode = apqBackward.popMin();
        visited.push_back(backwardNode);
        if(distForward[backwardNode] < INT_MAX){
            return visited.size();}
        if(distBackward[forwardNode] < INT_MAX){
            return visited.size();}

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                double forwardF = distForward[forwardEdge]+ estimate(forwardEdge, targetNode-1, potentials);



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

                double backwardF = distBackward[backwardEdge]+ estimate(backwardEdge, sourceNode-1, potentials);



                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }

            }
        }

    }

    return visited.size();
}

#endif //PRAKTIKUM_ALTBIDIRECTIONALSTOPPING_HH
