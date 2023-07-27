#ifndef PRAKTIKUM_ASTARBIDIRECTIONAL_HH
#define PRAKTIKUM_ASTARBIDIRECTIONAL_HH
#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>
double AStarBidirectional(Graph& myGraph, double& sourceNode, double& targetNode) {
    //Get two priority queues, visited vectors and dist arrays
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);

    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);

    // Init the minPath variable to infinity. If the current bestPath gets higher than this variable, we terminate
    double minPath = INT_MAX;
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();


        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                double forwardH = distance(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode-1]);

                double forwardF = distForward[forwardEdge] + forwardH;
                if(distForward[forwardEdge]+forwardWeight+distBackward[forwardEdge] < minPath){minPath=distForward[forwardEdge]+distBackward[forwardEdge];}
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

                double backwardH = distance(myGraph.nodes[backwardEdge], myGraph.nodes[sourceNode-1]);
                double backwardF = distBackward[backwardEdge] + backwardH;
                if(distForward[backwardEdge]+backwardWeight+distBackward[backwardEdge] < minPath){minPath=distForward[backwardEdge]+distBackward[backwardEdge];}

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
        if (minPath<=distForward[forwardNode]+ distance(myGraph.nodes[forwardNode], myGraph.nodes[backwardNode])+distBackward[backwardNode]) {
            return minPath;
        }
    }

    return minPath;
}

double AStarBidirectionalSaving(Graph& myGraph, double& sourceNode, double& targetNode, std::string exploredFileName) {
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> forwardPath(myGraph.nodes.size(), -1);
    std::vector<double> backwardPath(myGraph.nodes.size(), -1);
    double minPath = INT_MAX;
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;
    std::vector<double> forwardVisited;
    std::vector<double> backwardVisited;

    std::ofstream exploredForwardFile("experiments/"+exploredFileName+"_forward");
    std::ofstream exploredBackwardFile("experiments/"+exploredFileName+"_backward");
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();
        backwardVisited.push_back(backwardNode);
        forwardVisited.push_back(forwardNode);

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                forwardPath[forwardEdge] = forwardNode;

                double forwardH = distance(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode-1]);
                double forwardF = distForward[forwardEdge] + forwardH;
                if(distForward[forwardEdge]+distBackward[forwardEdge] < minPath){
                    minPath=distForward[forwardEdge]+distBackward[forwardEdge];
                }


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

                double backwardH = distance(myGraph.nodes[backwardEdge], myGraph.nodes[sourceNode-1]);
                double backwardF = distBackward[backwardEdge] + backwardH;
                if(distForward[backwardEdge]+distBackward[backwardEdge] < minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }


                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
        if (minPath<=distForward[forwardNode]+distBackward[backwardNode]){
            return minPath;
        }
    }

    for (auto forwardNode : forwardVisited){
        exploredForwardFile << myGraph.getNode(forwardNode).coordinateX << " " << myGraph.getNode(forwardNode).coordinateY << std::endl;  // Write explored node to the file
    }
    for (auto backwardNode : backwardVisited){
        exploredBackwardFile << myGraph.getNode(backwardNode).coordinateX << " " << myGraph.getNode(backwardNode).coordinateY << std::endl;  // Write explored node to the file

    }
    exploredForwardFile.close();
    exploredBackwardFile.close();

    return minPath;
}


double AStarBidirectionalSearchSpace(Graph& myGraph, double& sourceNode, double& targetNode) {
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<double> visited;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> forwardPath(myGraph.nodes.size(), -1);
    std::vector<double> backwardPath(myGraph.nodes.size(), -1);
    double minPath = INT_MAX;
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    visited.push_back(sourceNode);
    visited.push_back(targetNode);
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();


        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            visited.push_back(forwardEdge);
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                forwardPath[forwardEdge] = forwardNode;

                double forwardH = distance(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode-1]);
                double forwardF = distForward[forwardEdge] + forwardH;
                if(distForward[forwardEdge]+distBackward[forwardEdge] < minPath){
                    minPath=distForward[forwardEdge]+distBackward[forwardEdge];
                }


                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            visited.push_back(backwardEdge);
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                double backwardH = distance(myGraph.nodes[backwardEdge], myGraph.nodes[sourceNode-1]);
                double backwardF = distBackward[backwardEdge] + backwardH;
                if(distForward[backwardEdge]+distBackward[backwardEdge] < minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }


                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
        if (apqBackward.getMin().second >= minPath || apqForward.getMin().second >= minPath){
            return visited.size();
        }
    }





    return visited.size();
}


#endif