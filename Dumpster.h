//
// Created by alvar on 05/06/2023.
//
/*
#ifndef PRAKTIKUM_DUMPSTER_H
#define PRAKTIKUM_DUMPSTER_H

#include <cmath>
#include <unordered_set>
#include "Graph.hh"
#include "APQ.hh"
void BiBiALT (Graph& g, int source, int target, std::vector<std::vector<double>>& vec) {

    MinHeap Qf;
    MinHeap Qb;

    std::vector<double> df(g.getNodes().size(), INFINITY);
    std::vector<int> parentf(g.getNodes().size());

    std::vector<double> db(g.getNodes().size(), INFINITY);
    std::vector<int> parentb(g.getNodes().size());

    std::vector<bool> visitedT(g.getNodes().size(), false);
    std::vector<bool> visitedS(g.getNodes().size(), false);

    visitedS[source] = true;
    parentf[source] = source;
    df[source] = 0;
    Qf.insert(source,0);

    visitedT[target] = true;
    parentb[target] = target;
    db[target] = 0;
    Qb.insert(target,0);

    double mu = INFINITY;

    std::unordered_set<int> Sf;
    std::unordered_set<int> Sb;

    while (!Qf.isEmpty() && !Qb.isEmpty()){
        std::pair<int,double> uf = Qf.deleteMin();
        std::pair<int,double> ub = Qb.deleteMin();

        Sf.insert(uf.first);
        Sb.insert(ub.first);

        for (int v = 0; v < g.getNodes()[uf.first + 1].second - g.getNodes()[uf.first].second; v++) {

            int nodeV = g.getEdges()[g.getNodes()[uf.first].second+v].getId();
            double eucDist = euclideanDistance(g.getNodes()[uf.first].first, g.getNodes()[nodeV].first);

            if(df[uf.first] + eucDist< df[nodeV]){
                parentf[nodeV] = uf.first;
                df[nodeV] = df[uf.first] + eucDist;
                if(!visitedS[nodeV]) {
                    Qf.insert(nodeV, eucDist + df[uf.first] + getPotential(target, nodeV, vec));
                    visitedS[nodeV] = true;
                } else {
                    Qf.decreaseKey(nodeV, eucDist + df[uf.first] + getPotential(target, nodeV, vec));
                }
            }
            if(visitedS[nodeV] && df[uf.first] + eucDist + db[nodeV] < mu) {
                mu = df[uf.first] + eucDist + db[nodeV];
            }
        }

        for (int v = 0; v < g.getNodes()[ub.first + 1].second - g.getNodes()[ub.first].second; v++) {

            int nodeV = g.getEdges()[g.getNodes()[ub.first].second+v].getId();
            double eucDist = euclideanDistance(g.getNodes()[ub.first].first, g.getNodes()[nodeV].first);

            if(db[ub.first] + eucDist < db[nodeV]){
                parentb[nodeV] = ub.first;
                db[nodeV] = db[ub.first] + eucDist;
                if(!visitedT[nodeV]) {
                    Qb.insert(nodeV, eucDist + db[ub.first] + getPotential(source, nodeV, vec));
                    visitedT[nodeV] = true;
                } else {
                    Qb.decreaseKey(nodeV, eucDist + db[ub.first] + getPotential(source, nodeV, vec));
                }
            }
            if(visitedT[nodeV] && db[ub.first] + eucDist + df[nodeV] < mu) {
                mu = db[ub.first] + eucDist + df[nodeV];
            }
        }

        /*if (visitedT[uf.first] && visitedS[ub.first]) {
                std::cout<<"distance has been found with Bidirectional Dijkstra between " << source << " and "  << target << " it is " << mu <<std::endl;
                return;
        }*/

if(Qf.top().second >= mu || Qb.top().second >= mu) {
std::cout<<"distance has been found with Bidirectional Dijkstra between " << source << " and "  << target << " it is " << mu <<std::endl;
return;
}

 }
    return;
}
/*
double AStarBidirectionalSaving(Graph myGraph, double sourceNode, double targetNode, const std::string& exploredFileName, const std::string& pathFileName) {
    APQ forwardAPQ = APQ();
    APQ backwardAPQ = APQ();
    std::set<double> forwardVisited;
    std::set<double> backwardVisited;
    std::vector<double> forwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> backwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> forwardPath(myGraph.nodes.size(), -1);
    std::vector<double> backwardPath(myGraph.nodes.size(), -1);

    forwardAPQ.insertNode(sourceNode - 1, 0);
    backwardAPQ.insertNode(targetNode - 1, 0);
    forwardDist[sourceNode - 1] = 0;
    backwardDist[targetNode - 1] = 0;

    double bestPath = INT_MAX;

    std::ofstream exploredNodeFile(exploredFileName);
    std::ofstream pathFile(pathFileName);

    while (!forwardAPQ.isEmpty() && !backwardAPQ.isEmpty()) {
        double forwardCurrentNode = forwardAPQ.popMin();
        double backwardCurrentNode = backwardAPQ.popMin();

        forwardVisited.insert(forwardCurrentNode);
        backwardVisited.insert(backwardCurrentNode);



        double forwardStartEdge = (forwardCurrentNode > 0) ? myGraph.edgeStarts[forwardCurrentNode - 1] + 1 : 0;
        double forwardEndEdge = myGraph.edgeStarts[forwardCurrentNode];

        double backwardStartEdge = (backwardCurrentNode > 0) ? myGraph.edgeStarts[backwardCurrentNode - 1] + 1 : 0;
        double backwardEndEdge = myGraph.edgeStarts[backwardCurrentNode];

        // Forward search
        for (double forwardEdgeIndex = forwardStartEdge; forwardEdgeIndex <= forwardEndEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardCurrentNode], myGraph.nodes[forwardEdge]);

            if (forwardDist[forwardCurrentNode] + forwardWeight < forwardDist[forwardEdge]) {
                forwardDist[forwardEdge] = forwardDist[forwardCurrentNode] + forwardWeight;
                forwardPath[forwardEdge] = forwardCurrentNode;

                if (backwardVisited.count(forwardEdge) > 0) {
                    double pathCost = forwardDist[forwardEdge] + backwardDist[forwardEdge];
                    if (pathCost < bestPath) {
                        bestPath = pathCost;
                    }
                }

                double forwardH = distance(myGraph.nodes[forwardEdge], myGraph.nodes[targetNode - 1]);
                double forwardF = forwardDist[forwardEdge] + forwardH;

                if (forwardAPQ.contains(forwardEdge)) {
                    forwardAPQ.decreaseKey(forwardEdge, forwardF);
                } else {
                    forwardAPQ.insertNode(forwardEdge, forwardF);
                }
            }
        }

        // Backward search
        for (double backwardEdgeIndex = backwardStartEdge; backwardEdgeIndex <= backwardEndEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardCurrentNode], myGraph.nodes[backwardEdge]);

            if (backwardDist[backwardCurrentNode] + backwardWeight < backwardDist[backwardEdge]) {
                backwardDist[backwardEdge] = backwardDist[backwardCurrentNode] + backwardWeight;
                backwardPath[backwardEdge] = backwardCurrentNode;

                if (forwardVisited.count(backwardEdge) > 0) {
                    double pathCost = backwardDist[backwardEdge] + forwardDist[backwardEdge];
                    if (pathCost < bestPath) {
                        bestPath = pathCost;
                    }
                }

                double backwardH = distance(myGraph.nodes[backwardEdge], myGraph.nodes[sourceNode - 1]);
                double backwardF = backwardDist[backwardEdge] + backwardH;

                if (backwardAPQ.contains(backwardEdge)) {
                    backwardAPQ.decreaseKey(backwardEdge, backwardF);
                } else {
                    backwardAPQ.insertNode(backwardEdge, backwardF);
                }
            }
        }
        //Check if shortest path already found
        if (bestPath != INT_MAX) {
            break;
        }
    }

    std::cout << "All available edges relaxed" << std::endl;
    if (bestPath == INT_MAX) {
        return -1;
    }

    double currentNode = targetNode - 1;
    while (currentNode != sourceNode - 1) {
        pathFile << myGraph.nodes[currentNode].coordinateX << " " << myGraph.nodes[currentNode].coordinateY << "\n";
        currentNode = forwardPath[currentNode];
    }
    pathFile << myGraph.nodes[sourceNode - 1].coordinateX << " " << myGraph.nodes[sourceNode - 1].coordinateY << "\n";

    exploredNodeFile.close();
    pathFile.close();

    return bestPath;
}*/