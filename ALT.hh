// Created by alvaro on 29/05/2023.
//

#ifndef UNTITLED_ALT_HH
#define UNTITLED_ALT_HH

#include "Graph.hh"
#include <climits>
#include "APQ.hh"
#include "Dijkstra.hh"
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>

double estimate(double source, double target, const std::vector<std::vector<double>>& potentials,
                std::vector<double> usefulLandmarks) {
    double potential=0;
    //TODO: Inverted!!
    for (int i : usefulLandmarks) {
        double potentialPlus = potentials[i][target] - potentials[i][source];
        double potentialMinus = potentials[i][source] - potentials[i][target];
        if(potential < potentialMinus or potential < potentialPlus){
            potential = std::max(potentialMinus, potentialPlus);
        }
    }


    return potential;
}
void calculateEndPoints(const Node& A, const Node& B, Node& End1, Node& End2) {
    Node AB = { B.coordinateX - A.coordinateX, B.coordinateY - A.coordinateY };
    double length = sqrt(AB.coordinateX * AB.coordinateX + AB.coordinateY * AB.coordinateY);
    double angle_degrees = 15.0;
    double angle_radians = angle_degrees * M_PI / 180.0;
    double angle_offset_degrees = 7.5;
    double angle_offset_radians = angle_offset_degrees * M_PI / 180.0;

    // Calculate End1
    double end1_angle = atan2(AB.coordinateY, AB.coordinateX) - angle_radians;
    double end1_x = B.coordinateX + length * cos(end1_angle);
    double end1_y = B.coordinateY + length * sin(end1_angle);
    End1 = { end1_x, end1_y };

    // Calculate End2
    double end2_angle = atan2(AB.coordinateY, AB.coordinateX) + angle_radians;
    double end2_x = B.coordinateX + length * cos(end2_angle);
    double end2_y = B.coordinateY + length * sin(end2_angle);
    End2 = { end2_x, end2_y };
}

bool isInsideTriangle(const Node& A, const Node& End1, const Node& End2, const Node& D) {
    Node AB = { End1.coordinateX - A.coordinateX, End1.coordinateY - A.coordinateY };
    Node AD = { D.coordinateX - A.coordinateX, D.coordinateY - A.coordinateY };
    double cross1 = AB.coordinateX * AD.coordinateY - AB.coordinateY * AD.coordinateX;

    Node BC = { End2.coordinateX - End1.coordinateX, End2.coordinateY - End1.coordinateY };
    Node BD = { D.coordinateX - End1.coordinateX, D.coordinateY - End1.coordinateY };
    double cross2 = BC.coordinateX * BD.coordinateY - BC.coordinateY * BD.coordinateX;

    Node CA = { A.coordinateX - End2.coordinateX, A.coordinateY - End2.coordinateY };
    Node CD = { D.coordinateX - End2.coordinateX, D.coordinateY - End2.coordinateY };
    double cross3 = CA.coordinateX * CD.coordinateY - CA.coordinateY * CD.coordinateX;

    return cross1 >= 0 && cross2 >= 0 && cross3 >= 0;
}
std::vector<double> findUsefulLandmarks(Graph& myGraph, const double& sourceNode, const double& targetNode,
                                        const std::vector<double>& landmarks) {
    std::vector<double> usefulLandmarks;
    Node source = myGraph.nodes[sourceNode-1];
    Node target = myGraph.nodes[targetNode-1];
    Node vertex1;
    Node vertex2;
    calculateEndPoints(source, target, vertex1, vertex2);
    for (int i = 0; i<landmarks.size();i++){
        if (isInsideTriangle(source, vertex1, vertex2, myGraph.nodes[landmarks[i]-1])){
            usefulLandmarks.push_back(i);
        }
        if (usefulLandmarks.size()>3){
            return usefulLandmarks;
        }
    }
    if (usefulLandmarks.empty()){
        for (int i = 0; i<landmarks.size(); i++){
            usefulLandmarks.push_back(i);
        }
    }
    return usefulLandmarks;
}


double ALT(Graph& myGraph, double& sourceNode, double& targetNode, const std::vector<std::vector<double>>& potentials,
           std::vector<double>& landmarks) {
    std::vector<double> usefulLandmarks = findUsefulLandmarks(myGraph, sourceNode, targetNode, landmarks);

    APQ apq = APQ();
    std::set<double> visited;
    std::vector<double> dist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityDist(myGraph.nodes.size(), INT_MAX);


    apq.insertNode(sourceNode - 1, 0);
    dist[sourceNode - 1] = 0;
    priorityDist[sourceNode - 1] = 0;

    while (!apq.isEmpty()) {
        double currentNode = apq.popMin();
        visited.insert(currentNode);


        double startEdge = (currentNode > 0) ? myGraph.edgeStarts[currentNode - 1] + 1 : 0;
        double endEdge = myGraph.edgeStarts[currentNode];

        for (double edgeIndex = startEdge; edgeIndex <= endEdge; edgeIndex++) {
            double edge = myGraph.edges[edgeIndex] - 1;
            double weight = distance(myGraph.nodes[currentNode], myGraph.nodes[edge]);

            if (dist[currentNode]+ weight < dist[edge]) {
                dist[edge] = dist[currentNode] + weight;
                priorityDist[edge] = priorityDist[currentNode] - estimate(currentNode, targetNode - 1, potentials, usefulLandmarks) + weight + estimate(edge, targetNode - 1, potentials, usefulLandmarks);
                if (edge == targetNode-1){
                    return dist[edge];
                }
                double f = priorityDist[edge];
                if (apq.contains(edge)) {
                    apq.decreaseKey(edge, f);
                } else {
                    apq.insertNode(edge, f);
                }
            }
        }
    }

    std::cout << "All available edges relaxed" << std::endl;
    if (dist[targetNode - 1] == INT_MAX) {
        return -1;
    }
    return dist[targetNode - 1];
}





#endif
