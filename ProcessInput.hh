//
// Created by alvar on 11/06/2023.
//

#ifndef PRAKTIKUM_PROCESSINPUT_HH
#define PRAKTIKUM_PROCESSINPUT_HH

#include <chrono>
#include <random>
#include "Graph.hh"
#include "Dijkstra.hh"
#include "AStar.hh"
#include "ALT-Landmarks.hh"
#include "ALT.hh"
#include "ALT_Saving.h"
#include "ALT_Bidirectional.hh"
#include "HighwayHierarchiesStar.hh"

void callDijkstra(Graph myGraph, double sourceNode, double targetNode){
    // Call Dijkstra and outDegree timer
    auto start = std::chrono::high_resolution_clock::now();
    double Dis = Dijkstra(myGraph, sourceNode, targetNode);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the distance and the time
    std::cout<<"Distance from Node "<<sourceNode<<" to Node "<< targetNode<<" is: "<<Dis<<std::endl;
    std::cout<<"Dijkstra took: "<<time.count()<<"s"<<std::endl;
}

void callAStar(Graph myGraph, double sourceNode, double targetNode){
    auto start = std::chrono::high_resolution_clock::now();
    double Dis = AStar(myGraph, sourceNode, targetNode);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the distance and the time
    std::cout<<"A* from Node "<<sourceNode<<" to Node "<< targetNode<<" is: "<<Dis<<std::endl;
    std::cout<<"A* took: "<<time.count()<<"s"<<std::endl;
}

void callALTMaxDegree(Graph myGraph, double sourceNode, double targetNode, int numLandmarks){
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> Landmarks = selectLandmarks(myGraph, numLandmarks);
    std::unordered_map<int, double> landmarkDistances = computeLandmarkDistances(myGraph, Landmarks, numLandmarks);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the time
    std::cout<<"ALT Avoid preprocessing took:"<<time.count()<<"s"<<std::endl;

    //ALT query with MaxDegree (avoiding) Landmarks
    start = std::chrono::high_resolution_clock::now();
    double Dis = ALT(myGraph, sourceNode, targetNode, landmarkDistances);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"ALT MaxDegree from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<Dis<<std::endl;
    std::cout<<"ALT MaxDegree query took: "<<time.count()<<std::endl;
    //Save the results
    ALTSaving(myGraph, sourceNode, targetNode, landmarkDistances, "Max_explored_nodes.txt", "Max_path.txt");
}

void callALTAvoid(Graph myGraph, double sourceNode, double targetNode, int numLandmarks){
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> Landmarks = avoidLandmarkSelection(myGraph, numLandmarks);
    std::unordered_map<int, double> landmarkDistances = computeLandmarkDistances(myGraph, Landmarks, numLandmarks);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the time
    std::cout<<"ALT Avoid preprocessing took:"<<time.count()<<"s"<<std::endl;

    //ALT query with MaxDegree (avoiding) Landmarks
    start = std::chrono::high_resolution_clock::now();
    double Dis = ALT(myGraph, sourceNode, targetNode, landmarkDistances);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"ALT Avoid from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<Dis<<std::endl;
    std::cout<<"ALT Avoid query took: "<<time.count()<<std::endl;
    //Save the results
    ALTSaving(myGraph, sourceNode, targetNode, landmarkDistances, "Avoid_nodes.txt", "Avoid_path.txt");
}
std::vector<double> callComputeFurthestLandmarks(Graph myGraph, int numLandmarks, std::string graph){
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> Landmarks = computeFurthestLandmarks(myGraph, numLandmarks, "landmarksbelgium");
    std::unordered_map<int, double>landmarkDistances = computeLandmarkDistances(myGraph, Landmarks, numLandmarks);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the time
    std::cout<<"ALT Furthest preprocessing took:"<<time.count()<<"s"<<std::endl;
    return Landmarks;
}

std::vector<double> loadLandmarks(std::string graph){
    std::vector<double> vec;
    std::ifstream file("landmarks"+graph+".txt");
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << "landmarks.txt" << std::endl;
        return vec;
    }

    double value;
    while (file >> value) {
        vec.push_back(value);
    }

    file.close();
    return vec;
}
void callALTFurthest(Graph myGraph, double sourceNode, double targetNode, int numLandmarks, std::vector<double> Landmarks){
    std::unordered_map<int, double>landmarkDistances = computeLandmarkDistances(myGraph, Landmarks, numLandmarks);

    //ALT Furthest Query
    auto start = std::chrono::high_resolution_clock::now();
    double Dis = ALT(myGraph, sourceNode, targetNode, landmarkDistances);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the distance and the time
    std::cout<<"ALT Furthest from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<Dis<<std::endl;
    std::cout<<"ALT Furthest query took: "<<time.count()<<std::endl;
    //Save the exploration and path
    ALTSaving(myGraph, sourceNode, targetNode, landmarkDistances, "Furthest_explored_nodes.txt", "Furthest_path.txt");

}

void callALTBiAvoid(Graph myGraph, double sourceNode, double targetNode, int numLandmarks){
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> Landmarks = avoidLandmarkSelection(myGraph, numLandmarks);
    std::unordered_map<int, double> landmarkDistances = computeLandmarkDistances(myGraph, Landmarks, numLandmarks);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the time
    std::cout<<"ALT Avoid preprocessing took:"<<time.count()<<"s"<<std::endl;

    //ALT query with MaxDegree (avoiding) Landmarks
    start = std::chrono::high_resolution_clock::now();
    double Dis = ALTBidirectional(myGraph, sourceNode, targetNode, landmarkDistances);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"ALT Bidirectional Avoid from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<Dis<<std::endl;
    std::cout<<"ALT Bidirectional Avoid query took: "<<time.count()<<std::endl;
    //Save the results
    ALTSaving(myGraph, sourceNode, targetNode, landmarkDistances, "BiAvoid_explored_nodes.txt", "BiAvoid_path.txt");
}

void callHHS(Graph myGraph, double sourceNode, double targetNode, int numLandmarks){
    auto start = std::chrono::high_resolution_clock::now();
    preprocessHighwayHierarchiesStar(myGraph);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the time
    std::cout<<"HH* preprocessing took:"<<time.count()<<"s"<<std::endl;

    //ALT query with MaxDegree (avoiding) Landmarks
    start = std::chrono::high_resolution_clock::now();
    double HHS = highwayHierarchiesStar(myGraph, sourceNode, targetNode);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"HHS from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<HHS<<std::endl;
    std::cout<<"HHS query took: "<<time.count()<<std::endl;
}
#endif //PRAKTIKUM_PROCESSINPUT_HH
