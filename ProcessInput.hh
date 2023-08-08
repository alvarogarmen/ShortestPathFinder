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
#include "AStarBidirectional.hh"
#include "DijkstraBidirectional.hh"
#include "ALTBidirectionalStopping.hh"
#include "AStarBidirectionalStopping.hh"

std::vector<double> loadLandmarks(std::string filename){
    std::vector<double> vec;
    std::ifstream file("experiments/"+filename);
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
    DijkstraSaving(myGraph, sourceNode, targetNode);

    start = std::chrono::high_resolution_clock::now();
    Dis = BidirectionalDijkstra(myGraph, sourceNode, targetNode);
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    std::cout<<"Bidirectional Dijkstra from Node "<<sourceNode<<" to Node "<<targetNode<<" is: "<<Dis<<std::endl;
    std::cout<<"Bidirectional Dijkstra took: "<<time.count()<<"s"<<std::endl;
    BidirectionalDijkstraSaving(myGraph, sourceNode, targetNode, "DijkstraBidi_explored");
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
    AStarSaving(myGraph, sourceNode, targetNode);

    start=std::chrono::high_resolution_clock::now();
    Dis = AStarBidirectional(myGraph, sourceNode, targetNode);
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    std::cout<<"A* Bidirectional: "<<Dis<<std::endl;
    std::cout<<"Took: "<<time.count()<<"s"<<std::endl;
    AStarBidirectionalSaving(myGraph, sourceNode, targetNode, "AStar_Bidirectional_explored");

    start=std::chrono::high_resolution_clock::now();
    Dis = AStarBidirectionalStopping(myGraph, sourceNode, targetNode);
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    std::cout<<"Variant: "<<Dis<<std::endl;
    std::cout<<"Took: "<<time.count()<<"s"<<std::endl;
}


/*void callALTMaxDegree(Graph myGraph, double sourceNode, double targetNode, int numLandmarks){
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
}*/

void callALTAvoid(Graph myGraph, double sourceNode, double targetNode, int numLandmarks, int newLandmarks, std::string filename, int secureBidirectional){
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> potentials;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time;
    if(newLandmarks==0){
        std::vector<double> landmarks=loadLandmarks("Landmarks_Avoid_"+filename.substr(13,3));
        potentials = precomputePotentialsEuclidian(myGraph, landmarks);
        std::cout<<"Copy: "<<potentials.size()<<" "<<potentials[0].size()<<std::endl;
        std::cout<<"Precomputation Complete"<<std::endl;

        // Start measuring the execution time
        start = std::chrono::high_resolution_clock::now();
        // Run ALT using the precomputed potentials
        double shortestDistance = ALT(myGraph, sourceNode, targetNode, potentials);

        // Stop measuring the execution time
        end = std::chrono::high_resolution_clock::now();

        // Compute the duration in milliseconds
        time = end - start;

        // Print the shortest distance and execution time
        std::cout<<"ALT Avoid from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<shortestDistance<<std::endl;
        std::cout << "ALT Avoid took: " << time.count()<< "s" << std::endl;
        ALTSaving(myGraph, sourceNode, targetNode, potentials, "ALTAvoid_explored");
        std::cout<<"Saved"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        double ALTBI = ALTBidirectional(myGraph, sourceNode, targetNode, potentials);
        end = std::chrono::high_resolution_clock::now();
        time = end - start;
        std::cout<<"Avoid ALTBI: "<<ALTBI<<std::endl;
        std::cout<<"Took :"<<time.count()<<std::endl;
        ALTBidirectionalSaving(myGraph, sourceNode, targetNode, "ALTAvoidBidirectional", potentials);
    }
    else if(newLandmarks==1){
        std::cout<<"No loading"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> landmarksAvoid =avoidLandmarkSelection(myGraph, numLandmarks, "Landmarks_Avoid_"+filename.substr(13,3));
        for (double land : landmarksAvoid){
            std::cout<<land<<std::endl;
        }

        std::vector<std::vector<double>> potentials = precomputePotentialsEuclidian(myGraph, landmarksAvoid);
        std::cout<<"Copy: "<<potentials.size()<<" "<<potentials[0].size()<<std::endl;
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        std::cout<<"ALT Avoid preprocessing took: "<<time.count()<<std::endl;
        std::cout<<"Precomputation Complete"<<std::endl;

        // Start measuring the execution time
        start = std::chrono::high_resolution_clock::now();
        // Run ALT using the precomputed potentials
        double shortestDistance = ALT(myGraph, sourceNode, targetNode, potentials);

        // Stop measuring the execution time
        end = std::chrono::high_resolution_clock::now();

        // Compute the duration in milliseconds
        time = end - start;

        // Print the shortest distance and execution time
        std::cout<<"ALT Avoid from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<shortestDistance<<std::endl;
        std::cout << "ALT Avoid took: " << time.count()<< "s" << std::endl;
        ALTSaving(myGraph, sourceNode, targetNode, potentials, "ALTAvoid_explored");
        std::cout<<"Saved"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        double ALTBI = ALTBidirectional(myGraph, sourceNode, targetNode, potentials);
        end = std::chrono::high_resolution_clock::now();
        time = end - start;
        std::cout<<"Avoid ALTBI: "<<ALTBI<<std::endl;
        std::cout<<"Took :"<<time.count()<<std::endl;
        ALTBidirectionalSaving(myGraph, sourceNode, targetNode, "ALTAvoidBidirectional", potentials);
    }

}



void callALTFarthest(Graph myGraph, double sourceNode, double targetNode, int numLandmarks, int newLandmarks, std::string filename){
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> potentials;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time;
    if(newLandmarks==0){
        std::vector<double> landmarks=loadLandmarks("Landmarks_Farthest_"+filename.substr(13,3));
        potentials = precomputePotentialsEuclidian(myGraph, landmarks);
        std::cout<<"Copy: "<<potentials.size()<<" "<<potentials[0].size()<<std::endl;
        std::cout<<"Precomputation Complete"<<std::endl;

        // Start measuring the execution time
        start = std::chrono::high_resolution_clock::now();
        // Run ALT using the precomputed potentials
        double shortestDistance = ALT(myGraph, sourceNode, targetNode, potentials);

        // Stop measuring the execution time
        end = std::chrono::high_resolution_clock::now();

        // Compute the duration in milliseconds
        time = end - start;

        // Print the shortest distance and execution time
        std::cout<<"ALT Farthest from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<shortestDistance<<std::endl;
        std::cout << "ALT Farthest took: " << time.count()<< "s" << std::endl;
        ALTSaving(myGraph, sourceNode, targetNode, potentials, "ALTFarthest_explored");
        std::cout<<"Saved"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        double ALTBI = ALTBidirectional(myGraph, sourceNode, targetNode, potentials);
        end = std::chrono::high_resolution_clock::now();
        time = end - start;
        std::cout<<"Farthest ALTBI: "<<ALTBI<<std::endl;
        std::cout<<"Took :"<<time.count()<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        ALTBidirectionalStopping(myGraph, sourceNode, targetNode, potentials);
        end = std::chrono::high_resolution_clock::now();
        time = end - start;
        std::cout<<"Variant took :"<<time.count()<<std::endl;

        std::cout<<"Variant:"<<ALTBidirectionalStopping(myGraph, sourceNode, targetNode, potentials)<<std::endl;

        ALTBidirectionalSaving(myGraph, sourceNode, targetNode, "ALTAFarthestBidirectional", potentials);
    }
    else if(newLandmarks==1){
        std::cout<<"No loading"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> landmarksFarthest =farthestLandmarkSelection(myGraph, numLandmarks, "Landmarks_Farthest_"+
                                                                         filename.substr(13,3)+"_"+std::to_string(numLandmarks));

        std::vector<std::vector<double>> potentials = precomputePotentialsEuclidian(myGraph, landmarksFarthest);
        std::cout<<"Copy: "<<potentials.size()<<" "<<potentials[0].size()<<std::endl;
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        std::cout<<"ALT Farthest preprocessing took: "<<time.count()<<std::endl;
        std::cout<<"Precomputation Complete"<<std::endl;

        // Start measuring the execution time
        start = std::chrono::high_resolution_clock::now();
        // Run ALT using the precomputed potentials
        double shortestDistance = ALT(myGraph, sourceNode, targetNode, potentials);

        // Stop measuring the execution time
        end = std::chrono::high_resolution_clock::now();

        // Compute the duration in milliseconds
        time = end - start;

        // Print the shortest distance and execution time
        std::cout<<"ALT Farthest from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<shortestDistance<<std::endl;
        std::cout << "ALT Farthest took: " << time.count()<< "s" << std::endl;
        ALTSaving(myGraph, sourceNode, targetNode, potentials, "ALTFarthest_explored");
        std::cout<<"Saved"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        double ALTBI = ALTBidirectional(myGraph, sourceNode, targetNode, potentials);
        end = std::chrono::high_resolution_clock::now();
        time = end - start;
        std::cout<<"Farthest ALTBI: "<<ALTBI<<std::endl;
        std::cout<<"Took :"<<time.count()<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        ALTBidirectionalStopping(myGraph, sourceNode, targetNode, potentials);
        end = std::chrono::high_resolution_clock::now();
        time = end - start;
        std::cout<<"Variant took :"<<time.count()<<std::endl;

        std::cout<<"Variant:"<<ALTBidirectionalStopping(myGraph, sourceNode, targetNode, potentials)<<std::endl;
        ALTBidirectionalSaving(myGraph, sourceNode, targetNode, "ALTAFarthestBidirectional", potentials);
    }

}
void plotALTFarthest(Graph myGraph, double sourceNode, double targetNode, int numLandmarks, int newLandmarks, std::string filename, int secureBidirectional){
    std::vector<double> landmarks=loadLandmarks("Landmarks_Farthest_"+filename.substr(13,3));
    std::vector<std::vector<double>> potentials = precomputePotentialsEuclidian(myGraph, landmarks);

    ALTSaving(myGraph, sourceNode, targetNode, potentials, "ALTFarthest_explored");

    ALTBidirectionalSaving(myGraph, sourceNode, targetNode, "ALTAFarthestBidirectional", potentials);
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
