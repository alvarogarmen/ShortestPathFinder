//
// Created by alvar on 09/07/2023.
//

#ifndef PRAKTIKUM_EXPERIMENTS_HH
#define PRAKTIKUM_EXPERIMENTS_HH

#include <string>
#include <random>
#include <chrono>
#include "Graph.hh"
#include "readgraph.hh"
#include "Dijkstra.hh"
#include "AStar.hh"
#include "DijkstraBidirectional.hh"
#include "ALT.hh"
#include "ALT_Bidirectional.hh"
#include "ALT-Landmarks.hh"
#include "ProcessInput.hh"

std::vector<std::pair<int, int>> generatePoints(Graph myGraph, int numPoints, std::string filename){

    std::vector<std::pair<int, int>> pairs;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(1, myGraph.nodes.size());
    std::string numPairs = std::to_string(numPoints);
    std::ofstream outputFile(filename.substr(13, 3)+"_"+numPairs);

    for (int i = 0; i < numPoints; ++i) {
        int randomInt1 = distribution(gen);
        int randomInt2 = distribution(gen);
        pairs.push_back(std::make_pair(randomInt1, randomInt2));
        outputFile<<randomInt1 << " " << randomInt2 << std::endl;
    }
    outputFile.close();
    return pairs;
}

std::vector<std::pair<int, int>> readPairsFromFile(const std::string& filename) {
    std::vector<std::pair<int, int>> pairs;

    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Error opening the file " << filename << std::endl;
        return pairs;
    }

    int value1, value2;
    while (inputFile >> value1 >> value2) {
        pairs.push_back(std::make_pair(value1, value2));
    }

    inputFile.close();
    return pairs;
}

void callExperimentDijkstra(Graph myGraph, std::vector<std::pair<int, int>> Points, std::string filename){
    std::ofstream dijkstraTimes("dijkstraTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream dijkstraSearch("dijkstraSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (std::pair experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        Dijkstra(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        dijkstraSearch<<DijkstraSearchSpace(myGraph, sourceNode, targetNode);
        dijkstraTimes<<time.count()<<std::endl;
    }

    //Now Bidirectional
    std::ofstream bidiDijkstraTimes("dijkstraBidiTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream bidiDijkstraSearch("dijkstraBidiSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (std::pair experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        BidirectionalDijkstra(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        bidiDijkstraSearch<<BidirectionalDijkstraSearchSpace(myGraph, sourceNode, targetNode);
        bidiDijkstraTimes<<time.count()<<std::endl;
    }

}

void callExperimentAStar(Graph myGraph, std::vector<std::pair<int, int>> Points, std::string filename){
    std::ofstream aStarTimes("aStarTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream aStarSearch("aStarSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (std::pair experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        AStar(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        aStarTimes<<time.count()<<std::endl;
        aStarSearch<<AStarSearchSpace(myGraph, sourceNode, targetNode);
    }

    //Now Bidirectional
    std::ofstream bidiAStarTimes("aStarBidiTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream bidiAStarSearch("aStarBidiSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (std::pair experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        AStarBidirectional(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        bidiAStarSearch<<AStarBidirectionalSearchSpace(myGraph, sourceNode, targetNode);
        bidiAStarTimes<<time.count()<<std::endl;
    }

}


void callExperimentALT(Graph myGraph, std::vector<std::pair<int, int>> Points, std::string filename, int numLandmarks){
    std::ofstream preTimes("ALTFarthestPreTimes_"+filename.substr(13,3));
    for (int i = 2; i<=numLandmarks;){
        auto proStart = std::chrono::high_resolution_clock::now();
        std::vector<double> landmarksFarthest =farthestLandmarkSelection(myGraph, i, "Landmarks_Farthest_"+filename.substr(13,3));
        auto proEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> proTime = proEnd-proStart;
        preTimes<<proTime.count()<<" ";
        proStart = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> potentials = precomputePotentialsEuclidian(myGraph, landmarksFarthest);
        proEnd =  std::chrono::high_resolution_clock::now();
        proTime = proEnd - proStart;
        preTimes<<proTime.count()<<std::endl;
        writePotentialsToFile(potentials, "Potentials_Farthest_"+filename.substr(13,3)+std::to_string(potentials.size()));
        i+=2;
    }
    int i = 2;
    while(i<=numLandmarks) {
        std::vector<std::vector<double>> potentials = loadPotentials(
                "Potentials_Farthest" + filename.substr(13, 3) + std::to_string(i));
        std::ofstream ALTFarthestTimes(
                "ALTFarthestTimes_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));
        std::ofstream ALTFarthestSearch("ALTFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));
        for (std::pair experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            auto start = std::chrono::high_resolution_clock::now();
            ALT(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time = end - start;
            ALTFarthestTimes << time.count() << std::endl;
            ALTFarthestSearch << ALTSearchSpace(myGraph, sourceNode, targetNode, potentials);
        }
        i+=2;
    }
    //Now Bidirectional
    i=2;
    while(i<=numLandmarks) {
        std::vector<std::vector<double>> potentials = loadPotentials(
                "Potentials_Farthest" + filename.substr(13, 3) + std::to_string(i));
        std::ofstream ALTBidiFarthestTimes(
                "ALTBidiFarthestTimes_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));
        std::ofstream ALTBidiFarthestSearch("ALTBidiFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                                std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));
        for (std::pair experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            auto start = std::chrono::high_resolution_clock::now();
            ALTBidirectional(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time = end - start;
            ALTBidiFarthestTimes << time.count() << std::endl;
            ALTBidiFarthestSearch << ALTBidirectionalSearchSpace(myGraph, sourceNode, targetNode, potentials);
        }
        i+=2;
    }


}

void loadLandmarksExperiment(std::string filename, int numlandmarks){

}
#endif //PRAKTIKUM_EXPERIMENTS_HH