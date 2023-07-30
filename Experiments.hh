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

#include "ALT-Landmarks.hh"
#include "ProcessInput.hh"

std::vector<std::pair<int, int>> generatePoints(Graph myGraph, int numPoints, std::string filename){

    std::vector<std::pair<int, int>> pairs;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(1, myGraph.nodes.size());
    std::string numPairs = std::to_string(numPoints);
    std::ofstream outputFile("experiments/"+filename.substr(13, 3)+"_"+numPairs);

    for (int i = 0; i < numPoints; ++i) {
        int randomInt1 = distribution(gen);
        int randomInt2 = distribution(gen);
        pairs.push_back(std::make_pair(randomInt1, randomInt2));
        outputFile<<randomInt1 << " " << randomInt2 << std::endl;

    }
    outputFile.close();
    return pairs;
}

std::vector<std::pair<int, int>> loadPoints(const std::string& filename) {
    std::vector<std::pair<int, int>> pairs;

    std::ifstream inputFile("experiments/"+filename);
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
    std::ofstream dijkstraRanks("experiments/dijkstraRanks_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream dijkstraTimes("experiments/dijkstraTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream dijkstraSearch("experiments/dijkstraSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        Dijkstra(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        dijkstraRanks<<DijkstraRank(myGraph, sourceNode, targetNode)<<std::endl;
        dijkstraSearch<<DijkstraSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;
        dijkstraTimes<<time.count()<<std::endl;
    }
    dijkstraRanks.close();
    dijkstraTimes.close();
    dijkstraSearch.close();

    //Now Bidirectional
    std::ofstream bidiDijkstraTimes("experiments/dijkstraBidiTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream bidiDijkstraSearch("experiments/dijkstraBidiSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        BidirectionalDijkstra(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        bidiDijkstraSearch<<BidirectionalDijkstraSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;
        bidiDijkstraTimes<<time.count()<<std::endl;
    }
    bidiDijkstraSearch.close();
    bidiDijkstraTimes.close();

}

void callExperimentAStar(Graph myGraph, std::vector<std::pair<int, int>> Points, std::string filename){
    std::ofstream aStarTimes("experiments/aStarTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream aStarSearch("experiments/aStarSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        AStar(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        aStarTimes<<time.count()<<std::endl;
        aStarSearch<<AStarSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;
    }
    aStarSearch.close();
    aStarTimes.close();

    //Now Bidirectional
    std::ofstream bidiAStarTimes("experiments/aStarBidiTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream bidiAStarSearch("experiments/aStarBidiSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        AStarBidirectional(myGraph, sourceNode, targetNode);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        bidiAStarTimes<<time.count()<<std::endl;
        bidiAStarSearch<<AStarBidirectionalSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;

    }
    bidiAStarSearch.close();
    bidiAStarTimes.close();

}


void callExperimentALT(Graph myGraph, std::vector<std::pair<int, int>> Points, std::string filename,
                       int numLandmarks, int newLandmarks){
    std::ofstream preTimes("experiments/ALTFarthestPreTimes_"+filename.substr(13,3));
    if (newLandmarks!=0) {
        for (int i = 2; i <= numLandmarks;) {
            auto proStart = std::chrono::high_resolution_clock::now();
            std::vector<double> landmarksFarthest = farthestLandmarkSelection(myGraph, i, "Landmarks_Farthest_" +
                                                                                          filename.substr(13, 3)+
                                                                                          std::to_string(i));
            auto proEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> proTime = proEnd - proStart;
            preTimes << proTime.count() << " ";
            proStart = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<double>> potentials = precomputePotentialsEuclidian(myGraph, landmarksFarthest);
            proEnd = std::chrono::high_resolution_clock::now();
            proTime = proEnd - proStart;
            preTimes << proTime.count() << std::endl;
            writePotentialsToFile(potentials,
                                  "Potentials_Farthest_" + filename.substr(13, 3) + std::to_string(potentials.size()));
            i *= 2;
        }
        preTimes.close();
    }
    int i = 2;
    while(i<=numLandmarks) {
        std::vector<std::vector<double>> potentials = loadPotentials(
                "Potentials_Farthest_" + filename.substr(13, 3) + std::to_string(i));
        std::cout<<"Num Landmarks: "<<potentials.size()<<" "<<potentials[1].size()<<std::endl;
        std::vector<double> landmarks = loadLandmarks(
                "Landmarks_Farthest_" + filename.substr(13, 3) + std::to_string(i));
        //TODO: Fix writing and loading landmarks
        std::cout<<"Num Landmarks: "<<landmarks.size()<<std::endl;
        std::ofstream ALTFarthestTimes(
                "experiments/ALTFarthestTimes_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

        std::ofstream ALTFarthestSearch("experiments/ALTFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));
        for (auto experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            auto start = std::chrono::high_resolution_clock::now();
            ALT(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time = end - start;
            ALTFarthestTimes << time.count() << std::endl;
            ALTFarthestSearch << ALTSearchSpace(myGraph, sourceNode, targetNode, potentials)<<std::endl;

        }
        i=i*2;
        ALTFarthestSearch.close();
        ALTFarthestTimes.close();
    }
    //Now Bidirectional
    i=2;
    while(i<=numLandmarks) {
        std::vector<std::vector<double>> potentials = loadPotentials(
                "Potentials_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        std::vector<double> landmarks = loadLandmarks(
                "Landmarks_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        std::ofstream ALTBidiFarthestTimes(
                "experiments/ALTBidiFarthestTimes_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

        std::ofstream ALTBidiFarthestSearch("experiments/ALTBidiFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                                std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

        for (auto experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            auto start = std::chrono::high_resolution_clock::now();
            ALTBidirectional(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time = end - start;
            ALTBidiFarthestTimes << time.count() << std::endl;
            ALTBidiFarthestSearch << ALTBidirectionalSearchSpace(myGraph, sourceNode, targetNode, potentials)<<std::endl;
        }
        i=i*2;
        ALTBidiFarthestTimes.close();
        ALTBidiFarthestSearch.close();
    }

    //Now Bidirectional with new breaking condition (0.001% error)
    i=2;
    while(i<=numLandmarks) {
        std::vector<std::vector<double>> potentials = loadPotentials(
                "Potentials_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        std::vector<double> landmarks = loadLandmarks(
                "Landmarks_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        std::ofstream ALTBidiStopFarthestTimes(
                "experiments/ALTBidiStopFarthestTimes_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

        std::ofstream ALTBidiStopFarthestSearch("experiments/ALTBidiStopFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                                            std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

        for (auto experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            auto start = std::chrono::high_resolution_clock::now();
            ALTBidirectionalStopping(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time = end - start;
            ALTBidiStopFarthestTimes << time.count() << std::endl;
            ALTBidiStopFarthestSearch << ALTBidirectionalStoppingSearchSpace(myGraph, sourceNode, targetNode, potentials)<<std::endl;
        }
        i=i*2;
        ALTBidiStopFarthestTimes.close();
        ALTBidiStopFarthestSearch.close();
    }


}


#endif //PRAKTIKUM_EXPERIMENTS_HH
