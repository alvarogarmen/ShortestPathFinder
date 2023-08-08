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
#include "AStarBidirectionalStopping.hh"

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

void callRanks(Graph& myGraph, std::vector<std::pair<int, int>> Points, std::string filename){
    std::ofstream dijkstraRanks("experimentsDone/NewdijkstraRanks_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        dijkstraRanks<<DijkstraRank(myGraph, sourceNode, targetNode)<<std::endl;
    }
    dijkstraRanks.close();
}
void callExperimentDijkstra(Graph myGraph, std::vector<std::pair<int, int>> Points, std::string filename){
    //std::ofstream dijkstraRanks("experiments/dijkstraRanks_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    //std::ofstream dijkstraTimes("experiments/dijkstraTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream dijkstraSearch("experiments/dijkstraSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        //auto start = std::chrono::high_resolution_clock::now();
        //Dijkstra(myGraph, sourceNode, targetNode);
        //Stop the clock
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> time = end - start;
        //dijkstraRanks<<DijkstraRank(myGraph, sourceNode, targetNode)<<std::endl;
        dijkstraSearch<<DijkstraSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;
        //dijkstraTimes<<time.count()<<std::endl;
    }
    //dijkstraRanks.close();
    //dijkstraTimes.close();
    dijkstraSearch.close();

    //Now Bidirectional
    //std::ofstream bidiDijkstraTimes("experiments/dijkstraBidiTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream bidiDijkstraSearch("experiments/dijkstraBidiSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        //auto start = std::chrono::high_resolution_clock::now();
        //BidirectionalDijkstra(myGraph, sourceNode, targetNode);
        //Stop the clock
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> time = end - start;
        bidiDijkstraSearch<<BidirectionalDijkstraSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;
        //bidiDijkstraTimes<<time.count()<<std::endl;
    }
    bidiDijkstraSearch.close();
    //bidiDijkstraTimes.close();

}

void callExperimentAStar(Graph myGraph, std::vector<std::pair<int, int>> Points, std::string filename){
    //std::ofstream aStarTimes("experiments/aStarTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream aStarSearch("experiments/aStarSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
       // auto start = std::chrono::high_resolution_clock::now();
       // AStar(myGraph, sourceNode, targetNode);
        //Stop the clock
       // auto end = std::chrono::high_resolution_clock::now();
       // std::chrono::duration<double> time = end - start;
        //aStarTimes<<time.count()<<std::endl;
        aStarSearch<<AStarSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;
    }
    aStarSearch.close();
    //aStarTimes.close();

    //Now Bidirectional
    //std::ofstream bidiAStarTimes("experiments/aStarBidiTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream bidiAStarSearch("experiments/aStarBidiSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        //auto start = std::chrono::high_resolution_clock::now();
        //AStarBidirectional(myGraph, sourceNode, targetNode);
        //Stop the clock
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> time = end - start;
        //bidiAStarTimes<<time.count()<<std::endl;
        bidiAStarSearch<<AStarBidirectionalSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;

    }
    bidiAStarSearch.close();
    //bidiAStarTimes.close();
    //Now Bidirectional Stopping
    //std::ofstream bidiAStarStopTimes("experiments/aStarBidiStopTimes_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    std::ofstream bidiAStarStopSearch("experiments/aStarBidiStopSearchSpace_"+filename.substr(13,3)+"_"+std::to_string(Points.size()));
    for (auto experiment : Points){
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        //auto start = std::chrono::high_resolution_clock::now();
        //AStarBidirectionalStopping(myGraph, sourceNode, targetNode);
        //Stop the clock
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> time = end - start;
        //bidiAStarStopTimes<<time.count()<<std::endl;
        bidiAStarStopSearch<<AStarBidirectionalStoppingSearchSpace(myGraph, sourceNode, targetNode)<<std::endl;

    }
    bidiAStarStopSearch.close();
    //bidiAStarStopTimes.close();

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
        std::cout<<"Num Landmarks: "<<landmarks.size()<<std::endl;
        //std::ofstream ALTFarthestTimes(
        //        "experiments/ALTFarthestTimes_" + filename.substr(13, 3) + "_" +
        //        std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

        std::ofstream ALTFarthestSearch("experiments/ALTFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));
        for (auto experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            //auto start = std::chrono::high_resolution_clock::now();
            //ALT(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            //auto end = std::chrono::high_resolution_clock::now();
            //std::chrono::duration<double> time = end - start;
            //ALTFarthestTimes << time.count() << std::endl;
            ALTFarthestSearch << ALTSearchSpace(myGraph, sourceNode, targetNode, potentials)<<std::endl;

        }
        i=i*2;
        ALTFarthestSearch.close();
        //ALTFarthestTimes.close();
    }
    //Now Bidirectional
    i=2;
    while(i<=numLandmarks) {
        std::vector<std::vector<double>> potentials = loadPotentials(
                "Potentials_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        std::vector<double> landmarks = loadLandmarks(
                "Landmarks_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        //std::ofstream ALTBidiFarthestTimes(
        //        "experiments/ALTBidiFarthestTimes_" + filename.substr(13, 3) + "_" +
        //        std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

        std::ofstream ALTBidiFarthestSearch("experiments/ALTBidiFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                                std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

        for (auto experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            //auto start = std::chrono::high_resolution_clock::now();
            //ALTBidirectional(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            //auto end = std::chrono::high_resolution_clock::now();
            //std::chrono::duration<double> time = end - start;
            //ALTBidiFarthestTimes << time.count() << std::endl;
            ALTBidiFarthestSearch << ALTBidirectionalSearchSpace(myGraph, sourceNode, targetNode, potentials)<<std::endl;
        }
        i=i*2;
        //ALTBidiFarthestTimes.close();
        ALTBidiFarthestSearch.close();
    }

    //Now Bidirectional with new breaking condition (0.001% error)
    i=2;
    while(i<=numLandmarks) {
        std::vector<std::vector<double>> potentials = loadPotentials(
                "Potentials_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        std::vector<double> landmarks = loadLandmarks(
                "Landmarks_Farthest_" + filename.substr(13, 3) + std::to_string(i));

        //std::ofstream ALTBidiStopFarthestTimes(
        //        "experiments/ALTBidiStopFarthestTimes_" + filename.substr(13, 3) + "_" +
        //        std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

        std::ofstream ALTBidiStopFarthestSearch("experiments/ALTBidiStopFarthestSearchSpace_" + filename.substr(13, 3) + "_" +
                                            std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

        for (auto experiment: Points) {
            double sourceNode = experiment.first;
            double targetNode = experiment.second;
            //auto start = std::chrono::high_resolution_clock::now();
            //ALTBidirectionalStopping(myGraph, sourceNode, targetNode, potentials);
            //Stop the clock
            //auto end = std::chrono::high_resolution_clock::now();
            //std::chrono::duration<double> time = end - start;
            //ALTBidiStopFarthestTimes << time.count() << std::endl;
            ALTBidiStopFarthestSearch << ALTBidirectionalStoppingSearchSpace(myGraph, sourceNode, targetNode, potentials)<<std::endl;
        }
        i=i*2;
        //ALTBidiStopFarthestTimes.close();
        ALTBidiStopFarthestSearch.close();
    }


}

void callExperimentALTUseful(Graph& myGraph, std::vector<std::pair<int, int>>& Points, std::string filename){
    int i = 64;
    std::vector<std::vector<double>> potentials = loadPotentials(
            "Potentials_Farthest_" + filename.substr(13, 3) + std::to_string(i));

    std::vector<double> landmarks = loadLandmarks(
            "Landmarks_Farthest_" + filename.substr(13, 3) + std::to_string(i));

    /*std::ofstream ALTTriangleTimes(
            "experiments/ALTTriangleTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));
*/
    std::ofstream ALTTriangleSearch("experiments/ALTTriangleSearchSpace_" + filename.substr(13, 3) + "_" +
                                        std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> trianglePotentials = findUsefulPotentials(myGraph, sourceNode, targetNode, landmarks
                                                            , potentials);
        //ALT(myGraph, sourceNode, targetNode, trianglePotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        //ALTTriangleTimes << time.count() << std::endl;
        ALTTriangleSearch << ALTSearchSpace(myGraph, sourceNode, targetNode, trianglePotentials) << std::endl;
    }
    //ALTTriangleTimes.close();
    ALTTriangleSearch.close();

    //Now Top Average Potentials
    std::ofstream ALTTopMeanTimes(
            "experiments/ALTTopMeanTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTTopMeanSearch("experiments/ALTTopMeanSearchSpace_" + filename.substr(13, 3) + "_" +
                                    std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> topMeanPotentials = getTopKRows(potentials, 4);
        ALT(myGraph, sourceNode, targetNode, topMeanPotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTTopMeanTimes << time.count() << std::endl;
        ALTTopMeanSearch << ALTSearchSpace(myGraph, sourceNode, targetNode, topMeanPotentials) << std::endl;
    }
    ALTTopMeanTimes.close();
    ALTTopMeanSearch.close();

    //Now Closest Landmarks

    std::ofstream ALTClosestTimes(
            "experiments/ALTClosestTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTClosestSearch("experiments/ALTClosestSearchSpace_" + filename.substr(13, 3) + "_" +
                                    std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> closestPotentials = getClosestLandmarkPotentials(landmarks, potentials, targetNode,
                                                                                          myGraph, 4);
        ALT(myGraph, sourceNode, targetNode, closestPotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTClosestTimes << time.count() << std::endl;
        ALTClosestSearch << ALTSearchSpace(myGraph, sourceNode, targetNode, closestPotentials) << std::endl;
    }
    ALTClosestTimes.close();
    ALTClosestSearch.close();


    //Now Bidirectional

    //Triangle
    std::ofstream ALTBidiTriangleTimes(
            "experiments/ALTBidiTriangleTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTBidiTriangleSearch("experiments/ALTBidiTriangleSearchSpace_" + filename.substr(13, 3) + "_" +
                                    std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> trianglePotentials = findUsefulPotentials(myGraph, sourceNode, targetNode, landmarks
                , potentials);
        ALTBidirectional(myGraph, sourceNode, targetNode, trianglePotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTBidiTriangleTimes << time.count() << std::endl;
        ALTBidiTriangleSearch << ALTBidirectionalSearchSpace(myGraph, sourceNode, targetNode, trianglePotentials) << std::endl;
    }
    ALTBidiTriangleTimes.close();
    ALTBidiTriangleSearch.close();

    //Now Top Average Potentials
    std::ofstream ALTBidiTopMeanTimes(
            "experiments/ALTBidiTopMeanTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTBidiTopMeanSearch("experiments/ALTBidiTopMeanSearchSpace_" + filename.substr(13, 3) + "_" +
                                   std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> topMeanPotentials = getTopKRows(potentials, 4);
        ALTBidirectional(myGraph, sourceNode, targetNode, topMeanPotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTBidiTriangleTimes << time.count() << std::endl;
        ALTBidiTriangleSearch << ALTBidirectionalSearchSpace(myGraph, sourceNode, targetNode, topMeanPotentials) << std::endl;
    }
    ALTBidiTopMeanTimes.close();
    ALTBidiTopMeanSearch.close();

    //Now Closest Landmarks

    std::ofstream ALTBidiClosestTimes(
            "experiments/ALTBidiClosestTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTBidiClosestSearch("experiments/ALTBidiClosestSearchSpace_" + filename.substr(13, 3) + "_" +
                                   std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> closestPotentials = getClosestLandmarkPotentials(landmarks, potentials, targetNode,
                                                                                          myGraph, 4);
        ALTBidirectional(myGraph, sourceNode, targetNode, closestPotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTBidiClosestTimes << time.count() << std::endl;
        ALTBidiClosestSearch << ALTBidirectionalSearchSpace(myGraph, sourceNode, targetNode, closestPotentials) << std::endl;
    }
    ALTBidiClosestTimes.close();
    ALTBidiClosestSearch.close();

    //Now Stopping
    //Triangle
    std::ofstream ALTBidiStopTriangleTimes(
            "experiments/ALTBidiStopTriangleTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTBidiStopTriangleSearch("experiments/ALTBidiStopTriangleSearchSpace_" + filename.substr(13, 3) + "_" +
                                        std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> trianglePotentials = findUsefulPotentials(myGraph, sourceNode, targetNode, landmarks
                , potentials);
        ALTBidirectionalStopping(myGraph, sourceNode, targetNode, trianglePotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTBidiStopTriangleTimes << time.count() << std::endl;
        ALTBidiStopTriangleSearch << ALTBidirectionalStoppingSearchSpace(myGraph, sourceNode, targetNode, trianglePotentials) << std::endl;
    }
    ALTBidiStopTriangleTimes.close();
    ALTBidiStopTriangleSearch.close();

    //Now Top Average Potentials
    std::ofstream ALTBidiStopTopMeanTimes(
            "experiments/ALTBidiStopTopMeanTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTBidiStopTopMeanSearch("experiments/ALTBidiStopTopMeanSearchSpace_" + filename.substr(13, 3) + "_" +
                                       std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> topMeanPotentials = getTopKRows(potentials, 4);
        ALTBidirectionalStopping(myGraph, sourceNode, targetNode, topMeanPotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTBidiStopTriangleTimes << time.count() << std::endl;
        ALTBidiStopTriangleSearch << ALTBidirectionalStoppingSearchSpace(myGraph, sourceNode, targetNode, topMeanPotentials) << std::endl;
    }
    ALTBidiStopTopMeanTimes.close();
    ALTBidiStopTopMeanSearch.close();


    //Closest
    std::ofstream ALTBidiStopClosestTimes(
            "experiments/ALTBidiStopClosestTimes_" + filename.substr(13, 3) + "_" +
            std::to_string(Points.size()) + "_numLandmarks_" + std::to_string(i));

    std::ofstream ALTBidiStopClosestSearch("experiments/ALTBidiStopClosestSearchSpace_" + filename.substr(13, 3) + "_" +
                                       std::to_string(Points.size())+ "_numLandmarks_" + std::to_string(i));

    for (auto experiment: Points) {
        double sourceNode = experiment.first;
        double targetNode = experiment.second;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> closestPotentials = getClosestLandmarkPotentials(landmarks, potentials, targetNode,
                                                                                          myGraph, 4);
        ALTBidirectionalStopping(myGraph, sourceNode, targetNode, closestPotentials);
        //Stop the clock
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        ALTBidiStopClosestTimes << time.count() << std::endl;
        ALTBidiStopClosestSearch << ALTBidirectionalStoppingSearchSpace(myGraph, sourceNode, targetNode, closestPotentials) << std::endl;
    }
    ALTBidiStopClosestTimes.close();
    ALTBidiStopClosestSearch.close();

}
#endif //PRAKTIKUM_EXPERIMENTS_HH
