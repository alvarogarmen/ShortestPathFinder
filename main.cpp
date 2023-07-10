//
// Created by alvar on 26/04/2023.
//

#include <string>
#include "argtable3.h"
#include <chrono>
#include "readgraph.hh"
#include "Dijkstra.hh"
#include "AStar.hh"

#include "Experiments.hh"
#include <iostream>


void processInput(double sourceNode, double targetNode, std::string graph, int numLandmarks, int newLandmarks, int newPoints, int numPoints) {        //This will call Dijkstra and read the files
    // Give out the input
    std::cout << "Source Node: " << sourceNode << std::endl;
    std::cout << "Target Node: " << targetNode << std::endl;
    std::cout<<graph<<std::endl;
    std::cout << "Graph1: " << graph << std::endl;
    // Call the read-functions
    Graph myGraph = readCoordFile(graph);
    std::cout<<"Reading successful"<<std::endl;
    readOtherFile(graph, myGraph);
    std::cout<<"Reading successful"<<std::endl;
    std::vector<std::pair<int, int>> points;
    callDijkstra(myGraph, sourceNode, targetNode);
    callALTFarthest(myGraph, sourceNode, targetNode, numLandmarks, newLandmarks, graph);
    if (newPoints!=0){
        points = generatePoints(myGraph, numPoints, graph);
        std::cout<<"New Points generated!"<<std::endl;
    }
    else{
        points = loadPoints(graph.substr(13, 3)+"_"+std::to_string(numPoints));
        std::cout<<"Points loaded!"<<std::endl;
    }
    callExperimentDijkstra(myGraph, points, graph);
    std::cout<<"Dijkstra Experiments successful!"<<std::endl;
    callExperimentAStar(myGraph, points, graph);
    std::cout<<"AStar Experiments successful!"<<std::endl;
    callExperimentALT(myGraph, points, graph, numLandmarks, newLandmarks);
    std::cout<<"ALT Experiments successful!"<<std::endl;



}

void processInputForPlot(double sourceNode, double targetNode, std::string graph){
    // Give out the input
    std::cout << "Source Node: " << sourceNode << std::endl;
    std::cout << "Target Node: " << targetNode << std::endl;
    std::cout<<graph<<std::endl;
    std::cout << "Graph1: " << graph << std::endl;
    // Call the read-functions
    Graph myGraph = readCoordFile(graph);
    std::cout<<"Reading successful"<<std::endl;
    readOtherFile(graph, myGraph);
    std::cout<<"Reading successful"<<std::endl;
    // Call Dijkstra and outDegree timer
    auto start = std::chrono::high_resolution_clock::now();
    double Dis = DijkstraSaving(myGraph, sourceNode, targetNode);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the distance and the time
    std::cout<<"Distance from Node "<<sourceNode<<" to Node "<< targetNode<<" is: "<<Dis<<std::endl;
    std::cout<<"Dijkstra took: "<<time.count()<<"s"<<std::endl;

    // Call Dijkstra and outDegree timer
    start = std::chrono::high_resolution_clock::now();
    Dis = AStarSaving(myGraph, sourceNode, targetNode);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"Distance from Node "<<sourceNode<<" to Node "<< targetNode<<" is: "<<Dis<<std::endl;
    std::cout<<"A* took: "<<time.count()<<"s"<<std::endl;
}
int main(int argc, char* argv[]) {
    //Argtable 3
    double sourceNode;
    double targetNode;
    int numLandmarks;
    int newLandmarks;
    int newPoints;
    int numPoints;
    std::string graph;

    struct arg_dbl* sourceArg = arg_dbl0(NULL, "source", "<double>", "source node");
    struct arg_dbl* targetArg = arg_dbl0(NULL, "target", "<double>", "target node");
    struct arg_int* landmarkArg = arg_int0(NULL, "landmarks", "<integer>", "number of landmarks");
    struct arg_int* newLandmarkArg = arg_int0(NULL, "newLandmarks?", "<integer>", "0 for old landmarks");
    struct arg_str* graphArg = arg_str0(NULL, "graph", "<string>", "file with the graph, keep the .graph!");
    struct arg_int* newPointsArg = arg_int0(NULL, "newPoints?", "<integer>", "generate new Points?");
    struct arg_int* numPointsArg = arg_int0(NULL, "numNewPoints?", "<integer>", "how many new Points?");


    struct arg_end* end = arg_end(20); // Define the end marker for the argtable array

    void* argtable[] = { sourceArg, targetArg, landmarkArg, newLandmarkArg, graphArg, newPointsArg, numPointsArg, end };

    const char* progname = "myprogram";
    int nerrors = arg_parse(argc, argv, argtable);

    if (nerrors > 0) {
        arg_print_errors(stdout, end, progname);
        arg_print_syntax(stdout, argtable, "\n");
        return 1; // Handle parsing errors
    }

    if (sourceArg->count > 0) {
        sourceNode = sourceArg->dval[0];
    }

    if (targetArg->count > 0) {
        targetNode = targetArg->dval[0];
    }
    if (landmarkArg->count > 0){
        numLandmarks = landmarkArg->ival[0];
    }
    if(newLandmarkArg->count >0){
        newLandmarks = newLandmarkArg->ival[0];
    }
    if (graphArg->count > 0) {
        graph = graphArg->sval[0];
    }
    if (numPointsArg->count > 0) {
        numPoints = numPointsArg->ival[0];
    }
    if (newPointsArg->count > 0) {
        newPoints = newPointsArg->ival[0];
    }
    // Call out functions
    processInput(sourceNode, targetNode, graph, numLandmarks, newLandmarks, newPoints, numPoints);

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}