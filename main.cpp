//
// Created by alvar on 26/04/2023.
//

#include <string>
#include "argtable3.h"
#include <chrono>
#include "readgraph.hh"
#include "Dijkstra.hh"
#include "AStar.hh"
#include "ALT-Landmarks.hh"
#include "ALT.hh"
#include "Dijkstra_Saving.hh"
#include "AStar_Saving.hh"
#include "ALT_Saving.h"
#include "ProcessInput.hh"
#include <iostream>


void processInput(double sourceNode, double targetNode, std::string graph, int numLandmarks, int newLandmarks) {        //This will call Dijkstra and read the files
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

    //Dijkstra
    callDijkstra(myGraph, sourceNode, targetNode);

    //A*
    callAStar(myGraph, sourceNode, targetNode);
    //callAStarBidirectional(myGraph, sourceNode, targetNode);

    //Saving A* and Dijkstra
    DijkstraSaving(myGraph, sourceNode, targetNode);
    AStarSaving(myGraph, sourceNode, targetNode);
    std::cout<<"HERE"<<std::endl;
    //callALTFurthest(myGraph, sourceNode, targetNode);
    //Bidirectional A*
    //callAStarBidirectional(myGraph, sourceNode, targetNode);
    //ALT Avoid
    callALTAvoid(myGraph, sourceNode, targetNode, numLandmarks);
    std::cout<<"STRANGE"<<std::endl;

    //ALT Bidirectional Avoid
    //callALTBiAvoid(myGraph, sourceNode, targetNode, numLandmarks);
    //ALT with MaxDegree (avoiding) Landmarks
   //callALTMaxDegree(myGraph, sourceNode, targetNode, numLandmarks);

    //callALTFurthest(myGraph, sourceNode, targetNode, numLandmarks);

    //callHHS(myGraph, sourceNode, targetNode, numLandmarks);
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
    std::string graph;

    struct arg_dbl* sourceArg = arg_dbl0(NULL, "source", "<double>", "source node");
    struct arg_dbl* targetArg = arg_dbl0(NULL, "target", "<double>", "target node");
    struct arg_int* landmarkArg = arg_int0(NULL, "landmarks", "<integer>", "number of landmarks");
    struct arg_int* newLandmarkArg = arg_int0(NULL, "newLandmarks?", "<integer>", "0 for old landmarks");
    struct arg_str* graphArg = arg_str0(NULL, "graph", "<string>", "file with the graph, keep the .graph!");

    struct arg_end* end = arg_end(20); // Define the end marker for the argtable array

    void* argtable[] = { sourceArg, targetArg, landmarkArg, newLandmarkArg, graphArg, end };

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
    // Call out functions
    processInput(sourceNode, targetNode, graph, numLandmarks, newLandmarks);

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}