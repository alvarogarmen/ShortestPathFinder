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
#include <iostream>


void processInput(double sourceNode, double targetNode, std::string graph) {        //This will call Dijkstra and read the files
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
    double Dis = Dijkstra(myGraph, sourceNode, targetNode);
    //Stop the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    // Give out the distance and the time
    std::cout<<"Distance from Node "<<sourceNode<<" to Node "<< targetNode<<" is: "<<Dis<<std::endl;
    std::cout<<"Dijkstra took: "<<time.count()<<"s"<<std::endl;
    //A*
    start = std::chrono::high_resolution_clock::now();
    Dis = AStar(myGraph, sourceNode, targetNode);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"A* from Node "<<sourceNode<<" to Node "<< targetNode<<" is: "<<Dis<<std::endl;
    std::cout<<"A* took: "<<time.count()<<"s"<<std::endl;

    //Saving A* and Dijkstra
    DijkstraSaving(myGraph, sourceNode, targetNode);
    AStarSaving(myGraph, sourceNode, targetNode);

    //ALT with MaxDegree (avoiding) Landmarks
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> Landmarks = selectLandmarks(myGraph, 32);
    std::unordered_map<int, double> landmarkDistances = computeLandmarkDistances(myGraph, Landmarks, 32);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the time
    std::cout<<"ALT Avoid preprocessing took:"<<time.count()<<"s"<<std::endl;

    //ALT query with MaxDegree (avoiding) Landmarks
    start = std::chrono::high_resolution_clock::now();
    Dis = ALT(myGraph, sourceNode, targetNode, landmarkDistances);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"ALT Avoid from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<Dis<<std::endl;
    std::cout<<"ALT Avoid query took: "<<time.count()<<std::endl;
    //Save the exploration and path
    ALTSaving(myGraph, sourceNode, targetNode, landmarkDistances, "Max_explored_nodes.txt", "Max_path.txt");


    //ALT Furthest Landmarks
    start = std::chrono::high_resolution_clock::now();
    Landmarks = computeFurthestLandmarks(myGraph, 32);
    landmarkDistances = computeLandmarkDistances(myGraph, Landmarks, 32);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the time
    std::cout<<"ALT Furthest preprocessing took:"<<time.count()<<"s"<<std::endl;

    //ALT Furthest Query
    start = std::chrono::high_resolution_clock::now();
    Dis = ALT(myGraph, sourceNode, targetNode, landmarkDistances);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"ALT Furthest from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<Dis<<std::endl;
    std::cout<<"ALT Furthest query took: "<<time.count()<<std::endl;
    //Save the exploration and path
    ALTSaving(myGraph, sourceNode, targetNode, landmarkDistances, "Furthest_explored_nodes.txt", "Furthest_path.txt");


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
    std::string graph;

    struct arg_dbl* sourceArg = arg_dbl0(NULL, "source", "<double>", "source node");
    struct arg_dbl* targetArg = arg_dbl0(NULL, "target", "<double>", "target node");
    struct arg_str* graphArg = arg_str0(NULL, "graph", "<string>", "file with the graph, keep the .graph!");

    struct arg_end* end = arg_end(20); // Define the end marker for the argtable array

    void* argtable[] = { sourceArg, targetArg, graphArg, end };

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

    if (graphArg->count > 0) {
        graph = graphArg->sval[0];
    }
    // Call out functions
    processInput(sourceNode, targetNode, graph);

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}