//
// Created by alvar on 26/04/2023.
//

#ifndef UNTITLED_READGRAPH_HH
#define UNTITLED_READGRAPH_HH

#include <vector>
#include <string>
#include "Graph.hh"
#include <iostream>
#include <fstream>
#include <sstream>


Graph readCoordFile(std::string Coordfilename){
    int nodeId = 0;
    Graph myGraph;
    std::ifstream infile(Coordfilename+".xyz"); // open the input file
    std::string line;
    while (std::getline(infile, line)) { // read the file line by line
        std::istringstream iss(line);

        double node, coordx, coordy, coordz;

        iss >> coordx >> coordy >> coordz; // extract coordinates from line
        // store x and y for the node and push it to the graph
        Node myNode(coordx, coordy, nodeId);
        myGraph.insertNode(myNode);
        std::cout << "Node " << nodeId <<":" << "\n";
        std::cout << "X Coordinate: " << coordx << "\n";
        std::cout << "Y Coordinate: " << coordy << "\n";
        std::cout << "Z Coordinate: " << coordz << "\n";
        nodeId++;
    }
    for (int i=0; i<myGraph.nodes.size();i++){
        std::cout << myGraph.getNode(i).nodeId<<", ";
    }
    return myGraph;

}

void readOtherFile(std::string filename, Graph& myGraph) {
    std::ifstream infile(filename);
    std::string newline;
    std::vector<double> helper;
    double nodeId;
    double count = 0;
    double index = 0;

    while (std::getline(infile, newline)) { // read each line of the input file
        std::istringstream iss(newline); // create an input string stream for the line
        while (iss >> nodeId) { // read each value on the line
            if (count < 1) {
                double nodeCount, edgeCount;
                iss >> nodeCount >> edgeCount;
                myGraph.nodeCount = nodeCount;
                myGraph.edgeCount = edgeCount;
                count = 1;
            }
            else {
                myGraph.insertEdge(nodeId);
                helper.push_back(nodeId);

            }
        }
        //set the value for the node to be able to find its edges
        if(count==1) {
            if(helper.size()!=0) {
                myGraph.getNode(index-1).setStart(helper.size() - 1);
                myGraph.edgeStarts.push_back(helper.size()-1);
            }

            index++;
        }

    }

}



#endif //UNTITLED_READGRAPH_HH
