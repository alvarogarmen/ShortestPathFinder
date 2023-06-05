//
// Created by alvar on 05/06/2023.
//

#ifndef PRAKTIKUM_DUMPSTER_H
#define PRAKTIKUM_DUMPSTER_H

/*
/ALT with random Landmarks Preprocessing
    start = std::chrono::high_resolution_clock::now();
    std::unordered_map<int, double> landmarkDistancesRandom = computeLandmarkDistancesRandom(myGraph, 8);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the time
    std::cout<<"ALT Random preprocessing took:"<<time.count()<<"s"<<std::endl;

    //ALT query with random Landmarks
    start = std::chrono::high_resolution_clock::now();
    Dis = ALT(myGraph, sourceNode, targetNode, landmarkDistancesRandom);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the distance and the time
    std::cout<<"ALT Random from Node "<<sourceNode<<" to Node "<<targetNode<<" is "<<Dis<<std::endl;
    std::cout<<"ALT Random query took: "<<time.count()<<std::endl;

    //ALT with furthest Landmarks Preprocessing
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> landmarkDistancesFurthest = computeFurthestLandmarks(myGraph, 8);
    //Stop the clock
    end = std::chrono::high_resolution_clock::now();
    time = end - start;
    // Give out the time
    std::cout<<"ALT Furthest preprocessing took:"<<time.count()<<"s"<<std::endl;
*/
#endif //PRAKTIKUM_DUMPSTER_H
