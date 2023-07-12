//
// Created by alvar on 12/07/2023.
//

#ifndef PRAKTIKUM_ALTERNATIVE_HH
#define PRAKTIKUM_ALTERNATIVE_HH

#include "Dijkstra.hh"
double calculateAverageDistanceReduction(const std::vector<std::vector<double>>& distances) {
    int totalPairs = 0;
    double totalReduction = 0.0;

    int numVertices = distances.size();
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            if (i != j) {
                totalReduction += (distances[i][j] - distances[i][j]) / static_cast<double>(distances[i][j]);
                totalPairs++;
            }
        }
    }

    return totalReduction / totalPairs;
}

// Function to find the best landmark
int findBestLandmark(const std::vector<Node>& vertices, const std::vector<double>& edges) {
    int numVertices = vertices.size();
    std::vector<std::vector<double>> distances(numVertices, std::vector<double>(numVertices, std::numeric_limits<int>::max()));

    // Calculate initial distances using your graph traversal algorithm (e.g., Dijkstra's algorithm)
    // Update the distances matrix accordingly

    double bestReduction = 0.0;
    int bestLandmark = -1;

    for (int i = 0; i < numVertices; ++i) {
        double reduction = calculateAverageDistanceReduction(distances);

        if (reduction > bestReduction) {
            bestReduction = reduction;
            bestLandmark = i;
        }
    }

    return bestLandmark;
}
#endif //PRAKTIKUM_ALTERNATIVE_HH
