#include <bitset>

double ALTBidirectional(Graph& myGraph, double& sourceNode, double& targetNode,
                        std::vector<std::vector<double>>& potentials) {
    //Get two priority queues, visited vectors and dist arrays
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);

    //bool meeting = false;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityForDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityBackDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        if (apqForward.getMin().second < apqBackward.getMin().second){
            double forwardNode = apqForward.popMin();

            double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
            double endForwardEdge = myGraph.edgeStarts[forwardNode];

            for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
                double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
                double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);


                if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                    distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                    priorityForDist[forwardEdge] = priorityForDist[forwardNode] -
                                                   estimate(forwardNode, targetNode - 1, potentials) +
                                                   forwardWeight +
                                                   estimate(forwardEdge, targetNode - 1, potentials);

                    double forwardF = priorityForDist[forwardEdge];

                    if (apqForward.contains(forwardEdge)) {
                        apqForward.decreaseKey(forwardEdge, forwardF);
                    } else {
                        apqForward.insertNode(forwardEdge, forwardF);
                    }
                    if(distBackward[forwardEdge]<INT_MAX || forwardEdge==targetNode-1){
                        return distForward[forwardEdge]+distBackward[forwardEdge];
                    }
                }
            }
        }

        else {
            double backwardNode = apqBackward.popMin();

            double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
            double endBackwardEdge = myGraph.edgeStarts[backwardNode];

            for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
                double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
                double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

                if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                    distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;
                    priorityBackDist[backwardEdge] = priorityForDist[backwardNode] -
                                                     estimate(backwardNode, sourceNode - 1, potentials) +
                                                     backwardWeight +
                                                     estimate(backwardEdge, sourceNode - 1, potentials);

                    double backwardF = priorityForDist[backwardEdge];


                    if (apqBackward.contains(backwardEdge)) {
                        apqBackward.decreaseKey(backwardEdge, backwardF);
                    } else {
                        apqBackward.insertNode(backwardEdge, backwardF);
                    }
                    if (distForward[backwardEdge]<INT_MAX || backwardEdge==sourceNode-1){
                        return distForward[backwardEdge]+distBackward[backwardEdge];
                    }
                }
            }

        }

    }

    return -1;
}

double ALTBidirectionalSaving(Graph& myGraph, double& sourceNode, double& targetNode, std::string filename,
                              std::vector<std::vector<double>>& potentials) {
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<double> forwardVisited;
    std::vector<double> backwardVisited;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityForDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityBackDist(myGraph.nodes.size(), INT_MAX);
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    std::ofstream exploredForwardFile("experiments/"+filename+"_Forward_explored");
    std::ofstream exploredBackwardFile("experiments/"+filename+"_Backward_explored");
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        if (apqForward.getMin().second <= apqBackward.getMin().second){
            double forwardNode = apqForward.popMin();
            forwardVisited.push_back(forwardNode);
            double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
            double endForwardEdge = myGraph.edgeStarts[forwardNode];
            for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
                double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
                double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);
                if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                    distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                    if(distBackward[forwardEdge]<INT_MAX){
                        break;
                    }
                    priorityForDist[forwardEdge] = priorityForDist[forwardNode] -
                                                   estimate(forwardNode, targetNode - 1, potentials) +
                                                   forwardWeight +
                                                   estimate(forwardEdge, targetNode - 1, potentials);

                    double forwardF = priorityForDist[forwardEdge];


                    if (apqForward.contains(forwardEdge)) {
                        apqForward.decreaseKey(forwardEdge, forwardF);
                    } else {
                        apqForward.insertNode(forwardEdge, forwardF);
                    }
                }
            }
        }
        else{
            double backwardNode = apqBackward.popMin();

            backwardVisited.push_back(backwardNode);


            double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
            double endBackwardEdge = myGraph.edgeStarts[backwardNode];


            for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
                double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
                double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

                if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                    distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;
                    if(distForward[backwardEdge]){
                        break;
                    }

                    priorityForDist[backwardEdge] = priorityForDist[backwardNode] -
                                                    estimate(backwardNode, sourceNode - 1, potentials) +
                                                    backwardWeight +
                                                    estimate(backwardEdge, sourceNode - 1, potentials);

                    double backwardF = priorityForDist[backwardEdge];

                    if (apqBackward.contains(backwardEdge)) {
                        apqBackward.decreaseKey(backwardEdge, backwardF);
                    } else {
                        apqBackward.insertNode(backwardEdge, backwardF);
                    }

                }
            }
        }


    }

    for(auto forwardNode : forwardVisited){
        exploredForwardFile << myGraph.getNode(forwardNode).coordinateX << " " << myGraph.getNode(forwardNode).coordinateY << std::endl;  // Write explored node to the file
    }
    for(auto backwardNode : backwardVisited){
        exploredBackwardFile << myGraph.getNode(backwardNode).coordinateX << " " << myGraph.getNode(backwardNode).coordinateY << std::endl;  // Write explored node to the file
    }

    exploredForwardFile.close();
    exploredBackwardFile.close();

    return 0;
}

double ALTBidirectionalSearchSpace(Graph& myGraph, double& sourceNode, double& targetNode,
                                   std::vector<std::vector<double>>& potentials) {

    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<double> priorityForDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityBackDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> visited;
    double minPath = INT_MAX;
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    visited.push_back(sourceNode);
    visited.push_back(targetNode);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            visited.push_back(forwardEdge);
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                priorityForDist[forwardEdge] = priorityForDist[forwardNode] -
                                               estimate(forwardNode, targetNode - 1, potentials) +
                                               forwardWeight +
                                               estimate(forwardEdge, targetNode - 1, potentials);

                double forwardF = priorityForDist[forwardEdge];
                if(forwardF+distBackward[forwardEdge] < minPath){
                    minPath=distBackward[forwardEdge]+distForward[forwardEdge];
                }


                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            visited.push_back(backwardEdge);
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                priorityForDist[backwardEdge] = priorityForDist[forwardNode] -
                                                estimate(forwardNode, sourceNode - 1, potentials) +
                                                backwardWeight +
                                                estimate(backwardEdge, sourceNode - 1, potentials);

                double backwardF = priorityForDist[backwardEdge];
                if(backwardF+distForward[backwardEdge] < minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }


                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
        if (apqBackward.getMin().second >= minPath || apqForward.getMin().second >= minPath){
            return visited.size();
        }
    }

    return visited.size();
}

double ALTBidirectionalUseful(Graph& myGraph, double& sourceNode, double& targetNode,
                        std::vector<std::vector<double>>& potentials, std::vector<double>& landmarks) {
    std::vector<double> usefulLandmarksForward = findUsefulLandmarks(myGraph, sourceNode, targetNode, landmarks);
    std::vector<double> usefulLandmarksBackward = findUsefulLandmarks(myGraph, targetNode, sourceNode, landmarks);
    //Get two priority queues, visited vectors and dist arrays
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<bool> visitedForward(myGraph.nodes.size(), false);

    std::vector<bool> visitedBackward(myGraph.nodes.size(), false);
    //bool meeting = false;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityForDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityBackDist(myGraph.nodes.size(), INT_MAX);

    // Init the minPath variable to infinity. If the current bestPath gets higher than this variable, we terminate
    double minPath = INT_MAX;
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;
    visitedForward[sourceNode-1]= true;
    visitedBackward[targetNode-1]= true;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();



        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);
            visitedForward[forwardEdge]= true;


            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                priorityForDist[forwardEdge] = priorityForDist[forwardNode] -
                                               estimateUseful(forwardNode, targetNode - 1, potentials,usefulLandmarksForward) +
                                               forwardWeight +
                                               estimateUseful(forwardEdge, targetNode - 1, potentials, usefulLandmarksForward);

                double forwardF = priorityForDist[forwardEdge];
                if(distForward[forwardEdge]+distBackward[forwardEdge] < minPath){
                    minPath=distForward[forwardEdge]+distBackward[forwardEdge];
                }
                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
                if(visitedBackward[forwardEdge]){
                    return minPath;
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);
            visitedBackward[backwardNode]=true;

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;
                priorityBackDist[backwardEdge] = priorityForDist[backwardNode] -
                                                 estimateUseful(backwardNode, sourceNode - 1, potentials,usefulLandmarksBackward) +
                                                 backwardWeight +
                                                 estimateUseful(backwardEdge, sourceNode - 1, potentials, usefulLandmarksBackward);

                double backwardF = priorityForDist[backwardEdge];
                if(distForward[backwardEdge]+distBackward[backwardEdge] < minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
                if (visitedForward[backwardEdge]){
                    return minPath;
                }
            }
        }

    }

    return minPath;
}

double ALTBidirectionalUsefulSaving(Graph& myGraph, double& sourceNode, double& targetNode, std::string filename,
                              std::vector<std::vector<double>>& potentials, std::vector<double>& landmarks) {
    std::vector<double> usefulLandmarksForward = findUsefulLandmarks(myGraph, sourceNode, targetNode, landmarks);
    std::vector<double> usefulLandmarksBackward = findUsefulLandmarks(myGraph, targetNode, sourceNode, landmarks);
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<double> forwardVisited;
    std::vector<double> backwardVisited;
    std::vector<bool> visitedForward(myGraph.nodes.size(), false);
    std::vector<bool> visitedBackward(myGraph.nodes.size(), false);
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityForDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityBackDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> forwardPath(myGraph.nodes.size(), -1);
    std::vector<double> backwardPath(myGraph.nodes.size(), -1);
    double minPath = INT_MAX;
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;
    visitedForward[sourceNode-1] = true;
    visitedBackward[targetNode-1] = true;

    std::ofstream exploredForwardFile("experiments/"+filename+"_Forward_explored");
    std::ofstream exploredBackwardFile("experiments/"+filename+"_Backward_explored");
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        forwardVisited.push_back(forwardNode);
        backwardVisited.push_back(backwardNode);

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            visitedForward[forwardEdge] = true;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                forwardPath[forwardEdge] = forwardNode;

                priorityForDist[forwardEdge] = priorityForDist[forwardNode] -
                                               estimateUseful(forwardNode, targetNode - 1, potentials,usefulLandmarksForward) +
                                               forwardWeight +
                                               estimateUseful(forwardEdge, targetNode - 1, potentials, usefulLandmarksForward);

                double forwardF = priorityForDist[forwardEdge];
                if(forwardF+distBackward[forwardEdge] < minPath){
                    minPath=distBackward[forwardEdge]+distForward[forwardEdge];
                }


                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
                if (visitedBackward[forwardEdge]){
                    break;
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            visitedBackward[backwardEdge]=true;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                priorityBackDist[backwardEdge] = priorityForDist[backwardNode] -
                                                 estimateUseful(backwardNode, sourceNode - 1, potentials,usefulLandmarksBackward) +
                                                 backwardWeight +
                                                 estimateUseful(backwardEdge, sourceNode - 1, potentials, usefulLandmarksBackward);

                double backwardF = priorityForDist[backwardEdge];
                if(backwardF+distForward[backwardEdge] < minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }


                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
        if (apqBackward.getMin().second >= minPath || apqForward.getMin().second >= minPath){
            break;
        }
    }

    for(auto forwardNode : forwardVisited){
        exploredForwardFile << myGraph.getNode(forwardNode).coordinateX << " " << myGraph.getNode(forwardNode).coordinateY << std::endl;  // Write explored node to the file
    }
    for(auto backwardNode : backwardVisited){
        exploredBackwardFile << myGraph.getNode(backwardNode).coordinateX << " " << myGraph.getNode(backwardNode).coordinateY << std::endl;  // Write explored node to the file
    }

    exploredForwardFile.close();
    exploredBackwardFile.close();

    return minPath;
}

double ALTBidirectionalUsefulSearchSpace(Graph& myGraph, double& sourceNode, double& targetNode,
                                   std::vector<std::vector<double>>& potentials, std::vector<double>& landmarks) {
    std::vector<double> usefulLandmarksForward = findUsefulLandmarks(myGraph, sourceNode, targetNode, landmarks);
    std::vector<double> usefulLandmarksBackward = findUsefulLandmarks(myGraph, targetNode, sourceNode, landmarks);

    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);

    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> visited;
    double minPath = INT_MAX;
    std::vector<double> priorityForDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> priorityBackDist(myGraph.nodes.size(), INT_MAX);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    visited.push_back(sourceNode);
    visited.push_back(targetNode);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {

        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            visited.push_back(forwardEdge);
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                priorityForDist[forwardEdge] = priorityForDist[forwardNode] -
                                               estimateUseful(forwardNode, targetNode - 1, potentials,usefulLandmarksForward) +
                                               forwardWeight +
                                               estimateUseful(forwardEdge, targetNode - 1, potentials, usefulLandmarksForward);

                double forwardF = priorityForDist[forwardEdge];
                if(forwardF+distBackward[forwardEdge] < minPath){
                    minPath=distBackward[forwardEdge]+distForward[forwardEdge];
                }


                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            visited.push_back(backwardEdge);
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                priorityBackDist[backwardEdge] = priorityForDist[backwardNode] -
                                                 estimateUseful(backwardNode, sourceNode - 1, potentials,usefulLandmarksBackward) +
                                                 backwardWeight +
                                                 estimateUseful(backwardEdge, sourceNode - 1, potentials, usefulLandmarksBackward);

                double backwardF = priorityForDist[backwardEdge];
                if(backwardF+distForward[backwardEdge] < minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }


                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
        if (apqBackward.getMin().second >= minPath || apqForward.getMin().second >= minPath){
            return visited.size();
        }
    }

    return visited.size();
}