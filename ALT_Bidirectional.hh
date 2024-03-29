double ALTBidirectional(const Graph& myGraph, double sourceNode, double targetNode, const std::vector<std::vector<double>>& potentials) {
    APQ forwardAPQ = APQ(myGraph.nodeCount);
    APQ backwardAPQ = APQ(myGraph.nodeCount);
    std::vector<double> forwardDist(myGraph.nodes.size(), INT_MAX);
    std::vector<double> backwardDist(myGraph.nodes.size(), INT_MAX);

    forwardAPQ.insertNode(sourceNode - 1, 0);
    backwardAPQ.insertNode(targetNode - 1, 0);
    forwardDist[sourceNode - 1] = 0;
    backwardDist[targetNode - 1] = 0;

    double meetingNode = -1;
    double bestPath = INT_MAX;

    while (!forwardAPQ.isEmpty() && !backwardAPQ.isEmpty()) {
        // Perform a forward search step
        double forwardNode = forwardAPQ.popMin();

        double forwardStartEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double forwardEndEdge = myGraph.edgeStarts[forwardNode];

        for (double edgeIndex = forwardStartEdge; edgeIndex <= forwardEndEdge; edgeIndex++) {
            double forwardEdge = myGraph.edges[edgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (forwardDist[forwardNode] + forwardWeight < forwardDist[forwardEdge]) {
                forwardDist[forwardEdge] = forwardDist[forwardNode] + forwardWeight;
                double forwardF = forwardDist[forwardEdge] + estimate(forwardEdge, targetNode-1, potentials);
                if (forwardAPQ.contains(forwardEdge)) {
                    forwardAPQ.decreaseKey(forwardEdge, forwardF);
                } else if(forwardF < bestPath) {
                    forwardAPQ.insertNode(forwardEdge, forwardF);
                }

                // Check if the forward search reaches the backward search frontier
                if(backwardDist[forwardEdge]!=INT_MAX && forwardDist[forwardEdge]+backwardDist[forwardEdge]<bestPath){
                    bestPath=forwardDist[forwardEdge]+backwardDist[forwardEdge];
                }

            }
        }

        // Perform a backward search step
        double backwardNode = backwardAPQ.popMin();

        double backwardStartEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double backwardEndEdge = myGraph.edgeStarts[backwardNode];

        for (double edgeIndex = backwardStartEdge; edgeIndex <= backwardEndEdge; edgeIndex++) {
            double backwardEdge = myGraph.edges[edgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);
            if (backwardDist[backwardNode] + backwardWeight < backwardDist[backwardEdge]) {
                backwardDist[backwardEdge] = backwardDist[backwardNode] + backwardWeight;
                double backwardF = backwardDist[backwardEdge] + estimate(backwardEdge, sourceNode-1, potentials);
                if (backwardAPQ.contains(backwardEdge)) {
                    backwardAPQ.decreaseKey(backwardEdge, backwardF);
                } else if(backwardF < bestPath){
                    backwardAPQ.insertNode(backwardEdge, backwardF);
                }

                // Check if the backward search reaches the forward search frontier
                if(forwardDist[backwardEdge]!=INT_MAX && forwardDist[backwardEdge]+backwardDist[backwardEdge]<bestPath){
                    bestPath=forwardDist[backwardEdge]+backwardDist[backwardEdge];
                }

            }
        }
        if (backwardAPQ.isEmpty()||forwardAPQ.isEmpty()){
            return bestPath;
        }
        if (backwardAPQ.getMin().second >= bestPath || forwardAPQ.getMin().second >= bestPath){
            return bestPath;
        }
    }

    return bestPath;
}

/*double ALTBidirectional(Graph& myGraph, double& sourceNode, double& targetNode,
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
    double minPath = INT_MAX;
    bool meeting = false;

    //TODO: add poped array so that nodes dont get reinserted
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.popMin();
        double backwardNode = apqBackward.popMin();

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);


            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                double forwardF = distForward[forwardEdge]+
                                  estimate(forwardEdge, targetNode - 1, potentials);

                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
                if(distBackward[forwardEdge]+distForward[forwardEdge]<minPath){
                    minPath=distForward[forwardEdge]+distBackward[forwardEdge];
                }
                if(forwardEdge==backwardNode){meeting = true;}
            }
        }


        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                double backwardF = distBackward[backwardEdge]+
                                   estimate(backwardEdge, sourceNode - 1, potentials);


                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
                if(distBackward[backwardEdge]+distForward[backwardEdge]<minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }
                if(backwardEdge==forwardNode){
                    meeting=true;
                }
            }
        }
        if (forwardNode==targetNode-1||backwardNode==sourceNode-1){
            return minPath;
        }
        if (meeting){
            return minPath;
        }
    }

    return -1;
}
*/
double ALTBidirectionalSaving(Graph& myGraph, double& sourceNode, double& targetNode, std::string filename,
                              std::vector<std::vector<double>>& potentials) {
    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
    std::vector<double> forwardVisited;
    std::vector<double> backwardVisited;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;
    double bestPath = INT_MAX;

    std::ofstream exploredForwardFile("experiments/"+filename+"_Forward_explored");
    std::ofstream exploredBackwardFile("experiments/"+filename+"_Backward_explored");
    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.popMin();
        forwardVisited.push_back(forwardNode);
        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];
        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);
            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                double forwardF = distForward[forwardEdge]+estimate(forwardEdge, targetNode - 1, potentials);

                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
                if(distBackward[forwardEdge]!=INT_MAX && distForward[forwardEdge]+distBackward[forwardEdge]<bestPath){
                    bestPath=distForward[forwardEdge]+distBackward[forwardEdge];
                }
            }
        }


        double backwardNode = apqBackward.popMin();

        backwardVisited.push_back(backwardNode);


        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];


        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight  < distBackward[backwardEdge] ) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                double backwardF = distBackward[backwardEdge]+estimate(backwardEdge, sourceNode - 1, potentials);

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
                if(distForward[backwardEdge]!=INT_MAX && distForward[backwardEdge]+distBackward[backwardEdge]<bestPath){
                    bestPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }

            }

        }
        if(apqBackward.isEmpty()||apqForward.isEmpty()){
            break;
        }
        if (apqBackward.getMin().second >= bestPath || apqForward.getMin().second >= bestPath){
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

    return 0;
}

double ALTBidirectionalSearchSpace(Graph& myGraph, double& sourceNode, double& targetNode,
                                   std::vector<std::vector<double>>& potentials) {

    APQ apqForward = APQ(myGraph.nodeCount);
    APQ apqBackward = APQ(myGraph.nodeCount);
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
        visited.push_back(forwardNode);

        double backwardNode = apqBackward.popMin();
        visited.push_back(backwardNode);


        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;
                double forwardF = distForward[forwardEdge]+ estimate(forwardEdge, targetNode-1, potentials);
                if(distBackward[forwardEdge]!=INT_MAX && distForward[forwardEdge]+distBackward[forwardEdge]<minPath){
                    minPath=distForward[forwardEdge]+distBackward[forwardEdge];
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

                double backwardF = distBackward[backwardEdge]+ estimate(backwardEdge, sourceNode-1, potentials);
                if(distForward[backwardEdge]!=INT_MAX && distForward[backwardEdge]+distBackward[backwardEdge]<minPath){
                    minPath=distForward[backwardEdge]+distBackward[backwardEdge];
                }

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
        if(apqBackward.isEmpty()||apqForward.isEmpty()){
            return visited.size();
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
        if(apqForward.isEmpty()||apqBackward.isEmpty()){
            break;
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
        if(apqBackward.isEmpty()||apqForward.isEmpty()){
            return visited.size();
        }
        if (apqBackward.getMin().second >= minPath || apqForward.getMin().second >= minPath){
            return visited.size();
        }
    }

    return visited.size();
}