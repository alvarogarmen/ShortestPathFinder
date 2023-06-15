double ALTBidirectional(Graph myGraph, double sourceNode, double targetNode, const std::vector<std::vector<double>>& potentials) {
    APQ apqForward = APQ();  // Priority queue for forward search
    APQ apqBackward = APQ(); // Priority queue for backward search
    std::set<double> visitedForward;  // Set of visited nodes for forward search
    std::set<double> visitedBackward; // Set of visited nodes for backward search
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);  // Distances from source for forward search
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX); // Distances from target for backward search

    // Initialize the search
    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    double meetingNode = -1; // Node where forward and backward searches meet

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        // Perform one step of forward search
        double currentForwardNode = apqForward.getMin().first;
        visitedForward.insert(currentForwardNode);
        apqForward.popMin();

        double startForwardEdge = (currentForwardNode > 0) ? myGraph.edgeStarts[currentForwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[currentForwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[currentForwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[currentForwardNode] + forwardWeight < distForward[forwardEdge]) {
                distForward[forwardEdge] = distForward[currentForwardNode] + forwardWeight;

                if (visitedBackward.find(forwardEdge) != visitedBackward.end()) {
                    // Node visited in both searches, we have a meeting node
                    meetingNode = forwardEdge;
                    break;
                }

                double forwardF = distForward[forwardEdge] + potentials[forwardEdge][targetNode - 1];

                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        if (meetingNode != -1)
            break;

        // Perform one step of backward search
        double currentBackwardNode = apqBackward.getMin().first;
        visitedBackward.insert(currentBackwardNode);
        apqBackward.popMin();

        double startBackwardEdge = (currentBackwardNode > 0) ? myGraph.edgeStarts[currentBackwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[currentBackwardNode];

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[currentBackwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[currentBackwardNode] + backwardWeight < distBackward[backwardEdge]) {
                distBackward[backwardEdge] = distBackward[currentBackwardNode] + backwardWeight;

                if (visitedForward.find(backwardEdge) != visitedForward.end()) {
                    // Node visited in both searches, we have a meeting node
                    meetingNode = backwardEdge;
                    break;
                }

                double backwardF = distBackward[backwardEdge] + potentials[backwardEdge][sourceNode - 1];

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }

        if (meetingNode != -1)
            break;
    }

    if (meetingNode == -1) {
        // No path found
        std::cout << "No path found" << std::endl;
        return -1;
    }

    // Calculate the shortest path distance
    double shortestPathDistance = INT_MAX;
    for (double node : visitedForward) {
        if (visitedBackward.find(node) != visitedBackward.end()) {
            double pathDistance = distForward[node] + distBackward[node];
            if (pathDistance < shortestPathDistance)
                shortestPathDistance = pathDistance;
        }
    }

    return shortestPathDistance;
}
