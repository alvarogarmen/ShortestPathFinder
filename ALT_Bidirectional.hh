double ALTBidirectional(Graph& myGraph, double& sourceNode, double& targetNode, const std::vector<std::vector<double>>& potentials) {
    APQ apqForward = APQ();
    APQ apqBackward = APQ();
    std::set<double> visitedForward;
    std::set<double> visitedBackward;
    std::vector<double> distForward(myGraph.nodes.size(), INT_MAX);
    std::vector<double> distBackward(myGraph.nodes.size(), INT_MAX);

    apqForward.insertNode(sourceNode - 1, 0);
    apqBackward.insertNode(targetNode - 1, 0);
    distForward[sourceNode - 1] = 0;
    distBackward[targetNode - 1] = 0;

    double bestPath = INT_MAX;
    double meetingNode = -1;

    while (!apqForward.isEmpty() && !apqBackward.isEmpty()) {
        double forwardNode = apqForward.getMin().first;
        double backwardNode = apqBackward.getMin().first;

        visitedForward.insert(forwardNode);
        visitedBackward.insert(backwardNode);

        apqForward.popMin();
        apqBackward.popMin();

        if (distForward[forwardNode] + distBackward[forwardNode] < bestPath) {
            meetingNode = forwardNode;
            break;
        }

        if (distForward[forwardNode] + distBackward[forwardNode] < bestPath) {
            bestPath = distForward[forwardNode] + distBackward[forwardNode];
            meetingNode = forwardNode;
        }

        double startForwardEdge = (forwardNode > 0) ? myGraph.edgeStarts[forwardNode - 1] + 1 : 0;
        double endForwardEdge = myGraph.edgeStarts[forwardNode];

        double startBackwardEdge = (backwardNode > 0) ? myGraph.edgeStarts[backwardNode - 1] + 1 : 0;
        double endBackwardEdge = myGraph.edgeStarts[backwardNode];

        for (double forwardEdgeIndex = startForwardEdge; forwardEdgeIndex <= endForwardEdge; forwardEdgeIndex++) {
            double forwardEdge = myGraph.edges[forwardEdgeIndex] - 1;
            double forwardWeight = distance(myGraph.nodes[forwardNode], myGraph.nodes[forwardEdge]);

            if (distForward[forwardNode] + forwardWeight + estimate(forwardNode, forwardEdge, potentials) < distForward[forwardEdge] + estimate(forwardNode, forwardEdge, potentials)) {
                distForward[forwardEdge] = distForward[forwardNode] + forwardWeight;

                double forwardH = estimate(forwardEdge, targetNode, potentials);
                double forwardF = distForward[forwardEdge] + forwardH;

                if (apqForward.contains(forwardEdge)) {
                    apqForward.decreaseKey(forwardEdge, forwardF);
                } else {
                    apqForward.insertNode(forwardEdge, forwardF);
                }
            }
        }

        for (double backwardEdgeIndex = startBackwardEdge; backwardEdgeIndex <= endBackwardEdge; backwardEdgeIndex++) {
            double backwardEdge = myGraph.edges[backwardEdgeIndex] - 1;
            double backwardWeight = distance(myGraph.nodes[backwardNode], myGraph.nodes[backwardEdge]);

            if (distBackward[backwardNode] + backwardWeight + estimate(backwardEdge, backwardNode, potentials) < distBackward[backwardEdge] + estimate(backwardEdge, backwardNode, potentials)) {
                distBackward[backwardEdge] = distBackward[backwardNode] + backwardWeight;

                double backwardH = estimate(sourceNode, backwardEdge, potentials);
                double backwardF = distBackward[backwardEdge] + backwardH;

                if (apqBackward.contains(backwardEdge)) {
                    apqBackward.decreaseKey(backwardEdge, backwardF);
                } else {
                    apqBackward.insertNode(backwardEdge, backwardF);
                }
            }
        }
    }

    if (meetingNode != -1) {
        double shortestPath = distForward[meetingNode] + distBackward[meetingNode];
        if (shortestPath < bestPath)
            bestPath = shortestPath;
    }

    if (bestPath == INT_MAX) {
        return -1;
    }

    return bestPath;
}
