//
// Created by alvar on 05/06/2023.
//
/*
#ifndef PRAKTIKUM_DUMPSTER_H
#define PRAKTIKUM_DUMPSTER_H

#include <cmath>
#include <unordered_set>
#include "Graph.hh"
#include "APQ.hh"

void BiBiALT (Graph& g, int source, int target, std::vector<std::vector<double>>& vec) {

    APQ Qf;
    APQ Qb;

    std::vector<double> df(g.getNodes().size(), INFINITY);
    std::vector<int> parentf(g.getNodes().size());

    std::vector<double> db(g.getNodes().size(), INFINITY);
    std::vector<int> parentb(g.getNodes().size());

    std::vector<bool> visitedT(g.getNodes().size(), false);
    std::vector<bool> visitedS(g.getNodes().size(), false);

    visitedS[source] = true;
    parentf[source] = source;
    df[source] = 0;
    Qf.insertNode(source,0);

    visitedT[target] = true;
    parentb[target] = target;
    db[target] = 0;
    Qb.insertNode(target,0);

    double mu = INFINITY;

    std::unordered_set<int> Sf;
    std::unordered_set<int> Sb;

    while (!Qf.isEmpty() && !Qb.isEmpty()){
        std::pair<int,double> uf = Qf.getMin();
        std::pair<int,double> ub = Qb.getMin();
        Qf.popMin();
        Qb.popMin();

        Sf.insert(uf.first);
        Sb.insert(ub.first);

        for (int v = 0; v < g.getNodes()[uf.first + 1].second - g.getNodes()[uf.first].second; v++) {

            int nodeV = g.getEdges()[g.getNodes()[uf.first].second+v].getId();
            double eucDist = euclideanDistance(g.getNodes()[uf.first].first, g.getNodes()[nodeV].first);

            if(df[uf.first] + eucDist + getPotential(target, uf.first, vec) < df[nodeV] + getPotential(target, uf.first, vec)){
                parentf[nodeV] = uf.first;
                df[nodeV] = df[uf.first] + eucDist;
                if(!visitedT[nodeV]) {
                    Qf.insert(nodeV, eucDist + df[uf.first] + getPotential(target, nodeV, vec));
                    visitedT[nodeV] = true;
                } else {
                    Qf.decreaseKey(nodeV, eucDist + df[uf.first] + getPotential(target, nodeV, vec));
                }
            }
            if(CheckInSet(Sb, nodeV) && df[uf.first] + eucDist + db[nodeV] < mu) {
                mu = df[uf.first] + eucDist + db[nodeV];
            }
        }

        for (int v = 0; v < g.getNodes()[ub.first + 1].second - g.getNodes()[ub.first].second; v++) {

            int nodeV = g.getEdges()[g.getNodes()[ub.first].second+v].getId();
            double eucDist = euclideanDistance(g.getNodes()[ub.first].first, g.getNodes()[nodeV].first);

            if(db[ub.first] + eucDist + getPotential(source, ub.first, vec) < db[nodeV] + getPotential(source, ub.first, vec)){
                parentb[nodeV] = ub.first;
                db[nodeV] = db[ub.first] + eucDist;
                if(!visitedS[nodeV]) {
                    Qb.insert(nodeV, eucDist + db[ub.first] + getPotential(source, nodeV, vec));
                    visitedS[nodeV] = true;
                } else {
                    Qb.decreaseKey(nodeV, eucDist + db[ub.first] + getPotential(source, nodeV, vec));
                }
            }
            if(CheckInSet(Sf, nodeV) && db[ub.first] + eucDist + df[nodeV] < mu) {
                mu = db[ub.first] + eucDist + df[nodeV];
            }
        }



        if (CheckInSet(Sf, ub.first) && CheckInSet(Sb, uf.first)) {
            std::cout<<"distance has been found with Bidirectional Dijkstra between " << source << " and "  << target << " it is " << mu <<std::endl;
            return;
        }

        /*if(Qf.top().second >= mu || Qb.top().second >= mu) {
            std::cout<<"distance has been found with Bidirectional Dijkstra between " << source << " and "  << target << " it is " << mu <<std::endl;
            return;
        }

        }
    return;
}
#endif //PRAKTIKUM_DUMPSTER_H
*/
