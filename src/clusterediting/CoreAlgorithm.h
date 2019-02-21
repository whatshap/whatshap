#ifndef SRC_COREALGORITHM_H
#define SRC_COREALGORITHM_H

#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include <set>

#include "ClusterEditingSolutionLight.h"
#include "InducedCostHeuristic.h"
#include "StaticSparseGraph.h"
#include "Globals.h"

class CoreAlgorithm{

    public:

    CoreAlgorithm(
            StaticSparseGraph& graph,
            bool bundleEdges
    )
    :_graph(graph), bundleEdges(bundleEdges)
    {};

    ClusterEditingSolutionLight run();

    /**
     * Attempts a "clean" interrupt of the solving process by stopping CPLEX and setting a kill-flag which is checked throughout the process
     */
    void cancel();

    private:
        StaticSparseGraph& _graph;
        bool bundleEdges;
};

#endif /* SRC_COREALGORITHM_H */
