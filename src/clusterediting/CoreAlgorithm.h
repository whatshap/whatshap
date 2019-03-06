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
    CoreAlgorithm(DynamicSparseGraph& graph, bool bundleEdges) :
    graph(graph), 
    bundleEdges(bundleEdges)
    {};

    /**
     * Start the algorithm on the provided graph and settings
     */
    ClusterEditingSolutionLight run();

private:
    DynamicSparseGraph& graph;
    bool bundleEdges;
};

#endif /* SRC_COREALGORITHM_H */
