#ifndef CLUSTEREDITINGSOLVER_H
#define CLUSTEREDITINGSOLVER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include <set>

#include "clustereditingsolution.h"
#include "inducedcostheuristic.h"
#include "staticsparsegraph.h"

/**
 * Central solver for cluster editing instances. Uses a InducedCostHeuristic to
 * determine a low-cost transformation of the given graph into a clique graph.
 */
class ClusterEditingSolver {

public:
    /**
     * Creates a solver object using the provided weighted graph.
     * 
     * @param param_graph The graph instance to solve.
     * @param bundleEdges Changes how edge merges are handled. Merging an edge
     *      also means that the two end nodes must be merged, which induces
     *      even more edge merges. In the heuristic, nodes are not merged explicitly,
     *      it can just keep track of edges, which actually would be merged. If
     *      this option is enabled, these groups of edges are treated as bundles,
     *      with shared and combined induced costs. If disabled, each edge is treated
     *      individually.
     */
    ClusterEditingSolver(DynamicSparseGraph& graph, bool bundleEdges) :
		graph(graph), 
		bundleEdges(bundleEdges)
    {};

    /**
     * Starts the solving process on the provided instance.
     * 
     * @return An object containing the clusters (with node ids) and the total editing
     *      cost.
     */
    ClusterEditingSolution run();

private:
    DynamicSparseGraph& graph;
    bool bundleEdges;
};

#endif

