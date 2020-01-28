#include "clustereditingsolver.h"

ClusterEditingSolution ClusterEditingSolver::run() {
    StaticSparseGraph sGraph(graph);
    InducedCostHeuristic instance(sGraph, bundleEdges);
    ClusterEditingSolution solution = instance.solve();

    return solution;
}
