#include "clustereditingsolver.h"

ClusterEditingSolution ClusterEditingSolver::run() {
    StaticSparseGraph sGraph(m);
    InducedCostHeuristic instance(sGraph, bundleEdges);
    ClusterEditingSolution solution = instance.solve();

    return solution;
}
