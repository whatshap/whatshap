#ifndef INDUCEDCOSTHEURISTICLIGHT_H
#define INDUCEDCOSTHEURISTICLIGHT_H

#include "EdgeHeap.h"
#include "ClusterEditingSolutionLight.h"

class InducedCostHeuristic {

public:
    InducedCostHeuristic(StaticSparseGraph& param_graph, bool param_bundleEdges);
    ClusterEditingSolutionLight solve();

private:    
    void init();
    void choosePermanentEdge(const DynamicSparseGraph::Edge eIcf);
    void chooseForbiddenEdge(const DynamicSparseGraph::Edge eIcp);
    bool resolvePermanentForbidden();
    void setForbidden(const DynamicSparseGraph::Edge e);
    void setPermanent(const DynamicSparseGraph::Edge e);
    void printHeuristicProgress();
    
    /**
    * Updates icf and icp for the edge uw under the assumption that edge uv will be set to forbidden.
    */
    void updateTripleForbiddenUW(const DynamicSparseGraph::EdgeWeight uv, const DynamicSparseGraph::Edge uw, const DynamicSparseGraph::EdgeWeight vw);

    /**
    * Updates icf and icp for the edge uw under the assumption that edge uv will be set to permanent.
    */
    void updateTriplePermanentUW(const DynamicSparseGraph::EdgeWeight uv, const DynamicSparseGraph::Edge uw, const DynamicSparseGraph::EdgeWeight vw);
    
    bool bundleEdges;
    StaticSparseGraph graph;
    EdgeHeap edgeHeap;
    DynamicSparseGraph::EdgeWeight totalCost;
    uint64_t totalEdges;
};

#endif // INDUCEDCOSTHEURISTICLIGHT_H
