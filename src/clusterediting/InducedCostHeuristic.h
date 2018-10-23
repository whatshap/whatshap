#ifndef INDUCEDCOSTHEURISTICLIGHT_H
#define INDUCEDCOSTHEURISTICLIGHT_H

#include "LightCompleteGraph.h"
#include "EdgeHeap.h"
#include "ClusterEditingSolutionLight.h"

namespace ysk {

class InducedCostHeuristic {
  

public:
    InducedCostHeuristic(LightCompleteGraph& param_graph, bool param_pruneZeroEdges);
    ClusterEditingSolutionLight solve();

private:    
    void init();
    bool resolvePermanentForbidden();
    void setForbidden(const LightCompleteGraph::Edge e);
    void setPermanent(const LightCompleteGraph::Edge e);
    
    /**
    * Updates icf and icp for the edge uw under the assumption that edge uv will be set to forbidden.
    */
    void updateTripleForbiddenUW(const LightCompleteGraph::EdgeWeight uv, const LightCompleteGraph::Edge uw, const LightCompleteGraph::EdgeWeight vw);

    /**
    * Updates icf and icp for the edge uw under the assumption that edge uv will be set to permanent.
    */
    void updateTriplePermanentUW(const LightCompleteGraph::EdgeWeight uv, const LightCompleteGraph::Edge uw, const LightCompleteGraph::EdgeWeight vw);
    
    bool pruneZeroEdges;
    LightCompleteGraph graph;
    EdgeHeap edgeHeap;
    LightCompleteGraph::EdgeWeight totalCost;
};

} // namespace ysk

#endif // INDUCEDCOSTHEURISTICLIGHT_H
