#ifndef INDUCEDCOSTHEURISTIC_H
#define INDUCEDCOSTHEURISTIC_H

#include <cstdint>
#include "edgeheap.h"
#include "clustereditingsolution.h"

/**
 * A heuristic for large cluster editing instances. The heuristic works in iterations,
 * where in each iteration an edge of the graph is either set to permanent or forbidden.
 * Permanent means, the two end nodes will stay in one cluster, forbidden means, they
 * will stay in seperate clusters.
 * 
 * If provided with an unsolvable instance, it will return an empty clustering with
 * infinite editing costs.
 */
class InducedCostHeuristic {

public:
    /**
     * Creates a problem instance using the given graph.
     * 
     * @param param_graph The graph instance to solve.
     * @param bundleEdges Changes how edge merges are handled. Merging an edge
     *      also means that the two end nodes must be merged, which induces
     *      even more edge merges. In the heuristic, nodes are not merged explicitly,
     *      it can just keep track of edges, which actually would be merged. If
     *      this option is enabled, these groups of edges are treated as bundles,
     *      with shared and combined induced costs. If disabled, each edge is treated
     *      individually.
     * 
     */
    InducedCostHeuristic(StaticSparseGraph& param_graph, bool param_bundleEdges);
    
    /**
     * Starts the solving process.
     * 
     * @return An object containing the clusters (with node ids) and the total editing
     *      cost.
     */
    ClusterEditingSolution solve();

private:    
    /**
     * Invoked, when the heuristic decides the make the specified edge permanent.
     * All implications are resolved then.
     */
    void choosePermanentEdge(const StaticSparseGraph::Edge eIcf);

    /**
     * Invoked, when the heuristic decides the make the specified edge forbidden.
     * All implications are resolved then.
     */
    void chooseForbiddenEdge(const StaticSparseGraph::Edge eIcp);
    
    /**
     * Initial run on the instance to fill in all present implications. If nodes u and v
     * are connected by a permanent edge and v and w are connected by a forbidden edge,
     * the edge uw must also be forbidden.
     * 
     * @return false if instance contained contradictions, true otherwise
     */
    bool resolvePermanentForbidden();
    
    /**
     * Sets an edge to forbidden by actually manipulating the graph. Updates the induced
     * costs accordingly. This step does not resolve implications and dependencies
     * and must therefore be called by chooseForbiddenEdge or choosePermanentEdge.
     */
    void setForbidden(const StaticSparseGraph::Edge e);
    
    /**
     * Sets an edge to permanent by actually manipulating the graph. Updates the induced
     * costs accordingly. This step does not resolve implications and dependencies
     * and must therefore be called by chooseForbiddenEdge or choosePermanentEdge.
     */
    void setPermanent(const StaticSparseGraph::Edge e);
    
    /**
     * Updates icf and icp for the edge uw under the assumption that edge uv will be set to forbidden.
     */
    void updateTripleForbiddenUW(const StaticSparseGraph::EdgeWeight uv, const StaticSparseGraph::Edge uw, const StaticSparseGraph::EdgeWeight vw);

    /**
     * Updates icf and icp for the edge uw under the assumption that edge uv will be set to permanent.
     */
    void updateTriplePermanentUW(const StaticSparseGraph::EdgeWeight uv, const StaticSparseGraph::Edge uw, const StaticSparseGraph::EdgeWeight vw);
    
    bool bundleEdges;
    StaticSparseGraph graph;
    EdgeHeap edgeHeap;
    StaticSparseGraph::EdgeWeight totalCost;
    uint64_t totalEdges;
};

#endif
