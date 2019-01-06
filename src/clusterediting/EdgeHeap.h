#ifndef EDGEHEAP_H
#define EDGEHEAP_H

#include <vector>
#include <set>
#include "StaticSparseGraph.h"

namespace ysk {

class EdgeHeap {
public:

    /**
    * Constructs a new instance using the provided graph to precompute all icf and icp values.
    */
    EdgeHeap(StaticSparseGraph& param_graph, bool param_pruneZeroEdges);
  
    /**
     * Initializes the induced costs for all edges. May take quite long!
     */
    void initInducedCosts();
    /**
    * Returns the edge with the highest icf.
    */
    StaticSparseGraph::Edge getMaxIcfEdge() const;
    /**
    * Returns the edge with the highest icp.
    */
    StaticSparseGraph::Edge getMaxIcpEdge() const;
    /**
    * Returns the icf of the provided edge.
    */
    StaticSparseGraph::EdgeWeight getIcf(const StaticSparseGraph::Edge e) const;
    /**
    * Returns the icp of the provided edge.
    */
    StaticSparseGraph::EdgeWeight getIcp(const StaticSparseGraph::Edge e) const;
    /**
    * Adds the provided value to the icf of the given edge.
    */
    void increaseIcf(const StaticSparseGraph::Edge e, const StaticSparseGraph::EdgeWeight w);
    /**
    * Adds the provided value to the icp of the given edge.
    */
    void increaseIcp(const StaticSparseGraph::Edge e, const StaticSparseGraph::EdgeWeight w);
    /**
    * Removes the specified edge.
    */
    void removeEdge(const StaticSparseGraph::Edge e);
    /**
        * Computes the induced cost for the the triple uvw, if uv is set to forbidden
        */
    StaticSparseGraph::EdgeWeight getIcf(const StaticSparseGraph::EdgeWeight uw, const StaticSparseGraph::EdgeWeight vw) const;
    /**
        * Computes the induced cost for the the triple uvw, if uv is set to permanent
        */
    StaticSparseGraph::EdgeWeight getIcp(const StaticSparseGraph::EdgeWeight uw, const StaticSparseGraph::EdgeWeight vw) const;
    
    int numUnprocessed() const;

private:
    /**
    * Ensures that the heap structure of the given heap stays intact after the icf/icp value of an edge has been modified. 
    * Provided are the id of the modified edge, new and old value, an index (which maps edge ids to their position in the heap)
    * and a score vector (which maps an edge id to either its icf or icp).
    */
    void updateHeap(std::vector<StaticSparseGraph::RankId>& heap, const StaticSparseGraph::RankId e, const StaticSparseGraph::EdgeWeight change, 
            std::vector<StaticSparseGraph::RankId>& index, const std::vector<StaticSparseGraph::EdgeWeight>& score);
    
    bool pruneZeroEdges;
    StaticSparseGraph& graph;
    int unprocessed;
    std::vector<StaticSparseGraph::Edge> edges;             // internal/rank id -> edge
    std::vector<StaticSparseGraph::EdgeWeight> icf;		    // edge rank id -> icf of edge (zero edges have no icf)
    std::vector<StaticSparseGraph::EdgeWeight> icp;		    // edge rank id -> icp of edge (zero edges have no icp)
    std::vector<StaticSparseGraph::RankId> forb_rank2edge;	// max-heap over rank edge ids, sorted by icf
    std::vector<StaticSparseGraph::RankId> perm_rank2edge;	// max-heap over rank edge ids, sorted by icp
    std::vector<StaticSparseGraph::RankId> edge2forb_rank;	// rank edge id -> position in icf-heap
    std::vector<StaticSparseGraph::RankId> edge2perm_rank;	// rank edge id -> position in icp-heap
};

} // namespace ysk;

#endif // EDGEHEAP_H

