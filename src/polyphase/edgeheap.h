#ifndef EDGEHEAP_H
#define EDGEHEAP_H

#include <vector>
#include <set>
#include "staticsparsegraph.h"
    
/**
 * Data structure to manage the induced costs of a graph's edges.
 */
class EdgeHeap {
public:

    /**
    * Constructs a new instance using the provided graph to precompute all icf and icp values.
    */
    EdgeHeap(StaticSparseGraph& param_graph);
  
    /**
     * Initializes the induced costs for all edges. May take quite long!
     */
    void initInducedCosts();
    
    /**
    * Returns the edge with the highest icf.
    */
    DynamicSparseGraph::Edge getMaxIcfEdge() const;
    
    /**
    * Returns the edge with the highest icp.
    */
    DynamicSparseGraph::Edge getMaxIcpEdge() const;
    
    /**
    * Returns the icf of the provided edge.
    */
    DynamicSparseGraph::EdgeWeight getIcf(const DynamicSparseGraph::Edge e) const;
    
    /**
    * Returns the icp of the provided edge.
    */
    DynamicSparseGraph::EdgeWeight getIcp(const DynamicSparseGraph::Edge e) const;
    
    /**
    * Adds the provided value to the icf of the given edge.
    */
    void increaseIcf(const DynamicSparseGraph::Edge e, const DynamicSparseGraph::EdgeWeight w);
    
    /**
    * Adds the provided value to the icp of the given edge.
    */
    void increaseIcp(const DynamicSparseGraph::Edge e, const DynamicSparseGraph::EdgeWeight w);
    
    /**
     * Bundles two edges together. This indicates that if either of the edges is set to permanent or forbidden, the other
     * must be as well. One of the edges receives the induced costs for both edges, while the other is removed from the heap.
     */
    void mergeEdges(const DynamicSparseGraph::Edge e1, const DynamicSparseGraph::Edge e2);
    
    /**
    * Removes the specified edge.
    */
    void removeEdge(const DynamicSparseGraph::Edge e);
    
    /**
     * Computes the induced cost for the the triple uvw, if uv is set to forbidden
     */
    inline DynamicSparseGraph::EdgeWeight getIcf(const DynamicSparseGraph::EdgeWeight uw, const DynamicSparseGraph::EdgeWeight vw) const  {
        if (uw > 0 && vw > 0) {
            // if both other edges present, remove the cheapest of both
            return std::min(uw, vw); 
        } else {
            return 0;
        }
    }
    
    /*
     * Computes the induced cost for the the triple uvw, if uv is set to permanent
     */
    inline DynamicSparseGraph::EdgeWeight getIcp(const DynamicSparseGraph::EdgeWeight uw, const DynamicSparseGraph::EdgeWeight vw) const {
        if (uw < 0 && vw > 0) {
            return std::min(vw, -uw); 	// either add uw or remove vw
        } else if (uw > 0 && vw < 0) {
            return std::min(-vw, uw); 	// either add vw or remove uw
        } else {
            return 0;
        }
    }
    
    uint64_t numUnprocessed() const;

private:
    /**
     * Removes the edge with the specified rank
     */
    
    void removeEdge(const DynamicSparseGraph::RankId rId);
    /**
     * Ensures that the heap structure of the given heap stays intact after the icf/icp value of an edge has been modified. 
     * Provided are the id of the modified edge, new and old value, an index (which maps edge ids to their position in the heap)
     * and a score vector (which maps an edge id to either its icf or icp).
     */
    
    void updateHeap(std::vector<DynamicSparseGraph::RankId>& heap, const DynamicSparseGraph::RankId e, const DynamicSparseGraph::EdgeWeight change, 
            std::vector<DynamicSparseGraph::RankId>& index, const std::vector<DynamicSparseGraph::EdgeWeight>& score);
    
    StaticSparseGraph& graph;
    uint64_t unprocessed;
    std::vector<DynamicSparseGraph::Edge> edges;                         // internal/rank id -> edge
    std::vector<DynamicSparseGraph::EdgeWeight> icf;                     // edge rank id -> icf of edge (zero edges have no icf)
    std::vector<DynamicSparseGraph::EdgeWeight> icp;                     // edge rank id -> icp of edge (zero edges have no icp)
    std::vector<DynamicSparseGraph::RankId> forb_rank2edge;              // max-heap over rank edge ids, sorted by icf
    std::vector<DynamicSparseGraph::RankId> perm_rank2edge;              // max-heap over rank edge ids, sorted by icp
    std::vector<DynamicSparseGraph::RankId> edge2forb_rank;              // rank edge id -> position in icf-heap
    std::vector<DynamicSparseGraph::RankId> edge2perm_rank;              // rank edge id -> position in icp-heap
    std::vector<DynamicSparseGraph::RankId> edgeToBundle;                // edge rank id -> representant of edges bunch
    std::vector<std::vector<DynamicSparseGraph::RankId>> edgeBundles;    // edge bunch representant -> set of edges belonging to this bunch
};

#endif

