#ifndef STATICSPARSEGRAPH_H
#define STATICSPARSEGRAPH_H

#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <cmath>
#include <set>
#include <unordered_set>
#include <limits>
#include "dynamicsparsegraph.h"

/**
 * This class implements a pseudo-complete graph with floating point edge weights. Edges with weight zero are not
 * stored explicitly. If a zero edge is queried, the graph will detect that this edge is not stored and return
 * weight zero. The graph is static, i.e. it is not possible to modify edge weights. The only exception is, that
 * non-forbidden, non-permanent edges can be made forbidden or permanent ONCE. For non-zero edges, the stored
 * edge weight is modified accordingly. For zero edges, we use no storage, so whether a formally zero-edge is
 * forbidden or permanent is derived from implicit information:
 * Each node has a clique id. Nodes, which form a clique (connected by permanent edges) have the same clique id.
 * If an edge is set to permanent, then the clique ids are adjusted accordingly. Therefore these clique ids are
 * not stable and not accessable from outside. An edge, which is not stored can be detected as permanent, if the
 * clique ids of its end vertices is equal. This also means that zero-edges can be set to permanent implicitly,
 * if its end vertices end up in the same clique.
 * To keep track of forbidden edges, for each clique the graph stores a set of "forbidden" cliques. This relation
 * is symmetric and is basically a pair-relation between cliques to indicate which cliques cannot be merged
 * anymore. Again, for former zero-edges the graph answers a weight query be checking whether the clique ids
 * of its end vertices are in this forbidden-relation.
 */
class StaticSparseGraph {

public:    
    /**
     * Creates a hard copy of the provided graph.
     */
    StaticSparseGraph(StaticSparseGraph& other);
    
    /**
     * Creates a copy of the provided graph.
     */
    StaticSparseGraph(DynamicSparseGraph& other);
    
    /**
     * Creates a graph, which only contains a subset of nodes of the provided graph.
     */
    StaticSparseGraph(DynamicSparseGraph& other, std::vector<DynamicSparseGraph::NodeId>& nodes);

    /**
     * Returns the weight of an edge.
     */
    DynamicSparseGraph::EdgeWeight getWeight(const DynamicSparseGraph::Edge e);

    /**
     * Returns the weight of an edge by providing its rank id.
     */
    DynamicSparseGraph::EdgeWeight getWeight(const DynamicSparseGraph::RankId r);
    
    /**
     * Returns whether the given edge is permanent.
     */
    bool isPermanent(const DynamicSparseGraph::Edge e);
    
    /**
     * Returns whether the given edge is forbidden.
     */
    bool isForbidden(const DynamicSparseGraph::Edge e);
    
    /**
     * Makes the specified edge e=(u,v) permanent. Any implications raising from this must be handled by the caller.
     */
    void setPermanent(const DynamicSparseGraph::Edge e);
    
    /**
     * Makes the specified edge e=(u,v) forbidden. Any implications raising from this must be handled by the caller.
     */
    void setForbidden(const DynamicSparseGraph::Edge e);
    
    /**
     * Makes the specified edge e=(u,v) permanent. Any implications raising from this must be handled by the caller.
     */
    void setPermanent(const DynamicSparseGraph::Edge e, const DynamicSparseGraph::RankId r);
    
    /**
     * Makes the specified edge e=(u,v) forbidden. Any implications raising from this must be handled by the caller.
     */
    void setForbidden(const DynamicSparseGraph::Edge e, const DynamicSparseGraph::RankId r);
    
    /**
    * Returns the number of nodes in the graph.
    */
    uint64_t numNodes() const;
    
    /**
     * Returns the number of edges in the graph.
     */
    uint64_t numEdges() const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a perment edge, including v itself.
     */
    const std::vector<DynamicSparseGraph::NodeId>& getCliqueOf(const DynamicSparseGraph::NodeId v) const;
    
    /**
     * For a node v, returns all forbidden nodes, with which v cannot be connected.
     */
    const std::vector<DynamicSparseGraph::NodeId> getForbiddenNeighbors(const DynamicSparseGraph::NodeId v) const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a perment edge, including v itself.
     */
    DynamicSparseGraph::NodeId getCliqueIdOf(const DynamicSparseGraph::NodeId v) const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a non infinite, non-zero edge.
     */
    const std::vector<DynamicSparseGraph::NodeId>& getUnprunedNeighbours(const DynamicSparseGraph::NodeId v) const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a non-zero edge.
     */
    const std::vector<DynamicSparseGraph::NodeId>& getNonZeroNeighbours(const DynamicSparseGraph::NodeId v) const;
    
    /**
     * Returns an edge's index in the rank data structure. For non-existing edges, zero is returned.
     */
    DynamicSparseGraph::RankId findIndex(const DynamicSparseGraph::Edge e) const;
    
    /**
     * Returns an edge's index in the rank data structure. For non-existing edges, zero is returned.
     */
    DynamicSparseGraph::RankId findIndex(const DynamicSparseGraph::EdgeId id) const;

private:
	// masks and algorithm to compute popcounts on 64bit words
	const uint64_t m1  = 0x5555555555555555;
	const uint64_t m2  = 0x3333333333333333;
	const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
	const uint64_t h01 = 0x0101010101010101;
	inline uint64_t popcount(uint64_t bitv) const {
		// copied from Wikipedia (https://en.wikipedia.org/wiki/Hamming_weight)
		bitv -= (bitv >> 1) & m1;
		bitv = (bitv & m2) + ((bitv >> 2) & m2);
		bitv = (bitv + (bitv >> 4)) & m4;
		return (bitv * h01) >> 56;
		//return _mm_popcnt_u64(bitv); //TODO: can be used instead as a speedup, but requires <x86intrin.h> to be included
	}
	
    // used for sparse and fast storage of edge weights
    uint64_t size;
    std::vector<uint64_t> rank1;                // size = 8 bytes per 4096 possible edges (= size*(size-1)/1024)
    std::vector<uint64_t> offset1;              // size = rank1
    std::vector<uint64_t> rank2;                // size: best case = 8 bytes per 64 existing edges, worst case = 8 bytes per existing edge
    std::vector<uint64_t> offset2;              // size = rank 2
    std::vector<DynamicSparseGraph::EdgeWeight> weightv;            // size = 8 bytes per existing edge
    
    
    // additional information for efficient iteration over non-zero non-infinity neighbours
    std::vector<std::vector<DynamicSparseGraph::NodeId>> unprunedNeighbours;
    
    // additional information for efficient iteration over non-zero neighbours
    std::vector<std::vector<DynamicSparseGraph::NodeId>> nonzeroNeighbours;
    
    // additional information for detecting, whether a (zero-)edge is acutally permanent/forbidden due to implication
    std::vector<DynamicSparseGraph::NodeId> cliqueOfNode;           // clique id of every node
    std::vector<std::vector<DynamicSparseGraph::NodeId>> cliques;   // all nodes, which belong to a certain cliqueOf
    std::vector<std::unordered_set<DynamicSparseGraph::NodeId>> forbidden;    // for each cluster, a set of forbidden clusters
    
    /**
     * Inserts all added edges into the static data structure of this graph.
     */
    void compile(DynamicSparseGraph& dg, const std::vector<DynamicSparseGraph::NodeId>& nodes);
    
    /**
     * Refreshes interal data about edges. Necessary for consistency.
     */
    void refreshEdgeMetaData(const DynamicSparseGraph::Edge e, const DynamicSparseGraph::EdgeWeight oldW, const DynamicSparseGraph::EdgeWeight newW);
    
    /**
     * Removes a specific node id from the vector.
     */
    bool removeFromVector(std::vector<DynamicSparseGraph::NodeId>& vec, DynamicSparseGraph::NodeId v);
    
    std::string int2bin(uint64_t u) {
        std::string s(64, '0');
        uint64_t mask = (1UL << 63); // fill in values right-to-left
        for (int i = 0; i < 64; i++, mask >>= 1)
            if ((u & mask) != 0)
                s[i] = '1';
        return s;
    }
};

#endif
