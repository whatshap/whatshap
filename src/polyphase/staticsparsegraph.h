#ifndef STATICSPARSEGRAPH_H
#define STATICSPARSEGRAPH_H

#include <cstdint>
#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <cmath>
#include <set>
#include <unordered_set>
#include <limits>
#include "trianglesparsematrix.h"

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
    typedef uint32_t NodeId;
    typedef uint64_t EdgeId;
    typedef uint32_t RankId;
    typedef float EdgeWeight;
    
    static inline uint64_t popcount(uint64_t bitv) {
		// copied from Wikipedia (https://en.wikipedia.org/wiki/Hamming_weight)
		bitv -= (bitv >> 1) & m1;
		bitv = (bitv & m2) + ((bitv >> 2) & m2);
		bitv = (bitv + (bitv >> 4)) & m4;
		return (bitv * h01) >> 56;
		//return _mm_popcnt_u64(bitv); //TODO: can be used instead as a speedup, but requires <x86intrin.h> to be included
	}
    
    /**
    * Compact data structure to represent an edge. It consists of two node indices.
    */
    struct Edge {
        NodeId u;
        NodeId v;
        
        Edge(NodeId pu, NodeId pv) {
            int ordered = pu < pv;
            u = ordered*pu + (1-ordered)*pv;
            v = ordered*pv + (1-ordered)*pu;
        }
        
        Edge(EdgeId i) {
            u = std::ceil(std::sqrt(2*(i+1)+0.25) - 0.5);
            v = (NodeId)(i - (uint64_t)u * (uint64_t)(u-1) / 2);
        }
        
        Edge() : u(0), v(1) {};
        
        /**
        * Returns the id of this edge for a triangle adjacency matrix representation.
        */
        EdgeId id() const {
            return ((uint64_t)v)*((uint64_t)(v-1))/2 + (uint64_t)u;
        }
        
        bool operator==(const Edge& other) const {
            return u == other.u && v == other.v;
        }
    };
    
    static const EdgeWeight Forbidden;
    static const EdgeWeight Permanent;
    static const Edge InvalidEdge;
    static const EdgeId InvalidEdgeId;
    static const NodeId InvalidNodeId;
    
    /**
     * Creates a hard copy of the provided graph.
     */
    StaticSparseGraph(StaticSparseGraph& other);
    
    /**
     * Creates a hard copy of the provided graph.
     */
    StaticSparseGraph(TriangleSparseMatrix& other);

    /**
     * Returns the weight of an edge.
     */
    EdgeWeight getWeight(const Edge e);

    /**
     * Returns the weight of an edge by providing its rank id.
     */
    EdgeWeight getWeight(const RankId r);
    
    /**
     * Returns whether the given edge is permanent.
     */
    bool isPermanent(const Edge e);
    
    /**
     * Returns whether the given edge is forbidden.
     */
    bool isForbidden(const Edge e);
    
    /**
     * Makes the specified edge e=(u,v) permanent. Any implications raising from this must be handled by the caller.
     */
    void setPermanent(const Edge e);
    
    /**
     * Makes the specified edge e=(u,v) forbidden. Any implications raising from this must be handled by the caller.
     */
    void setForbidden(const Edge e);
    
    /**
     * Makes the specified edge e=(u,v) permanent. Any implications raising from this must be handled by the caller.
     */
    void setPermanent(const Edge e, const RankId r);
    
    /**
     * Makes the specified edge e=(u,v) forbidden. Any implications raising from this must be handled by the caller.
     */
    void setForbidden(const Edge e, const RankId r);
    
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
    const std::vector<NodeId>& getCliqueOf(const NodeId v) const;
    
    /**
     * For a node v, returns all forbidden nodes, with which v cannot be connected.
     */
    const std::vector<NodeId> getForbiddenNeighbors(const NodeId v) const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a perment edge, including v itself.
     */
    NodeId getCliqueIdOf(const NodeId v) const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a non infinite, non-zero edge.
     */
    const std::vector<NodeId>& getUnprunedNeighbours(const NodeId v) const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a non-zero edge.
     */
    const std::vector<NodeId>& getNonZeroNeighbours(const NodeId v) const;
    
    /**
     * Returns an edge's index in the rank data structure. For non-existing edges, zero is returned.
     */
    RankId findIndex(const Edge e) const;
    
    /**
     * Returns an edge's index in the rank data structure. For non-existing edges, zero is returned.
     */
    RankId findIndex(const EdgeId id) const;

private:
	// masks and algorithm to compute popcounts on 64bit words
	static const uint64_t m1  = 0x5555555555555555;
	static const uint64_t m2  = 0x3333333333333333;
	static const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
	static const uint64_t h01 = 0x0101010101010101;
	
    // used for sparse and fast storage of edge weights
    uint64_t size;
    std::vector<uint64_t> rank1;                // size = 8 bytes per 4096 possible edges (= size*(size-1)/1024)
    std::vector<uint64_t> offset1;              // size = rank1
    std::vector<uint64_t> rank2;                // size: best case = 8 bytes per 64 existing edges, worst case = 8 bytes per existing edge
    std::vector<uint64_t> offset2;              // size = rank 2
    std::vector<EdgeWeight> weightv;            // size = 8 bytes per existing edge
    
    
    // additional information for efficient iteration over non-zero non-infinity neighbours
    std::vector<std::vector<NodeId>> unprunedNeighbours;
    
    // additional information for efficient iteration over non-zero neighbours
    std::vector<std::vector<NodeId>> nonzeroNeighbours;
    
    // additional information for detecting, whether a (zero-)edge is acutally permanent/forbidden due to implication
    std::vector<NodeId> cliqueOfNode;           // clique id of every node
    std::vector<std::vector<NodeId>> cliques;   // all nodes, which belong to a certain cliqueOf
    std::vector<std::unordered_set<NodeId>> forbidden;    // for each cluster, a set of forbidden clusters
    
    /**
     * Inserts all added edges into the static data structure of this graph.
     */
    void compile(TriangleSparseMatrix& m);
    
    /**
     * Refreshes interal data about edges. Necessary for consistency.
     */
    void refreshEdgeMetaData(const Edge e, const EdgeWeight oldW, const EdgeWeight newW);
    
    /**
     * Removes a specific node id from the vector.
     */
    bool removeFromVector(std::vector<NodeId>& vec, NodeId v);
};

#endif
