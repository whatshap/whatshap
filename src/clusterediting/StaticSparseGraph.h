#ifndef STATICSPARSEGRAPH_H
#define STATICSPARSEGRAPH_H

#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <cmath>
#include <set>
#include <limits>
#include "Globals.h"

namespace ysk {

class StaticSparseGraph {

public:
    typedef uint32_t NodeId;
    typedef uint64_t EdgeId;
    typedef uint64_t RankId;
    typedef double EdgeWeight;
    
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
            v = i - u * (u-1) / 2;
        }
        
        Edge() : u(0), v(1) {};
        
        /**
        * Returns the id of this edge for a triangle adjacency matrix representation.
        */
        EdgeId id() const {
            return v*(v-1)/2 + u;
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
     * Constructs a new empty graph with the specified number of nodes. Initializes all (non-)edges with weight 0.
     */
    StaticSparseGraph(uint32_t numNodes);
    
    /**
     * Creates a hard copy of the provided graph.
     */
    StaticSparseGraph(StaticSparseGraph& other);
    
    /**
     * Clears all edges from the graph and resizes it to the new number of nodes.
     */
    void clearAndResize(const uint32_t newSize);
    
    /**
     * Inserts all added edges into the static data structure of this graph.
     */
    void compile();
    
    /**
     * Adds the given edge to the graph. Effect will only be visible after the graph has been compiled.
     */
    void addEdge(const Edge e, const EdgeWeight w);
    
    /**
     * Adds the given edge to the graph. Effect will only be visible after the graph has been compiled.
     */
    void addEdge(const NodeId v, const NodeId u, const EdgeWeight w);
    
    /**
     * Returns the weight of an edge.
     */
    EdgeWeight getWeight(const Edge e);
    
    /**
     * Returns the weight of an edge by providing its rank id.
     */
    EdgeWeight getWeight(const RankId r);
    
    /**
     * Makes the specified edge e=(u,v) permanent. If u and v are connected to other nodes via permanent edges, then all edges
     * between u's clique and v's clique are made permanent as well.
     */
    void setPermanent(const Edge e);
    
    /**
     * Makes the specified edge e=(u,v) forbidden. If u and v are connected to other nodes via permanent edges, then all edges
     * between u's clique and v's clique are made forbidden as well.
     */
    void setForbidden(const Edge e);
    
    /**
     * Modifies the weight of an edge, given the Edge.
     */
    void setWeight(const Edge e, const EdgeWeight w);

    /**
     * Modifies the weight of an edge, given the Node IDs.
     */
    void setWeight(const NodeId v, const NodeId u, const EdgeWeight w);
    
    /**
    * Returns the number of nodes in the graph.
    */
    unsigned int numNodes() const;
    
    /**
     * Returns the number of edges in the graph.
     */
    unsigned long numEdges() const;
    
    /**
     * For a node v, returns all adjacent nodes, which are connected to v via a perment edge, including v itself.
     */
    const std::vector<NodeId>& getCliqueOf(const NodeId v) const;
    
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
    // used for inefficient popcount
    const uint64_t m1  = 0x5555555555555555;
    const uint64_t m2  = 0x3333333333333333;
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
    const uint64_t h01 = 0x0101010101010101;
    
    // used for dynamic addition of edges
    uint64_t size;
    std::map<EdgeId, EdgeWeight> weights;
    std::vector<std::vector<NodeId>> neighbours;
    
    // used for sparse and fast storage of edge weights
    bool compiled;
    std::vector<uint64_t> rank1;                // size = 8 bytes per 4096 possible edges (= size*(size-1)/1024)
    std::vector<uint64_t> offset1;              // size = rank1
    std::vector<uint64_t> rank2;                // size: best case = 8 bytes per 64 existing edges, worst case = 8 bytes per existing edge
    std::vector<uint64_t> offset2;              // size = rank 2
    std::vector<EdgeWeight> weightv;            // size = 8 bytes per existing edge
    
    // additional information for detecting, whether a (zero-)edge is acutally permanent/forbidden due to implication
    std::vector<NodeId> cliqueOfNode;           // clique id of every node
    std::vector<std::vector<NodeId>> cliques;   // all nodes, which belong to a certain cliqueOf
    std::vector<std::set<NodeId>> forbidden;    // for each cluster, a set of forbidden clusters
    
    // additional information for efficient iteration over non-zero non-infinity neighbours
    std::vector<std::vector<NodeId>> unprunedNeighbours;
    
    // additional information for efficient iteration over non-zero neighbours
    std::vector<std::vector<NodeId>> nonzeroNeighbours;
    
    /**
     * Refreshes interal data about edges. Necessary for consistency.
     */
    void refreshEdgeMetaData(const Edge e, const EdgeWeight oldW, const EdgeWeight newW);
    
    /**
     * Removes a specific node id from the vector.
     */
    bool removeFromVector(std::vector<NodeId>& vec, NodeId v);
    
    /**
     * Returns the number of set bits in a 64bit-word.
     */
    uint64_t popcount(uint64_t bitv) const;
    
    std::string int2bin(uint64_t u) {
        std::string s(64, '0');
        uint64_t mask = (1UL << 63); // fill in values right-to-left
        for (int i = 0; i < 64; i++, mask >>= 1)
            if ((u & mask) != 0)
                s[i] = '1';
        return s;
    }
};

} //namespace ysk

#endif // STATICSPARSEGRAPH_H
