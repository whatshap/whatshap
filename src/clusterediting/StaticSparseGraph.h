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
#include "Globals.h"
#include "DynamicSparseGraph.h"

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
     * Makes the specified edge e=(u,v) permanent. If u and v are connected to other nodes via permanent edges, then all edges
     * between u's clique and v's clique are made permanent as well.
     */
    void setPermanent(const DynamicSparseGraph::Edge e);
    
    /**
     * Makes the specified edge e=(u,v) forbidden. If u and v are connected to other nodes via permanent edges, then all edges
     * between u's clique and v's clique are made forbidden as well.
     */
    void setForbidden(const DynamicSparseGraph::Edge e);
    
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

#endif // STATICSPARSEGRAPH_H
