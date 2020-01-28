#ifndef DYNAMICSPARSEGRAPH_H
#define DYNAMICSPARSEGRAPH_H

#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <cmath>
#include <set>
#include <limits>

/**
 * A sparse graph data structure, which is based on adjacency lists and a hashtable
 * to map edge ids to weights.
 */
class DynamicSparseGraph {

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
     * Constructs a new empty graph with no nodes or edges.
     */
    DynamicSparseGraph();
    
    /**
     * Constructs a new graph with the specified number of nodes and no edges (non-existing edges implicitly have a weight of 0)
     */
    DynamicSparseGraph(uint32_t numNodes);
    
    /**
     * Creates a hard copy of the provided graph.
     */
    DynamicSparseGraph(DynamicSparseGraph& other);
    
    /**
     * Clears all edges from the graph and resizes it to the new number of nodes.
     */
    void clearAndResize(const uint32_t newSize);
    
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
    EdgeWeight getWeight(const Edge e) const;
    
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
     * Returns a reference to the internal adjacency lists. Node i's adjacency is stored in the i-th list of the reference.
     */
    const std::vector<NodeId>& getNeighbours(NodeId u) const;
    
    /**
     * Returns the positively connected components of the graph as a vector over a vector of node indices.
     */
    std::vector<std::vector<NodeId>> getPositiveComponentes() const;
     

private:
    // used for dynamic addition of edges
    uint64_t size;
    std::map<EdgeId, EdgeWeight> weights;
    std::vector<std::vector<NodeId>> neighbours;
};

#endif
