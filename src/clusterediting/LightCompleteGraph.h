#ifndef LIGHTCOMPLETEGRAPH_H
#define LIGHTCOMPLETEGRAPH_H

#include <vector>
#include <iostream>
#include <list>
#include <cmath>
#include <set>
#include <limits>

namespace ysk {

class LightCompleteGraph {
  
// template<typename T>
// using vector = std::vector<T>;

public:
    typedef uint32_t NodeId;
    typedef uint32_t EdgeId;
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
    LightCompleteGraph(uint32_t numNodes, bool param_pruneZeroEdges);

    //   /**
    //    * Constructs a new graph using the provided lemon graph.
    //    */
    //   LightCompleteGraph(WorkingCopyInstance& inst, bool param_pruneZeroEdges);
    
    /**
    * Creates a hard copy of the provided graph.
    */
    LightCompleteGraph(LightCompleteGraph& other);
    
    /**
    * Clears all edges from the graph and resizes it to the new number of nodes.
    */
    void clearAndResize(const uint32_t newSize);
    
    /**
    * Returns the weight of an edge.
    */
    double getWeight(const Edge e) const;
    
    /**
    * Modifies the weight of an edge, given the Edge.
    */
    void setWeight(const Edge e, const EdgeWeight w);

    /**
    * Modifies the weight of an edge, given the Node IDs.
    */
    void setWeight(NodeId v, NodeId u, const EdgeWeight w);
    
    /**
    * Returns the number of nodes in the graph.
    */
    unsigned int numNodes() const;
    
    /**
   template<typename T>
using vector = std::vector<T>; * Returns the number of edges in the graph.
    */
    unsigned long numEdges() const;
    
    /**
    * For a node v, returns all adjacent nodes, which are connected to v via a perment edge, including v itself.
    */
    const std::vector<NodeId>& getCliqueOf(const NodeId v) const;
    
    /**
    * For a node v, returns all adjacent nodes, which are connected to v via a real valued, non-zero edge.
    */
    const std::vector<NodeId>& getUnprunedNeighbours(const NodeId v) const;

    /**
    * Contracts the specified edge, merging the two end nodes together. This will change the internal addresses of the nodes,
    * so the other two methods for addressing (getInternalId and getOriginalIds) should be used.
    */
    //   void contract(const Edge e);
    
    /**
    * For a node id of the original graph, returns the id of the internal node, which contains the given node.
    */
    //   NodeId getInternalId(const NodeId v) const;
    
    /**
    * For an internal node id, returns the set of original node ids, which are contained in the internal node. The multiplicity
    * might be greater than 1, if edges have been contracted on this graph.
    */
    //   vector<NodeId> getOriginalIds(const NodeId v) const;
    
    private:
    unsigned int size;
    std::vector<double> weights;
    bool pruneZeroEdges;
    //   vector<NodeId> origToCompr;
    //   vector<vector<NodeId>> comprToOrig;
    std::vector<std::vector<NodeId>> cliqueOf;
    std::vector<std::vector<NodeId>> unprunedNeighbours;
    /**
    * Removes a specific node id from the vector.
    */
    bool removeFromVector(std::vector<NodeId> vec, NodeId v);
};

} //namespace ysk

#endif // LIGHTCOMPLETEGRAPH_H
