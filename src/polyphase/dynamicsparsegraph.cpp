#include "dynamicsparsegraph.h"
#include <algorithm>
#include <queue>

const DynamicSparseGraph::EdgeWeight DynamicSparseGraph::Forbidden = -std::numeric_limits<EdgeWeight>::infinity();
const DynamicSparseGraph::EdgeWeight DynamicSparseGraph::Permanent = std::numeric_limits<EdgeWeight>::infinity();
const DynamicSparseGraph::Edge DynamicSparseGraph::InvalidEdge = {std::numeric_limits<NodeId>::max(), std::numeric_limits<NodeId>::max()};
const DynamicSparseGraph::EdgeId DynamicSparseGraph::InvalidEdgeId = -1;
const DynamicSparseGraph::NodeId DynamicSparseGraph::InvalidNodeId = -1;

using Edge = DynamicSparseGraph::Edge;
using EdgeWeight = DynamicSparseGraph::EdgeWeight;
using EdgeId = DynamicSparseGraph::EdgeId;
using RankId = DynamicSparseGraph::RankId;
using NodeId = DynamicSparseGraph::NodeId;

DynamicSparseGraph::DynamicSparseGraph() :
    size(0),
    neighbours(size, std::vector<NodeId>(0))
{
    
}

DynamicSparseGraph::DynamicSparseGraph(uint32_t numNodes) :
    size(numNodes),
    neighbours(size, std::vector<NodeId>(0))
{
    
}

DynamicSparseGraph::DynamicSparseGraph(DynamicSparseGraph& other) :
    size(other.size),
    weights(other.weights),
    neighbours(other.neighbours)
{
    
}

void DynamicSparseGraph::clearAndResize(const uint32_t newSize) {
    // clear old stuff
    neighbours.clear();
    
    // resize
    size = newSize;
    neighbours.clear();
    weights.clear();
    neighbours.resize(size, std::vector<NodeId>(0));
}

void DynamicSparseGraph::addEdge(const Edge e, const EdgeWeight w) {
    if (w != 0.0) {
        if (e.v >= size) {
            neighbours.resize(e.v+1, std::vector<NodeId>(0));
            size = e.v+1;
        }
        neighbours[e.v].push_back(e.u);
        weights[e.id()] = w;
    }
}

void DynamicSparseGraph::addEdge(const NodeId v, const NodeId u, const EdgeWeight w) {
    addEdge(Edge(v, u), w);
}

EdgeWeight DynamicSparseGraph::getWeight(const Edge e) const {
    std::map<EdgeId, EdgeWeight>::const_iterator it = weights.find(e.id());
    if (it == weights.end()) {
        return 0;
    } else {
        return it->second;
    }
}

void DynamicSparseGraph::setWeight(const Edge e, const EdgeWeight w) {
    std::map<EdgeId, EdgeWeight>::iterator it = weights.find(e.id());
    if (it == weights.end()) {
        std::cout<<"setWeight tried to modify non existing edge ("<<e.u<<","<<e.v<<") to "<<w<<std::endl;
        return;
    } else {
        it->second = w;
    }
}

void DynamicSparseGraph::setWeight(NodeId v, NodeId u, const EdgeWeight w){
    Edge e(v,u);
    setWeight(e,w);
}

unsigned int DynamicSparseGraph::numNodes() const {
    return size;
}

unsigned long DynamicSparseGraph::numEdges() const {
    return weights.size() - 1;
}

const std::vector<NodeId>& DynamicSparseGraph::getNeighbours(NodeId u) const {
    return neighbours[u];
}

std::vector<std::vector<NodeId>> DynamicSparseGraph::getPositiveComponentes() const {
    // create reverse neighbours temporarily
    std::vector<std::vector<NodeId>> revNeighbours(size, std::vector<NodeId>(0));
    for (NodeId u = 0; u < size; u++) {
        for (NodeId v : neighbours[u]) {
            revNeighbours[v].push_back(u);
        }
    }
    
    std::vector<std::vector<NodeId>> components;
    std::vector<int> componentOfNode(size, -1);
    for (NodeId u = size - 1; u < size; u--) {
        // add component if not explored yet
        int c = componentOfNode[u];
        if (c >= 0) {
            continue;
        }
		c = components.size();
		componentOfNode[u] = c;
		components.push_back(std::vector<NodeId>(0));
        
        // add all connected nodes to same component
        std::queue<NodeId> remaining;
        remaining.push(u);
        while (!remaining.empty()) {
            NodeId current = remaining.front();
            remaining.pop();
            components[c].push_back(current);
            for (NodeId v : neighbours[current]) {
                if (componentOfNode[v] == -1 && getWeight(Edge(current, v)) > 0.0) {
                    componentOfNode[v] = c;
                    remaining.push(v);
                }
            }
            for (NodeId v : revNeighbours[current]) {
                if (componentOfNode[v] == -1 && getWeight(Edge(current, v)) > 0.0) {
                    componentOfNode[v] = c;
                    remaining.push(v);
                }
            }
        }
        std::sort(components[c].begin(), components[c].end());
    }
    
    // check if result correct
    for (NodeId u = 0; u < size; u++) {
        if ( componentOfNode[u] == -1)
            std::cout<<"Component of node "<<u<<" was -1!"<<std::endl;
    }
    uint32_t sum = 0;
    for (uint32_t i = 0; i < components.size(); i++) {
        sum += components[i].size();
    }
    if (sum != size)
        std::cout<<"Sum of components was "<<sum<<" instead of "<<size<<"!"<<std::endl;
    
    return components;
}
