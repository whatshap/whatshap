#include "LightCompleteGraph.h"

using namespace std;

namespace ysk {

const LightCompleteGraph::EdgeWeight LightCompleteGraph::Forbidden = -std::numeric_limits< EdgeWeight >::infinity();
const LightCompleteGraph::EdgeWeight LightCompleteGraph::Permanent = std::numeric_limits< EdgeWeight >::infinity();
const LightCompleteGraph::Edge LightCompleteGraph::InvalidEdge = {std::numeric_limits<NodeId>::max(), std::numeric_limits<NodeId>::max()};
const LightCompleteGraph::EdgeId LightCompleteGraph::InvalidEdgeId = -1;
const LightCompleteGraph::NodeId LightCompleteGraph::InvalidNodeId = -1;

LightCompleteGraph::LightCompleteGraph(uint32_t numNodes, bool param_pruneZeroEdges) :
    size(numNodes),
    weights(size*(size-1)/2),
    pruneZeroEdges(param_pruneZeroEdges),
    cliqueOf(size),
    unprunedNeighbours(size, vector<NodeId>(0)) {
}

LightCompleteGraph::LightCompleteGraph(LightCompleteGraph& other) :
    size(other.size),
    weights(other.weights),
    pruneZeroEdges(other.pruneZeroEdges),
    //   origToCompr(other.origToCompr),
    //   comprToOrig(other.comprToOrig),
    cliqueOf(other.cliqueOf),
    unprunedNeighbours(other.unprunedNeighbours)
{
}

void LightCompleteGraph::clearAndResize(const uint32_t newSize) {
    // clear old stuff
    std::fill(weights.begin(), weights.end(), 0.0);
    for (unsigned int i = 0; i < size; i++) {
        cliqueOf[i].clear();
        unprunedNeighbours[i].clear();
    }
    // resize
    size = newSize;
    weights.resize(size*(size-1)/2, 0.0);
    cliqueOf.resize(size, vector<NodeId>(0));
    unprunedNeighbours.resize(size, vector<NodeId>(0));
}

double LightCompleteGraph::getWeight(const Edge e) const {
    return weights[e.id()];
}

void LightCompleteGraph::setWeight(const Edge e, const EdgeWeight w) {
    if (weights[e.id()] != Permanent && w == Permanent) {
        cliqueOf[e.u].push_back(e.v);
        cliqueOf[e.v].push_back(e.u);
    } else if (weights[e.id()] == Permanent && w != Permanent) {
        if (!removeFromVector(cliqueOf[e.u], e.v))
            std::cout<<"Error: Permanent neighbour not found"<<std::endl;
        if (!removeFromVector(cliqueOf[e.v], e.u))
            std::cout<<"Error: Permanent neighbour not found"<<std::endl;
    }
    if ((weights[e.id()] == Forbidden || weights[e.id()] == Permanent || (weights[e.id()] == 0.0 && pruneZeroEdges)) && ((w != 0.0 || !pruneZeroEdges) && w != Forbidden && w != Permanent)) {
        unprunedNeighbours[e.u].push_back(e.v);
        unprunedNeighbours[e.v].push_back(e.u);
    } else if (weights[e.id()] != Forbidden && weights[e.id()] != Permanent && (weights[e.id()] != 0.0 || !pruneZeroEdges) && ((w == 0.0 && pruneZeroEdges) || w == Forbidden || w == Permanent)) {
        if (!removeFromVector(unprunedNeighbours[e.u], e.v))
            std::cout<<"Error: Non-zero real neighbour not found"<<std::endl;
        if (!removeFromVector(unprunedNeighbours[e.v], e.u))
            std::cout<<"Error: Non-zero real neighbour not found"<<std::endl;
    }
    weights[e.id()] = w;
}

void LightCompleteGraph::setWeight(const NodeId v, const NodeId u, const EdgeWeight w){
    Edge e(v,u);
    setWeight(e,w);
}

unsigned int LightCompleteGraph::numNodes() const {
    return size;
}

unsigned long LightCompleteGraph::numEdges() const {
    return weights.size();
}

const std::vector<LightCompleteGraph::NodeId>& LightCompleteGraph::getCliqueOf(const LightCompleteGraph::NodeId v) const {
    return cliqueOf[v];
}

const std::vector<LightCompleteGraph::NodeId>& LightCompleteGraph::getUnprunedNeighbours(const LightCompleteGraph::NodeId v) const {
    return unprunedNeighbours[v];
}

// void LightCompleteGraph::contract(const LightCompleteGraph::Edge e) {
//   NodeId u = e.u;
//   NodeId v = e.v;
//   NodeId last = numNodes() - 1;
//   
//   for (NodeId w = 0; w < numNodes(); w++) {
//     if (w == u || w == v)
//       continue;
//     Edge uw(u, w);
//     Edge vw(v, w);
//     Edge wl(w, last);
//     EdgeWeight weight_V = getWeight(vw);
//     // combine edge weights of u and v
//     setWeight(uw, getWeight(uw) + weight_V);
//     // swap v with last node, by assigning edge weights of last node to v
//     setWeight(wl, weight_V);
//   }
//   
//   // link changes
//   origToCompr[last] = v;
//   origToCompr[v] = u;
//   comprToOrig[u].push_back(v);
//   
//   // delete last internal node
//   size--;
//   weights.resize(size*(size-1)/2);
// }
/*
LightCompleteGraph::NodeId LightCompleteGraph::getInternalId(const LightCompleteGraph::NodeId v) const {
  return origToCompr[v];
}

LightCompleteGraph::vector< uint32_t > LightCompleteGraph::getOriginalIds(const LightCompleteGraph::NodeId v) const {
  return comprToOrig[v];
}*/

bool LightCompleteGraph::removeFromVector(std::vector<NodeId> vec, LightCompleteGraph::NodeId v) {
    for (unsigned int i = 0; i < vec.size(); i++) {
        if (vec[i] == v) {
            vec.erase(vec.begin()+i);
            return true;
        }
    }
    return false;
}


} // namespace ysk
