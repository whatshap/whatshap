#include "StaticSparseGraph.h"
#include <x86intrin.h>
#include <bitset>
#include <algorithm>

using namespace std;

namespace ysk {

const StaticSparseGraph::EdgeWeight StaticSparseGraph::Forbidden = -std::numeric_limits< EdgeWeight >::infinity();
const StaticSparseGraph::EdgeWeight StaticSparseGraph::Permanent = std::numeric_limits< EdgeWeight >::infinity();
const StaticSparseGraph::Edge StaticSparseGraph::InvalidEdge = {std::numeric_limits<NodeId>::max(), std::numeric_limits<NodeId>::max()};
const StaticSparseGraph::EdgeId StaticSparseGraph::InvalidEdgeId = -1;
const StaticSparseGraph::NodeId StaticSparseGraph::InvalidNodeId = -1;

using Edge = StaticSparseGraph::Edge;
using EdgeWeight = StaticSparseGraph::EdgeWeight;
using EdgeId = StaticSparseGraph::EdgeId;
using RankId = StaticSparseGraph::RankId;
using NodeId = StaticSparseGraph::NodeId;

StaticSparseGraph::StaticSparseGraph(uint32_t numNodes) :
    size(numNodes),
    neighbours(size, vector<NodeId>(0)),
    compiled(false),
    rank1((size*(size-1)/2) / 4096 + 1),
    offset1((size*(size-1)/2) / 4096 + 1),
    rank2(0),
    offset2(0),
    weightv(0),
    cliqueOfNode(size),
    cliques(size, vector<NodeId>(0)),
    forbidden(size),
    unprunedNeighbours(size, vector<NodeId>(0)),
    nonzeroNeighbours(size, vector<NodeId>(0)){
    for (NodeId u = 0; u < size; u++) {
        cliques[u].push_back(u);
        cliqueOfNode[u] = u;
    }
}

StaticSparseGraph::StaticSparseGraph(StaticSparseGraph& other) :
    size(other.size),
    weights(other.weights),
    neighbours(other.neighbours),
    compiled(other.compiled),
    rank1(other.rank1),
    offset1(other.offset1),
    rank2(other.rank2),
    offset2(other.offset2),
    weightv(other.weightv),
    cliqueOfNode(other.cliqueOfNode),
    cliques(other.cliques),
    forbidden(other.forbidden),
    unprunedNeighbours(other.unprunedNeighbours),
    nonzeroNeighbours(other.nonzeroNeighbours)
{
    compile();
}

void StaticSparseGraph::clearAndResize(const uint32_t newSize) {
    // clear old stuff
    for (unsigned int i = 0; i < size; i++) {
        unprunedNeighbours[i].clear();
        nonzeroNeighbours[i].clear();
        cliques[i].clear();
        forbidden[i].clear();
    }
    neighbours.clear();
    unprunedNeighbours.clear();
    nonzeroNeighbours.clear();
    cliqueOfNode.clear();
    cliques.clear();
    forbidden.clear();
    
    // resize
    size = newSize;
    neighbours.clear();
    weights.clear();
    rank1.resize((size*(size-1)/2 - 1) / 4096 + 1, 0UL);
    offset1.resize((size*(size-1)/2 - 1) / 4096 + 1, 0UL);
    rank2.clear();
    offset2.clear();
    weightv.clear();
    neighbours.resize(size, vector<NodeId>(0));
    unprunedNeighbours.resize(size, vector<NodeId>(0));
    nonzeroNeighbours.resize(size, vector<NodeId>(0));
    cliqueOfNode.resize(size, 0);
    cliques.resize(size, vector<NodeId>(0));
    forbidden.resize(size);
    for (NodeId i = 0; i < size; i++) {
        cliqueOfNode[i] = i;
        cliques[i].push_back(i);
    }
    compiled = false;
}

void StaticSparseGraph::compile() {
    if (compiled)
        return;
    if (verbosity >= 3) {
        std::cout<<"Compiling graph with "<<size<<" nodes and "<<weights.size()<<" edges ..."<<std::endl;
    }
    // initialize weight vector with a zero for non-registered edges
    weightv.clear();
    weightv.push_back(0.0);
    rank1.clear();
    rank1.resize((size*(size-1)/2 - 1) / 4096 + 1, 0UL);
    offset1.clear();
    offset1.resize((size*(size-1)/2 - 1) / 4096 + 1, 0UL);
    rank2.clear();
    offset2.clear();
    
    unprunedNeighbours.clear();
    nonzeroNeighbours.clear();
    cliqueOfNode.clear();
    cliques.clear();
    forbidden.clear();
    unprunedNeighbours.resize(size, vector<NodeId>(0));
    nonzeroNeighbours.resize(size, vector<NodeId>(0));
    cliqueOfNode.resize(size, 0);
    cliques.resize(size, vector<NodeId>(0));
    forbidden.resize(size);
    for (NodeId i = 0; i < size; i++) {
        cliqueOfNode[i] = i;
        cliques[i].push_back(i);
    }
    
    compiled = true;
    
    // iterate over all edges
    for (NodeId i = 0; i < size; i++) {
        for (NodeId j : neighbours[i]) {
            // add edge to weight vector
            EdgeId id = Edge(i,j).id();
             
            // insert entry into rank structure
            uint64_t block1 = id / 4096UL;
            uint64_t block2 = (id/64UL) % 64UL;
            uint64_t bitv = rank1[block1] >> (63 - block2);
            
            // check if new block in rank1 has been reached
            if (rank1[block1] == 0UL) {
                offset1[block1] = rank2.size();
            }
            
            // check if rank1 already has a one at this position
            if ((bitv & 1UL) == 0) {
                // set bit in rank1 to one
                rank1[block1] |= (1UL << (63 - block2));
                bitv |= 1UL;
                // create new block for rank2
                rank2.push_back(0UL);
                offset2.push_back(weightv.size());
            }
            
            block2 = offset1[block1] + popcount(bitv) - 1; //only count ones BEFORE current position
            uint64_t block3 = id % 64UL;
            bitv = rank2[block2] >> (63 - block3);
            
            // insert edge
            if((bitv & 1UL) != 0) {
                std::cout<<"Assertion violated (Edge already inserted): "<<i<<" "<<j<<std::endl;
            }
            rank2[block2] |= (1UL << (63 - block3));
            bitv |= 1UL;
            
            if(offset2[block2] + popcount(bitv) - 1 != weightv.size()) {
                std::cout<<"Assertion violated (Weight vector incorrect size): "<<i<<" "<<j<<" "<<(offset2[block2] + popcount(bitv) - 1)<<" "<<(weightv.size())<<std::endl;
            }
            weightv.push_back(weights[id]);
            
            refreshEdgeMetaData(Edge(i,j), 0.0, weights[id]);
            
            EdgeWeight checkWeight = getWeight(Edge(i,j));
            if(weights[id] != checkWeight) {
                std::cout<<"Assertion violated (Get != Set): "<<i<<" "<<j<<" "<<(weights[id])<<" "<<(getWeight(Edge(i,j)))<<std::endl;
            }
            
            if (verbosity >= 4)
                std::cout<<"Compiled edge ("<<i<<","<<j<<") with weight "<<weights[id]<<std::endl;
        }
    }
    
//     std::cout<<size<<std::endl;
//     for (NodeId i = 0; i < size; i++) {
//         std::cout<<"read_"<<i<<std::endl;
//     }
//     for (NodeId i = 0; i < size; i++) {
//         if (i+1 < size)
//             std::cout<<getWeight(Edge(i, i+1));
//         for (NodeId j = i+2; j < size; j++) {
//             std::cout<<" "<<getWeight(Edge(i, j));
//         }
//         std::cout<<std::endl;
//     }
    
    weightv.shrink_to_fit();
    compiled = true;
}

void StaticSparseGraph::addEdge(const Edge e, const EdgeWeight w) {
    if (w != 0.0) {
        neighbours[e.v].push_back(e.u);
        weights[e.id()] = w;
        compiled = false;
    }
}

void StaticSparseGraph::addEdge(const NodeId v, const NodeId u, const EdgeWeight w) {
    addEdge(Edge(v, u), w);
}

EdgeWeight StaticSparseGraph::getWeight(const Edge e) {
    if (!compiled) {
        std::cout<<"Cannot get weight from uncompiled graph!"<<std::endl;
        compile();
    }
    
    RankId r = findIndex(e);
    if (r > 0) {
        return weightv[r];
    } else {
        NodeId cu = cliqueOfNode[e.u];
        NodeId cv = cliqueOfNode[e.v];
        if (cu == cv) {
            return Permanent;
        } else if (forbidden[cu].size() < forbidden[cv].size() && forbidden[cu].find(cv) != forbidden[cu].end()) {
            return Forbidden;
        } else if (forbidden[cv].find(cu) != forbidden[cv].end()) {
            return Forbidden;
        }
    }
    
    return 0.0;
}

EdgeWeight StaticSparseGraph::getWeight(const RankId r) {
    if (!compiled) {
        std::cout<<"Cannot get weight from uncompiled graph!"<<std::endl;
        compile();
    }
    
    return weightv[r];
}

void StaticSparseGraph::setPermanent(const Edge e) {
    if (!compiled) {
        compile();
    }
    
    /* Zero edges will not explicitly set to permanent. Since zero edges are only modified by implication,
     * the meta is being taken care of other edges, so nothing to do here*/
    RankId rankIndex = findIndex(e);
    if (rankIndex == 0) {
        return;
    }
    
    if (forbidden[cliqueOfNode[e.u]].find(cliqueOfNode[e.v]) != forbidden[cliqueOfNode[e.u]].end()) {
        std::cout<<"Making forbidden edge permanent ("<<e.u<<", "<<e.v<<")."<<std::endl;
        return;
    }
    refreshEdgeMetaData(e, weightv[rankIndex], Permanent);
    weightv[rankIndex] = Permanent;
}

void StaticSparseGraph::setForbidden(const Edge e) {
    if (!compiled) {
        compile();
    }
    
    /* Zero edges will not explicitly set to forbidden. Since zero edges are only modified by implication,
     * the meta is being taken care of other edges, so nothing to do here*/
    RankId rankIndex = findIndex(e);
    if (rankIndex == 0) {
        return;
    }
    
    NodeId cu = cliqueOfNode[e.u];
    NodeId cv = cliqueOfNode[e.v];
    
    if (cu == cv) {
        std::cout<<"Making permanent edge forbidden ("<<e.u<<", "<<e.v<<")."<<std::endl;
//         std::exit(EXIT_FAILURE);
        return;
    }
    
    refreshEdgeMetaData(e, weightv[rankIndex], Forbidden);
    weightv[rankIndex] = Forbidden;
}

void StaticSparseGraph::setWeight(const Edge e, const EdgeWeight w) {
    std::map<EdgeId, EdgeWeight>::iterator it = weights.find(e.id());
    if (it == weights.end()) {
        std::cout<<"setWeight tried to modify non existing edge ("<<e.u<<","<<e.v<<") to "<<w<<std::endl;
        return;
    } else {
        if (it->second == Permanent && w != Permanent) {
            std::cout<<"setWeight tried to make permanent edge unpermanent ("<<e.u<<","<<e.v<<") to"<<w<<std::endl;
            return;
        }
        if (it->second == Forbidden && w != Forbidden) {
            std::cout<<"setWeight tried to make forbidden edge unforbidden ("<<e.u<<","<<e.v<<") to"<<w<<std::endl;
            return;
        }
        it->second = w;
    }
    compiled = false;
}

void StaticSparseGraph::setWeight(NodeId v, NodeId u, const EdgeWeight w){
    Edge e(v,u);
    setWeight(e,w);
}

unsigned int StaticSparseGraph::numNodes() const {
    return size;
}

unsigned long StaticSparseGraph::numEdges() const {
    return weightv.size() - 1;
}

const std::vector<NodeId>& StaticSparseGraph::getCliqueOf(const NodeId v) const {
    return cliques[cliqueOfNode[v]];
}

NodeId StaticSparseGraph::getCliqueIdOf(const NodeId v) const {
    return cliqueOfNode[v];
}

const std::vector<NodeId>& StaticSparseGraph::getUnprunedNeighbours(const NodeId v) const {
    return unprunedNeighbours[v];
}

const std::vector<NodeId>& StaticSparseGraph::getNonZeroNeighbours(const NodeId v) const {
    return nonzeroNeighbours[v];
}

void StaticSparseGraph::refreshEdgeMetaData(const Edge e, const EdgeWeight oldW, const EdgeWeight newW) {
    if (oldW != Permanent && newW == Permanent) {
        NodeId merged, discarded; // merge smaller cluster into greater cluster        
        if (cliques[cliqueOfNode[e.u]].size() < cliques[cliqueOfNode[e.v]].size()) {
            merged = cliqueOfNode[e.v];
            discarded = cliqueOfNode[e.u];
        } else {
            merged = cliqueOfNode[e.u];
            discarded = cliqueOfNode[e.v];
        }
        if (merged != discarded) {
            // move nodes from discarded to merged cluster
            for (NodeId d : cliques[discarded]) {
                cliqueOfNode[d] = merged;
                cliques[merged].push_back(d);
            }
            cliques[discarded].clear();
            
            // copy forbidden connections to merged cluster and update references
            for (NodeId f : forbidden[discarded]) {
                forbidden[merged].insert(f);
                forbidden[f].insert(merged);
                forbidden[f].erase(discarded);
            }
            forbidden[discarded].clear();
            
            if (cliqueOfNode[e.u] != cliqueOfNode[e.v]) {
                std::cout<<"Error 1000 "<<cliqueOfNode[e.u]<<" != "<<cliqueOfNode[e.v]<<std::endl;
            }
        }
    } else if (oldW != Forbidden && newW == Forbidden) {
        NodeId cu = cliqueOfNode[e.u];
        NodeId cv = cliqueOfNode[e.v];
        if (cu != cv) {
            // mark cluster pair as forbidden
            forbidden[cu].insert(cv);
            forbidden[cv].insert(cu);
            
            if (forbidden[cu].find(cv) == forbidden[cu].end()) {
                std::cout<<"Error 1002"<<std::endl;
            }
            if (forbidden[cv].find(cu) == forbidden[cv].end()) {
                std::cout<<"Error 1003"<<std::endl;
            }
        }
    }
    if ((oldW == Forbidden || oldW == Permanent || (oldW == 0.0)) && ((newW != 0.0) && newW != Forbidden && newW != Permanent)) {
        unprunedNeighbours[e.u].push_back(e.v);
        unprunedNeighbours[e.v].push_back(e.u);
    } else if (oldW != Forbidden && oldW != Permanent && (oldW != 0.0) && ((newW == 0.0) || newW == Forbidden || newW == Permanent)) {
        if (!removeFromVector(unprunedNeighbours[e.u], e.v))
            std::cout<<"Error: Non-zero real neighbour "<<e.v<<" of "<<e.u<<" not found. Weight was set from "<<oldW<<" to "<<newW<<std::endl;
        if (!removeFromVector(unprunedNeighbours[e.v], e.u))
            std::cout<<"Error: Non-zero real neighbour "<<e.u<<" of "<<e.v<<" not found. Weight was set from "<<oldW<<" to "<<newW<<std::endl;
        if (std::find(unprunedNeighbours[e.u].begin(), unprunedNeighbours[e.u].end(), e.v) != unprunedNeighbours[e.u].end())
            std::cout<<"Removed unpruned neighbour "<<e.v<<" from "<<e.u<<", but is still there"<<std::endl;
        if (std::find(unprunedNeighbours[e.v].begin(), unprunedNeighbours[e.v].end(), e.u) != unprunedNeighbours[e.v].end())
            std::cout<<"Removed unpruned neighbour "<<e.u<<" from "<<e.v<<", but is still there"<<std::endl;
    }
    if (oldW == 0.0 && newW != 0.0) {
        nonzeroNeighbours[e.u].push_back(e.v);
        nonzeroNeighbours[e.v].push_back(e.u);
    } else if (oldW != 0.0 && newW == 0.0) {
        if (!removeFromVector(nonzeroNeighbours[e.u], e.v))
            std::cout<<"Error: Non-zero neighbour "<<e.v<<" of "<<e.u<<" not found. Weight was set from "<<oldW<<" to "<<newW<<std::endl;
        if (!removeFromVector(nonzeroNeighbours[e.v], e.u))
            std::cout<<"Error: Non-zero neighbour "<<e.u<<" of "<<e.v<<" not found. Weight was set from "<<oldW<<" to "<<newW<<std::endl;
    }
}

bool StaticSparseGraph::removeFromVector(std::vector<NodeId>& vec, StaticSparseGraph::NodeId v) {
    bool found = false;
    for (unsigned int i = 0; i < vec.size(); i++) {
        if (vec[i] == v) {
            vec[i] = vec.back();
            vec.pop_back();
            found = true;
            break;
        }
    }
    return found;
}

RankId StaticSparseGraph::findIndex(const Edge e) const {
    return findIndex(e.id());
}

RankId StaticSparseGraph::findIndex(const EdgeId id) const {
    u_int64_t block1 = id / 4096;
    u_int64_t block2 = (id/64) % 64;
    u_int64_t bitv = rank1[block1] >> (63 - block2);
    
    // check if corresponding bit in rank block is unset
    if ((bitv & 1UL) == 0) {
        return 0;
    }
    
    block2 = offset1[block1] + popcount(bitv) - 1;
    u_int64_t block3 = id % 64;
    bitv = rank2[block2] >> (63 - block3);
    
    if ((bitv & 1UL) == 0) {
        return 0;
    }
    
    return offset2[block2] + popcount(bitv) - 1;
}

uint64_t StaticSparseGraph::popcount(uint64_t bitv) const {
        // copied from Wikipedia (https://en.wikipedia.org/wiki/Hamming_weight)
        bitv -= (bitv >> 1) & m1;
        bitv = (bitv & m2) + ((bitv >> 2) & m2);
        bitv = (bitv + (bitv >> 4)) & m4;
        return (bitv * h01) >> 56;
        // return _mm_popcnt_u64(bitv)
    }

} // namespace ysk
