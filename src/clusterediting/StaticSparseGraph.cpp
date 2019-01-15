#include "StaticSparseGraph.h"
#include <x86intrin.h>
#include <bitset>

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
    //neighbours(size, vector<NodeId>(0)),
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
//     std::cout<<"Clear and resize graph with size "<<newSize<<std::endl;
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
//     EdgeId critical = 0;
    
    // iterate over all edges
    for (NodeId i = 0; i < size; i++) {
        for (NodeId j : neighbours[i]) {
            // add edge to weight vector
            EdgeId id = Edge(i,j).id();
//             std::cout<<"Adding edge "<<i<<" "<<j<<" ("<<id<<") ============================================================"<<std::endl;
             
            // insert entry into rank structure
            uint64_t block1 = id / 4096UL;
            uint64_t block2 = (id/64UL) % 64UL;
            uint64_t bitv = rank1[block1] >> (63 - block2);
            
//             if (block1 == 165207) {
//                 if (critical != 0) {
//                     std::cout<<"Critical edge with block1==165207 was "<<id<<std::endl;
//                 }
//                 critical = id;
//             }
            
//             std::cout<<"rank1["<<block1<<"] = \t"<<int2bin(rank1[block1])<<" "<<rank1[block1]<<std::endl;
            
            // check if new block in rank1 has been reached
            if (rank1[block1] == 0UL) {
                offset1[block1] = rank2.size();
            }
            
//             std::cout<<"offset1["<<block1<<"] = \t"<<offset1[block1]<<std::endl;
//             std::cout<<"block2 = \t"<<block2<<std::endl;
            
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
            
//             std::cout<<"rank2["<<block2<<"] = \t"<<int2bin(rank2[block2])<<" "<<rank2[block2]<<std::endl;
//             std::cout<<"offset2["<<block2<<"] = \t"<<offset2[block2]<<std::endl;
//             std::cout<<"block3 = \t"<<block3<<std::endl;
            
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
            
            if(weights[id] != getWeight(Edge(i,j))) {
                std::cout<<"Assertion violated (Get != Set): "<<i<<" "<<j<<" "<<(weights[id])<<" "<<(getWeight(Edge(i,j)))<<std::endl;
            }
            
            refreshEdgeMetaData(Edge(i,j), 0.0, weights[id]);
        }
    }
    
    weightv.shrink_to_fit();
    compiled = true;
}

void StaticSparseGraph::addEdge(const Edge e, const EdgeWeight w) {
    if (w != 0.0) {
        neighbours[e.v].push_back(e.u);
        weights[e.id()] = w;
//      std::cout<<"addEdge "<<e.u<<" "<<e.v<<" "<<e.id()<<" / "<<cliqueOfNode[e.u]<<" "<<cliqueOfNode[e.v]<<std::endl;
        compiled = false;
    }
}

void StaticSparseGraph::addEdge(const NodeId v, const NodeId u, const EdgeWeight w) {
    addEdge(Edge(v, u), w);
}

EdgeWeight StaticSparseGraph::getWeight(const Edge e) {
    
    NodeId cu = cliqueOfNode[e.u];
    NodeId cv = cliqueOfNode[e.v];
    
    // Important: Check clique-information first!
    if (cu == cv) {
        //std::cout<<"cu = "<<cu<<" cv = "<<cv<<std::endl;
        return Permanent;
    } else if (forbidden[cu].size() < forbidden[cv].size() && forbidden[cu].find(cv) != forbidden[cu].end()) {
        return Forbidden;
    } else if (forbidden[cv].find(cu) != forbidden[cv].end()) {
        return Forbidden;
    } else {
        if (compiled)
            return weightv[findIndex(e)];
        else {
            std::map<EdgeId, EdgeWeight>::iterator it = weights.find(e.id());
            if (it == weights.end())
                return 0.0;
            else
                return it->second;
        }
    }
}

void StaticSparseGraph::setPermanent(const Edge e) {
    if (!compiled) {
        compile();
    }
    if (cliqueOfNode[e.u] == cliqueOfNode[e.v]) {
//         std::cout<<"Making permanent edge permanent again ("<<e.u<<", "<<e.v<<")."<<std::endl;
        return;
    }
    
    if (forbidden[cliqueOfNode[e.u]].find(cliqueOfNode[e.v]) != forbidden[cliqueOfNode[e.u]].end()) {
        std::cout<<"Making forbidden edge permanent ("<<e.u<<", "<<e.v<<")."<<std::endl;
        return;
    }
    
    uint64_t rankIndex = findIndex(e);
    if (rankIndex == 0) {
        std::cout<<"Making zero edge permanent explicitly ("<<e.u<<", "<<e.v<<")."<<std::endl;
        refreshEdgeMetaData(e, 0.0, Permanent);
        return;
    } else {
        refreshEdgeMetaData(e, weightv[rankIndex], Permanent);
        weightv[rankIndex] = Permanent;
    }
}

void StaticSparseGraph::setForbidden(const Edge e) {
    if (!compiled) {
        compile();
    }
    NodeId cu = cliqueOfNode[e.u];
    NodeId cv = cliqueOfNode[e.v];
    
    if (cu == cv) {
        std::cout<<"Making forbidden edge permanent ("<<e.u<<", "<<e.v<<")."<<std::endl;
        return;
    }
    
    if (forbidden[cu].find(cv) != forbidden[cu].end()) {
//         std::cout<<"Making forbidden edge forbidden again ("<<e.u<<", "<<e.v<<")."<<std::endl;
        return;
    }
    
    uint64_t rankIndex = findIndex(e);
    if (rankIndex == 0) {
        std::cout<<"Making zero edge forbidden explicitly ("<<e.u<<", "<<e.v<<")."<<std::endl;
        refreshEdgeMetaData(e, 0.0, Forbidden);
        return;
    } else {
        refreshEdgeMetaData(e, weightv[rankIndex], Forbidden);
        weightv[rankIndex] = Forbidden;
    }
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

const std::vector<StaticSparseGraph::NodeId>& StaticSparseGraph::getCliqueOf(const StaticSparseGraph::NodeId v) const {
    return cliques[cliqueOfNode[v]];
}

const std::vector<StaticSparseGraph::NodeId>& StaticSparseGraph::getUnprunedNeighbours(const StaticSparseGraph::NodeId v) const {
    return unprunedNeighbours[v];
}

const std::vector<StaticSparseGraph::NodeId>& StaticSparseGraph::getNonZeroNeighbours(const StaticSparseGraph::NodeId v) const {
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
//         std::cout<<"Merging clusters "<<merged<<" and "<<discarded<<std::endl;
        
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
            std::cout<<"Error 1000"<<std::endl;
        }
    } else if (oldW == Permanent && newW != Permanent) {
        NodeId cu = cliqueOfNode[e.u];
        NodeId cv = cliqueOfNode[e.v];
        // mark cluster pair as forbidden
        forbidden[cu].insert(cv);
        forbidden[cv].insert(cu);
        
        if (cliqueOfNode[e.u] == cliqueOfNode[e.v]) {
            std::cout<<"Error 1001"<<std::endl;
        }
        if (forbidden[cu].find(cv) == forbidden[cu].end()) {
            std::cout<<"Error 1002"<<std::endl;
        }
        if (forbidden[cv].find(cu) == forbidden[cv].end()) {
            std::cout<<"Error 1003"<<std::endl;
        }
    }
    if ((oldW == Forbidden || oldW == Permanent || (oldW == 0.0)) && ((newW != 0.0) && newW != Forbidden && newW != Permanent)) {
        unprunedNeighbours[e.u].push_back(e.v);
        unprunedNeighbours[e.v].push_back(e.u);
    } else if (oldW != Forbidden && oldW != Permanent && (oldW != 0.0) && ((newW == 0.0) || newW == Forbidden || newW == Permanent)) {
        if (!removeFromVector(unprunedNeighbours[e.u], e.v))
            std::cout<<"Error: Non-zero real neighbour not found"<<std::endl;
        if (!removeFromVector(unprunedNeighbours[e.v], e.u))
            std::cout<<"Error: Non-zero real neighbour not found"<<std::endl;
    }
    if (oldW == 0.0 && newW != 0.0) {
        nonzeroNeighbours[e.u].push_back(e.v);
        nonzeroNeighbours[e.v].push_back(e.u);
    } else if (oldW != 0.0 && newW == 0.0) {
        if (!removeFromVector(nonzeroNeighbours[e.u], e.v))
            std::cout<<"Error: Non-zero real neighbour not found"<<std::endl;
        if (!removeFromVector(nonzeroNeighbours[e.v], e.u))
            std::cout<<"Error: Non-zero real neighbour not found"<<std::endl;
    }
}

bool StaticSparseGraph::removeFromVector(std::vector<NodeId> vec, StaticSparseGraph::NodeId v) {
    for (unsigned int i = 0; i < vec.size(); i++) {
        if (vec[i] == v) {
            vec.erase(vec.begin()+i);
            return true;
        }
    }
    return false;
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
//         std::cout<<"findIndex("<<Edge(id).u<<","<<Edge(id).v<<") = "<<0<<" size of weights = "<<weights.size()<<std::endl;
        return 0;
    }
    
    block2 = offset1[block1] + popcount(bitv) - 1;
    u_int64_t block3 = id % 64;
    bitv = rank2[block2] >> (63 - block3);
    
    if ((bitv & 1UL) == 0) {
//         std::cout<<"findIndex("<<Edge(id).u<<","<<Edge(id).v<<") = "<<0<<" size of weights = "<<weights.size()<<std::endl;
        return 0;
    }
    
    block3 = offset2[block2] + popcount(bitv) - 1;
//     std::cout<<"findIndex("<<Edge(id).u<<","<<Edge(id).v<<") = "<<block3<<" size of weights = "<<weights.size()<<std::endl;
    return block3;
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
