#include "staticsparsegraph.h"
#include <algorithm>

const StaticSparseGraph::EdgeWeight StaticSparseGraph::Forbidden = -std::numeric_limits<EdgeWeight>::infinity();
const StaticSparseGraph::EdgeWeight StaticSparseGraph::Permanent = std::numeric_limits<EdgeWeight>::infinity();
const StaticSparseGraph::Edge StaticSparseGraph::InvalidEdge = {std::numeric_limits<NodeId>::max(), std::numeric_limits<NodeId>::max()};
const StaticSparseGraph::EdgeId StaticSparseGraph::InvalidEdgeId = -1;
const StaticSparseGraph::NodeId StaticSparseGraph::InvalidNodeId = -1;

using Edge = StaticSparseGraph::Edge;
using EdgeWeight = StaticSparseGraph::EdgeWeight;
using EdgeId = StaticSparseGraph::EdgeId;
using RankId = StaticSparseGraph::RankId;
using NodeId = StaticSparseGraph::NodeId;

StaticSparseGraph::StaticSparseGraph(StaticSparseGraph& other) :
    size(other.size),
    rank1(other.rank1),
    offset1(other.offset1),
    rank2(other.rank2),
    offset2(other.offset2),
    weightv(other.weightv),
    unprunedNeighbours(other.unprunedNeighbours),
    nonzeroNeighbours(other.nonzeroNeighbours),
    cliqueOfNode(other.cliqueOfNode),
    cliques(other.cliques),
    forbidden(other.forbidden) {}

StaticSparseGraph::StaticSparseGraph(TriangleSparseMatrix& m) :
    size((uint64_t)m.getMaxDim()),
    rank1(std::max((int64_t)0, (int64_t)(size*(size-1)/2 - 1) / 4096 + 1), 0UL),
    offset1(std::max((int64_t)0, (int64_t)(size*(size-1)/2 - 1) / 4096 + 1), 0UL),
    rank2(0),
    offset2(0),
    weightv(0),
    unprunedNeighbours(size, std::vector<NodeId>(0)),
    nonzeroNeighbours(size, std::vector<NodeId>(0)),
    cliqueOfNode(size, 0),
    cliques(size, std::vector<NodeId>(0)),
    forbidden(size)
{
    std::vector<NodeId> nodes(size, 0);
    for (NodeId i = 0; i < size; i++) {
        nodes[i] = i;
        cliqueOfNode[i] = i;
        cliques[i].push_back(i);
    }
    compile(m);
}

void StaticSparseGraph::compile(TriangleSparseMatrix& m) {
    weightv.push_back(0.0);
    
    // iterate over all sorted edges
    for (EdgeId id : m.getIndices()) {
        NodeId u = std::ceil(std::sqrt(2*(id+1)+0.25) - 0.5);
        NodeId v = (NodeId)(id - (uint64_t)u * (uint64_t)(u-1) / 2);
        EdgeWeight w = m.get(u, v);
        
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
            std::cout<<"Assertion violated (Edge already inserted): "<<u<<" "<<v<<std::endl;
        }
        rank2[block2] |= (1UL << (63 - block3));
        bitv |= 1UL;
        
        if(offset2[block2] + popcount(bitv) - 1 != weightv.size()) {
            std::cout<<"Assertion violated (Weight vector incorrect size): "<<u<<" "<<v<<" "<<(offset2[block2] + popcount(bitv) - 1)<<" "<<(weightv.size())<<std::endl;
        }
        weightv.push_back(w);
        if (w == StaticSparseGraph::Forbidden)
            setForbidden(Edge(u, v), weightv.size()-1);
        else if (w == StaticSparseGraph::Permanent)
            setPermanent(Edge(u, v), weightv.size()-1);
        
        refreshEdgeMetaData(Edge(u,v), 0.0, w);
        
        EdgeWeight checkWeight = getWeight(Edge(u,v));
        if(w != checkWeight) {
            std::cout<<"Assertion violated (Get != Set): "<<u<<" "<<v<<" "<<w<<" "<<(getWeight(Edge(u,v)))<<std::endl;
        }
    }
    
    weightv.shrink_to_fit();
}

EdgeWeight StaticSparseGraph::getWeight(const Edge e) {
    return weightv[findIndex(e)];
}

EdgeWeight StaticSparseGraph::getWeight(const RankId r) {
    return weightv[r];
}

bool StaticSparseGraph::isPermanent(const Edge e) {
    NodeId cu = cliqueOfNode[e.u];
    NodeId cv = cliqueOfNode[e.v];
    return cu == cv;
}

bool StaticSparseGraph::isForbidden(const Edge e) {
    NodeId cu = cliqueOfNode[e.u];
    NodeId cv = cliqueOfNode[e.v];
    return forbidden[cu].find(cv) != forbidden[cu].end();
}

void StaticSparseGraph::setPermanent(const Edge e) {
    /* Zero edges will not explicitly set to permanent. Since zero edges are only modified by implication,
     * the meta is being taken care of other edges, so nothing to do here*/
    RankId rankIndex = findIndex(e);
    if (rankIndex > 0) {
        setPermanent(e, rankIndex);
    }
}

void StaticSparseGraph::setPermanent(const Edge e, RankId r) {
    if (forbidden[cliqueOfNode[e.u]].find(cliqueOfNode[e.v]) != forbidden[cliqueOfNode[e.u]].end()) {
        std::cout<<"Making forbidden edge permanent ("<<e.u<<", "<<e.v<<")."<<std::endl;
        return;
    }
    
    NodeId cu = cliqueOfNode[e.u];
    NodeId cv = cliqueOfNode[e.v];
    NodeId merged, discarded; // merge smaller cluster into greater cluster        
    if (cliques[cu].size() < cliques[cv].size()) {
        merged = cv;
        discarded = cu;
    } else {
        merged = cu;
        discarded = cv;
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
    refreshEdgeMetaData(e, weightv[r], StaticSparseGraph::Permanent);
    if (r > 0)
        weightv[r] = StaticSparseGraph::Permanent;
}

void StaticSparseGraph::setForbidden(const Edge e) {
    /* Zero edges will not explicitly set to forbidden. Since zero edges are only modified by implication,
     * the meta is being taken care of other edges, so nothing to do here*/
    RankId rankIndex = findIndex(e);
    if (rankIndex > 0) {
        setForbidden(e, rankIndex);
    }
}

void StaticSparseGraph::setForbidden(const Edge e, RankId r) {
    NodeId cu = cliqueOfNode[e.u];
    NodeId cv = cliqueOfNode[e.v];
    
    if (cu == cv) {
        std::cout<<"Making permanent edge forbidden ("<<e.u<<", "<<e.v<<")."<<std::endl;
        return;
    } else {
        // mark cluster pair as forbidden
        forbidden[cu].insert(cv);
        forbidden[cv].insert(cu);
    }
    refreshEdgeMetaData(e, weightv[r], StaticSparseGraph::Forbidden);
    if (r > 0)
        weightv[r] = StaticSparseGraph::Forbidden;
}

uint64_t StaticSparseGraph::numNodes() const {
    return size;
}

uint64_t StaticSparseGraph::numEdges() const {
    return weightv.size() - 1;
}

const std::vector<NodeId>& StaticSparseGraph::getCliqueOf(const NodeId v) const {
    return cliques[cliqueOfNode[v]];
}

const std::vector<NodeId> StaticSparseGraph::getForbiddenNeighbors(const NodeId v) const {
    std::vector<NodeId> f;
    NodeId cv = cliqueOfNode[v];
    for (NodeId forbiddenClique : forbidden[cv]) {
        for (NodeId forbiddenNode : cliques[forbiddenClique]) {
            f.push_back(forbiddenNode);
        }
    }
    return f;
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
    bool oldPruned = (oldW == StaticSparseGraph::Forbidden) | (oldW == StaticSparseGraph::Permanent) | (oldW == 0.0);
    bool newPruned = (newW == StaticSparseGraph::Forbidden) | (newW == StaticSparseGraph::Permanent) | (newW == 0.0);
    bool oldZero = (oldW == 0.0);
    bool newZero = (newW == 0.0);
    if (oldPruned & !newPruned) {
        unprunedNeighbours[e.u].push_back(e.v);
        unprunedNeighbours[e.v].push_back(e.u);
    } else if (!oldPruned & newPruned) {
        removeFromVector(unprunedNeighbours[e.u], e.v);
        removeFromVector(unprunedNeighbours[e.v], e.u);
    }
    if (oldZero & !newZero) {
        nonzeroNeighbours[e.u].push_back(e.v);
        nonzeroNeighbours[e.v].push_back(e.u);
    } else if (!oldZero & newZero) {
        removeFromVector(nonzeroNeighbours[e.u], e.v);
        removeFromVector(nonzeroNeighbours[e.v], e.u);
    }
}

bool StaticSparseGraph::removeFromVector(std::vector<NodeId>& vec, NodeId v) {
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
    u_int64_t block2 = (id / 64) % 64;
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
