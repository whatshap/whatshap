#include "StaticSparseGraph.h"
#include <algorithm>

using Edge = DynamicSparseGraph::Edge;
using EdgeWeight = DynamicSparseGraph::EdgeWeight;
using EdgeId = DynamicSparseGraph::EdgeId;
using RankId = DynamicSparseGraph::RankId;
using NodeId = DynamicSparseGraph::NodeId;

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

StaticSparseGraph::StaticSparseGraph(DynamicSparseGraph& other) :
    size(other.numNodes()),
    rank1(std::max(0, ((int32_t)size*((int32_t)size-1)/2 - 1) / 4096 + 1), 0UL),
    offset1(std::max(0, ((int32_t)size*((int32_t)size-1)/2 - 1) / 4096 + 1), 0UL),
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
    compile(other, nodes);
}

StaticSparseGraph::StaticSparseGraph(DynamicSparseGraph& other, std::vector<NodeId>& nodes) :
    size(other.numNodes()),
    rank1((size*(size-1)/2 - 1) / 4096 + 1, 0UL),
    offset1((size*(size-1)/2 - 1) / 4096 + 1, 0UL),
    rank2(0),
    offset2(0),
    weightv(0),
    unprunedNeighbours(size, std::vector<NodeId>(0)),
    nonzeroNeighbours(size, std::vector<NodeId>(0)),
    cliqueOfNode(size, 0),
    cliques(size, std::vector<NodeId>(0)),
    forbidden(size)
{
    for (NodeId i = 0; i < size; i++) {
        cliqueOfNode[i] = i;
        cliques[i].push_back(i);
    }
    compile(other, nodes);
}

void StaticSparseGraph::compile(DynamicSparseGraph& dg, const std::vector<NodeId >& nodes) {
    if (verbosity >= 3) {
        std::cout<<"Compiling graph with "<<nodes.size()<<" nodes ..."<<std::endl;
    }
    weightv.push_back(0.0);
    
    // create mapping from global to local node id
    std::map<NodeId, NodeId> globalToLocal;
    for (uint32_t i = 0; i < nodes.size(); i++) {
        globalToLocal[nodes[i]] = i;
    }
        
    // iterate over all edges
    for (NodeId i : nodes) {
        for (NodeId j : dg.getNeighbours(i)) {
            // check whether (i,j) is in induced graph and whether edge is relevant
            if (globalToLocal.find(j) == globalToLocal.end())
                continue;
            EdgeWeight w = dg.getWeight(Edge(i,j));
            if (w == 0.0)
                continue;
            
            // u and v are indices of this new graph, while i and j are the indices of the source graph
            NodeId u = globalToLocal[i];
            NodeId v = globalToLocal[j];
            
            // add edge to weight vector
            EdgeId id = Edge(u, v).id();
             
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
            
            refreshEdgeMetaData(Edge(u,v), 0.0, w);
            
            EdgeWeight checkWeight = getWeight(Edge(u,v));
            if(w != checkWeight) {
                std::cout<<"Assertion violated (Get != Set): "<<u<<" "<<v<<" "<<w<<" "<<(getWeight(Edge(u,v)))<<std::endl;
            }
            
            if (verbosity >= 4)
                std::cout<<"Compiled edge ("<<u<<","<<v<<") with weight "<<w<<std::endl;
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
    
//     std::cout<<size<<std::endl;
//     for (NodeId i = 0; i < size; i++) {
//         std::cout<<"read_"<<i<<std::endl;
//     }
//     for (NodeId i = 0; i < size; i++) {
//		 for (NodeId j : getNonZeroNeighbours(i)) {
//			 if (i < j) {
//				 std::cout<<i<<" "<<j<<" "<<getWeight(Edge(i,j))<<std::endl;
//			 }
//		 }
//     }
    
    weightv.shrink_to_fit();
}

EdgeWeight StaticSparseGraph::getWeight(const Edge e) {
    RankId r = findIndex(e);
    if (r > 0) {
        return weightv[r];
    } else {
        NodeId cu = cliqueOfNode[e.u];
        NodeId cv = cliqueOfNode[e.v];
        if (cu == cv) {
            return DynamicSparseGraph::Permanent;
        } else if (forbidden[cu].size() * forbidden[cv].size() == 0) {
            return 0.0;
        } else if (forbidden[cu].size() < forbidden[cv].size() && forbidden[cu].find(cv) != forbidden[cu].end()) {
            return DynamicSparseGraph::Forbidden;
        } else if (forbidden[cv].find(cu) != forbidden[cv].end()) {
            return DynamicSparseGraph::Forbidden;
        }
    }
    return 0.0;
}

EdgeWeight StaticSparseGraph::getWeight(const RankId r) {
    return weightv[r];
}

void StaticSparseGraph::setPermanent(const Edge e) {
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
    refreshEdgeMetaData(e, weightv[rankIndex], DynamicSparseGraph::Permanent);
    weightv[rankIndex] = DynamicSparseGraph::Permanent;
}

void StaticSparseGraph::setForbidden(const Edge e) {
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
        return;
    }
    
    refreshEdgeMetaData(e, weightv[rankIndex], DynamicSparseGraph::Forbidden);
    weightv[rankIndex] = DynamicSparseGraph::Forbidden;
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
    if (oldW != DynamicSparseGraph::Permanent && newW == DynamicSparseGraph::Permanent) {
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
    } else if (oldW != DynamicSparseGraph::Forbidden && newW == DynamicSparseGraph::Forbidden) {
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
    if ((oldW == DynamicSparseGraph::Forbidden || oldW == DynamicSparseGraph::Permanent || (oldW == 0.0)) && ((newW != 0.0) && newW != DynamicSparseGraph::Forbidden && newW != DynamicSparseGraph::Permanent)) {
        unprunedNeighbours[e.u].push_back(e.v);
        unprunedNeighbours[e.v].push_back(e.u);
    } else if (oldW != DynamicSparseGraph::Forbidden && oldW != DynamicSparseGraph::Permanent && (oldW != 0.0) && ((newW == 0.0) || newW == DynamicSparseGraph::Forbidden || newW == DynamicSparseGraph::Permanent)) {
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
