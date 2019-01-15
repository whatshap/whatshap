#include "EdgeHeap.h"
#include <cmath>
#include <algorithm>

namespace ysk {
  
using Edge = StaticSparseGraph::Edge;
using EdgeWeight = StaticSparseGraph::EdgeWeight;
using EdgeId = StaticSparseGraph::EdgeId;
using RankId = StaticSparseGraph::RankId;
using NodeId = StaticSparseGraph::NodeId;

EdgeHeap::EdgeHeap(StaticSparseGraph& param_graph, bool param_pruneZeroEdges) :
    pruneZeroEdges(param_pruneZeroEdges),
    graph(param_graph),
    unprocessed(0),
    edges(1+param_graph.numEdges(), StaticSparseGraph::InvalidEdge),
    icf(1+param_graph.numEdges(), StaticSparseGraph::Forbidden),
    icp(1+param_graph.numEdges(), StaticSparseGraph::Forbidden),
    edge2forb_rank(2+2*param_graph.numEdges(), 0),
    edge2perm_rank(2+2*param_graph.numEdges(), 0)
{
    //initInducedCosts();
}

void EdgeHeap::initInducedCosts() {
    if (verbosity >= 1)
        std::cout<<"Compute induced costs.."<<std::endl;
        
    // compute array: edge -> icf/icp
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        if (verbosity >= 1 && u % 100 == 0)
            std::cout<<"Completed "<<(((2UL*graph.numNodes()-(uint64_t)u+1UL)*(uint64_t)u*50UL)/(((uint64_t)graph.numNodes()*((uint64_t)graph.numNodes()-1UL))/2UL))<<"%\r"<<std::flush;
        for (NodeId v : graph.getNonZeroNeighbours(u)) {
            if (v < u)
                continue;
            
            // iterate over all edges uv
            Edge uv(u,v);
            EdgeId id = uv.id();
            RankId rId = graph.findIndex(id);
            
            // Zero edges have no icp/icf
            if (graph.findIndex(id) == 0) {
                continue;
            } else {
                edges[rId] = uv;
            }
            
            EdgeWeight w_uv = graph.getWeight(uv);
            if (w_uv == 0.0)
                std::cout<<"Zero edge has reached icf/icp computation."<<std::endl;

            if (w_uv == 0.0 || w_uv == StaticSparseGraph::Forbidden || w_uv == StaticSparseGraph::Permanent) {
                continue;
            } else {
                icf[rId] = 0.0;
                icp[rId] = 0.0;
                unprocessed++;
            }
            
            // costs for the edge uv itself
            if (w_uv >= 0) {	
                icf[rId] += w_uv;	// costs for removing uv
            } else {
                icp[rId] += -w_uv;	// costs for adding uv
            }
            
            // look at all triangles uvw containing uv. Triangles with a zero edge can be ignored
            std::vector<NodeId> w_vec;
            std::set_intersection(graph.getNonZeroNeighbours(u).begin(), graph.getNonZeroNeighbours(u).end(), 
                                  graph.getNonZeroNeighbours(v).begin(), graph.getNonZeroNeighbours(v).end(), back_inserter(w_vec));

            for (NodeId w : w_vec) {
                Edge uw(u,w);
                Edge vw(v,w);
                EdgeWeight w_uw = graph.getWeight(uw);
                EdgeWeight w_vw = graph.getWeight(vw);
                icf[rId] += getIcf(w_uw, w_vw);
                icp[rId] += getIcp(w_uw, w_vw);
            }
//             std::cout<<"Rank "<<rId<<" = ("<<u<<","<<v<<")"<<" icp="<<icp[rId]<<" icf="<<icf[rId]<<std::endl;
        }
    }
    
    for (unsigned int i = 0; i < icf.size(); i++){
        if(std::isnan(icf[graph.findIndex(i)])) {
            std::cout<<"NaN! in icf"<<std::endl;
            break;
        }
        if(std::isnan(icp[graph.findIndex(i)])) {
            std::cout<<"NaN! in icp"<<std::endl;
            break;
        }
    }
    
    // sort edges by icf and icp values
    for (RankId id = 0; id < icf.size(); id++) {
        forb_rank2edge.push_back(id);
        perm_rank2edge.push_back(id);
    }
    
    std::sort(forb_rank2edge.begin(), forb_rank2edge.end(), [this] (const EdgeId& a, const EdgeId& b) { return icf[a] > icf[b]; });
    std::sort(perm_rank2edge.begin(), perm_rank2edge.end(), [this] (const EdgeId& a, const EdgeId& b) { return icp[a] > icp[b]; });
    
//     NodeId over = graph.numNodes()*(graph.numNodes()-1)/2+1;
//     for (EdgeId i = 0; i < icf.size(); i++) {
//         icf.push_back(StaticSparseGraph::Forbidden);
//         icp.push_back(StaticSparseGraph::Forbidden);
//         forb_rank2edge.push_back(over+i);
//         perm_rank2edge.push_back(over+i);
//     }
    
    // save index in sorted vectors for each edge
    for (RankId i = 0; i < icf.size(); i++) {
        edge2forb_rank[forb_rank2edge[i]] = i;
        edge2perm_rank[perm_rank2edge[i]] = i;
    }
    
//     std::cout<<"icf = ";
//     for (const auto& i: icf)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
//     std::cout<<"forb_rank2edge = ";
//     for (const auto& i: forb_rank2edge)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
//     std::cout<<"icp = ";
//     for (const auto& i: icp)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
//     std::cout<<"perm_rank2edge = ";
//     for (const auto& i: perm_rank2edge)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
}

Edge EdgeHeap::getMaxIcfEdge() const {
    RankId ei = forb_rank2edge[0];
    if (forb_rank2edge.size() <= 1) {
        // only rank 0 entry left
//         std::cout<<"only 1 edge left on icf heap ("<<edges[ei].u<<","<<edges[ei].v<<")"<<std::endl;
        return StaticSparseGraph::InvalidEdge;
    }
    if (icf[ei] < 0) {
//         std::cout<<"negative icf costs on top of icf heap ("<<edges[ei].u<<","<<edges[ei].v<<") ("<<ei<<") : "<<icf[ei]<<std::endl;
        return StaticSparseGraph::InvalidEdge;
    }
//     std::cout<<"Max icf edge = ("<<ei<<") = ("<<edges[ei].u<<","<<edges[ei].v<<") weight ("<<icf[ei]<<")"<<std::endl;
    return edges[ei];
//     NodeId u = std::ceil(std::sqrt(2*(ei+1)+0.25) - 0.5);
//     NodeId v = ei - u * (u-1) / 2;
//     return Edge(u, v);
}

Edge EdgeHeap::getMaxIcpEdge() const {
    RankId ei = perm_rank2edge[0];
    if (perm_rank2edge.size() <= 1) {
        // only rank 0 entry left
//         std::cout<<"only 1 edge left on icf heap ("<<edges[ei].u<<","<<edges[ei].v<<")"<<std::endl;
        return StaticSparseGraph::InvalidEdge;
    }
    if (icp[ei] < 0) {
//         std::cout<<"negative icf costs on top of icp heap ("<<edges[ei].u<<","<<edges[ei].v<<") : "<<icp[ei]<<std::endl;
        return StaticSparseGraph::InvalidEdge;
    }
//     std::cout<<"Max icp edge = ("<<ei<<") = ("<<edges[ei].u<<","<<edges[ei].v<<") weight ("<<icp[ei]<<")"<<std::endl;
    return edges[ei];
//     NodeId u = std::ceil(std::sqrt(2*(ei+1)+0.25) - 0.5);
//     NodeId v = ei - u * (u-1) / 2;
//     return Edge(u, v);
}

EdgeWeight EdgeHeap::getIcf(const Edge e) const {
    return icf[graph.findIndex(e)];
}

EdgeWeight EdgeHeap::getIcp(const Edge e) const {
    return icp[graph.findIndex(e)];
}

void EdgeHeap::increaseIcf(const Edge e, const EdgeWeight w) {
    RankId rId = graph.findIndex(e);
    if (rId == 0)
        std::cout<<"increaseIcf called on zero edge"<<std::endl;
    if (rId > 0 && w != 0 && icf[rId] >= 0) {
//         std::cout<<"Increase icf on "<<rId<<" ("<<e.u<<","<<e.v<<") from "<<icf[rId]<<" by "<<w<<std::endl;
        icf[rId] += w;
        icf[rId] = std::max(icf[rId], 0.0);
        updateHeap(forb_rank2edge, rId, w, edge2forb_rank, icf);
    }
}

void EdgeHeap::increaseIcp(const Edge e, const EdgeWeight w) {
    RankId rId = graph.findIndex(e);
    if (rId == 0)
        std::cout<<"increaseIcp called on zero edge"<<std::endl;
    if (rId > 0 && w != 0 && icp[rId] >= 0) {
//         std::cout<<"Increase icp on "<<rId<<" ("<<e.u<<","<<e.v<<") from "<<icp[rId]<<" by "<<w<<std::endl;
        icp[rId] += w;
        icp[rId] = std::max(icp[rId], 0.0);
        updateHeap(perm_rank2edge, rId, w, edge2perm_rank, icp);
    }
}

void EdgeHeap::removeEdge(const Edge e) {
    RankId rId = graph.findIndex(e);
    if (rId == 0)
        std::cout<<"removeEdge called on zero edge"<<std::endl;
//     else
//         std::cout<<"Removing edge ("<<e.u<<","<<e.v<<") from heap ("<<rId<<")"<<std::endl;
    if (rId > 0 && icf[rId] != StaticSparseGraph::Forbidden && icp[rId] != StaticSparseGraph::Forbidden) {
        icf[rId] = StaticSparseGraph::Forbidden;
        icp[rId] = StaticSparseGraph::Forbidden;
        updateHeap(forb_rank2edge, rId, StaticSparseGraph::Forbidden, edge2forb_rank, icf);
        updateHeap(perm_rank2edge, rId, StaticSparseGraph::Forbidden, edge2perm_rank, icp);
        unprocessed--;
    }
    
//     std::cout<<"icf = ";
//     for (const auto& i: icf)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
//     std::cout<<"forb_rank2edge = ";
//     for (const auto& i: forb_rank2edge)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
//     std::cout<<"icp = ";
//     for (const auto& i: icp)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
//     std::cout<<"perm_rank2edge = ";
//     for (const auto& i: perm_rank2edge)
//         std::cout << i << ' ';
//     std::cout<<std::endl;
}

StaticSparseGraph::EdgeWeight EdgeHeap::getIcf(const EdgeWeight uw, const EdgeWeight vw) const {
    if (uw > 0 && vw > 0) {
        // if both other edges present, remove the cheapest of both
        return std::min(uw, vw); 
    } else {
        return 0;
    }
}

StaticSparseGraph::EdgeWeight EdgeHeap::getIcp(const EdgeWeight uw, const EdgeWeight vw) const {
    if (uw < 0 && vw > 0) {
        return std::min(vw, -uw); 	// either add uw or remove vw
    } else if (uw > 0 && vw < 0) {
        return std::min(-vw, uw); 	// either add vw or remove uw
    } else {
        return 0;
    }
}

int EdgeHeap::numUnprocessed() const {
    return unprocessed;
}

void EdgeHeap::updateHeap(std::vector<RankId>& heap, const RankId e, const EdgeWeight change, std::vector<RankId>& index, const std::vector<EdgeWeight>& score) {
//     std::cout<<"   Score of "<<e<<": "<<score[e]<<" -> "<<(score[e]+change)<<std::endl;
    unsigned int pos = index[e];
    /*
     * index arithemetic for zero based array: parent = (index-1)/2, children = 2*index+1 and 2*index+2
     */
    if (change > 0) {
        // value increased -> move edge upwards in heap
        while(pos > 0 && score[heap[(pos-1)/2]] < score[heap[pos]]) {
            // swap pos and pos/2
            std::swap(heap[pos], heap[(pos-1)/2]);
            index[heap[pos]] = pos;
            index[heap[(pos-1)/2]] = (pos-1)/2;
//             std::cout<<"   Move "<<e<<" upwards in heap ("<<pos<<" -> "<<((pos-1)/2)<<")"<<std::endl;
            pos = (pos-1)/2;
        }
    } else {
        // value decreased -> move edge downwards in heap
        while((2*pos+1 < heap.size() && score[heap[pos]] < score[heap[2*pos+1]])
            || (2*pos+2 < heap.size() && score[heap[pos]] < score[heap[2*pos+2]]) ) {
            if (2*pos+2 < heap.size() && score[heap[2*pos+1]] < score[heap[2*pos+2]]) {
                // element 2*pos+2 exists and is larger than 2*pos+1 -> swap pos with 2*pos+2
                std::swap(heap[pos], heap[2*pos+2]);
                index[heap[pos]] = pos;
                index[heap[pos*2+2]] = pos*2+2;
//                 std::cout<<"   Move "<<e<<" downwards in heap ("<<pos<<" -> "<<(2*pos+2)<<")"<<std::endl;
                pos = 2*pos+2;
            } else {
                // else swap with 2*pos+1
                std::swap(heap[pos], heap[2*pos+1]);
                index[heap[pos]] = pos;
                index[heap[pos*2+1]] = pos*2+1;
//                 std::cout<<"   Move "<<e<<" downwards in heap ("<<pos<<" -> "<<(2*pos+1)<<")"<<std::endl;
                pos = 2*pos+1;
            }
        }
    }
}

} // namespace ysk

