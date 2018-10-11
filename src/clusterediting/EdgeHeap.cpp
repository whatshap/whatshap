#include "EdgeHeap.h"
#include <cmath>
#include <algorithm>

namespace ysk {
  
using Edge = LightCompleteGraph::Edge;
using EdgeWeight = LightCompleteGraph::EdgeWeight;
using EdgeId = LightCompleteGraph::EdgeId;
using NodeId = LightCompleteGraph::NodeId;

EdgeHeap::EdgeHeap(LightCompleteGraph& param_graph, bool param_pruneZeroEdges) :
    pruneZeroEdges(param_pruneZeroEdges),
    graph(param_graph),
    unprocessed(0),
    icf(param_graph.numEdges(), 0.0),
    icp(param_graph.numEdges(), 0.0),
    edge2forb_rank(2*param_graph.numEdges(), 0),
    edge2perm_rank(2*param_graph.numEdges(), 0)
{
    initInducedCosts();
}

void EdgeHeap::initInducedCosts() {
    if (verbosity >= 1)
        std::cout<<"Compute induced cost."<<std::endl;
    // preprocessing: sorted vector of non-zero neighbours for each node
    std::vector<std::vector<NodeId>> nonZeroNeighbours;
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        std::vector<NodeId> n;
        for (NodeId v = 0; v < graph.numNodes(); v++) {
        if (u != v && graph.getWeight(Edge(u,v)) != 0.0) {
        n.push_back(v);
        }
        }
        nonZeroNeighbours.push_back(n);
    }
    
    // compute array: edge -> icf/icp
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        if (verbosity >= 1)
            std::cout<<"Completed "<<(((2*graph.numNodes()-u)*(u+1)/2)*100/graph.numEdges())<<"%\r"<<std::flush;
        for (NodeId v = u + 1; v < graph.numNodes(); v++) {
            // iterate over all edges uv
            Edge uv(u,v);
            EdgeWeight w_uv = graph.getWeight(uv);
            EdgeId id = uv.id();
            
            if (w_uv == 0.0 || w_uv == LightCompleteGraph::Forbidden || w_uv == LightCompleteGraph::Permanent) {
                icf[id] = LightCompleteGraph::Forbidden;
                icp[id] = LightCompleteGraph::Forbidden;
                continue;
            } else {
                unprocessed++;
            }
            
            // costs for the edge uv itself
            if (w_uv >= 0) {	
                icf[id] += w_uv;	// costs for removing uv
            } else {
                icp[id] += -w_uv;	// costs for adding uv
            }
            
            // look at all triangles uvw containing uv. Triangles with a zero edge can be ignored
            std::vector<NodeId> w_vec;
            std::set_intersection(nonZeroNeighbours[u].begin(), nonZeroNeighbours[u].end(), nonZeroNeighbours[v].begin(), nonZeroNeighbours[v].end(), back_inserter(w_vec));

            for (NodeId w : w_vec) {
                Edge uw(u,w);
                Edge vw(v,w);
                EdgeWeight w_uw = graph.getWeight(uw);
                EdgeWeight w_vw = graph.getWeight(vw);
                icf[id] += getIcf(w_uw, w_vw);
                icp[id] += getIcp(w_uw, w_vw);
            }
        }
    }
    
    for (unsigned int i = 0; i < icf.size(); i++){
        if(std::isnan(icf[i])) {
        std::cout<<"NaN! in icf"<<std::endl;
        break;
        }
        if(std::isnan(icp[i])) {
        std::cout<<"NaN! in icp"<<std::endl;
        break;
        }
    }
    
    // sort edges by icf and icp values
    for (EdgeId i = 0; i < graph.numEdges(); i++) {
        forb_rank2edge.push_back(i);
        perm_rank2edge.push_back(i);
    }
    
    std::sort(forb_rank2edge.begin(), forb_rank2edge.end(), [this] (const EdgeId& a, const EdgeId& b) { return icf[a] > icf[b]; });
    std::sort(perm_rank2edge.begin(), perm_rank2edge.end(), [this] (const EdgeId& a, const EdgeId& b) { return icp[a] > icp[b]; });
    
    for (EdgeId i = 0; i < graph.numEdges(); i++) {
        icf.push_back(LightCompleteGraph::Forbidden);
        icp.push_back(LightCompleteGraph::Forbidden);
        forb_rank2edge.push_back(graph.numEdges()+i);
        perm_rank2edge.push_back(graph.numEdges()+i);
    }
    
    // save index in sorted vectors for each edge
    for (EdgeId i = 0; i < 2*graph.numEdges(); i++) {
        edge2forb_rank[forb_rank2edge[i]] = i;
        edge2perm_rank[perm_rank2edge[i]] = i;
    }
}

Edge EdgeHeap::getMaxIcfEdge() const {
    if (forb_rank2edge.size() <= 0) {
        return LightCompleteGraph::InvalidEdge;
    }
    EdgeId ei = forb_rank2edge[0];
    if (icf[ei] < 0) {
        return LightCompleteGraph::InvalidEdge;
    }
    NodeId u = std::ceil(std::sqrt(2*(ei+1)+0.25) - 0.5);
    NodeId v = ei - u * (u-1) / 2;
    return Edge(u, v);
}

Edge EdgeHeap::getMaxIcpEdge() const {
    if (perm_rank2edge.size() <= 0) {
        return LightCompleteGraph::InvalidEdge;
    }
    EdgeId ei = perm_rank2edge[0];
    if (icp[ei] < 0) {
        return LightCompleteGraph::InvalidEdge;
    }
    NodeId u = std::ceil(std::sqrt(2*(ei+1)+0.25) - 0.5);
    NodeId v = ei - u * (u-1) / 2;
    return Edge(u, v);
}

EdgeWeight EdgeHeap::getIcf(const Edge e) const {
    return icf[e.id()];
}

EdgeWeight EdgeHeap::getIcp(const Edge e) const {
    return icp[e.id()];
}

void EdgeHeap::increaseIcf(const Edge e, const EdgeWeight w) {
    EdgeId id = e.id();
    if (w != 0 && icf[id] >= 0) {
        icf[id] += w;
        icf[id] = std::max(icf[id], 0.0);
        updateHeap(forb_rank2edge, id, w, edge2forb_rank, icf);
    }
}

void EdgeHeap::increaseIcp(const Edge e, const EdgeWeight w) {
    EdgeId id = e.id();
    if (w != 0 && icp[id] >= 0) {
        icp[id] += w;
        icp[id] = std::max(icp[id], 0.0);
        updateHeap(perm_rank2edge, id, w, edge2perm_rank, icp);
    }
}

void EdgeHeap::removeEdge(const Edge e) {
    EdgeId id = e.id();
    if (icf[id] != LightCompleteGraph::Forbidden && icp[id] != LightCompleteGraph::Forbidden) {
        icf[id] = LightCompleteGraph::Forbidden;
        icp[id] = LightCompleteGraph::Forbidden;
        updateHeap(forb_rank2edge, id, LightCompleteGraph::Forbidden, edge2forb_rank, icf);
        updateHeap(perm_rank2edge, id, LightCompleteGraph::Forbidden, edge2perm_rank, icp);
        unprocessed--;
    }
}

LightCompleteGraph::EdgeWeight EdgeHeap::getIcf(const EdgeWeight uw, const EdgeWeight vw) const {
    if (uw > 0 && vw > 0) {
        // if both other edges present, remove the cheapest of both
        return std::min(uw, vw); 
    } else {
        return 0;
    }
}

LightCompleteGraph::EdgeWeight EdgeHeap::getIcp(const EdgeWeight uw, const EdgeWeight vw) const {
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

void EdgeHeap::updateHeap(std::vector<EdgeId>& heap, const EdgeId e, const EdgeWeight change, std::vector<EdgeId>& index, const std::vector<EdgeWeight>& score) {
    unsigned int pos = index[e];
    if (change > 0) {
        // value increased -> move edge upwards in heap
        while(score[heap[pos/2]] < score[heap[pos]]) {
            // swap pos and pos/2
            std::swap(heap[pos], heap[pos/2]);
            index[heap[pos]] = pos;
            index[heap[pos/2]] = pos/2;
            pos = pos/2;
        }
    } else {
    // value decreased -> move edge downwards in heap
    while((2*pos < heap.size() && score[heap[pos]] < score[heap[2*pos]])
        || (2*pos+1 < heap.size() && score[heap[pos]] < score[heap[2*pos+1]]) ) {
        if (2*pos+1 < heap.size() && score[heap[2*pos]] < score[heap[2*pos+1]]) {
            // element 2*pos+1 exists and is larger than 2*pos -> swap pos with 2*pos+1
            std::swap(heap[pos], heap[2*pos+1]);
            index[heap[pos]] = pos;
            index[heap[pos*2+1]] = pos*2+1;
            pos = 2*pos+1;
        } else {
            // else swap with 2*pos
            std::swap(heap[pos], heap[2*pos]);
            index[heap[pos]] = pos;
            index[heap[pos*2]] = pos*2;
            pos = 2*pos;
        }
    }
}
//   bool wrong = false;
//   unsigned int where = 0;
//   for (unsigned int i = 0; i < heap.size()/2; i++) {
//     if (score[heap[i]] < score[heap[2*i]]) {
//       wrong = true;
//       where = i;
//     }
//     if (2*i+1 < heap.size() && score[heap[i]] < score[heap[2*i+1]]) {
//       wrong = true;
//       where = i;
//     }
//   }
//   if (wrong) {
//     std::cout<<"Error in heap ("<<where<<") "<<score[heap[where/2]]<<" >= "<<score[heap[where]]<<" >= "<<score[heap[2*where]]<<" & "<<score[heap[2*where+1]]<<std::endl;
//   }
}

} // namespace ysk
