#include "edgeheap.h"
#include <cmath>
#include <algorithm>
  
using Edge = StaticSparseGraph::Edge;
using EdgeWeight = StaticSparseGraph::EdgeWeight;
using EdgeId = StaticSparseGraph::EdgeId;
using RankId = StaticSparseGraph::RankId;
using NodeId = StaticSparseGraph::NodeId;

EdgeHeap::EdgeHeap(StaticSparseGraph& param_graph) :
    graph(param_graph),
    unprocessed(0),
    edges(1+param_graph.numEdges(), StaticSparseGraph::InvalidEdge),
    icf(1+param_graph.numEdges(), StaticSparseGraph::Forbidden),
    icp(1+param_graph.numEdges(), StaticSparseGraph::Forbidden),
    edge2forb_rank(1+param_graph.numEdges(), 0),
    edge2perm_rank(1+param_graph.numEdges(), 0),
    edgeToBundle(1+param_graph.numEdges(), 0),
    edgeBundles(1+param_graph.numEdges(), std::vector<RankId>(0))
{}

void EdgeHeap::initInducedCosts() {
    // for each node: sorted list of all non-zero neighbours with higher index
    std::vector<std::vector<NodeId>> higherNeighbour(graph.numNodes(), std::vector<NodeId>(0));
    // compute initial icp/icf values
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        for (NodeId v : graph.getNonZeroNeighbours(u)) {
            if (v < u)
                continue;

            // iterate over all edges uv with u < v
            Edge uv(u,v);
            RankId rId = graph.findIndex(uv);
            edges[rId] = uv;
            icf[rId] = 0.0;
            icp[rId] = 0.0;

            EdgeWeight w_uv = graph.getWeight(rId);
            if (w_uv == 0.0 || w_uv == StaticSparseGraph::Forbidden || w_uv == StaticSparseGraph::Permanent)
                continue;

            unprocessed++;
            higherNeighbour[u].push_back(v);

            // costs for the edge uv itself
            if (w_uv >= 0)
                icf[rId] += w_uv;	// costs for removing uv
            else
                icp[rId] += -w_uv;	// costs for adding uv
        }
        std::sort(higherNeighbour[u].begin(), higherNeighbour[u].end());
    }
    // compute contributions of common neighbors for each edge (u,v)
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        for (uint32_t i = 0; i < higherNeighbour[u].size(); i++) {
            NodeId v = higherNeighbour[u][i];
            Edge uv(u,v);
            RankId rUV = graph.findIndex(uv);
            EdgeWeight w_uv = graph.getWeight(rUV);
            
            // look at all triangles uvw containing uv. Triangles with a zero edge can be ignored
            std::vector<NodeId> w_vec;
            std::set_intersection(higherNeighbour[u].begin() + i, higherNeighbour[u].end(),
                                  higherNeighbour[v].begin(), higherNeighbour[v].end(),
                                  back_inserter(w_vec));

            for (NodeId w : w_vec) {
                // to avoid too many getWeight calls, increase icf/icp for each edge in the triple, but process triple only once
                Edge uw(u,w);
                Edge vw(v,w);
                RankId rUW = graph.findIndex(uw);
                RankId rVW = graph.findIndex(vw);
                EdgeWeight w_uw = graph.getWeight(rUW);
                EdgeWeight w_vw = graph.getWeight(rVW);
                icf[rUV] += getIcf(w_uw, w_vw);
                icp[rUV] += getIcp(w_uw, w_vw);
                icf[rUW] += getIcf(w_uv, w_vw);
                icp[rUW] += getIcp(w_uv, w_vw);
                icf[rVW] += getIcf(w_uv, w_uw);
                icp[rVW] += getIcp(w_uv, w_uw);
            }
        }
    }
    
    // sort edges by icf and icp values
    for (RankId id = 0; id < icf.size(); id++) {
        forb_rank2edge.push_back(id);
        perm_rank2edge.push_back(id);
    }
    
    std::sort(forb_rank2edge.begin(), forb_rank2edge.end(), [this] (const EdgeId& a, const EdgeId& b) { return icf[a] > icf[b]; });
    std::sort(perm_rank2edge.begin(), perm_rank2edge.end(), [this] (const EdgeId& a, const EdgeId& b) { return icp[a] > icp[b]; });
    
    // save index in sorted vectors for each edge
    for (RankId i = 0; i < icf.size(); i++) {
        edge2forb_rank[forb_rank2edge[i]] = i;
        edge2perm_rank[perm_rank2edge[i]] = i;
    }
    
    // initialize edge bundles
    for (RankId id = 0; id < icf.size(); id++) {
        edgeToBundle[id] = id;
        edgeBundles[id].push_back(id);
    }
}

Edge EdgeHeap::getMaxIcfEdge() const {
    RankId ei = forb_rank2edge[0];
    if (forb_rank2edge.size() <= 1) {
        // only rank 0 entry left
        return StaticSparseGraph::InvalidEdge;
    }
    if (icf[ei] < 0) {
        return StaticSparseGraph::InvalidEdge;
    }
    return edges[ei];
}

Edge EdgeHeap::getMaxIcpEdge() const {
    RankId ei = perm_rank2edge[0];
    if (perm_rank2edge.size() <= 1) {
        // only rank 0 entry left
        return StaticSparseGraph::InvalidEdge;
    }
    if (icp[ei] < 0) {
        return StaticSparseGraph::InvalidEdge;
    }
    return edges[ei];
}

EdgeWeight EdgeHeap::getIcf(const Edge e) const {
    return icf[edgeToBundle[graph.findIndex(e)]];
}

EdgeWeight EdgeHeap::getIcp(const Edge e) const {
    return icp[edgeToBundle[graph.findIndex(e)]];
}

void EdgeHeap::increaseIcf(const Edge e, const EdgeWeight w) {
    RankId rId = graph.findIndex(e);
    if (rId > 0 && icf[edgeToBundle[rId]] >= 0) {
        RankId eb = edgeToBundle[rId];
        icf[eb] = std::max(icf[eb] + w, 0.0f);
        updateHeap(forb_rank2edge, eb, w, edge2forb_rank, icf);
    }
}

void EdgeHeap::increaseIcp(const Edge e, const EdgeWeight w) {
    RankId rId = graph.findIndex(e);
    if (rId > 0 && icp[edgeToBundle[rId]] >= 0) {
        RankId eb = edgeToBundle[rId];
        icp[eb] = std::max(icp[eb] + w, 0.0f);
        updateHeap(perm_rank2edge, eb, w, edge2perm_rank, icp);
    }
}

void EdgeHeap::mergeEdges(const Edge e1, const Edge e2) {
    RankId r1 = graph.findIndex(e1);
    RankId r2 = graph.findIndex(e2);
    if ((r1 & r2) == 0)
        return;
    RankId eb1 = edgeToBundle[r1];
    RankId eb2 = edgeToBundle[r2];
    if (eb1 == eb2)
        return;
    
    if (edgeBundles[eb1].size() > edgeBundles[eb2].size()) {
        for (RankId toDelete : edgeBundles[eb2]) {
            edgeBundles[eb1].push_back(toDelete);
            edgeToBundle[toDelete] = eb1;
        }
        edgeBundles[eb2].clear();
        if (icf[eb2] < 0.0f) {
            std::cout<<"Merged edge has negative icf"<<std::endl;
        } else {
            icf[eb1] += icf[eb2];
        }
        if (icp[eb2] < 0.0f) {
            std::cout<<"Merged edge has negative icp"<<std::endl;
        } else {
            icp[eb1] += icp[eb2];
        }
        removeEdge(eb2);
    } else {
        for (RankId toDelete : edgeBundles[eb1]) {
            edgeBundles[eb2].push_back(toDelete);
            edgeToBundle[toDelete] = eb2;
        }
        edgeBundles[eb1].clear();
        if (icf[eb1] < 0.0f) {
            std::cout<<"Merged edge has negative icf"<<std::endl;
        } else {
            icf[eb2] += icf[eb1];
        }
        if (icp[eb1] < 0.0f) {
            std::cout<<"Merged edge has negative icp"<<std::endl;
        } else {
            icp[eb2] += icp[eb1];
        }
        removeEdge(eb1);
    }
}

void EdgeHeap::removeEdge(const Edge e) {
    removeEdge(graph.findIndex(e));
}

void EdgeHeap::removeEdge(const RankId rId) {
    if (rId == 0) {
        return;
    }
    if (rId > 0 && icf[rId] != StaticSparseGraph::Forbidden && icp[rId] != StaticSparseGraph::Forbidden) {
        icf[rId] = StaticSparseGraph::Forbidden;
        icp[rId] = StaticSparseGraph::Forbidden;
        updateHeap(forb_rank2edge, rId, StaticSparseGraph::Forbidden, edge2forb_rank, icf);
        updateHeap(perm_rank2edge, rId, StaticSparseGraph::Forbidden, edge2perm_rank, icp);
        unprocessed--;
    }
}

uint64_t EdgeHeap::numUnprocessed() const {
    return unprocessed;
}

void EdgeHeap::updateHeap(std::vector<RankId>& heap, const RankId e, const EdgeWeight change, std::vector<RankId>& index, const std::vector<EdgeWeight>& score) {
    uint64_t pos = index[e];
    /*
     * index arithemetic for zero based array: parent = (index-1)/2, children = 2*index+1 and 2*index+2
     */
    if (change > 0) {
        // value increased -> move edge upwards in heap
        uint64_t parent = (pos-1)/2;
        while(pos > 0 && score[heap[parent]] < score[heap[pos]]) {
            // swap pos and pos/2
            std::swap(heap[pos], heap[parent]);
            index[heap[pos]] = pos;
            index[heap[parent]] = parent;
            pos = parent;
            parent = (pos-1)/2;
        }
    } else {
        uint64_t lChild = pos + (pos + 1) * ((2*pos + 1) < heap.size());
        uint64_t rChild = pos + (pos + 2) * ((2*pos + 2) < heap.size());
        uint64_t next = lChild * (score[heap[rChild]] <= score[heap[lChild]]) + rChild * (score[heap[lChild]] < score[heap[rChild]]);
        while (score[heap[pos]] < score[heap[next]]) {
            std::swap(heap[pos], heap[next]);
            index[heap[pos]] = pos;
            index[heap[next]] = next;
            pos = next;
            lChild = pos + (pos + 1) * ((2*pos + 1) < heap.size());
            rChild = pos + (pos + 2) * ((2*pos + 2) < heap.size());
            next = lChild * (score[heap[rChild]] <= score[heap[lChild]]) + rChild * (score[heap[lChild]] < score[heap[rChild]]);
        }
    }
}
