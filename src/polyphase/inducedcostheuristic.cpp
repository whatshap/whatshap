#include "inducedcostheuristic.h"
#include <queue>
#include <algorithm>
#include <unordered_map>
  
using Edge = StaticSparseGraph::Edge;
using EdgeWeight = StaticSparseGraph::EdgeWeight;
using EdgeId = StaticSparseGraph::EdgeId;
using NodeId = StaticSparseGraph::NodeId;
using RankId = StaticSparseGraph::RankId;

InducedCostHeuristic::InducedCostHeuristic(StaticSparseGraph& param_graph, bool param_bundleEdges) :
    bundleEdges(param_bundleEdges),
    graph(param_graph),
    edgeHeap(graph),
    totalCost(0.0),
    totalEdges(0)
{
    /* Preprocessing: Find all forbidden and permanent edges, which are already in the graph. These
     * edges may either imply other edges to be permanent or forbidden or they might cause
     * contradictions, i.e. the cost the make it a clique graph is infinity.*/
    edgeHeap.initInducedCosts();
	if (!resolvePermanentForbidden()) {
        totalCost = std::numeric_limits<EdgeWeight>::infinity();
    }
    totalEdges = edgeHeap.numUnprocessed();
}

ClusterEditingSolution InducedCostHeuristic::solve() {
    // check if instance is solvable at all
    if (totalCost == std::numeric_limits<EdgeWeight>::infinity()) {
        // if resolving permanent and forbidden edges lead to contradiction, cost are infinte here, thus cancel algorithm
        std::cout<<"Instance is infeasible!" <<std::endl;
        ClusterEditingSolution sol;
        return sol;
    }
  
    /* In each iteration, extract edge with highest induced cost (either for becoming permanent or forbidden).
     * This edge will be set to permanent or forbidden, depending on what is cheaper according to icf and icp */
    for (uint64_t i = 0; i < graph.numEdges() + 1; i++) {
        Edge eIcf = edgeHeap.getMaxIcfEdge();
        Edge eIcp = edgeHeap.getMaxIcpEdge();
        
        // if edge heap returns an invalid edge, we know the heap is empty and all edges are processed
        if (eIcf == StaticSparseGraph::InvalidEdge || eIcp == StaticSparseGraph::InvalidEdge) {
            break;
        }
        
        // determine the induced costs for thw two edges
        EdgeWeight mIcf = edgeHeap.getIcf(eIcf);
        EdgeWeight mIcp = edgeHeap.getIcp(eIcp);
        
        if (mIcf >= mIcp) {
            // forbidden cost are the highest, set the corresponding edge to permanent
            choosePermanentEdge(eIcf);
        } else {
            // permanent cost are the highest, set the correspondong edge to forbidden
            chooseForbiddenEdge(eIcp);
        }
    }

    /* Construct the clusters, by finding group of nodes, which are connected via a permanent edge.
     * Assuming the correctness of the heuristic above, there should not be three nodes u, v, w such that
     * (u, v) and (v, w) is permanent, but (u, w) is not. Zero edges, which have not been set to either
     * permanent or forbidden, are considered forbidden, i.e. not part of any clique.*/
    std::vector<std::vector<NodeId>> clusters;
    std::vector<int> clusterOfNode(graph.numNodes(), -1);
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        // if u is not in a cluster yet, create a new cluster with u in it
        if (clusterOfNode[u] == -1) {
            int c = clusters.size();
            clusterOfNode[u] = c;
            clusters.push_back(std::vector<NodeId>(1, u));
            for (NodeId v : graph.getCliqueOf(u)) {
                if (u == v)
                    continue;
                clusterOfNode[v] = c;
                clusters[c].push_back(v);
            }
        }
    }
    
    // sort node ids in each cluster in ascending order
    for (std::vector<NodeId>& cluster : clusters) {
        std::sort(cluster.begin(), cluster.end());
    }
    return ClusterEditingSolution(totalCost, clusters);
}

void InducedCostHeuristic::choosePermanentEdge(const StaticSparseGraph::Edge eIcf) {
    /* We cannot just set the edge eIcf=(u,v) to permanent, because we have to handle implications of this step as well.
     * According to the heuristic, u and v must be merged into one node. However, we do not do this here, but instead
     * make sure that u and v are handled as a clique:
     * Node u and v might already be permanently connected to other nodes. If u and v are connected, the cliques of u and v
     * must be pairwise connected, too. For non-zero edges, we could let this be handled by the heuristic,
     * but for zero edges, we must do this.*/
    std::vector<Edge> implications;
    
    /* All pairs from uClique and vClique are found. We must be careful not to add eIcf to our list (we already have it) and
     * to not add edges, which are already permanent. We delay the actual modification of the edge until we know, which edges
     * must become permanent, because the weight for zero-edges implicitly changes, if other edges change. This would disturb
     * the search process.*/
    //TODO: Find reason, why some edges are already permanent here.
    for (NodeId x : graph.getCliqueOf(eIcf.u)) {
        for (NodeId y : graph.getCliqueOf(eIcf.v)) {
            Edge e = Edge(x,y);
            if (x == y || graph.findIndex(e) == 0 || (x == eIcf.u && y == eIcf.v)) {
                continue;
            }
            implications.push_back(e);
        }
    }
    
    /* The cliques we are connecting here might already be forbidden to other nodes/cliques. So, we need a second list of
     * implications to collect edges, which must be set to forbidden afterwards.*/
    std::vector<Edge> implicationsForbidden;
    for (NodeId f : graph.getForbiddenNeighbors(eIcf.u)) {
        if (graph.isForbidden(Edge(eIcf.v, f)))
            continue;
        for (NodeId x : graph.getCliqueOf(eIcf.v)) {
            Edge e = Edge(f,x);
            if (graph.findIndex(e) != 0 && !graph.isForbidden(e)) {
                implicationsForbidden.push_back(e);
            }
        }
    }
    for (NodeId f : graph.getForbiddenNeighbors(eIcf.v)) {
        if (graph.isForbidden(Edge(eIcf.u, f)))
            continue;
        for (NodeId x : graph.getCliqueOf(eIcf.u)) {
            Edge e = Edge(f,x);
            if (graph.findIndex(e) != 0 && !graph.isForbidden(e)) {
                implicationsForbidden.push_back(e);
            }
        }
    }

    // First, modify eIcf ...
    setPermanent(eIcf);
    edgeHeap.removeEdge(eIcf);
    
    // ... then all implications ...
    for (Edge e : implications) {
        setPermanent(e);
        edgeHeap.removeEdge(e);
   }
    
    // ... and all forbidden implications
    for (Edge e : implicationsForbidden) {
        setForbidden(e);
        edgeHeap.removeEdge(e);
    }
    
    if (bundleEdges) {
        /* Setting an edge to permanent must make u and v (and their cliques) to act as single node. Specifically, for every neighbor
         * of the clique, there must be uniform induced costs for making the connecting edge forbidden or permanent. To accomplish this
         * the edge heap organizes edges in bundles. At first every edge is its own bundle. If two nodes u and v are merged, all edges,
         * which to the same node w (w != u and w != v and w not in the same clique as u and v) are bundled together.*/
        std::unordered_map<NodeId, Edge> cliqueToRepresentative;
        /* Idea: We iterate over all outgoing edges from the combined clique of u and v. If we reach another cluster, we have not seen
         * before, we memorize the outgoing edge as representative edge for this cluster. When we find another edge leading to a cluster,
         * we have already seen, we merge this edge with the representative of this cluster.*/
        
        std::vector<NodeId> uClique(graph.getCliqueOf(eIcf.u));
        std::vector<NodeId> vClique(graph.getCliqueOf(eIcf.v));
        uClique.insert(uClique.end(), vClique.begin(), vClique.end());
        for (NodeId x : uClique) {
            for (NodeId xn : graph.getUnprunedNeighbours(x)) {
                // this edge should not be inside the current cluster, as all internal edges should be permanent by now
                Edge ex(x, xn);
                NodeId cxn = graph.getCliqueIdOf(xn);
                
                if (std::find(uClique.begin(), uClique.end(), xn) != uClique.end()) {
                    continue;
                }
                if (graph.findIndex(ex) == 0) {
                    continue;
                }
                // if new cluster is "discovered", set edge as representative, otherwise bundle with present representative
                if (cliqueToRepresentative.find(cxn) == cliqueToRepresentative.end()) {
                    cliqueToRepresentative[cxn] = ex;
                } else {
                    edgeHeap.mergeEdges(ex, cliqueToRepresentative[cxn]);
                }
            }
        }
    }
}

void InducedCostHeuristic::chooseForbiddenEdge(const StaticSparseGraph::Edge eIcp) {
    /* We cannot just set the edge eIcf=(u,v) to forbidden, because we have to handle implications of this step as well.
     * Node u and v might already permanently connected to other nodes. If we decide to not u and v into one clique, then
     * all other pair of nodes in the same clique as u and v must be a forbidden pair. For non-zero edges, we could let 
     * this be handled by the heuristic, but for zero edges, we must do this.*/
    std::vector<Edge> implications;
    
    /* All pairs from uClique and vClique are found. We must be careful not to add eIcp to our list (we already have it) and
     * to not add edges, which are already forbidden. We delay the actual modification of the edge until we know, which edges
     * must become forbidden, because the weight for zero-edges implicitly changes, if other edges change. This would disturb
     * the search process.*/
    //TODO: Find reason, why some edges are already forbidden here.
    for (NodeId x : graph.getCliqueOf(eIcp.u)) {
        for (NodeId y : graph.getCliqueOf(eIcp.v)) {
            Edge e = Edge(x,y);
            if (x == y || graph.findIndex(e) == 0 || (x == eIcp.u && y == eIcp.v)) {
                continue;
            }
            implications.push_back(e);
        }
    }

    // First, modify eIcp ...
    setForbidden(eIcp);
    edgeHeap.removeEdge(eIcp);
    
    // ... then all implications
    for (Edge e : implications) {
        setForbidden(e);
        edgeHeap.removeEdge(e);
    }
}


bool InducedCostHeuristic::resolvePermanentForbidden() {
    // make cliques by connecting all nodes with inf path between them
    std::vector<bool> processed(graph.numNodes(), false);
    std::vector<std::vector<NodeId>> cliques;
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        if (processed[u]) {
            continue;
        }
        std::vector<NodeId> clique;
        std::queue<NodeId> remaining;
        remaining.push(u);
        processed[u] = true;
        while (!remaining.empty()) {
            NodeId current = remaining.front();
            remaining.pop();
            clique.push_back(current);
            for (NodeId v : graph.getCliqueOf(current)) {
                if (!processed[v]) {
                    remaining.push(v);
                    processed[v] = true;
                }
            }
        }
        cliques.push_back(clique);
        for (uint32_t i = 0; i < clique.size(); i++) {
            for (uint32_t j = i + 1; j < clique.size(); j++) {
                Edge e (clique[i], clique[j]);
                EdgeWeight w = graph.getWeight(e);
                if (w == StaticSparseGraph::Forbidden)
                    return false;
                else {
                    edgeHeap.removeEdge(e);
                    if (w < 0.0)
                        totalCost -= w;
                    if (w != StaticSparseGraph::Permanent)
                        graph.setPermanent(e);
                }
            }
        }
    }
    
    // disconnect all cliques which have a forbidden edge between them
    for (unsigned int k = 0; k < cliques.size(); k++) {
        for (unsigned int l = k + 1; l < cliques.size(); l++) {
            // search for forbidden edge between
            bool found = false;
            for (NodeId u : cliques[k]) {
                if (found) break;
                for (NodeId v : cliques[l]) {
                    if (u == v)
                        continue;
                    if (graph.getWeight(Edge(u, v)) == StaticSparseGraph::Forbidden) {
                        found = true;
                        break;
                    }
                }
            }
            // make all edges forbidden, if one forbidden edge was found
            if (found) {
                for (NodeId u : cliques[k]) {
                    for (NodeId v : cliques[l]) {
                        Edge e(u, v);
                        edgeHeap.removeEdge(e);
                        EdgeWeight w = graph.getWeight(e);
                        if (w > 0.0)
                            totalCost += w;
                        if (w != StaticSparseGraph::Forbidden)
                            graph.setForbidden(e);
                    }
                }
            }
        }
    }
    return true;
}

void InducedCostHeuristic::setForbidden(const Edge e) {
    // this has to be called to update ic, even if the edge already is forbidden
    NodeId u = e.u;
    NodeId v = e.v;
    RankId id = graph.findIndex(e);
    
    /* If the edge was a zero edge in the original graph, it might have been implicitly set to
     * permanent or forbidden without updating the ic. Therefore we assume the weight to be 0 here.*/
    EdgeWeight uv = graph.getWeight(id);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        RankId r = graph.findIndex(vw);
        if (r > 0)
            updateTripleForbiddenUW(uv, uw, graph.getWeight(r));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        RankId r = graph.findIndex(uw);
        if (r > 0)
            updateTripleForbiddenUW(uv, vw, graph.getWeight(r));
    }
    if (uv > 0)
        totalCost += uv;
    if (id > 0)
        graph.setForbidden(e, id);
}

void InducedCostHeuristic::setPermanent(const Edge e) {
    // this has to be called to update ic, even if the edge already is permanent
    NodeId u = e.u;
    NodeId v = e.v;
    RankId id = graph.findIndex(e);
    
    /* If the edge was a zero edge in the original graph, it might have been implicitly set to
     * permanent or forbidden without updating the ic. Therefore we assume the weight to be 0 here.*/
    EdgeWeight uv = graph.getWeight(id);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        RankId r = graph.findIndex(vw);
        if (r > 0)
            updateTriplePermanentUW(uv, uw, graph.getWeight(r));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        RankId r = graph.findIndex(uw);
        if (graph.findIndex(r) > 0)
            updateTriplePermanentUW(uv, vw, graph.getWeight(r));
    }
    if (uv < 0)
        totalCost -= uv;
    if (id > 0)
        graph.setPermanent(e, id);
}

void InducedCostHeuristic::updateTripleForbiddenUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
    EdgeWeight icfChange = 0.0 - edgeHeap.getIcf(uv, vw);
    EdgeWeight icpChange = std::max(0.0f, vw) - edgeHeap.getIcp(uv, vw);
    if (icfChange != 0)
        edgeHeap.increaseIcf(uw, icfChange);
    if (icpChange != 0)
        edgeHeap.increaseIcp(uw, icpChange);
}

void InducedCostHeuristic::updateTriplePermanentUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
    EdgeWeight icfChange = std::max(0.0f, vw) - edgeHeap.getIcf(uv, vw);
    EdgeWeight icpChange = std::max(0.0f, -vw) - edgeHeap.getIcp(uv, vw);
    if (icfChange != 0)
        edgeHeap.increaseIcf(uw, icfChange);
    if (icpChange != 0)
        edgeHeap.increaseIcp(uw, icpChange);
}
