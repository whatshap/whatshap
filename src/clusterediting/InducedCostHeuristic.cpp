#include "InducedCostHeuristic.h"
#include <queue>
#include <algorithm>
#include <unordered_map>
  
using Edge = DynamicSparseGraph::Edge;
using EdgeWeight = DynamicSparseGraph::EdgeWeight;
using EdgeId = DynamicSparseGraph::EdgeId;
using NodeId = DynamicSparseGraph::NodeId;
using RankId = DynamicSparseGraph::RankId;

InducedCostHeuristic::InducedCostHeuristic(StaticSparseGraph& param_graph, bool param_bundleEdges) :
    bundleEdges(param_bundleEdges),
    graph(param_graph),
    edgeHeap(graph),
    totalCost(0.0)
{
    if (!resolvePermanentForbidden()) {
        totalCost = std::numeric_limits<EdgeWeight>::infinity();
    }
    edgeHeap.initInducedCosts();
}

ClusterEditingSolutionLight InducedCostHeuristic::solve() {
    // execute algorithm
    if (verbosity >= 1) {
        if (verbosity == 1)
            std::cout<<"Running heuristic.. \r" <<std::flush;
        else
            std::cout<<"Running heuristic.. " <<std::endl;
    }
    if (totalCost == std::numeric_limits<EdgeWeight>::infinity()) {
        // if resolving permanent and forbidden edges lead to contradiction, cost are infinte here, thus cancel algorithm
        std::cout<<"Instance is infeasible!" <<std::endl;
        ClusterEditingSolutionLight sol;
        return sol;
    }
  
    int totalEdges = edgeHeap.numUnprocessed();
    for (uint64_t i = 0; i < graph.numEdges() + 1; i++) {
        Edge eIcf = edgeHeap.getMaxIcfEdge();
        Edge eIcp = edgeHeap.getMaxIcpEdge();
        if (eIcf == DynamicSparseGraph::InvalidEdge || eIcp == DynamicSparseGraph::InvalidEdge) {
            break;
        }
        EdgeWeight mIcf = edgeHeap.getIcf(eIcf);
        EdgeWeight mIcp = edgeHeap.getIcp(eIcp);
        
        if (mIcf >= mIcp) {
            // set eIcf to permanent
            if (verbosity >= 5) {
                std::cout<<"Setting edge ("<<eIcf.u<<","<<eIcf.v<<") to permanent."<<std::endl;
            }

            // get edges, which will become permanent due to implication. This must be done before the actual edge is
            // modified, as this would implicitly modify the weights in the graph as well
            std::vector<Edge> implications;
            std::vector<NodeId> uClique(graph.getCliqueOf(eIcf.u));
            std::vector<NodeId> vClique(graph.getCliqueOf(eIcf.v));
            if (verbosity >= 5) {
                std::cout<<"Clique of "<<eIcf.u<<": ";
                for (const auto& i: uClique)
                    std::cout << i << ' ';
                std::cout<<std::endl;
                std::cout<<"Clique of "<<eIcf.v<<": ";
                for (const auto& i: vClique)
                    std::cout << i << ' ';
                std::cout<<std::endl;
            }
            
            for (NodeId x : uClique) {
                for (NodeId y : vClique) {
                    Edge e = Edge(x,y);
                    if (x == y || graph.getWeight(e) == DynamicSparseGraph::Permanent || (x == eIcf.u && y == eIcf.v)) {
                        if (verbosity >= 5) {
                            std::cout<<"Making ("<<x<<","<<y<<") silently not permanent due to implication."<<std::endl;
                        }
                        continue;
                    }
                    if (verbosity >= 5) {
                        std::cout<<"Making ("<<x<<","<<y<<") permanent due to implication."<<std::endl;
                    }
                    implications.push_back(e);
                }
            }

            // set actual edge first
            setPermanent(eIcf);
            edgeHeap.removeEdge(eIcf);
            
            // then all implications
            for (Edge e : implications) {
                setPermanent(e);
                edgeHeap.removeEdge(e);
                if (verbosity >= 1 && edgeHeap.numUnprocessed() % 1000 == 0) {
                    std::cout<<"Running heuristic.. "<<((totalEdges - edgeHeap.numUnprocessed())*100 / totalEdges)<<"%\r"<<std::flush;
                }
            }
            
            if (bundleEdges) {
                /* All outgoing edges of the current clique, which lead to the same clique, must be bundled together. If one of these
                edges is made permanent or forbidden, then all other edges must be as well. Therefore not the induced costs for all
                edges must be taken into account, not only the cost of the chosen edge.*/
                NodeId cu = graph.getCliqueIdOf(eIcf.u);
                if (verbosity >= 4)
                    std::cout<<"Contracting nodes of cluster id ("<<cu<<")."<<std::endl;
                std::unordered_map<NodeId, Edge> cliqueToRepresentative;
                uClique.insert(uClique.end(), vClique.begin(), vClique.end());
                for (NodeId x : uClique) {
                    for (NodeId xn : graph.getUnprunedNeighbours(x)) {
                        // this edge should not be inside the current cluster, as all internal edges should be permanent by now
                        Edge ex(x, xn);
                        NodeId cxn = graph.getCliqueIdOf(xn);
                        
                        if (std::find(uClique.begin(), uClique.end(), xn) != uClique.end()) {
                            if (verbosity >= 5)
                                std::cout<<"Observed edge ("<<x<<","<<xn<<") was inside the cluster!"<<std::endl;
                            continue;
                        }
                        if (graph.findIndex(ex) == 0) {
                            std::cout<<"Observed edge ("<<x<<","<<xn<<") was pruned edge with weight "<<graph.getWeight(ex)<<std::endl;
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
        } else {
            // set eIcp fo forbidden
            if (verbosity >= 5) {
                std::cout<<"Setting edge ("<<eIcp.u<<","<<eIcp.v<<") to forbidden."<<std::endl;
            }

            // get edges, which will become permanent due to implication. This must be done before the actual edge is
            // modified, as this would implicitly modify the weights in the graph as well
            std::vector<Edge> implications;
            std::vector<NodeId> uClique(graph.getCliqueOf(eIcp.u));
            std::vector<NodeId> vClique(graph.getCliqueOf(eIcp.v));
            for (NodeId x : uClique) {
                for (NodeId y : vClique) {
                    Edge e = Edge(x,y);
                    if (x == y || graph.getWeight(e) == DynamicSparseGraph::Forbidden || (x == eIcp.u && y == eIcp.v)) {
                        if (verbosity >= 5) {
                            std::cout<<"Making ("<<x<<","<<y<<") silently not forbidden due to implication."<<std::endl;
                        }
                        continue;
                    }
                    if (verbosity >= 5) {
                        std::cout<<"Making ("<<x<<","<<y<<") forbidden due to implication."<<std::endl;
                    }
                    implications.push_back(e);
                }
            }

            // set actual edge first
            setForbidden(eIcp);
            edgeHeap.removeEdge(eIcp);
            
            // then all implications
            for (Edge e : implications) {
                setForbidden(e);
                edgeHeap.removeEdge(e);
                if (verbosity >= 1 && edgeHeap.numUnprocessed() % 1000 == 0) {
                    std::cout<<"Running heuristic.. "<<((totalEdges - edgeHeap.numUnprocessed())*100 / totalEdges)<<"%\r"<<std::flush;
                }
            }
        }
        if (verbosity >= 1 && edgeHeap.numUnprocessed() % 100 == 0) {
            std::cout<<"Running heuristic.. "<<((totalEdges - edgeHeap.numUnprocessed())*100 / totalEdges)<<"%\r"<<std::flush;
        }
    }

    if (verbosity >= 1)
        std::cout<<"Running heuristic.. 100%   "<<std::endl;
  
    // calculate clustering
    if (verbosity >= 1) {
        if (verbosity == 1)
            std::cout<<"Constructing result.. \r"<<std::flush;
        else
            std::cout<<"Constructing result.. "<<std::endl;
    }
    std::vector<std::vector<NodeId>> clusters;
    std::vector<int> clusterOfNode(graph.numNodes(), -1);
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        // add cluster if not explored yet
        if (verbosity >= 4) {
            std::cout<<"Processing node "<<u<<std::endl;
        }
        if (verbosity >= 1 && verbosity <= 4) {
            std::cout<<"Constructing result.. "<<(u*100/graph.numNodes())<<"%\r"<<std::flush;
        }
        if (clusterOfNode[u] == -1) {
            int c = clusters.size();
            if (verbosity >= 4) {
                std::cout<<"Node "<<u<<" not in any cluster yet. Creating new cluster "<<c<<" for this"<<std::endl;
            }
            clusterOfNode[u] = c;
            clusters.push_back(std::vector<NodeId>(1, u));
            for (NodeId v = u+1; v < graph.numNodes(); v++) {
                if (graph.getWeight(Edge(u, v)) == DynamicSparseGraph::Permanent) {
                    clusterOfNode[v] = c;
                    if (verbosity >= 4) {
                        std::cout<<"Adding connected node "<<v<<" in same cluster."<<std::endl;
                    }
                    clusters[c].push_back(v);
                }
            }
        }
    }
    if (verbosity >= 1)
        std::cout<<"Constructing result.. 100%   "<<std::endl;
    return ClusterEditingSolutionLight(totalCost, clusters);
}

bool InducedCostHeuristic::resolvePermanentForbidden() {
    if (verbosity >= 1) {
        if (verbosity == 1)
            std::cout<<"Resolving forbidden and permanent edges.. \r" <<std::flush;
        else
            std::cout<<"Resolving forbidden and permanent edges.. " <<std::endl;
    }
    // make cliques by connecting all nodes with inf path between them
    std::vector<bool> processed(graph.numNodes(), false);
    std::vector<std::vector<NodeId>> cliques;
    std::vector<std::vector<NodeId>> moreThanOneCliques;
    for (NodeId u = 0; u < graph.numNodes(); u++) {
        if (processed[u])
            continue;
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
        if (clique.size() > 1) {
            moreThanOneCliques.push_back(clique);
        }
        for (NodeId x : clique) {
            for (NodeId y : clique) {
                if (x != y) {
                    Edge e (x,y);
                    EdgeWeight w = graph.getWeight(e);
                    if (w == DynamicSparseGraph::Forbidden)
                        return false;
                    else if (w != DynamicSparseGraph::Permanent) {
                        if (w < 0.0)
                            totalCost -= w;
                        graph.setPermanent(Edge(x,y));
                        if (verbosity >= 5) {
                            std::cout<<"Making ("<<x<<","<<y<<") permanent due to implication."<<std::endl;
                        }
                    }
                }
            }
        }
    }
    
    // disconnect all cliques which have a forbidden edge between them
    for (unsigned int k = 0; k < cliques.size(); k++) {
        if (verbosity >= 1 && k % 100 == 0) {
            std::cout<<"Resolving forbidden and permanent edges.. "<<(((2UL*cliques.size()-(uint64_t)k+1UL)*(uint64_t)k*100UL)/((uint64_t)cliques.size()*((uint64_t)cliques.size()-1UL)))<<"%\r"<<std::flush;
        }
        for (unsigned int l = 0; l < moreThanOneCliques.size(); l++) {
            // search for forbidden edge between
            bool found = false;
            for (NodeId u : cliques[k]) {
                if (found) break;
                for (NodeId v : moreThanOneCliques[l]) {
                    if (graph.getWeight(Edge(u, v)) == DynamicSparseGraph::Forbidden) {
                        found = true;
                        break;
                    }
                }
            }
            // make all edges forbidden, if one forbidden edge was found
            if (found) {
                for (NodeId u : cliques[k]) {
                    for (NodeId v : moreThanOneCliques[l]) {
                        Edge e(u,v);
                        if (graph.getWeight(e) != DynamicSparseGraph::Forbidden) {
                            graph.setForbidden(e);
                            if (verbosity >= 5) {
                                std::cout<<"Making ("<<u<<","<<v<<") forbidden due to implication."<<std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
    if (verbosity >= 1)
        std::cout<<"Resolving forbidden and permanent edges.. 100%   "<<std::endl;
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
//     EdgeWeight uv = id == 0 ? 0.0 : graph.getWeight(id);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        if (graph.findIndex(vw) > 0)
            updateTripleForbiddenUW(uv, uw, graph.getWeight(vw));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        if (graph.findIndex(vw) > 0)
            updateTripleForbiddenUW(uv, vw, graph.getWeight(uw));
    }
    if (uv > 0) {
        totalCost += uv;
        graph.setForbidden(e);
    } else if (uv < 0) {
        graph.setForbidden(e);
    }
}

void InducedCostHeuristic::setPermanent(const Edge e) {
    // this has to be called to update ic, even if the edge already is permanent
    NodeId u = e.u;
    NodeId v = e.v;
    RankId id = graph.findIndex(e);
    
    /* If the edge was a zero edge in the original graph, it might have been implicitly set to
     * permanent or forbidden without updating the ic. Therefore we assume the weight to be 0 here.*/
    EdgeWeight uv = graph.getWeight(id);
//     EdgeWeight uv = id == 0 ? 0.0 : graph.getWeight(id);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        if (graph.findIndex(vw) > 0)
            updateTriplePermanentUW(uv, uw, graph.getWeight(vw));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        if (graph.findIndex(vw) > 0)
            updateTriplePermanentUW(uv, vw, graph.getWeight(uw));
    }
    if (uv < 0) {
        totalCost -= uv;
        graph.setPermanent(e);
    } else if (uv > 0) {
        graph.setPermanent(e);
    }
}

void InducedCostHeuristic::updateTripleForbiddenUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
    EdgeWeight icf_old = edgeHeap.getIcf(uv, vw);
    EdgeWeight icf_new = 0.0;
    EdgeWeight icp_old = edgeHeap.getIcp(uv, vw);
    EdgeWeight icp_new = std::max(0.0, vw);
    if (icf_new != icf_old)
        edgeHeap.increaseIcf(uw, icf_new - icf_old);
    if (icp_new != icp_old)
        edgeHeap.increaseIcp(uw, icp_new - icp_old);
}

void InducedCostHeuristic::updateTriplePermanentUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
    EdgeWeight icf_old = edgeHeap.getIcf(uv, vw);
    EdgeWeight icf_new = std::max(0.0, vw);
    EdgeWeight icp_old = edgeHeap.getIcp(uv, vw);
    EdgeWeight icp_new = std::max(0.0, -vw);
    if (icf_new != icf_old)
        edgeHeap.increaseIcf(uw, icf_new - icf_old);
    if (icp_new != icp_old)
        edgeHeap.increaseIcp(uw, icp_new - icp_old);
}
