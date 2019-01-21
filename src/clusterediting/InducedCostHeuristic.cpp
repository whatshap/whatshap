#include "InducedCostHeuristic.h"
#include <queue>
#include <algorithm>
#include <unordered_map>

namespace ysk {
  
using Edge = StaticSparseGraph::Edge;
using EdgeWeight = StaticSparseGraph::EdgeWeight;
using EdgeId = StaticSparseGraph::EdgeId;
using NodeId = StaticSparseGraph::NodeId;
using RankId = StaticSparseGraph::RankId;

InducedCostHeuristic::InducedCostHeuristic(StaticSparseGraph& param_graph, bool param_pruneZeroEdges) :
    pruneZeroEdges(param_pruneZeroEdges),
    graph(param_graph),
    edgeHeap(graph, param_pruneZeroEdges),
    totalCost(0.0)
{
    if (!resolvePermanentForbidden()) {
        totalCost = std::numeric_limits<EdgeWeight>::infinity();
    }
    graph.compile();
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
    for (unsigned int i = 0; i < graph.numEdges() + 1; i++) {
        Edge eIcf = edgeHeap.getMaxIcfEdge();
        Edge eIcp = edgeHeap.getMaxIcpEdge();
        if (eIcf == StaticSparseGraph::InvalidEdge || eIcp == StaticSparseGraph::InvalidEdge) {
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
                    if (x == y || /*graph.findIndex(e) == 0 ||*/ graph.getWeight(e) == StaticSparseGraph::Permanent || (x == eIcf.u && y == eIcf.v)) {
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
            
//             /* contract all nodes in the new cluster by choosing one outgoing edge for each adjacent cluster and
//              * cumulating all weights and induced costs to these representative edges
//              * other edges get value zero.*/
//             NodeId cu = graph.getCliqueIdOf(eIcf.u);
//             if (verbosity >= 4)
//                 std::cout<<"Contracting nodes of cluster id ("<<cu<<")."<<std::endl;
//             std::unordered_map<NodeId, EdgeWeight> weightToCluster;
//             std::unordered_map<NodeId, EdgeWeight> positiveWeightToCluster;
//             std::unordered_map<NodeId, EdgeWeight> negativeWeightToCluster;
//             std::vector<Edge> removedEdges;
//             std::vector<Edge> representatives;
//             uClique.insert(uClique.end(), vClique.begin(), vClique.end());
//             for (NodeId x : uClique) {
//                 for (NodeId xn : graph.getNonZeroNeighbours(x)) {
//                     // this edge should not be inside the current cluster, as all internal edges should be permanent by now
//                     Edge ex(x, xn);
//                     NodeId cxn = graph.getCliqueIdOf(xn);
//                     EdgeWeight wxn = graph.getWeight(ex);
//                     if (verbosity >= 4 && std::find(uClique.begin(), uClique.end(), xn) != uClique.end()) {
//                         std::cout<<"Observed edge ("<<x<<","<<xn<<") was inside the cluster!"<<std::endl;
//                         continue;
//                     }
//                     if (graph.findIndex(ex) == 0 /*|| wxn == StaticSparseGraph::Forbidden || wxn == StaticSparseGraph::Permanent*/) {
//                         std::cout<<"Observed edge ("<<x<<","<<xn<<") was pruned edge with weight "<<wxn<<std::endl;
//                         continue;
//                     }
//                     // if new cluster is "discovered", add edge to representatives, otherwise to removedEdges
//                     if (weightToCluster.find(cxn) == weightToCluster.end()) {
//                         if (verbosity >= 5) {
//                             std::cout<<"Edge ("<<x<<","<<xn<<") connects to new cluster "<<cxn<<std::endl;
//                         }
//                         weightToCluster[cxn] = wxn;
//                         representatives.push_back(ex);
//                         if (wxn > 0)
//                             positiveWeightToCluster[cxn] = wxn;
//                         else
//                             negativeWeightToCluster[cxn] = -wxn;
//                     } else {
//                         if (verbosity >= 5) {
//                             std::cout<<"Edge ("<<x<<","<<xn<<") connects to already discovered cluster "<<cxn<<std::endl;
//                         }
//                         weightToCluster[cxn] += wxn;
//                         removedEdges.push_back(ex);
//                         if (wxn > 0)
//                             positiveWeightToCluster[cxn] += wxn;
//                         else
//                             negativeWeightToCluster[cxn] -= wxn;
//                     }
//                 }
//             }
//             
//             // remove all unused edges
//             for (Edge e : removedEdges) {
//                 setZero(e);
//                 edgeHeap.removeEdge(e);
//             }
//             
//             // accumulate weights and induced costs on representative edges
//             for (Edge e : representatives) {
//                 NodeId u = e.u;
//                 NodeId v = e.v;
//                 NodeId cxn = graph.getCliqueIdOf(u) == cu ? graph.getCliqueIdOf(v) : graph.getCliqueIdOf(u);
//                 if (std::isnan(weightToCluster[cxn]))
//                     std::cout<<"Permanent and forbidden edge was added together. Confict triple!"<<std::endl;
//                 else if (verbosity >= 4)
//                     std::cout<<"Changing weight of ("<<e.u<<","<<e.v<<") to "<<weightToCluster[cxn]<<std::endl;
//                 setToWeight(e, weightToCluster[cxn]);
//                 // recompute induced cost for this edge
//                 EdgeWeight icf = positiveWeightToCluster[cxn];
//                 EdgeWeight icp = negativeWeightToCluster[cxn];
//                 for (NodeId w : graph.getUnprunedNeighbours(u)) {
//                     if (w == v)
//                         continue;
//                     icp += edgeHeap.getIcp(graph.getWeight(Edge(u,w)), graph.getWeight(Edge(v,w)));
//                     icf += edgeHeap.getIcf(graph.getWeight(Edge(u,w)), graph.getWeight(Edge(v,w)));
//                 }                
//                 
//                 if (graph.findIndex(e) != 0) {
// //                     std::cout<<"Reset icf ("<<e.id()<<") with rank id "<<graph.findIndex(e)<<std::endl;
//                     edgeHeap.increaseIcf(e, icf - edgeHeap.getIcf(e));
// //                     std::cout<<"Reset icp ("<<e.id()<<") with rank id "<<graph.findIndex(e)<<std::endl;
//                     edgeHeap.increaseIcp(e, icp - edgeHeap.getIcp(e));
//                 } else  {
//                     std::cout<<"Zero edge was representative edge!"<<std::endl;
//                 }
//             }
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
                    if (x == y || /*graph.findIndex(e) == 0 ||*/ graph.getWeight(e) == StaticSparseGraph::Forbidden || (x == eIcp.u && y == eIcp.v)) {
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
                if (graph.getWeight(Edge(u, v)) == StaticSparseGraph::Permanent) {
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
                    if (w == StaticSparseGraph::Forbidden)
                        return false;
                    else if (w != StaticSparseGraph::Permanent) {
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
                    if (graph.getWeight(Edge(u, v)) == StaticSparseGraph::Forbidden) {
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
                        if (graph.getWeight(e) != StaticSparseGraph::Forbidden) {
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
    EdgeWeight uv = id == 0 ? 0.0 : graph.getWeight(e);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        updateTripleForbiddenUW(uv, uw, graph.getWeight(vw));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
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
    EdgeWeight uv = id == 0 ? 0.0 : graph.getWeight(e);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        updateTriplePermanentUW(uv, uw, graph.getWeight(vw));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        updateTriplePermanentUW(uv, vw, graph.getWeight(uw));
    }
    if (uv < 0) {
        totalCost -= uv;
        graph.setPermanent(e);
    } else if (uv > 0) {
        graph.setPermanent(e);
    }
}

void InducedCostHeuristic::setZero(const Edge e) {
    // this has to be called to update ic, even if the edge already is permanent
    NodeId u = e.u;
    NodeId v = e.v;
    RankId id = graph.findIndex(e);
    
    // If edge is already unpresent in graph, there is nothing to do
    if (id == 0)
        return;
    
    EdgeWeight uv = graph.getWeight(e);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        updateTripleZeroUW(uv, uw, graph.getWeight(vw));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        updateTripleZeroUW(uv, vw, graph.getWeight(uw));
    }
    graph.removeEdge(e);
}

void InducedCostHeuristic::setToWeight(const Edge e, const EdgeWeight newWeight) {
    // this has to be called to update ic, even if the edge already is permanent
    NodeId u = e.u;
    NodeId v = e.v;
    RankId id = graph.findIndex(e);
    
    /* If the edge was a zero edge in the original graph, it might have been implicitly set to
     * permanent or forbidden without updating the ic. Therefore we assume the weight to be 0 here.*/
    EdgeWeight uv = id == 0 ? 0.0 : graph.getWeight(e);
    
    for (NodeId w : graph.getUnprunedNeighbours(u)) {
        if (w == v)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        updateTripleCustomWeightUW(uv, newWeight, uw, graph.getWeight(vw));
    }
    for (NodeId w : graph.getUnprunedNeighbours(v)) {
        if (w == u)
            continue;
        Edge uw(u, w);
        Edge vw(v, w);
        updateTripleCustomWeightUW(uv, newWeight, vw, graph.getWeight(uw));
    }
    graph.setWeight(e, newWeight);
}

void InducedCostHeuristic::updateTripleForbiddenUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
    EdgeWeight icf_old = edgeHeap.getIcf(uv, vw);
    EdgeWeight icf_new = 0.0;
    EdgeWeight icp_old = edgeHeap.getIcp(uv, vw);
    EdgeWeight icp_new = std::max(0.0, vw);
    if (uv == StaticSparseGraph::Forbidden)
        std::cout<<"updateTripleForbidden on already forbidden edge"<<std::endl;
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
    if (uv == StaticSparseGraph::Permanent)
        std::cout<<"updateTriplePermanent on already permanent edge"<<std::endl;
    if (icf_new != icf_old)
        edgeHeap.increaseIcf(uw, icf_new - icf_old);
    if (icp_new != icp_old)
        edgeHeap.increaseIcp(uw, icp_new - icp_old);
}

void InducedCostHeuristic::updateTripleZeroUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
    EdgeWeight icf_old = edgeHeap.getIcf(uv, vw);
    EdgeWeight icp_old = edgeHeap.getIcp(uv, vw);
    edgeHeap.increaseIcf(uw, - icf_old);
    edgeHeap.increaseIcp(uw, - icp_old);
}

void InducedCostHeuristic::updateTripleCustomWeightUW(const EdgeWeight uv, const EdgeWeight uv_new, const Edge uw, const EdgeWeight vw) {
    EdgeWeight icf_old = edgeHeap.getIcf(uv, vw);
    EdgeWeight icf_new = edgeHeap.getIcf(uv_new, vw);
    EdgeWeight icp_old = edgeHeap.getIcp(uv, vw);
    EdgeWeight icp_new = edgeHeap.getIcp(uv_new, vw);
    if (icf_new != icf_old)
        edgeHeap.increaseIcf(uw, icf_new - icf_old);
    if (icp_new != icp_old)
        edgeHeap.increaseIcp(uw, icp_new - icp_old);
}

} // namespace ysk
