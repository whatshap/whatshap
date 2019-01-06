#include "InducedCostHeuristic.h"
#include <queue>

namespace ysk {
  
using Edge = StaticSparseGraph::Edge;
using EdgeWeight = StaticSparseGraph::EdgeWeight;
using EdgeId = StaticSparseGraph::EdgeId;
using NodeId = StaticSparseGraph::NodeId;

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
	if (verbosity >= 1)
		std::cout<<"Running heuristic." << "."<<std::endl;
    if (totalCost == std::numeric_limits<EdgeWeight>::infinity()) {
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
                    if (x == y || graph.getWeight(e) == 0 || graph.getWeight(e) == StaticSparseGraph::Permanent || (x == eIcf.u && y == eIcf.v))
						continue;
                    if (verbosity >= 5) {
                        std::cout<<"Making ("<<x<<","<<y<<") permanent due to implication."<<std::endl;
                    }
                    implications.push_back(e);
// 					setPermanent(e);
// 					edgeHeap.removeEdge(e);
				}
			}
			
			// set actual edge first
			setPermanent(eIcf);
			edgeHeap.removeEdge(eIcf);
            
            // then all implications
            for (Edge e : implications) {
                setPermanent(e);
                edgeHeap.removeEdge(e);
            }
			
		} else {
			// set eIcp fo forbidden
// 			setForbidden(eIcp);
// 			edgeHeap.removeEdge(eIcp);
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
					if (x == y || graph.getWeight(e) == 0 || graph.getWeight(e) == StaticSparseGraph::Forbidden || (x == eIcp.u && y == eIcp.v))
						continue;
                    implications.push_back(e);
// 					setForbidden(e);
// 					edgeHeap.removeEdge(e);
				}
			}
			
			// set actual edge first
			setForbidden(eIcp);
			edgeHeap.removeEdge(eIcp);
            
            // then all implications
            for (Edge e : implications) {
                setForbidden(e);
                edgeHeap.removeEdge(e);
            }
		}
		if (verbosity >= 1 && i % 100 == 0) {
			std::cout<<"Completed "<<((totalEdges - edgeHeap.numUnprocessed())*100 / totalEdges)<<"%\r"<<std::flush;
		}
	}
  
	// calculate clustering
	if (verbosity >= 1)
		std::cout<<"Constructing result.."<<std::endl;
	std::vector<std::vector<NodeId>> clusters;
	std::vector<int> clusterOfNode(graph.numNodes(), -1);
	for (NodeId u = 0; u < graph.numNodes(); u++) {
		// add cluster if not explored yet
        if (verbosity >= 4) {
            std::cout<<"Processing node "<<u<<std::endl;
        }
        if (verbosity >= 1 && verbosity <= 4) {
            std::cout<<"Completed "<<(u*100/graph.numNodes())<<"%\r"<<std::flush;
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
	return ClusterEditingSolutionLight(totalCost, clusters);
}

bool InducedCostHeuristic::resolvePermanentForbidden() {
    if (verbosity >= 1)
		std::cout<<"Resolving forbidden and permanent edges." << "."<<std::endl;
    // make cliques by connecting all nodes with inf path between them
    std::vector<bool> processed(graph.numNodes(), false);
    std::vector<std::vector<NodeId>> cliques;
    for (NodeId u = 0; u < graph.numNodes() - 1; u++) {
        if (processed[u])
            continue;
        std::vector<NodeId> clique;
        std::queue<NodeId> remaining;
        remaining.push(u);
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
        for (unsigned int l = k+1; l < cliques.size(); l++) {
            // search for forbidden edge between
            bool found = false;
            for (NodeId u : cliques[k]) {
                if (found) break;
                for (NodeId v : cliques[l]) {
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
    
    return true;
}

void InducedCostHeuristic::setForbidden(const Edge e) {
//     if (graph.getWeight(e) == StaticSparseGraph::Forbidden) {
// 		return;
// 	}
	NodeId u = e.u;
	NodeId v = e.v;
	EdgeWeight uv = graph.getWeight(e);
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
	if (graph.getWeight(e) > 0)
		totalCost += graph.getWeight(e);
    graph.setForbidden(e);
}

void InducedCostHeuristic::setPermanent(const Edge e) {
//     if (graph.getWeight(e) == StaticSparseGraph::Permanent) {
// 		return;
// 	}
	NodeId u = e.u;
	NodeId v = e.v;
	EdgeWeight uv = graph.getWeight(e);
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
	if (graph.getWeight(e) < 0)
		totalCost -= graph.getWeight(e);
    graph.setPermanent(e);
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

} // namespace ysk
