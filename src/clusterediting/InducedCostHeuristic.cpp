#include "InducedCostHeuristic.h"

namespace ysk {
  
using Edge = LightCompleteGraph::Edge;
using EdgeWeight = LightCompleteGraph::EdgeWeight;
using EdgeId = LightCompleteGraph::EdgeId;
using NodeId = LightCompleteGraph::NodeId;

InducedCostHeuristic::InducedCostHeuristic(LightCompleteGraph& param_graph, bool param_pruneZeroEdges) :
  pruneZeroEdges(param_pruneZeroEdges),
  graph(param_graph),
  edgeHeap(graph, param_pruneZeroEdges),
  totalCost(0.0)
{
}

ClusterEditingSolutionLight InducedCostHeuristic::solve() {
	// execute algorithm
	if (verbosity >= 1)
		std::cout<<"Running heuristic." << "."<<std::endl;
  
	int totalEdges = edgeHeap.numUnprocessed();
	for (unsigned int i = 0; i < graph.numEdges() + 1; i++) {
		Edge eIcf = edgeHeap.getMaxIcfEdge();
		Edge eIcp = edgeHeap.getMaxIcpEdge();
		if (eIcf == LightCompleteGraph::InvalidEdge || eIcp == LightCompleteGraph::InvalidEdge) {
		break;
		}
		EdgeWeight mIcf = edgeHeap.getIcf(eIcf);
		EdgeWeight mIcp = edgeHeap.getIcp(eIcp);
		if (mIcf >= mIcp) {
			// set eIcf to permanent
			setPermanent(eIcf);
			edgeHeap.removeEdge(eIcf);

			// resolve implications
			std::vector<NodeId> uClique(graph.getCliqueOf(eIcf.u));
			std::vector<NodeId> vClique(graph.getCliqueOf(eIcf.v));
			for (NodeId x : uClique) {
				for (NodeId y : vClique) {
					if (x == y)
						continue;
					Edge e = Edge(x,y);
					setPermanent(e);
					edgeHeap.removeEdge(e);
				}
			}
		} else {
			// set eIcp fo forbidden
			setForbidden(eIcp);
			edgeHeap.removeEdge(eIcp);
			
			// resolve implications
			std::vector<NodeId> uClique(graph.getCliqueOf(eIcp.u));
			std::vector<NodeId> vClique(graph.getCliqueOf(eIcp.v));
			for (NodeId x : uClique) {
				for (NodeId y : vClique) {
					if (x == y)
						continue;
					Edge e = Edge(x,y);
					setForbidden(e);
					edgeHeap.removeEdge(e);
				}
			}
		}
		if (verbosity >= 1 && i % 100 == 0) {
			std::cout<<"Completed "<<((totalEdges - edgeHeap.numUnprocessed())*100 / totalEdges)<<"%\r"<<std::flush;
		}
	}
  
	// calculate clustering
	if (verbosity >= 1)
		std::cout<<"Constructing result."<<std::endl;
	std::vector<std::vector<NodeId>> clusters;
	std::vector<int> clusterOfNode(graph.numNodes(), -1);
	for (NodeId u = 0; u < graph.numNodes(); u++) {
		// add cluster if not explored yet
		std::cout<<"Completed "<<(u*100/graph.numNodes())<<"%\r"<<std::flush;
		if (clusterOfNode[u] == -1) {
			int c = clusters.size();
			clusterOfNode[u] = c;
			clusters.push_back(std::vector<NodeId>(1, u));
			for (NodeId v = u+1; v < graph.numNodes(); v++) {
				if (graph.getWeight(Edge(u, v)) == LightCompleteGraph::Permanent) {
				clusterOfNode[v] = c;
				clusters[c].push_back(v);
				}
			}
		}
	}
	return ClusterEditingSolutionLight(totalCost, clusters);
}

void InducedCostHeuristic::setForbidden(const Edge e) {
	if (graph.getWeight(e) == LightCompleteGraph::Forbidden) {
		return;
	}
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
	graph.setWeight(e, LightCompleteGraph::Forbidden);
}

void InducedCostHeuristic::setPermanent(const Edge e) {
	if (graph.getWeight(e) == LightCompleteGraph::Permanent) {
		return;
	}
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
	graph.setWeight(e, LightCompleteGraph::Permanent);
}

void InducedCostHeuristic::updateTripleForbiddenUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
	EdgeWeight icf_old = edgeHeap.getIcf(uv, vw);
	EdgeWeight icf_new = 0.0;
	EdgeWeight icp_old = edgeHeap.getIcp(uv, vw);
	EdgeWeight icp_new = std::max(0.0, vw);
	edgeHeap.increaseIcf(uw, icf_new - icf_old);
	edgeHeap.increaseIcp(uw, icp_new - icp_old);
}

void InducedCostHeuristic::updateTriplePermanentUW(const EdgeWeight uv, const Edge uw, const EdgeWeight vw) {
	EdgeWeight icf_old = edgeHeap.getIcf(uv, vw);
	EdgeWeight icf_new = std::max(0.0, vw);
	EdgeWeight icp_old = edgeHeap.getIcp(uv, vw);
	EdgeWeight icp_new = std::max(0.0, -vw);
	edgeHeap.increaseIcf(uw, icf_new - icf_old);
	edgeHeap.increaseIcp(uw, icp_new - icp_old);
}

} // namespace ysk
