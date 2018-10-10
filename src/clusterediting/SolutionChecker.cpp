#include "SolutionChecker.h"

using namespace ysk;
using namespace std;
using NodeId = LightCompleteGraph::NodeId;
using Edge = LightCompleteGraph::Edge;
using EdgeWeight = LightCompleteGraph::EdgeWeight;

bool SolutionChecker::verifySolution(const LightCompleteGraph& graph, const ClusterEditingSolutionLight& solution){
	if (verbosity > 0) {
		cout << "Checking solution integrity ..." << endl;
	}
    
    EdgeWeight totalCost = 0;
    
    // check costs for inserted edges inside clusters
    for (unsigned int c = 0; c < solution.getNumClusters(); c++){
		vector<NodeId> cluster = solution.getCluster(c);
		// check pairwise edges whether they were negative in the original graph
		for (unsigned int m = 0; m < cluster.size()-1; m++) {
			for (unsigned int n = m+1; n < cluster.size(); n++) {
				EdgeWeight w = graph.getWeight(Edge(cluster[m], cluster[n]));
				if (w < 0)
					totalCost -= w;
			}
		}
	}
	
	// check costs for deleted edges between clusters
    for (unsigned int c1 = 0; c1 < solution.getNumClusters(); c1++){
		vector<NodeId> cluster1 = solution.getCluster(c1);
		for (unsigned int c2 = c1+1; c2 < solution.getNumClusters(); c2++){
			vector<NodeId> cluster2 = solution.getCluster(c2);
			// check pairwise edges whether they were positive in the original graph
			for (unsigned int m = 0; m < cluster1.size(); m++) {
				for (unsigned int n = 0; n < cluster2.size(); n++) {
					EdgeWeight w = graph.getWeight(Edge(cluster1[m], cluster2[n]));
					if (w > 0)
						totalCost += w;
				}
			}
		}
    }
    
    if (abs(solution.getTotalCost() - totalCost) > MARGIN_OF_ERROR){
		cout << "Could not verify the editing costs!" << endl;
		return false;
    } else if (verbosity > 0) {
		cout << "Solution validity verified!" << endl;
	}
    return true;
}
