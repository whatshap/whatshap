#include "clustereditingsolution.h"
  
using NodeId = StaticSparseGraph::NodeId;

ClusterEditingSolution::ClusterEditingSolution() :
    valid(false),
    totalCost(0.0),
    clusters(0)
{}

ClusterEditingSolution::ClusterEditingSolution(double pTotalCost, std::vector<std::vector<NodeId>> &pClusters) :
    valid(true),
    totalCost(pTotalCost),
    clusters(pClusters)
{}

bool ClusterEditingSolution::isValid() const {
    return valid;
}

double ClusterEditingSolution::getTotalCost() const {
  return totalCost;
}

unsigned int ClusterEditingSolution::getNumClusters() const {
  return clusters.size();
}

const std::vector<NodeId>& ClusterEditingSolution::getCluster(const unsigned int index) const {
  return clusters[index];
}
