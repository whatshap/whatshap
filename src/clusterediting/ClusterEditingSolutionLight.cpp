#include "ClusterEditingSolutionLight.h"

namespace ysk {
  
using NodeId = LightCompleteGraph::NodeId;

ClusterEditingSolutionLight::ClusterEditingSolutionLight() :
    valid(false),
    totalCost(0.0),
    clusters(0)
{}

ClusterEditingSolutionLight::ClusterEditingSolutionLight(double pTotalCost, std::vector<std::vector<NodeId>> &pClusters) :
    valid(true),
    totalCost(pTotalCost),
    clusters(pClusters)
{}

bool ClusterEditingSolutionLight::isValid() const {
    return valid;
}

double ClusterEditingSolutionLight::getTotalCost() const {
  return totalCost;
}

unsigned int ClusterEditingSolutionLight::getNumClusters() const {
  return clusters.size();
}

const std::vector<NodeId>& ClusterEditingSolutionLight::getCluster(const unsigned int index) const {
  return clusters[index];
}

} //namespace ysk
