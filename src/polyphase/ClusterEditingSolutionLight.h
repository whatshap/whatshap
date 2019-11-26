#ifndef CLUSTEREDITINGSOLUTIONLIGHT_H
#define CLUSTEREDITINGSOLUTIONLIGHT_H

#include "DynamicSparseGraph.h"
#include <vector>

class ClusterEditingSolutionLight 
{
public:
    ClusterEditingSolutionLight();
    ClusterEditingSolutionLight(double pTotalCost, std::vector<std::vector<DynamicSparseGraph::NodeId>> &pClusters);
    bool isValid() const;
    unsigned int getNumClusters() const;
    const std::vector<DynamicSparseGraph::NodeId>& getCluster(const unsigned int index) const;
    double getTotalCost() const;
  
private:
    bool valid;
    double totalCost;
    std::vector<std::vector<DynamicSparseGraph::NodeId>> clusters;
};

#endif // CLUSTEREDITINGSOLUTIONLIGHT_H
