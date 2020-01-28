#ifndef CLUSTEREDITINGSOLUTIONLIGHT_H
#define CLUSTEREDITINGSOLUTIONLIGHT_H

#include "dynamicsparsegraph.h"
#include <vector>

class ClusterEditingSolution {
public:
    ClusterEditingSolution();
    ClusterEditingSolution(double pTotalCost, std::vector<std::vector<DynamicSparseGraph::NodeId>> &pClusters);
    bool isValid() const;
    unsigned int getNumClusters() const;
    const std::vector<DynamicSparseGraph::NodeId>& getCluster(const unsigned int index) const;
    double getTotalCost() const;
  
private:
    bool valid;
    double totalCost;
    std::vector<std::vector<DynamicSparseGraph::NodeId>> clusters;
};

#endif
