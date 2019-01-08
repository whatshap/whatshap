#ifndef CLUSTEREDITINGSOLUTIONLIGHT_H
#define CLUSTEREDITINGSOLUTIONLIGHT_H

#include "StaticSparseGraph.h"
#include <vector>

namespace ysk {

class ClusterEditingSolutionLight 
{
public:
    ClusterEditingSolutionLight();
    ClusterEditingSolutionLight(double pTotalCost, std::vector<std::vector<StaticSparseGraph::NodeId>> &pClusters);
    bool isValid() const;
    unsigned int getNumClusters() const;
    const std::vector<StaticSparseGraph::NodeId>& getCluster(const unsigned int index) const;
    double getTotalCost() const;
  
private:
    bool valid;
    double totalCost;
    std::vector<std::vector<StaticSparseGraph::NodeId>> clusters;
};

} // namespace ysk

#endif // CLUSTEREDITINGSOLUTIONLIGHT_H
