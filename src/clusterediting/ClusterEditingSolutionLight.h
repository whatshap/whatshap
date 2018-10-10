#ifndef CLUSTEREDITINGSOLUTIONLIGHT_H
#define CLUSTEREDITINGSOLUTIONLIGHT_H

#include "LightCompleteGraph.h"
#include <vector>

namespace ysk {

class ClusterEditingSolutionLight 
{
public:
    ClusterEditingSolutionLight();
    ClusterEditingSolutionLight(double pTotalCost, std::vector<std::vector<LightCompleteGraph::NodeId>> &pClusters);
    bool isValid() const;
    unsigned int getNumClusters() const;
    const std::vector<LightCompleteGraph::NodeId>& getCluster(const unsigned int index) const;
    double getTotalCost() const;
  
private:
    bool valid;
    double totalCost;
    std::vector<std::vector<LightCompleteGraph::NodeId>> clusters;
};

} // namespace ysk

#endif // CLUSTEREDITINGSOLUTIONLIGHT_H
