#ifndef HAPLOTHREADER_H
#define HAPLOTHREADER_H

#include <cstdint>
#include "tuple.h"
#include "tupleconverter.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <iostream>
#include <string> 
#include <sstream>
#include <algorithm>
#include <limits>
    
typedef uint32_t Position;
typedef float ThreadScore;

struct TupleEntry {
    ThreadScore score;
    ClusterTuple pred;
    
    TupleEntry() :
    score(std::numeric_limits<ThreadScore>::infinity()),
    pred(ClusterTuple::INVALID_TUPLE) {}
    
    TupleEntry(const ThreadScore score, const ClusterTuple pred) :
    score(score),
    pred(pred){}
    
    bool operator==(const TupleEntry& other) const {
        return pred == other.pred && score == other.score;
    }
    
    bool operator!=(const TupleEntry& other) const {
        return pred != other.pred || score != other.score;
    }
};

typedef std::pair<ClusterTuple, TupleEntry> ColumnEntry;

/**
 * ===
 * 
 * A class to compute a threading of haplotypes through a set of clusters.
 * 
 * ==
 */

class HaploThreader {

public:
    /**
     * Sets up an instance with fixed ploidy, switch costs and affine switch costs.
     * 
     * @param ploidy The number of paths, which have to be threaded through the clusters
     * @param switchCost The factor how much a single cluster switches is penalized over a wrong copy number of a cluster (compared to its coverage)
     * @param affineSwitchCost Penalty for a position, in which a cluster switch occurs
     * @param maxClusterGap Maximum allowed number of consecutive zero-coverage positions until cluster is not considered anymore
     * @param rowLimit Keeps at most this number of cluster tuples as candidates for each position. 0 means no limit
     */
    HaploThreader (uint32_t ploidy, double switchCost, double affineSwitchCost, uint32_t maxClusterGap, uint32_t rowLimit);
    
    /**
     * Computes a number of paths (depending on the provided ploidy), which run through the provided clusters. For each variant the result
     * contains a tuple of cluster ids, which represent the paths.
     * 
     * @param blockStarts A list of positions, from which the phasing runs have to start
     * @param covMap A vector, which for every position contains the global cluster ids of the present clusters
     * @param coverage A vector, which for every position contains the relative coverage of all local clusters at this position
     * @param consensus For every position contains the consensus of all local clusters at this position
     * @param genotypes The genotype for every positon
     * 
     */
    std::vector<std::vector<GlobalClusterId>> computePaths (const std::vector<Position>& blockStarts,
                    const std::vector<std::vector<GlobalClusterId>>& covMap,
                    const std::vector<std::unordered_map<GlobalClusterId, std::unordered_map<uint32_t, uint32_t>>>& alleleDepths
                   ) const;
                  
    /**
     * Computes a number of paths (depending on the provided ploidy), which run through the provided clusters. For each variant the result
     * contains a tuple of cluster ids, which represent the paths.
     * 
     * @param start First position to phase
     * @param end First position not to phase
     * @param covMap A vector, which for every position contains the global cluster ids of the present clusters
     * @param coverage A vector, which for every position contains the relative coverage of all local clusters at this position
     * @param consensus For every position contains the consensus of all local clusters at this position
     * @param genotypes The genotype for every positon
     */
    std::vector<std::vector<GlobalClusterId>> computePaths (Position start, Position end,
                    const std::vector<std::vector<GlobalClusterId>>& covMap,
                    const std::vector<std::unordered_map<GlobalClusterId, std::unordered_map<uint32_t, uint32_t>>>& alleleDepths,
                    Position displayedEnd = 0
                   ) const;

private:
    uint32_t ploidy;
    double switchCost;
    double affineSwitchCost;
    uint32_t maxClusterGap;
    uint32_t rowLimit;
    
    /**
     * Computes the coverage cost of a tuple, considering the following coverage distribution. All cluster
     * indices are local.
     */
    ThreadScore getCoverageCost(ClusterTuple tuple,
                          const uint32_t coverage,
                          const std::vector<uint32_t>& clusterCoverage) const;
    
    /**
     * Computes the switch cost between one tuple and all permutations of another tuple. The tuples must have global cluster ids,
     * by which they must be sorted in ascending order.
     */
    ThreadScore getSwitchCostAllPerms(const std::vector<GlobalClusterId>& prevTuple, const std::vector<GlobalClusterId>& curTuple) const;

    std::vector<ClusterTuple> computeRelevantTuples (const std::vector<std::vector<uint32_t>>& clusterCoverage, Position pos) const;

    /**
     * Computes the position-wise coverage (total and per-cluster) based on the allele counts for each cluster and position.
     */
    void computeCoverage(const std::vector<std::unordered_map<GlobalClusterId, std::unordered_map<uint32_t, uint32_t>>>& alleleDepths,
                         const std::vector<std::vector<GlobalClusterId>>& covMap,
                         std::vector<uint32_t>& coverage,
                         std::vector<std::vector<uint32_t>>& clusterCoverage) const;
                         
    void printStartOfColumn(Position pos,
                            const std::vector<ClusterTuple>& relevantTuples,
                            const std::vector<GlobalClusterId>& gIds) const {
        std::cout<<"Position "<<pos<<": "<<(relevantTuples.size())<<" tuples on "<<gIds.size()<<" clusters (";
        for (GlobalClusterId gid : gIds)
            std::cout<<" "<<gid;
        std::cout<< " ). "<<std::endl;
    }
};

#endif
