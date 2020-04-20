#ifndef HAPLOTHREADER_H
#define HAPLOTHREADER_H

#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <string> 
#include <sstream>
#include <algorithm>
    
typedef uint32_t GlobalClusterId;
typedef uint32_t LocalClusterId;
typedef uint32_t Position;
typedef double Score;
typedef uint64_t TupleCode;

/**
 * Struct to represent cluster tuples. Each tuple is encoded by a 64bit unsigned integer. The first cluster is encoded by the lowest 5 bits, the second
 * cluster by the next highest 5 bits, etc. There is space for 12 clusters, each ranging from id 0 to 31. Each cluster is encoded as a local cluster id. 
 * Multiple ClusterTuples of the same variant/position are comparable, while they are not for different variants/positions. The local ids must be mapped 
 * back to global cluster ids. The ploidy of the tuple is not saved inside the tuple and must be provided for several methods.
 */

struct ClusterTuple {
    const static uint32_t BITS_PER_CLUSTER = 5;
    const static uint32_t MAX_CLUSTERS_PER_COLUMN = 1 << BITS_PER_CLUSTER; // must be 2^BITS_PER_CLUSTER
    const static uint32_t MAX_PLOIDY = 64/BITS_PER_CLUSTER;
    static constexpr uint64_t TUPLE_MASKS[12] = {31UL << 0, 31UL << 5, 31UL << 10, 31UL << 15, 
        31UL << 20, 31UL << 25, 31UL << 30, 31UL << 35, 
        31UL << 40, 31UL << 45, 31UL << 50, 31UL << 55};
    const static ClusterTuple INVALID_TUPLE;
    
    TupleCode tuple;
    
    /**
     * Initializes as invalid cluster
     */
    ClusterTuple() {
        tuple = -1;
    }
    
    /**
     * Initializes using a vector of local cluster ids.
     */
    ClusterTuple(const std::vector<LocalClusterId>& clusters) {
        tuple = 0;
        for (uint32_t i = 0; i < clusters.size(); i++) {
            tuple += (((TupleCode)clusters[i]) << (i*BITS_PER_CLUSTER));
        }
    }
    
    /**
     * Constructs a tuple from a numeric representation.
     */
    ClusterTuple(TupleCode other) {
        tuple = other;
    }
    
    /**
     * Sets the cluster at the given index.
     */
    void set(LocalClusterId c, uint32_t index) {
        tuple |= TUPLE_MASKS[index];
        tuple -= TUPLE_MASKS[index];
        tuple += (((TupleCode)c) << index*BITS_PER_CLUSTER);
    }
    
    /**
     * Rearranges the elements of the tuple using the provided permutation vector. E.g. the cluster at index i will be
     * set to the value that is at index perm[i] right now.
     */
    void permute(std::vector<uint32_t> perm) {
        uint64_t newTuple = 0;
        for (uint32_t i = 0; i < perm.size(); i++) {
            newTuple += ((uint64_t)get(perm[i])) << i*BITS_PER_CLUSTER;
        }
        tuple = newTuple;
    }
    
    /**
     * Counts how often a certain cluster id is present in the tuple.
     */
    uint32_t count(LocalClusterId id, uint32_t ploidy) const {
        uint32_t c = 0;
        for (uint32_t i = 0; i < ploidy; i++) {
            c += get(i) == id;
        }
        return c;
    }
    
    /**
     * Returns the cluster at the given index.
     */
    LocalClusterId get(uint32_t i) const {
        if (i >= MAX_PLOIDY) {
            std::cout<<"Accessed element "<<i<<" > MAX_PLOIDY of a tuple!"<<std::endl;
        }
        return (tuple & TUPLE_MASKS[i]) >> (i*BITS_PER_CLUSTER);
    }
    
    /**
     * Returns the numeric representation of the tuple.
     */
    TupleCode asNumber() const {
        return tuple;
    }
    
    TupleCode fingerprint(uint32_t ploidy) const {
        std::vector<LocalClusterId> clusters;
        for (uint32_t i = 0; i < ploidy; i++) {
            clusters.push_back(get(i));
        }
        std::sort(clusters.begin(), clusters.end());
        TupleCode fp = 0;
        for (uint32_t i = 0; i < ploidy; i++) {
            fp += (((TupleCode)clusters[i]) << (i*BITS_PER_CLUSTER));
        }
        return fp;
    }
    
    /**
     * Converts the tuple into a vector with global cluster ids. Each local cluster id is mapped to a global id using the provided vector.
     */
    std::vector<GlobalClusterId> asVector(uint32_t ploidy, std::vector<GlobalClusterId> globalIds) const {
        std::vector<LocalClusterId> clusters;
        for (uint32_t i = 0; i < ploidy; i++) {
            if (get(i) >= globalIds.size()) {
                std::cout<<"Stored local cluster id was higher than size of global id vector! "<<i<<" "<<get(i)<<" "<<globalIds.size()<<std::endl;
                std::vector<GlobalClusterId> empty;
                return empty;
            } else {
                clusters.push_back(globalIds[get(i)]);
            }
        }
        return clusters;
    }
    
    /**
     * Returns a strin representation of the tuple.
     */
    std::string asString(uint32_t ploidy, std::vector<GlobalClusterId> globalIds) const {
        std::stringstream s;
        s<<"[";
        for (uint32_t i = 0; i < ploidy; i++) {
            s<<globalIds[get(i)];
            if (i < ploidy - 1) {
                s<<", ";
            } 
        }
        s<<"]";
        return s.str();
    }
    
    /**
     * Returns a strin representation of the tuple.
     */
    std::string asString(uint32_t ploidy) const {
        std::stringstream s;
        s<<"[";
        for (uint32_t i = 0; i < ploidy; i++) {
            s<<get(i);
            if (i < ploidy - 1) {
                s<<", ";
            } 
        }
        s<<"]";
        return s.str();
    }
    
    bool operator==(const ClusterTuple& other) const {
        return tuple == other.tuple;
    }
    
    bool operator!=(const ClusterTuple& other) const {
        return tuple != other.tuple;
    }
};

namespace std {
    template <>
    struct hash<ClusterTuple> {
        size_t operator()(const ClusterTuple& t) const {
            return hash<TupleCode>()(t.tuple);
        }
    };
}

struct ClusterEntry {
    Score score;
    ClusterTuple pred;
    
    ClusterEntry() :
    score(std::numeric_limits<Score>::infinity()),
    pred(ClusterTuple::INVALID_TUPLE) {}
    
    ClusterEntry(const Score score, const ClusterTuple pred) :
    score(score),
    pred(pred){}
    
    bool operator==(const ClusterEntry& other) const {
        return pred == other.pred && score == other.score;
    }
    
    bool operator!=(const ClusterEntry& other) const {
        return pred != other.pred || score != other.score;
    }
};

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
     */
    HaploThreader (uint32_t ploidy, double switchCost, double affineSwitchCost, bool symmetryOptimization, uint32_t rowLimit);
    
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
                    const std::vector<std::vector<double>>& coverage, 
                    const std::vector<std::vector<uint32_t>>& consensus,
                    const std::vector<std::unordered_map<uint32_t, uint32_t>>& genotypes
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
                    const std::vector<std::vector<double>>& coverage, 
                    const std::vector<std::vector<uint32_t>>& consensus,
                    const std::vector<std::unordered_map<uint32_t, uint32_t>>& genotypes,
                    Position displayedEnd = 0
                   ) const;

private:
    uint32_t ploidy;
    double switchCost;
    double affineSwitchCost;
    bool symmetryOptimization;
    uint32_t rowLimit;
    
    /**
     * Computes the coverage cost of a tuple, considering the following coverage distribution. All cluster
     * indices are local.
     */
    Score getCoverageCost(ClusterTuple tuple, const std::vector<double>& coverage) const;

    /**
     * Computes the switch cost between one tuple and all permutations of another tuple. The tuples must have global cluster ids,
     * by which they must be sorted in ascending order.
     */
    Score getSwitchCostAllPerms(const std::vector<GlobalClusterId>& prevTuple, const std::vector<GlobalClusterId>& curTuple) const;
                        
    /**
     * Computes the switch cost between one tuple and all permutations of another tuple. The tuples must have global cluster ids,
     * by which they must be sorted in ascending order. Writes residual positions (positions in either tuple, which could not be
     * matched to the other one) into the provided vectors.
     */
    Score getSwitchCostAllPerms(const std::vector<GlobalClusterId>& prevTuple, const std::vector<GlobalClusterId>& curTuple,
                                std::vector<uint32_t>& residualPosPrev, std::vector<uint32_t>& residualPosCur) const;
    
    
    std::vector<ClusterTuple> computeGenotypeConformTuples (const std::vector<GlobalClusterId>& covMap,
                                                            const std::vector<uint32_t>& consensus, 
                                                            const std::unordered_map<uint32_t, uint32_t>& genotype) const;
                                                            
    /**
     * Computes all sets of clusters having are have a genotype, which is exactly distance away from the specified genotype.
     * Result is returned as a vector of tuples over the local index space, i.e. from 0 to size of clusters - 1. No tuple
     * is a permutation of another one (they have to be extended if required).
     */
    std::vector<ClusterTuple> getGenotypeConformTuples (const std::vector<GlobalClusterId>& clusters, const std::vector<uint32_t>& consensus, 
                                                        const std::unordered_map<uint32_t, uint32_t>& genotype) const;
};

#endif
