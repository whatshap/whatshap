#ifndef SWITCHFLIP_H
#define SWITCHFLIP_H

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <string> 
#include <sstream>
#include <algorithm>
#include <limits>

typedef uint32_t Position;
typedef uint32_t Haplotype;
typedef double Score;
typedef uint64_t PermutationCode;

/**
 * Struct to represent cluster tuples. Each tuple is encoded by a 64bit unsigned integer. The first cluster is encoded by the lowest 5 bits, the second
 * cluster by the next highest 5 bits, etc. There is space for 12 clusters, each ranging from id 0 to 31. Each cluster is encoded as a local cluster id. 
 * Multiple ClusterTuples of the same variant/position are comparable, while they are not for different variants/positions. The local ids must be mapped 
 * back to global cluster ids. The ploidy of the tuple is not saved inside the tuple and must be provided for several methods.
 */

struct Permutation {
    const static uint32_t BITS_PER_BLOCK = 4;
    const static uint32_t MAX_BLOCK_VALUE = 1 << BITS_PER_BLOCK; // must be 2^BITS_PER_CLUSTER
    const static uint32_t MAX_PLOIDY = 63/BITS_PER_BLOCK;
    static constexpr uint64_t TUPLE_MASKS[16] = {15UL << 0, 15UL << 4, 15UL << 8, 15UL << 12, 
        15UL << 16, 15UL << 20, 15UL << 24, 15UL << 28, 
        15UL << 32, 15UL << 36, 15UL << 40, 15UL << 44,
        15UL << 48, 15UL << 52, 15UL << 56, 15UL << 60
    };
    const static Permutation INVALID;
    
    PermutationCode code;
    
    /**
     * Initializes as invalid cluster
     */
    Permutation() {
        code = -1;
    }
    
    /**
     * Initializes using a vector of local cluster ids.
     */
    Permutation(Haplotype hapOrdering[], uint32_t ploidy) {
        code = 0;
        for (uint32_t i = 0; i < ploidy; i++) {
            code += (((PermutationCode)hapOrdering[i]) << (i*BITS_PER_BLOCK));
        }
    }
    
    /**
     * Initializes using a vector of local cluster ids.
     */
    Permutation(const std::vector<Haplotype>& hapOrdering) {
        code = 0;
        for (uint32_t i = 0; i < hapOrdering.size(); i++) {
            code += (((PermutationCode)hapOrdering[i]) << (i*BITS_PER_BLOCK));
        }
    }
    
    /**
     * Copy constructor
     */
    Permutation(const Permutation& other) {
        code = other.code;
    }
    
    /**
     * Constructs a tuple from a numeric representation.
     */
    Permutation(PermutationCode other) {
        code = other;
    }
    
    /**
     * Sets the cluster at the given index.
     */
    void set(Haplotype h, uint32_t index) {
        code |= TUPLE_MASKS[index];
        code -= TUPLE_MASKS[index];
        code += (((PermutationCode)h) << index*BITS_PER_BLOCK);
    }
    
    /**
     * Returns the cluster at the given index.
     */
    Haplotype get(uint32_t i) const {
        if (i >= MAX_PLOIDY) {
            std::cout<<"Accessed element "<<i<<" > MAX_PLOIDY of a tuple!"<<std::endl;
        }
        return (code & TUPLE_MASKS[i]) >> (i*BITS_PER_BLOCK);
    }
    
    /**
     * Returns the numeric representation of the tuple.
     */
    PermutationCode asNumber() const {
        return code;
    }
    
    /**
     * Converts the tuple into a vector with global cluster ids. Each local cluster id is mapped to a global id using the provided vector.
     */
    std::vector<Haplotype> asVector(uint32_t ploidy) const {
        std::vector<Haplotype> permutation;
        for (uint32_t i = 0; i < ploidy; i++) {
            permutation.push_back(get(i));
        }
        return permutation;
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
    
    bool operator==(const Permutation& other) const {
        return code == other.code;
    }
    
    bool operator!=(const Permutation& other) const {
        return code != other.code;
    }
};

namespace std {
    template <>
    struct hash<Permutation> {
        size_t operator()(const Permutation& t) const {
            return hash<PermutationCode>()(t.code);
        }
    };
}

struct PermutationEntry {
    Score score;
    Permutation pred;
    
    PermutationEntry() :
    score(std::numeric_limits<Score>::infinity()), pred(Permutation::INVALID) {}
    
    PermutationEntry(const Score score, const Permutation pred) :
    score(score),
    pred(pred){}
    
    PermutationEntry(const PermutationEntry& other) {
        pred = other.pred;
        score = other.score;
    }
    
    bool operator==(const PermutationEntry& other) const {
        return pred == other.pred && score == other.score;
    }
    
    bool operator!=(const PermutationEntry& other) const {
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

class SwitchFlipCalculator {

public:
    
    /**
     * Sets up an instance with fixed ploidy, switch costs and affine switch costs.
     * 
     * @param ploidy The number of paths, which have to be threaded through the clusters
     * @param switchCost The factor how much a single cluster switches is penalized over a wrong copy number of a cluster (compared to its coverage)
     * @param flipCost The factor how much a single cluster switches is penalized over a wrong copy number of a cluster (compared to its coverage)
     * @param affineSwitchCost Penalty for a position, in which a cluster switch occurs
     */
    SwitchFlipCalculator (uint32_t ploidy, double switchCost, double flipCost);
    
    /**
     * Computes a number of paths (depending on the provided ploidy), which run through the provided clusters. For each variant the result
     * contains a tuple of cluster ids, which represent the paths.
     * 
     * @param phasing0 A list of positions, from which the phasing runs have to start
     * @param phasing1
     * 
     */
    std::pair<Score, Score> compare (const std::vector<std::vector<uint32_t>>& phasing0,
                    const std::vector<std::vector<uint32_t>>& phasing1,
                    std::vector<uint32_t>& switchesInColumn,
                    std::vector<std::vector<uint32_t>>& flippedHapsInColumn,
                    std::vector<std::vector<uint32_t>>& permInColumn
                    ) const;

private:
    uint32_t ploidy;
    double switchCost;
    double flipCost;
    
    /**
     * Computes the coverage cost of a tuple, considering the following coverage distribution. All cluster
     * indices are local.
     */
    Score getNumFlips(Permutation p, const std::vector<uint32_t>& phase0, const std::vector<uint32_t>& phase1) const;
    
    /**
     * Computes the coverage cost of a tuple, considering the following coverage distribution. All cluster
     * indices are local.
     */
    std::vector<uint32_t> getFlippedHaps(Permutation p, const std::vector<uint32_t>& phase0, const std::vector<uint32_t>& phase1) const;
    
    /**
     * Computes the switch cost of the two tuples. All cluster indices are local. Uniform switch cost between clusters
     * are assumed
     */
    Score getNumSwitches(const Permutation p1, const Permutation p2) const;
    
    /**
     * Computes all permutations as tuples.
     */
    std::vector<Permutation> getPermutations() const;
};

#endif
