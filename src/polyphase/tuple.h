#ifndef TUPLE_H
#define TUPLE_H

#include <cstdint>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>

typedef uint32_t GlobalClusterId;
typedef uint32_t LocalClusterId;
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
     * Initialize by global ids using mapping from local to global ids.
     */
    ClusterTuple(const std::vector<GlobalClusterId>& clusters, const std::vector<GlobalClusterId>& covMap) {
        std::vector<LocalClusterId> v;
        for (uint32_t i = 0; i < clusters.size(); i++) {
            LocalClusterId c = MAX_CLUSTERS_PER_COLUMN+1;
            for (LocalClusterId j = 0; j < covMap.size(); j++) {
                if (covMap[j] == clusters[i]) {
                    c = j;
                    break;
                }
            }
            if (c > MAX_CLUSTERS_PER_COLUMN) {
                tuple = -1;
                return;
            } else {
                v.push_back(c);
            }
        }
        tuple = 0;
        for (uint32_t i = 0; i < v.size(); i++) {
            tuple += (((TupleCode)v[i]) << (i*BITS_PER_CLUSTER));
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

#endif
