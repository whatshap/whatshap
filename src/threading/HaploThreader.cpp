#include "HaploThreader.h"
#include <limits>
#include <algorithm>
#include <unordered_set>

constexpr uint64_t ClusterTuple::TUPLE_MASKS[];
const ClusterTuple ClusterTuple::INVALID_TUPLE = ClusterTuple((TupleCode)-1);

HaploThreader::HaploThreader (uint32_t ploidy, double switchCost, double affineSwitchCost) :
    ploidy(ploidy),
    switchCost(switchCost),
    affineSwitchCost(affineSwitchCost)
{
    for (uint32_t p = 0; p <= ploidy; p++) {
        std::vector<std::vector<uint32_t>> pPerms;
        std::vector<uint32_t> n;
        for (uint32_t i = 0; i < p; i++)
            n.push_back(i);
        do {
            pPerms.push_back(n);
        } while (std::next_permutation(n.begin(), n.end()));
        perms.push_back(pPerms);
    }
}

std::vector<std::vector<GlobalClusterId>> HaploThreader::computePaths (Position num_vars, 
                    std::vector<std::vector<GlobalClusterId>>& covMap,
                    std::vector<std::vector<double>>& coverage, 
                    std::vector<std::vector<uint32_t>>& consensus,
                    std::vector<uint32_t>& genotypes
                   ) const {
    
    // setup data structures
    std::unordered_map<ClusterTuple, Score> pred;
    std::vector<std::unordered_map<ClusterTuple, Score>> d;
    std::vector<std::unordered_map<ClusterTuple, ClusterTuple>> b;
    Position firstUnthreadedPosition = 0;
    std::vector<std::vector<GlobalClusterId>> path;
    std::vector<std::vector<ClusterTuple>> confSets;
    std::vector<std::vector<ClusterTuple>> confSetsExt;
    
    // aux data structures for symmetry pruning
    GlobalClusterId numClusters = 0;
    for (std::vector<GlobalClusterId>& clusters : covMap) {
        for (GlobalClusterId cluster : clusters) {
            numClusters = std::max(numClusters, cluster);
        }
    }
    std::vector<int32_t> mapPos(numClusters*ploidy, -1);
    std::vector<uint32_t> numUsed(numClusters, 0);

    // initialize first column
    confSets.push_back(computeGenotypeConformTuples(covMap[0], consensus[0], genotypes[0], false));
    confSetsExt.push_back(computeGenotypeConformTuples(covMap[0], consensus[0], genotypes[0], true));
    for (ClusterTuple t : confSets[0]) {
        for (uint32_t i = 0; i < ploidy; i++) {
            if (t.get(i) >= covMap[0].size()) {
                std::cout<<"Size error in pos "<<0<<std::endl;
            }
        }
    }
    for (ClusterTuple t : confSetsExt[0]) {
        for (uint32_t i = 0; i < ploidy; i++) {
            if (t.get(i) >= covMap[0].size()) {
                std::cout<<"Size ext error in pos "<<0<<std::endl;
            }
        }
    }
    std::unordered_map<ClusterTuple, Score> scoreColumn;
    std::unordered_map<ClusterTuple, ClusterTuple> btColumn;
    if (confSetsExt[0].size() == 0) {
        std::cout<<"First variant has no clusters!"<<std::endl;
        return path;
    }
    Score minimumInColumn = std::numeric_limits<Score>::infinity();
    ClusterTuple minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
    ClusterTuple minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;
    for (ClusterTuple t : confSetsExt[0]) {
        scoreColumn[t] = getCoverageCost(t, coverage[0]);
        btColumn[t] = ClusterTuple::INVALID_TUPLE;
        firstUnthreadedPosition = 1;
        if (scoreColumn[t] < minimumInColumn) {
            minimumInColumn = scoreColumn[t];
            minimumTupleInColumn = t;
        }
    }
    d.push_back(scoreColumn);
    b.push_back(btColumn);
    
//     std::cout<<"Best score in column 0 is "<<minimumInColumn<<"("<<getCoverageCost(minimumTupleInColumn, coverage[0])<<") by "<<(minimumTupleInColumn.asString(ploidy, covMap[0]))<<std::endl;
    
    // iterate over positions
    for (Position pos = 1; pos < num_vars; pos++) {
        // compute conform tuples
        confSets.push_back(computeGenotypeConformTuples(covMap[pos], consensus[pos], genotypes[pos], false));
        confSetsExt.push_back(computeGenotypeConformTuples(covMap[pos], consensus[pos], genotypes[pos], true));
        
        // delete old tuples to save memory
        if (pos >= 2) {
            confSets[pos-2].clear();
            confSetsExt[pos-2].clear();
        }
        
        // reset variables
        scoreColumn.clear();
        btColumn.clear();
        
        Score minimum = std::numeric_limits<Score>::infinity();
        ClusterTuple minimumPred = ClusterTuple::INVALID_TUPLE;
        bool minExists = false;
        minimumInColumn = std::numeric_limits<Score>::infinity();
        minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
        minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;
        
        // iterate over rows
        std::cout<<"Threading haplotypes through clusters .. ("<<pos<<"/"<<num_vars<<")\r"<<std::flush;
        for (ClusterTuple rowTuple : confSetsExt[pos]) {
            minimum = std::numeric_limits<Score>::infinity();
            minimumPred = ClusterTuple::INVALID_TUPLE;
            
            uint32_t extSize = 0;
            if (confSetsExt[pos-1].size() > 100000*ploidy*ploidy) {
                // optimize symmetry
                // set up position map
                for (uint32_t i = 0; i < ploidy; i++) {
                    GlobalClusterId c = covMap[pos][rowTuple.get(i)];
                    mapPos[c + numClusters*numUsed[c]] = i;
                    numUsed[c]++;
                }
                // clear array after usage
                for (uint32_t i = 0; i < ploidy; i++) {
                    numUsed[covMap[pos][rowTuple.get(i)]] = 0;
                }
                
                // iterate over non-permuted previous tuples and only generate necessary permutations from them
                for (ClusterTuple predTup : confSets[pos-1]) {
                    // generate extensions
                    for (ClusterTuple extTup : getNonSymmetricExtensions(predTup, covMap[pos-1], mapPos, numUsed, numClusters)) {
                        extSize++;
                        Score s = d[pos-1][predTup] + getSwitchCost(rowTuple, extTup, covMap[pos], covMap[pos-1]);
                        if (s < minimum) {
                            minExists = true;
                            minimum = s;
                            minimumPred = extTup;
                        }
                    }
                    
                    // clear array after usage
                    for (uint32_t i = 0; i < numClusters; i++) {
                        numUsed[i] = 0;
                    }
                }
                
                // clean up after usage of arrays
                for (uint32_t i = 0; i < ploidy; i++) {
                    GlobalClusterId c = covMap[pos][rowTuple.get(i)];
                    for (uint32_t j = 0; j < ploidy; j++)
                        mapPos[c + j*numClusters] = -1;
                }
                for (uint32_t i = 0; i < ploidy; i++) {
                    for (uint32_t j = 0; j < numClusters; j++)
                        mapPos[j + i*numClusters] = -1;
                }
//                 std::cout<<"Ext: "<<confSetsExt[pos-1].size()<<"  Norm: "<<confSets[pos-1].size()<<" -> "<<extSize<<std::endl;
            } else {
                // brute force
                for (ClusterTuple predTup : confSetsExt[pos-1]) {
                    Score s = d[pos-1][predTup] + getSwitchCost(rowTuple, predTup, covMap[pos], covMap[pos-1]);
                    if (s < minimum) {
                        minExists = true;
                        minimum = s;
                        minimumPred = predTup;
                    }
                }
            }
            
            // report best recursion
            if (minExists) {
                scoreColumn[rowTuple] = minimum + getCoverageCost(rowTuple, coverage[pos]);
                btColumn[rowTuple] = minimumPred;
            } else {
                scoreColumn[rowTuple] = getCoverageCost(rowTuple, coverage[pos]);
                btColumn[rowTuple] = ClusterTuple::INVALID_TUPLE;
            }
            firstUnthreadedPosition = pos+1;
            if (scoreColumn[rowTuple] < minimumInColumn) {
                minimumInColumn = scoreColumn[rowTuple];
                minimumTupleInColumn = rowTuple;
                minimumPredTupleInColumn = minimumPred;
            }
        }
        
        // write column into dp table(s)
        d.push_back(scoreColumn);
        b.push_back(btColumn);
        
//         std::cout<<"Best score in column "<<pos<<" is "<<minimumInColumn<<"("<<getCoverageCost(minimumTupleInColumn, coverage[pos])<<" + "<<(minimumInColumn-getCoverageCost(minimumTupleInColumn, coverage[pos]))<<") by "<<(minimumTupleInColumn.asString(ploidy, covMap[pos]))<<" from "<<(minimumPredTupleInColumn.asString(ploidy, covMap[pos-1]))<<std::endl;
    }
    
    std::cout<<"Threading haplotypes through clusters .. ("<<num_vars<<"/"<<num_vars<<")"<<std::endl;
    
    // backtracking start
//     std::cout<<"Start backtracking in column "<<(firstUnthreadedPosition-1)<<std::endl;
    
    ClusterTuple currentRow = ClusterTuple::INVALID_TUPLE;
    Score minimum = std::numeric_limits<Score>::infinity();
    for (std::pair<ClusterTuple, Score> entry : d[firstUnthreadedPosition-1]) {
        if (entry.second < minimum) {
            minimum = entry.second;
            currentRow = entry.first;
        }
    }
    if (currentRow == ClusterTuple::INVALID_TUPLE) {
        std::cout<<"No minimum in last threaded column!"<<std::endl;
    } else {
        if (currentRow.asVector(ploidy, covMap[firstUnthreadedPosition-1]).size() == 0)
            std::cout<<"Problem occured at position "<<(firstUnthreadedPosition-1)<<" in row "<<currentRow.asString(ploidy, covMap[firstUnthreadedPosition-1])<<std::endl;
        path.push_back(currentRow.asVector(ploidy, covMap[firstUnthreadedPosition-1]));
    }
    
//     std::cout<<"Optimal solution value = "<<(d[firstUnthreadedPosition-1][currentRow])<<std::endl;
    
    // backtracking iteration
    for (Position pos = firstUnthreadedPosition-1; pos > 0; pos--) {
        currentRow = b[pos][currentRow];
        if (currentRow.asVector(ploidy, covMap[pos-1]).size() == 0) {
            std::cout<<"Problem occured at position "<<(pos-1)<<" in row "<<currentRow.asString(ploidy, covMap[pos-1])<<std::endl;
            std::vector<GlobalClusterId> fallback;
            for (uint32_t i = 0; i < ploidy; i++)
                fallback.push_back(0);
            path.push_back(fallback);
        } else {
            path.push_back(currentRow.asVector(ploidy, covMap[pos-1]));
        }
    }
    
    // reverse as we constructed it back to front
    std::reverse(path.begin(), path.end());
    
    return path;
}

Score HaploThreader::getCoverageCost(ClusterTuple tuple, const std::vector<double>& coverage) const {
    // tuple contains local cluster ids, which have to be translated with covMap to get the global ids
    Score cost = 0.0;
    for (uint32_t i = 0; i < ploidy; i++) {
        double cov = coverage[tuple.get(i)];
        if (cov == 0) {
            return std::numeric_limits<double>::infinity();
        } else {
            uint32_t expCount = std::round(cov*(double)ploidy);
            if (tuple.count(tuple.get(i), ploidy) != expCount) {
                cost += 1.0;
            }
        }
    }
    return cost;
}

Score HaploThreader::getSwitchCost(const ClusterTuple tuple1, const ClusterTuple tuple2, 
                                   const std::vector<GlobalClusterId>& clusters1, const std::vector<GlobalClusterId>& clusters2) const {
    // tuple contains global cluster ids
    uint32_t switches = 0;
    for (uint32_t i = 0; i < ploidy; i++) {
        switches += clusters1[tuple1.get(i)] != clusters2[tuple2.get(i)];
    }
    return switches*switchCost;
}

std::vector<ClusterTuple> HaploThreader::computeGenotypeConformTuples (const std::vector<GlobalClusterId>& clusters,
                                                                       const std::vector<uint32_t>& consensus,
                                                                       uint32_t genotype, bool allowPermutations) const {
    std::vector<ClusterTuple> perfect = getGenotypeConformTuples (clusters, consensus, genotype, 0, allowPermutations);
    if (perfect.size() > 0) {
        return perfect;
    } else {
        std::vector<ClusterTuple> oneOff = getGenotypeConformTuples (clusters, consensus, genotype, 1, allowPermutations);
        if (oneOff.size() > 0) {
            return oneOff;
        } else {
            std::vector<ClusterTuple> noMatch = getCombinations(clusters.size(), allowPermutations);
            return noMatch;
        }
    }
}

std::vector<ClusterTuple> HaploThreader::getGenotypeConformTuples (const std::vector<GlobalClusterId>& clusters, const std::vector<uint32_t>& consensus, 
                                                                uint32_t genotype, uint32_t distance, bool allowPermutations) const {
    std::unordered_set<ClusterTuple> conformTuples;
    uint32_t maxElem = clusters.size();
    std::vector<LocalClusterId> curPerm(ploidy, 0);
    
    // enumerate all permutations of clusters
    while (curPerm[ploidy-1] < maxElem) {
        // check if genotype distance is correct
        uint32_t g = 0;
        for(uint32_t i = 0; i < ploidy; i++) {
            g += consensus[curPerm[i]];
        }
        
        if ((g - genotype == distance) | (genotype - g == distance)) {
            conformTuples.insert(ClusterTuple(curPerm));
        }
        
        // increment permutation iterator
        curPerm[0]++;
        for (uint32_t i = 1; i < ploidy; i++) {
            if (curPerm[i-1] >= maxElem) {
                curPerm[i]++;
                curPerm[i-1] = allowPermutations ? 0 : curPerm[i];
            }
        }
        if (!allowPermutations)
            for (uint32_t i = ploidy-1; i > 0; i--)
                if (curPerm[i-1] >= maxElem)
                    curPerm[i-1] = curPerm[i];
    }
    
    std::vector<ClusterTuple> uniqueTuples(conformTuples.begin(), conformTuples.end());
    return uniqueTuples;
}

std::vector<ClusterTuple> HaploThreader::getCombinations(uint32_t maxElem, bool allowPermutations) const {        
    std::unordered_set<ClusterTuple> tuples;
    std::vector<LocalClusterId> curPerm(ploidy, 0);
    
    // enumerate all permutations of clusters
    while (curPerm[ploidy-1] < maxElem) {
        tuples.insert(ClusterTuple(curPerm));
        
        // increment permutation iterator
        curPerm[0]++;
        for (uint32_t i = 1; i < ploidy; i++) {
            if (curPerm[i-1] >= maxElem) {
                curPerm[i]++;
                curPerm[i-1] = allowPermutations ? 0 : curPerm[i];
            }
        }
        if (!allowPermutations)
            for (uint32_t i = ploidy-1; i > 0; i--)
                if (curPerm[i-1] >= maxElem)
                    curPerm[i-1] = curPerm[i];
    }
    
    std::vector<ClusterTuple> uniqueTuples(tuples.begin(), tuples.end());
    return uniqueTuples;
}

std::vector<ClusterTuple> HaploThreader::getNonSymmetricExtensions(ClusterTuple tuple, const std::vector<GlobalClusterId>& clusters, 
                                                                   const std::vector<int32_t>& mapPos, std::vector<uint32_t>& numUsed, 
                                                                   uint32_t numClusters) const {
    std::unordered_set<ClusterTuple> s;
    
    ClusterTuple bluePrint(0);
    std::vector<bool> posUsed(ploidy, false);
    std::vector<uint32_t> openLocalClusters;
    std::vector<uint32_t> openPos;
    std::vector<GlobalClusterId> cs;
    
    // for every entry in given tuple, find position for blue print
    for (uint32_t i = 0; i < ploidy; i++) {
        GlobalClusterId c = clusters[tuple.get(i)];
        int32_t pos = mapPos[c + numClusters*numUsed[c]];
        if (pos >= 0) {
            if ((uint32_t)pos >= ploidy) {
                std::cout<<"Map pos was higher than ploidy!"<<std::endl;
            }
            bluePrint.set(pos, tuple.get(i));
            numUsed[c]++;
            cs.push_back(c);
            posUsed[pos] = true;
        } else {
            openLocalClusters.push_back(tuple.get(i));
        }
    }
    // clean numUsed array
    for (GlobalClusterId c : cs) {
        numUsed[c] = 0;
    }
    
    // collect all open positions in the blue print
    for (uint32_t i = 0; i < ploidy; i++) {
        if (!posUsed[i])
            openPos.push_back(i);
    }
    
    if (openPos.size() != openLocalClusters.size()) {
        std::cout<<"Unequal number of open positions and open clusters!"<<std::endl;
    }
    
    // permute only open clusters over open positions
//     std::vector<std::vector<uint32_t>> perms = getPermutationsVector(openPos.size());
//     std::cout<<perms.size()<<std::endl;
    for (std::vector<uint32_t> perm : perms[openPos.size()]) {
        ClusterTuple t(bluePrint);
        for (uint32_t i = 0; i < perm.size(); i++) {
            t.set(openPos[i], openLocalClusters[perm[i]]);
        }
        s.insert(t);
    }
    
    return std::vector<ClusterTuple>(s.begin(), s.end());
}

std::vector<std::vector<uint32_t>> HaploThreader::getPermutationsVector(uint32_t maxElem) const {
    std::vector<std::vector<uint32_t>> perms;
    std::vector<LocalClusterId> curPerm(ploidy, 0);
    
    // enumerate all permutations of clusters
    while (curPerm[ploidy-1] < maxElem) {
        perms.push_back(curPerm);
        
        // increment permutation iterator
        curPerm[0]++;
        for (uint32_t i = 1; i < ploidy; i++) {
            if (curPerm[i-1] >= maxElem) {
                curPerm[i]++;
                curPerm[i-1] = 0;
            }
        }
    }
    return perms;
}

std::vector<ClusterTuple> HaploThreader::extendAll(const std::vector<ClusterTuple>& combinations) const {
    std::vector<ClusterTuple> tuples;
    for (uint32_t i = 0; i < combinations.size(); i++) {
        for (std::vector<uint32_t> perm : getPermutationsVector(ploidy)) {
            ClusterTuple t(combinations[i]);
            t.permute(perm);
            tuples.push_back(t);
        }
    }
    
    std::unordered_set<TupleCode> s;
    for (ClusterTuple t : tuples) {
        s.insert(t.asNumber());
    }
    if (tuples.size() != s.size()) {
        std::cout<<"Extended size was "<<tuples.size()<<" but only "<<s.size()<<" are distinct"<<std::endl;
    }
    
    return tuples;
}
