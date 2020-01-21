#include "haplothreader.h"
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <random>

constexpr uint64_t ClusterTuple::TUPLE_MASKS[];
const ClusterTuple ClusterTuple::INVALID_TUPLE = ClusterTuple((TupleCode)-1);

HaploThreader::HaploThreader (uint32_t ploidy, double switchCost, double affineSwitchCost, bool symmetryOptimization, uint32_t rowLimit) :
    ploidy(ploidy),
    switchCost(switchCost),
    affineSwitchCost(affineSwitchCost),
    symmetryOptimization(symmetryOptimization),
    rowLimit(rowLimit)
{
}

std::vector<std::vector<GlobalClusterId>> HaploThreader::computePaths (const std::vector<Position>& blockStarts, 
                    const std::vector<std::vector<GlobalClusterId>>& covMap,
                    const std::vector<std::vector<double>>& coverage, 
                    const std::vector<std::vector<uint32_t>>& consensus,
                    const std::vector<std::unordered_map<uint32_t, uint32_t>>& genotypes,
                    const std::vector<std::vector<std::vector<double>>>& clusterDissim
                   ) const {
    Position numVars = covMap.size();
    std::vector<std::vector<GlobalClusterId>> path;
    for (uint32_t i = 0; i < blockStarts.size(); i++) {
        Position start = blockStarts[i];
        Position end = i == blockStarts.size()-1 ? numVars : blockStarts[i+1];
        if (end > start) {
            std::vector<std::vector<GlobalClusterId>> section = computePaths(blockStarts[i], end, covMap, coverage, consensus, genotypes, clusterDissim, numVars);
            for (auto tuple : section) {
                path.push_back(tuple);
            }
        }
    }
    return path;
}

std::vector<std::vector<GlobalClusterId>> HaploThreader::computePaths (Position start, Position end, 
                    const std::vector<std::vector<GlobalClusterId>>& covMap,
                    const std::vector<std::vector<double>>& coverage, 
                    const std::vector<std::vector<uint32_t>>& consensus,
                    const std::vector<std::unordered_map<uint32_t, uint32_t>>& genotypes,
                    const std::vector<std::vector<std::vector<double>>>& clusterDissim,
                    Position displayedEnd
                   ) const {
    
    // setup data structures
    std::vector<std::unordered_map<ClusterTuple, ClusterEntry>> m;
    Position firstUnthreadedPosition = start;
    std::vector<std::vector<GlobalClusterId>> path;

    // initialize first column
    if (displayedEnd == 0)
        displayedEnd = end;
    std::vector<ClusterTuple> confTuples = computeGenotypeConformTuples(covMap[start], consensus[start], genotypes[start], false);
    std::unordered_map<ClusterTuple, ClusterEntry> column;
    if (confTuples.size() == 0) {
        std::cout<<"First variant has no clusters!"<<std::endl;
        return path;
    }
    Score minimumInColumn = std::numeric_limits<Score>::infinity();
    ClusterTuple minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
    ClusterTuple minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;
    uint32_t allEntries = 0;
    uint32_t keptEntries = 0;
    for (ClusterTuple t : confTuples) {
        column[t] = ClusterEntry(getCoverageCost(t, coverage[start]), ClusterTuple::INVALID_TUPLE);
        firstUnthreadedPosition = start + 1;
        if (column[t].score < minimumInColumn) {
            minimumInColumn = column[t].score;
            minimumTupleInColumn = t;
        }
    }
    
    // cut down rows if parameter is set
    if (rowLimit > 0 && column.size() >= rowLimit) {
        std::vector<std::pair<ClusterTuple, ClusterEntry>> tuplePairs(column.begin(), column.end());
        std::sort(tuplePairs.begin(), tuplePairs.end(), [this] (const std::pair<ClusterTuple, ClusterEntry>& a, const std::pair<ClusterTuple, ClusterEntry>& b) { return a.second.score < b.second.score; });
        for (uint32_t i = rowLimit; i < tuplePairs.size(); i++) {
            column.erase(tuplePairs[i].first);
        }
    }
    
    m.push_back(std::unordered_map<ClusterTuple, ClusterEntry>(column.begin(), column.end()));
    allEntries += confTuples.size();
    keptEntries += column.size();
    
    // iterate over positions
    for (Position pos = start + 1; pos < end; pos++) {
        
        // reset variables
        confTuples.clear();
        column.clear();
        Score minimum = std::numeric_limits<Score>::infinity();
        Score minDissim = std::numeric_limits<Score>::infinity();
        ClusterTuple minimumPred = ClusterTuple::INVALID_TUPLE;
        bool minExists = false;
        minimumInColumn = std::numeric_limits<Score>::infinity();
        minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
        minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;
        
        // convert genotype map to vector
        uint32_t maxAllele = 0;
        for (std::pair<uint32_t, uint32_t> entry : genotypes[pos]) {
            maxAllele = std::max(maxAllele, entry.first);
        }
        std::vector<uint32_t> genotypeVec(maxAllele, 0);
        for (std::pair<uint32_t, uint32_t> entry : genotypes[pos]) {
            genotypeVec[entry.first] = entry.second;
        }
        
        // compute conform tuples
        confTuples = computeGenotypeConformTuples(covMap[pos], consensus[pos], genotypes[pos], true);
        
        // iterate over rows
        for (ClusterTuple rowTuple : confTuples) {
            minimum = std::numeric_limits<Score>::infinity();
            minimumPred = ClusterTuple::INVALID_TUPLE;
            
            // iterate over previous rows
            for (std::pair<ClusterTuple, ClusterEntry> predEntry : m[pos-1-start]) {
                Score s = predEntry.second.score + getSwitchCost(rowTuple, predEntry.first, covMap[pos], covMap[pos-1]);
                if (s <= minimum) {
                    Score d = predEntry.second.score + getSwitchDissimilarity(rowTuple, predEntry.first, clusterDissim[pos-1]);
                    if (s < minimum || d < minDissim) {
                        minExists = true;
                        minimum = s;
                        minDissim = d;
                        minimumPred = predEntry.first;
                    }
                }
            }
            
            Score coverageCost = getCoverageCost(rowTuple, coverage[pos]);
            
            // report best recursion
            if (minExists) {
                column[rowTuple] = ClusterEntry(minimum + coverageCost, minimumPred);
            } else {
                column[rowTuple] = ClusterEntry(coverageCost, ClusterTuple::INVALID_TUPLE);
            }
            firstUnthreadedPosition = pos+1;
            if (column[rowTuple].score < minimumInColumn) {
                minimumInColumn = column[rowTuple].score;
                minimumTupleInColumn = rowTuple;
                minimumPredTupleInColumn = minimumPred;
            }
        }
        
        if (symmetryOptimization) {
            // remove permutations among entries
            std::unordered_map<TupleCode, std::pair<ClusterTuple, Score>> bestPerms;
            std::vector<ClusterTuple> duplicates;
            for (std::pair<ClusterTuple, ClusterEntry> e : column) {
                TupleCode fp = e.first.fingerprint(ploidy);
                if (e.second.score < bestPerms[fp].second || bestPerms[fp].first == ClusterTuple::INVALID_TUPLE) {
                    bestPerms[fp] = std::pair<ClusterTuple, Score>(e.first, e.second.score);
                } else {
                    duplicates.push_back(e.first);
                }
            }
            for (ClusterTuple t : duplicates) {
                column.erase(t);
            }
            
            // remove non-profitable entries
            std::vector<ClusterTuple> profitableTuples;
            profitableTuples.push_back(minimumTupleInColumn);
            
            for (ClusterTuple t : confTuples) {
                if (t == minimumTupleInColumn)
                    continue;
                bool profitable = true;
                for (ClusterTuple p : profitableTuples) {
                    if (column[t].score >= column[p].score + getSwitchCost(t, p, covMap[pos], covMap[pos])) {
                        profitable = false;
                        break;
                    }
                }
                if (profitable) {
                    if (profitableTuples.size() < ploidy) {
                        profitableTuples.push_back(t);
                    }
                } else {
                    column.erase(t);
                }
            }
        }
        
        // cut down rows if parameter is set
        if (rowLimit > 0 && column.size() >= rowLimit) {
            std::vector<std::pair<ClusterTuple, ClusterEntry>> tuplePairs(column.begin(), column.end());
            std::sort(tuplePairs.begin(), tuplePairs.end(), [this] (const std::pair<ClusterTuple, ClusterEntry>& a, const std::pair<ClusterTuple, ClusterEntry>& b) { return a.second.score < b.second.score; });
            for (uint32_t i = rowLimit; i < tuplePairs.size(); i++) {
                column.erase(tuplePairs[i].first);
            }
        }
        
        // write column into dp table(s)
        m.push_back(std::unordered_map<ClusterTuple, ClusterEntry>(column.begin(), column.end()));
        allEntries += confTuples.size();
        keptEntries += column.size();
    }
    
    // backtracking start
    ClusterTuple currentRow = ClusterTuple::INVALID_TUPLE;
    Score minimum = std::numeric_limits<Score>::infinity();
    for (std::pair<ClusterTuple, ClusterEntry> entry : m[firstUnthreadedPosition-1-start]) {
        if (entry.second.score < minimum) {
            minimum = entry.second.score;
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
    
    // backtracking iteration
    for (Position pos = firstUnthreadedPosition-1; pos > start; pos--) {
        currentRow = m[pos-start][currentRow].pred;
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
            uint32_t realCount = tuple.count(tuple.get(i), ploidy);
            if (realCount != expCount) {
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
    return switches*switchCost + affineSwitchCost * (switches > 0);
}

Score HaploThreader::getSwitchDissimilarity(const ClusterTuple tuple1, const ClusterTuple tuple2,
                                   const std::vector<std::vector<double>>& clusterDissim
                                  ) const {
    // computes switch costs based on local cluster similarity
    double dissim = 0;
    for (uint32_t i = 0; i < ploidy; i++) {
        dissim += (1 + clusterDissim[tuple2.get(i)][tuple1.get(i)]) / 2;
    }
    return dissim;
}

std::vector<ClusterTuple> HaploThreader::computeGenotypeConformTuples (const std::vector<GlobalClusterId>& clusters,
                                                                       const std::vector<uint32_t>& consensus,
                                                                       const std::unordered_map<uint32_t, uint32_t>& genotype, bool allowPermutations) const {
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
                                                                const std::unordered_map<uint32_t, uint32_t>& genotype, uint32_t distance, bool allowPermutations) const {
    std::unordered_set<ClusterTuple> conformTuples;
    uint32_t maxElem = clusters.size();
    std::vector<LocalClusterId> curPerm(ploidy, 0);
    
    // convert genotype map to vector
    uint32_t maxAllele = 0;
    for (std::pair<uint32_t, uint32_t> entry : genotype) {
        maxAllele = std::max(maxAllele, entry.first);
    }
    std::vector<uint32_t> genotypeVec(maxAllele, 0);
    for (std::pair<uint32_t, uint32_t> entry : genotype) {
        genotypeVec[entry.first] = entry.second;
    }
    
    // enumerate all permutations of clusters
    while (curPerm[ploidy-1] < maxElem) {
        // check if genotype distance is correct
        std::vector<uint32_t> tupGenotypeVec(maxAllele, 0);
        for (uint32_t i = 0; i < ploidy; i++) {
            tupGenotypeVec[consensus[curPerm[i]]] += 1;
        }
        
        uint32_t difference = 0;
        for (uint32_t i = 0; i < maxAllele; i++) {
            difference += std::max(tupGenotypeVec[i] - genotypeVec[i], genotypeVec[i] - tupGenotypeVec[i]);
        }
        
        if (difference == distance) {
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

Score HaploThreader::getGenotypeDist(const ClusterTuple tuple, 
                                     const std::vector<uint32_t>& consensus,
                                     std::vector<uint32_t>& genotypeVec) const {
    std::vector<uint32_t> tupGenotypeVec(genotypeVec.size(), 0);
    for (uint32_t i = 0; i < ploidy; i++) {
        tupGenotypeVec[consensus[tuple.get(i)]] += 1;
    }
    
    uint32_t difference = 0;
    for (uint32_t i = 0; i < genotypeVec.size(); i++) {
        difference += std::max(tupGenotypeVec[i] - genotypeVec[i], genotypeVec[i] - tupGenotypeVec[i]);
    }
    
    return difference;
}
