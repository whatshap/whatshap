#include "haplothreader.h"
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <cassert>
#include "../genotype.h"
#include "../binomial.h"

constexpr uint64_t ClusterTuple::TUPLE_MASKS[];
const ClusterTuple ClusterTuple::INVALID_TUPLE = ClusterTuple((TupleCode)-1);

HaploThreader::HaploThreader (uint32_t ploidy, double switchCost, double affineSwitchCost, bool carryOverPreviousTuples, uint32_t rowLimit) :
    ploidy(ploidy),
    switchCost(switchCost),
    affineSwitchCost(affineSwitchCost),
    carryOverPreviousTuples(carryOverPreviousTuples),
    rowLimit(rowLimit)
{
}

std::vector<std::vector<GlobalClusterId>> HaploThreader::computePaths (const std::vector<Position>& blockStarts, 
                    const std::vector<std::vector<GlobalClusterId>>& covMap,
                    const std::vector<std::vector<std::unordered_map<uint32_t, uint32_t>>>& alleleDepths) const {
    Position numVars = covMap.size();
    std::vector<std::vector<GlobalClusterId>> path;
    for (uint32_t i = 0; i < blockStarts.size(); i++) {
        Position start = blockStarts[i];
        Position end = i == blockStarts.size()-1 ? numVars : blockStarts[i+1];
        if (end > start) {
            std::vector<std::vector<GlobalClusterId>> section = computePaths(blockStarts[i], end, covMap, alleleDepths, numVars);
            for (auto tuple : section) {
                path.push_back(tuple);
            }
        }
    }
    return path;
}

std::vector<std::vector<GlobalClusterId>> HaploThreader::computePaths (Position start, Position end, 
                    const std::vector<std::vector<GlobalClusterId>>& covMap,
                    const std::vector<std::vector<std::unordered_map<uint32_t, uint32_t>>>& alleleDepths,
                    Position displayedEnd) const {
    
    //  compute coverage and consensus based on allele depths
    std::vector<uint32_t> coverage(alleleDepths.size(), 0);
    std::vector<std::vector<uint32_t>> clusterCoverage(alleleDepths.size(), std::vector<uint32_t>());
    computeCoverage(alleleDepths, covMap, coverage, clusterCoverage);
    
    // the actual DP table with sparse columns
    std::vector<std::unordered_map<ClusterTuple, TupleEntry>> m;
    
    // data structure to store the final result
    std::vector<std::vector<GlobalClusterId>> path;

    // initialize first column
    if (displayedEnd == 0)
        displayedEnd = end;
    Position firstUnthreadedPosition = start;
    
    // auxiliary vector to store a vector for each tuple of a column. This vector contains the global cluster ids, sorted in ascending order
    std::unordered_map<ClusterTuple, std::vector<GlobalClusterId>> sortedGlobalTuples;
    
    // allocated space to store the current column, before it is stored in the DP table
    std::unordered_map<ClusterTuple, TupleEntry> column;
    
    /*
     * The basic idea of this algorithm is, that for every position we generate candidate tuples, which represent
     * the multiset of clusters, through which the haplotypes are threaded at this exact position. Therefore, for
     * every candidate we have to compute the coverage costs and the best predecessor tuple from the last column.
     * The best predecessor is the tuple, which minimizes the sum over its own total cost, plus the switch costs
     * to the candidate of the current column.
     * 
     * There is heavy symmetry optimization included in this algorithm. First, for a completely computed column
     * there is no point in having two tuples t1 and t2, which are permutations of each other. Since we consider
     * all genotype conform tuples in the next column anyways, we only need to store the better of the two t1 and
     * t2. This elimination of permutations is mainly done by the candidate generation and by the switch cost
     * functions:
     * 
     * 1. The candidate generator avoids permutations by construction.
     * 2. The advanced switch cost function is able to compute the minimal switch costs between a tuple t1 and all
     *    permutations of a tuple t2. If one candidate for the current column is processed, all of its
     *    permutations are actually processed as well. We only keep the permutation of t2 with the lowest switch
     *    cost, since all other permutations have equal coverage costs and can thus be discarded.
     * 
     * After a column is computed, it is optionally pruned by removing non-profitable tuples. Let t1 and t2 be
     * tuples in the current column. Then t2 is non-profitable if it holds that
     * 
     * total_cost(t2) >= total_cost(t1) + switch_cost(t1, t2)
     * 
     */
    
    std::vector<ClusterTuple> relevantTuples;
    std::unordered_set<ClusterTuple> enteringTuples;
    std::unordered_set<ClusterTuple> leftTuples;
    
    for (Position pos = start; pos < end; pos++) {
        // reset variables
        enteringTuples.clear();
        leftTuples.clear();
        relevantTuples.clear();
        column.clear();
        ThreadScore minimumSwitchCost = std::numeric_limits<ThreadScore>::infinity();
        ThreadScore optPredScore = std::numeric_limits<ThreadScore>::infinity();
        ClusterTuple optPredTuple = ClusterTuple::INVALID_TUPLE;
        Position offset = 1;
        
        // compute genotype conform tuples
        relevantTuples = computeRelevantTuples(clusterCoverage, pos);
        if (pos > start) {
            addCarriedTuples(relevantTuples, enteringTuples, leftTuples, covMap[pos], covMap[pos-1], m[pos-1]);
        } else {
            offset = 0;
            for (ClusterTuple t : relevantTuples)
                enteringTuples.insert(t);
        }
        
        // create translater between tuples of new and old position
        TupleConverter tc(covMap[pos - offset], covMap[pos], ploidy);
        
        // iterate over generated tuples
        printStartOfColumn(pos, relevantTuples, enteringTuples, covMap[pos]);
        if (relevantTuples.size()+enteringTuples.size() == 0) {
            std::cout<<"No tuples for this position. Aborting ..."<<std::endl;
            break;
        }
        
        for (ClusterTuple tuple : relevantTuples) {
            // variables to store best score and backtracking direction
            minimumSwitchCost = std::numeric_limits<ThreadScore>::infinity();
            
            // auxiliary data, is precomputed once here
            std::vector<GlobalClusterId> tupleGlobal = tuple.asVector(ploidy, covMap[pos]);
            std::sort(tupleGlobal.begin(), tupleGlobal.end());
            
            //TODO: For non-entering tuples, just compute the identical predecessor
            
            // compare each new tuple with every tuple from previous column
            if (pos > start) {
                for (std::pair<ClusterTuple, TupleEntry> predEntry : m[pos-1-start]) {
                    
                    // retrieve precomputed sorted vector over global ids and compute optimal switch cost
                    std::vector<GlobalClusterId> predTupleGlobal = sortedGlobalTuples[predEntry.first];
                    ThreadScore s = predEntry.second.score + getSwitchCostAllPerms(predTupleGlobal, tupleGlobal);
//                     std::cout<<"   switch [ ";
//                     for (auto& a: predTupleGlobal)
//                         std::cout<<a<<" ";
//                     std::cout<<"] -> [ ";
//                     for (auto& a: tupleGlobal)
//                         std::cout<<a<<" ";
//                     std::cout<<"] = "<<s<<std::endl;
                    
                    if (s < minimumSwitchCost) {
                        minimumSwitchCost = s;
                        optPredScore = predEntry.second.score;
                        optPredTuple = predEntry.first;
                    }
                }
            } else {
                optPredScore = 0.0;
                optPredTuple = tuple;
            }
            
            // in addition to best score over all predecessors, we need the best permutation of tuple to achieve this
            ClusterTuple bestPerm = tc.permuteAgainstOld(tuple, optPredTuple);
            ThreadScore coverageCost = getCoverageCost(tuple, coverage[pos], clusterCoverage[pos]);
//             ThreadScore coverageCost2 = getCoverageCost(bestPerm, coverage[pos], clusterCoverage[pos]);
//             if (coverageCost != coverageCost2)  {
//                 std::cout<<"   Inconsistent coverage cost: "<<coverageCost<<" != "<<coverageCost2<<std::endl;
//                 std::cout<<"   "<<tuple.asString(ploidy)<<std::endl;
//                 std::cout<<"   "<<bestPerm.asString(ploidy)<<std::endl;
//             }

            column[bestPerm] = TupleEntry(optPredScore + coverageCost, coverageCost, optPredTuple);
            firstUnthreadedPosition = pos+1;
//             std::cout<<"   "<<bestPerm.asString(ploidy, covMap[pos])<<" -> "<<column[bestPerm].score<<" via "<<optPredTuple.asString(ploidy, covMap[pos - offset])<<std::endl;
        }
        
        // precompute the sorted vectors with global cluster ids for this column (will be reused in next column)
        sortedGlobalTuples.clear();
        for (auto& a : column) {
            std::vector<GlobalClusterId> tupleGlobal = a.first.asVector(ploidy, covMap[pos]);
            std::sort(tupleGlobal.begin(), tupleGlobal.end());
            sortedGlobalTuples[a.first] = tupleGlobal;
        }
        
        // cut down rows if parameter is set
        if (rowLimit > 0 && column.size() >= rowLimit) {
            std::vector<std::pair<ClusterTuple, TupleEntry>> tuplePairs(column.begin(), column.end());
            std::sort(tuplePairs.begin(), tuplePairs.end(), [] (const std::pair<ClusterTuple, TupleEntry>& a, const std::pair<ClusterTuple, TupleEntry>& b) { return a.second.score < b.second.score; });
            for (uint32_t i = rowLimit; i < tuplePairs.size(); i++) {
                column.erase(tuplePairs[i].first);
            }
        }
        
        // write column into dp table
        m.push_back(std::unordered_map<ClusterTuple, TupleEntry>(column.begin(), column.end()));
    }
    
    // discard auxiliary data
    sortedGlobalTuples.clear();
    
    // backtracking start
    ClusterTuple currentRow = ClusterTuple::INVALID_TUPLE;
    ThreadScore minimum = std::numeric_limits<ThreadScore>::infinity();
    for (std::pair<ClusterTuple, TupleEntry> entry : m[firstUnthreadedPosition - 1 - start]) {
        if (entry.second.score < minimum) {
            minimum = entry.second.score;
            currentRow = entry.first;
        }
    }
    if (currentRow == ClusterTuple::INVALID_TUPLE) {
        std::cout<<"No minimum in last threaded column!"<<std::endl;
    } else {
        if (currentRow.asVector(ploidy, covMap[firstUnthreadedPosition - 1]).size() == 0)
            std::cout<<"Problem occured at position "<<(firstUnthreadedPosition - 1)<<" in row "<<currentRow.asString(ploidy, covMap[firstUnthreadedPosition-1])<<std::endl;
        path.push_back(currentRow.asVector(ploidy, covMap[firstUnthreadedPosition - 1]));
    }
    
    // backtracking iteration
    for (Position pos = firstUnthreadedPosition - 1; pos > start; pos--) {
        currentRow = m[pos-start][currentRow].pred;
        if (currentRow.asVector(ploidy, covMap[pos - 1]).size() == 0) {
            std::cout<<"Problem occured at position "<<(pos-1)<<" in row "<<currentRow.asString(ploidy, covMap[pos - 1])<<std::endl;
            std::vector<GlobalClusterId> fallback;
            for (uint32_t i = 0; i < ploidy; i++)
                fallback.push_back(0);
            path.push_back(fallback);
        } else {
            path.push_back(currentRow.asVector(ploidy, covMap[pos - 1]));
        }
    }
    
    // reverse as we constructed it back to front
    std::reverse(path.begin(), path.end());
    
    return path;
}

ThreadScore HaploThreader::getCoverageCost(ClusterTuple tuple, 
                                     const uint32_t coverage, 
                                     const std::vector<uint32_t>& clusterCoverage
                                    ) const {
    // tuple contains local cluster ids, which have to be translated with covMap to get the global ids
    double cost = 1.0;
    double opt = 1.0;
    double count = 0.0;
    uint32_t unthreadedReads = 0;
    
    std::vector<uint32_t> clustMult(clusterCoverage.size(), 0);
    for (uint32_t i = 0; i < ploidy; i++) {
        uint32_t cid = tuple.get(i);
        clustMult[cid] += 1;
        opt *= binom_pmf(coverage, ((double)coverage)/ploidy, 1.0/ploidy);
    }
    
    uint32_t cov = 0;
    for (uint32_t cid = 0; cid < clusterCoverage.size(); cid++) {
        cov += clusterCoverage [cid];
    }
    
    for (uint32_t cid = 0; cid < clusterCoverage.size(); cid++) {
        if (clustMult[cid] == 0) {
            unthreadedReads += clusterCoverage[cid];
        } else {
            double p = (0.95*clustMult[cid])/ploidy;
            cost *= binom_pmf(coverage, clusterCoverage[cid], p);
        }
        uint32_t expMult = (uint32_t)(0.5 + ploidy * ((ThreadScore)clusterCoverage[cid] / (ThreadScore)cov));
        if (expMult != clustMult[cid])
            count += 1.0;
    }
    
    cost *= binom_pmf(coverage, unthreadedReads, 0.05);

    return std::log(opt/cost);
//     return count;
}

ThreadScore HaploThreader::getSwitchCostAllPerms(const std::vector<GlobalClusterId>& prevTuple, const std::vector<GlobalClusterId>& curTuple) const {
    uint32_t pIdx = 0;
    uint32_t cIdx = 0;
    uint32_t switches = 0;
    // compare zig-zag-wise over sorted tuples
    while ((pIdx < ploidy) & (cIdx < ploidy)) {
        if (prevTuple[pIdx] == curTuple[cIdx]) {
            pIdx++; cIdx++;
        } else if (prevTuple[pIdx] < curTuple[cIdx]) {
            switches++;
            pIdx++;
        } else {
            cIdx++;
        }
    }
    switches += (ploidy - pIdx);
    
    return switchCost*switches + affineSwitchCost*(switches > 0);
}

std::vector<ClusterTuple> HaploThreader::computeRelevantTuples (const std::vector<std::vector<uint32_t>>& clusterCoverage, 
                                                                Position pos) const {
                                                                    
    uint32_t coverage = 0;
    for (uint32_t clustCov : clusterCoverage[pos])
        coverage += clustCov;
    
    std::vector<LocalClusterId> relevantClusters;
    std::vector<std::vector<uint32_t>> consensusLists;
    std::cout<<"   relevant ("<<coverage<<"):";
    for (uint32_t cid = 0; cid < clusterCoverage[pos].size(); cid++) {
        uint32_t smoothedCov = clusterCoverage[pos][cid];
        if (pos > 0)
            smoothedCov = std::max(clusterCoverage[pos - 1][cid], smoothedCov);
        if (pos < clusterCoverage.size() - 1)
            smoothedCov = std::max(clusterCoverage[pos + 1][cid], smoothedCov);
        if (4 * ploidy * smoothedCov >= coverage) {
            relevantClusters.push_back(cid);
            std::cout<<" ["<<cid<<","<<smoothedCov<<"]";
        } else {
            std::cout<<" ("<<cid<<","<<smoothedCov<<")";
        }
    }
    std::cout<<std::endl;
    
    // create vector of combinations
    std::vector<ClusterTuple> relevantTuples;
    
    uint32_t maxElem = relevantClusters.size();
    std::vector<uint32_t> v(ploidy, 0);
    
    /*
     * store vector, until maxElem-1 is in the last field, because then we must just
     * have stored the vector [maxElem-1, maxElem-1, ..., maxElem-1] and we are finished.
     */
    while (v[ploidy-1] < maxElem) {
        // translate to local cluster ids and store in combsOfAllele
        std::vector<uint32_t> c(ploidy, 0);
        for (uint32_t i = 0; i < ploidy; i++)
            c[i] = relevantClusters[v[i]];
        relevantTuples.push_back(ClusterTuple(c));
//         std::cout<<"      "<<ClusterTuple(c).asString(ploidy)<<std::endl;
        
        // increment like a counter
        v[0]++;
        
        // if element i-1 overflowed, increase element i by 1 (which then might also overflow and so on)
        for (uint32_t i = 1; i < ploidy; i++)
            if (v[i-1] >= maxElem)
                v[i]++;
        
        // any element i-1 which overflowed will be set to its minimum, i.e. the value of element i
        for (uint32_t i = ploidy-1; i > 0; i--)
            if (v[i-1] >= maxElem)
                v[i-1] = v[i];
    }
    
    return relevantTuples;
}

void HaploThreader::addCarriedTuples (std::vector<ClusterTuple>& relevantTuples,
                                      std::unordered_set<ClusterTuple>& enteringTuples,
                                      std::unordered_set<ClusterTuple>& leftTuples,
                                      const std::vector<GlobalClusterId>& curClusterIds,
                                      const std::vector<GlobalClusterId>& prevClusterIds,
                                      const std::unordered_map<ClusterTuple, TupleEntry>& prevColumn) const {

    std::unordered_set<uint32_t> matchedIndices;
    std::unordered_map<ClusterTuple, uint32_t> matchSet;
    for (uint32_t i = 0; i < relevantTuples.size(); i++) {
        matchSet[relevantTuples[i].fingerprint(ploidy)] = i;
    }
    
    // compute mapping betweem local cluster ids from previous to current position
    std::unordered_map<GlobalClusterId, LocalClusterId> idMap;
    std::unordered_map<LocalClusterId, LocalClusterId> directMap;
    for (LocalClusterId c = 0; c < curClusterIds.size(); c++)
        idMap[curClusterIds[c]] = c;
    for (LocalClusterId c = 0; c < prevClusterIds.size(); c++) {
        GlobalClusterId g = prevClusterIds[c];
        if (idMap.find(g) != idMap.end())
            directMap[c] = idMap[g];
    }

    // for all tuples from previous position: either add them to leftTuples or bring the currently used tuple into right order
    for (auto& a : prevColumn) {
        std::vector<LocalClusterId> v;
        for (uint32_t i = 0; i < ploidy; i++) {
            LocalClusterId c = a.first.get(i);
            if (idMap.find(c) != directMap.end())
                v.push_back(directMap[c]);
        }
        if (v.size() == ploidy) {
            ClusterTuple t(v);
            TupleCode fp = t.fingerprint(ploidy);
            if (matchSet.find(fp) != matchSet.end()) {
                relevantTuples[matchSet[fp]] = t;
                matchedIndices.insert(matchSet[fp]);
            } else {
                leftTuples.insert(a.first);
            }
        } else {
            leftTuples.insert(a.first);
        }
    }
    
    // for all current tuples not matched by a previous one: add to enteringTuples
    for (uint32_t i = 0; i < relevantTuples.size(); i++) {
        if (matchedIndices.find(i) != matchedIndices.end())
            enteringTuples.insert(relevantTuples[i]);
    }
}

void HaploThreader::computeCoverage (const std::vector<std::vector<std::unordered_map<uint32_t, uint32_t>>>& alleleDepths,
                                     const std::vector<std::vector<GlobalClusterId>>& covMap,
                                     std::vector<uint32_t>& coverage,
                                     std::vector<std::vector<uint32_t>>& clusterCoverage) const {
    for (Position pos = 0; pos < alleleDepths.size(); pos++) {
        uint32_t total = 0;
        for (uint32_t cid = 0; cid < alleleDepths[pos].size(); cid++) {
            uint32_t local = 0;
            for (auto& ad: alleleDepths[pos][cid])
                local += ad.second;
            total += local;
            clusterCoverage[pos].push_back(local);
        }
        for (uint32_t cid = alleleDepths[pos].size() + 1; cid < covMap[pos].size(); cid++)
            clusterCoverage[pos].push_back(0);
        coverage[pos] = total;
    }
}
