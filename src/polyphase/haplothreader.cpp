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

HaploThreader::HaploThreader (uint32_t ploidy, double switchCost, double affineSwitchCost, uint32_t maxClusterGap, uint32_t rowLimit) :
    ploidy(ploidy),
    switchCost(switchCost),
    affineSwitchCost(affineSwitchCost),
    maxClusterGap(maxClusterGap),
    rowLimit(rowLimit)
{
}

std::vector<std::vector<GlobalClusterId>> HaploThreader::computePaths (const std::vector<Position>& blockStarts, 
                    const std::vector<std::vector<GlobalClusterId>>& covMap,
                    const std::vector<std::unordered_map<GlobalClusterId, std::unordered_map<uint32_t, uint32_t>>>& alleleDepths) const {
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
                    const std::vector<std::unordered_map<GlobalClusterId, std::unordered_map<uint32_t, uint32_t>>>& alleleDepths,
                    Position displayedEnd) const {
    
    //  compute coverage and consensus based on allele depths
    std::vector<uint32_t> coverage(alleleDepths.size(), 0);
    std::vector<std::vector<uint32_t>> clusterCoverage(alleleDepths.size(), std::vector<uint32_t>());
    computeCoverage(alleleDepths, covMap, coverage, clusterCoverage);
    
    // the actual DP table with sparse columns
    std::vector<std::unordered_map<ClusterTuple, TupleEntry>> m;
    
    // initialize first column
    if (displayedEnd == 0)
        displayedEnd = end;
    Position firstUnthreadedPosition = start;
    
    // auxiliary vector to store a vector for each tuple of a column. This vector contains the global cluster ids, sorted in ascending order
    std::unordered_map<ClusterTuple, std::vector<GlobalClusterId>> sortedGlobalTuples;

    // auxiliary vector to sort previous column entries by score. map to efficiently find tuples from current column in previous one
    std::vector<ColumnEntry> scoreSortedPreds;
    std::unordered_map<TupleCode, ClusterTuple> fpToTuple;
    
    // allocated space to store the current column, before it is stored in the DP table
    std::unordered_map<ClusterTuple, TupleEntry> column;
    
    /*
     * The basic idea of this algorithm is, that for every position we generate candidate tuples, which represent
     * the multiset of clusters, through which the haplotypes are threaded at this exact position. The algorithm
     * is designed as a DP, where each column stands for one position and contains all candidate cluster multisets
     * (represented as cluster tuples). An entry m[t][i] contains the best threading ending in position i on
     * tuple t. The recursion is defined as
     *
     * m[0][t] = cov_cost(t)
     * m[i][t] = cov_cost(t) + {switch_cost(t') + m[i-1][t'] for all t' in m[i-1]}
     *
     * There is heavy symmetry optimization included in this algorithm. First, for a completely computed column
     * there is no point in having two tuples t1 and t2, which are permutations of each other. Since we can
     * reorder the tuples in the next column anyways, we only need to store the best permutations:
     * 
     * 1. The candidate generator avoids permutations by construction.
     * 2. The switch cost function computes the minimal switch cost for a tuple t and all permutations of a
     *    second tuple t' simultaneously. The minimizing permutation is reported and taken for the DP column.
     * 
     * For each column we only consider candidate tuples with coverage cost less than 4 times the best coverage
     * cost to not waste too much time on non-promising tuples.
     */
    
    std::vector<ClusterTuple> relevantTuples;
    
    for (Position pos = start; pos < end; pos++) {
        // reset variables
        relevantTuples.clear();
        fpToTuple.clear();
        scoreSortedPreds.clear();
        column.clear();
        Position offset = pos > start;
        
        // compute genotype conform tuples
        relevantTuples = computeRelevantTuples(clusterCoverage, pos);
        
        // create translater between tuples of new and old position
        TupleConverter tc(covMap[pos - offset], covMap[pos], ploidy);

        //printStartOfColumn(pos, relevantTuples, covMap[pos]);
        if (relevantTuples.size() == 0) {
            std::cout<<"No tuples for this position. Aborting ..."<<std::endl;
            break;
        }

        // precompute coverage costs for all tuples
        std::vector<ThreadScore> coverageCosts;
        ThreadScore minCovCost = std::numeric_limits<ThreadScore>::infinity();
        for (ClusterTuple tuple : relevantTuples) {
            ThreadScore coverageCost = getCoverageCost(tuple, coverage[pos], clusterCoverage[pos]);
            coverageCosts.push_back(coverageCost);
            if (coverageCost < minCovCost)
                minCovCost = coverageCost;
        }

        // precompute fingerprints for previous column and list sorted by score
        if (pos > start) {
            for (ColumnEntry predEntry : m[pos - 1 - start]) {
                fpToTuple[predEntry.first.fingerprint(ploidy)] = predEntry.first;
                scoreSortedPreds.push_back(predEntry);
            }
            std::sort(scoreSortedPreds.begin(), scoreSortedPreds.end(), [] (const ColumnEntry& a, const ColumnEntry& b) { return a.second.score < b.second.score; });
        }

        // compute full score for each tuple
        for (uint32_t tid = 0; tid < relevantTuples.size(); tid++) {
            // prune tuples with too high coverage cost
            ClusterTuple tuple = relevantTuples[tid];
            ThreadScore coverageCost = coverageCosts[tid];
            if (coverageCost > 30 + minCovCost)
                continue;

            // variables to store best score and backtracking direction
            ThreadScore optPredScore = std::numeric_limits<ThreadScore>::infinity();
            ClusterTuple optPredTuple = ClusterTuple::INVALID_TUPLE;
            
            // auxiliary data, is precomputed once here
            std::vector<GlobalClusterId> tupleGlobal = tuple.asVector(ploidy, covMap[pos]);
            std::sort(tupleGlobal.begin(), tupleGlobal.end());

            // use recursion over previous column for position >= 1
            if (fpToTuple.size() > 0) {
                // retrieve equivalent tuple from previous position and its score (if exists)
                TupleCode c = tc.convertNewToOld(tuple).fingerprint(ploidy);
                if (fpToTuple.find(c) != fpToTuple.end()) {
                    optPredTuple = fpToTuple[c];
                    if (m[pos - 1 -start].find(optPredTuple) != m[pos - 1 -start].end()) {
                        optPredScore = m[pos - 1 - start][optPredTuple].score;
                    }
                }

                // compare each new tuple with every tuple from previous column
                for (uint32_t pid = 0; pid < scoreSortedPreds.size(); pid++) {
                    ColumnEntry predEntry = scoreSortedPreds[pid];
                    // stop if costs cannot improve anymore
                    if (predEntry.second.score + switchCost + affineSwitchCost >= optPredScore)
                        break;

                    // retrieve precomputed sorted vector over global ids and compute optimal switch cost
                    std::vector<GlobalClusterId> predTupleGlobal = sortedGlobalTuples[predEntry.first];
                    ThreadScore s = predEntry.second.score + getSwitchCostAllPerms(predTupleGlobal, tupleGlobal);
                    if (s < optPredScore) {
                        optPredScore = s;
                        optPredTuple = predEntry.first;
                    }
                }
            } else {
                optPredScore = 0.0;
                optPredTuple = tuple;
            }

            // in addition to best score over all predecessors, we need the best permutation of tuple to achieve this
            ClusterTuple bestPerm = tc.permuteAgainstOld(tuple, optPredTuple);

            column[bestPerm] = TupleEntry(optPredScore + coverageCost, optPredTuple);
            firstUnthreadedPosition = pos+1;
        }
        
        // precompute the sorted vectors with global cluster ids for this column (will be reused in next column)
        sortedGlobalTuples.clear();
        for (auto& a : column) {
            std::vector<GlobalClusterId> tupleGlobal = a.first.asVector(ploidy, covMap[pos]);
            std::sort(tupleGlobal.begin(), tupleGlobal.end());
            sortedGlobalTuples[a.first] = tupleGlobal;
        }
        
        //std::cout<<"   Unpruned tuples: "<<column.size()<<std::endl;
        
        // cut down rows if parameter is set
        if (rowLimit > 0 && column.size() >= rowLimit) {
            std::vector<ColumnEntry> tuplePairs(column.begin(), column.end());
            std::sort(tuplePairs.begin(), tuplePairs.end(), [] (const ColumnEntry& a, const ColumnEntry& b) { return a.second.score < b.second.score; });
            for (uint32_t i = rowLimit; i < tuplePairs.size(); i++) {
                column.erase(tuplePairs[i].first);
            }
        }
        
        // write column into dp table
        m.push_back(std::unordered_map<ClusterTuple, TupleEntry>(column.begin(), column.end()));
    }
    
    // discard auxiliary data
    sortedGlobalTuples.clear();
    
    /*
     * ==================
     * backtracking start
     * ==================
     */

    // data structure to store the final result
    std::vector<std::vector<GlobalClusterId>> path;
    
    ClusterTuple currentRow = ClusterTuple::INVALID_TUPLE;
    ThreadScore minimum = std::numeric_limits<ThreadScore>::infinity();
    for (ColumnEntry entry : m[firstUnthreadedPosition - 1 - start]) {
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
    double llh = 1.0;
    uint32_t unthreadedReads = 0;
    
    std::vector<uint32_t> clustMult(clusterCoverage.size(), 0);
    for (uint32_t i = 0; i < ploidy; i++) {
        uint32_t cid = tuple.get(i);
        clustMult[cid] += 1;
    }
    
    uint32_t cov = 0;
    for (uint32_t cid = 0; cid < clusterCoverage.size(); cid++) {
        cov += clusterCoverage [cid];
    }
    
    for (uint32_t cid = 0; cid < clusterCoverage.size(); cid++) {
        if (clustMult[cid] == 0) {
            unthreadedReads += clusterCoverage[cid];
        } else {
            double p = (0.975*clustMult[cid])/ploidy;
            llh *= binom_pmf(coverage, clusterCoverage[cid], p);
        }
    }
    
    llh *= binom_pmf(coverage, unthreadedReads, 0.025);

    return -std::log(llh);
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
    //std::cout<<"Position "<<pos<<": relevant ("<<coverage<<"):";
    for (uint32_t cid = 0; cid < clusterCoverage[pos].size(); cid++) {
        uint32_t cov = clusterCoverage[pos][cid];
        if (4 * ploidy * cov >= coverage) {
            relevantClusters.push_back(cid);
            //std::cout<<" ["<<cid<<","<<cov<<"]";
        //} else {
        //    std::cout<<" ("<<cid<<","<<cov<<")";
        }
    }
    //std::cout<<std::endl;
    
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

void HaploThreader::computeCoverage (const std::vector<std::unordered_map<GlobalClusterId, std::unordered_map<uint32_t, uint32_t>>>& alleleDepths,
                                     const std::vector<std::vector<GlobalClusterId>>& covMap,
                                     std::vector<uint32_t>& coverage,
                                     std::vector<std::vector<uint32_t>>& clusterCoverage) const {
    Position numPos = alleleDepths.size();
    std::vector<std::unordered_map<GlobalClusterId, uint32_t>> clusterCoverageGlobal(numPos, std::unordered_map<GlobalClusterId, uint32_t>());
    for (Position pos = 0; pos < numPos; pos++) {
        for (uint32_t i = 0; i < covMap[pos].size(); i++) {
            uint32_t local = 0;
            GlobalClusterId cid = covMap[pos][i];
            for (auto& a: alleleDepths[pos].at(covMap[pos][i])) {
                local += a.second;
            }
            clusterCoverageGlobal[pos][cid] = local;
        }
    }
    for (Position pos = 0; pos < numPos; pos++) {
        uint32_t total = 0;
        for (uint32_t i = 0; i < covMap[pos].size(); i++) {
            GlobalClusterId cid = covMap[pos][i];
            uint32_t smoothedCov = 0;
            uint32_t numNonZero = 0;
            Position min = pos - maxClusterGap / 2;
            Position max = std::min(numPos - 1, pos + (maxClusterGap + 1) / 2);
            min *= min < max;
            for (Position p = min; p <= max; p++) {
                uint32_t cov = clusterCoverageGlobal[p][cid];
                if (cov > 0) {
                    smoothedCov += cov;
                    numNonZero++;
                }
            }
            if (numNonZero == 0)
                numNonZero = 1;
            clusterCoverage[pos].push_back(smoothedCov / numNonZero);
            total += clusterCoverage[pos][i];
        }
        coverage[pos] = total;
    }
}
