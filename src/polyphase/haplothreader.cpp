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
                    const std::vector<std::unordered_map<uint32_t, uint32_t>>& genotypes
                   ) const {
    Position numVars = covMap.size();
    std::vector<std::vector<GlobalClusterId>> path;
    for (uint32_t i = 0; i < blockStarts.size(); i++) {
        Position start = blockStarts[i];
        Position end = i == blockStarts.size()-1 ? numVars : blockStarts[i+1];
        if (end > start) {
            std::vector<std::vector<GlobalClusterId>> section = computePaths(blockStarts[i], end, covMap, coverage, consensus, genotypes, numVars);
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
                    Position displayedEnd
                   ) const {
    
    // the actual DP table with sparse columns
    std::vector<std::unordered_map<ClusterTuple, ClusterEntry>> m;
    
    // data structure to store the final result
    std::vector<std::vector<GlobalClusterId>> path;

    // initialize first column
    if (displayedEnd == 0)
        displayedEnd = end;
    Position firstUnthreadedPosition = start;
    
    /*
     * Compute the genotype conform tuples for the first column and quit if this set is empty.
     * Note that tuples in general only contain local cluster ids, which must be mapped by covMap[column_id]
     * in order to retrieve global cluster ids. The local ids make for a compact representation, but the global
     * ids are necessary to compare tuples from different column.
     */
    std::vector<ClusterTuple> confTuples = computeGenotypeConformTuples(covMap[start], consensus[start], genotypes[start]);
    if (confTuples.size() == 0) {
        std::cout<<"First variant has no clusters!"<<std::endl;
        return path;
    }
    
    // auxiliary vector to store the optimal permutations of the conform tuples
    std::vector<ClusterTuple> permedTuples;
    
    // auxiliary vector to store a vector for each tuple of a column. This vector contains the global cluster ids and is sorted
    std::unordered_map<ClusterTuple, std::vector<GlobalClusterId>> sortedGlobalTuples;
    
    // allocated space to store the current column, before it is stored in the DP table
    std::unordered_map<ClusterTuple, ClusterEntry> column;
    
    // allocated space to store the optima for the current column
    Score minimumInColumn = std::numeric_limits<Score>::infinity();
    ClusterTuple minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
    ClusterTuple minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;

    // fill first column by only using the coverage cost of each candidate tuple
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
    
    // store first column in DP table
    m.push_back(std::unordered_map<ClusterTuple, ClusterEntry>(column.begin(), column.end()));
    
    // precompute the sorted vector with global cluster ids for every entry in first column
    for (std::pair<ClusterTuple, ClusterEntry> predEntry : m[0]) {
        std::vector<GlobalClusterId> tupleGlobal = predEntry.first.asVector(ploidy, covMap[0]);
        std::sort(tupleGlobal.begin(), tupleGlobal.end());
        sortedGlobalTuples[predEntry.first] = tupleGlobal;
    }
    
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
    for (Position pos = start + 1; pos < end; pos++) {

        // reset variables
        confTuples.clear();
        permedTuples.clear();
        column.clear();
        Score minimum = std::numeric_limits<Score>::infinity();
        ClusterTuple minimumPred = ClusterTuple::INVALID_TUPLE;
        bool minExists = false;
        minimumInColumn = std::numeric_limits<Score>::infinity();
        minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
        minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;
        
        // compute genotype conform tuples
        confTuples = computeGenotypeConformTuples(covMap[pos], consensus[pos], genotypes[pos]);
        
        // iterate over generated tuples
        for (ClusterTuple rowTuple : confTuples) {
            // variables to store best score and backtracking direction
            minimum = std::numeric_limits<Score>::infinity();
            minimumPred = ClusterTuple::INVALID_TUPLE;
            
            // auxiliary data, is precomputed once here
            std::vector<GlobalClusterId> rowTupleGlobal = rowTuple.asVector(ploidy, covMap[pos]);
            std::sort(rowTupleGlobal.begin(), rowTupleGlobal.end());
            
            // this is the tuple into which the rowTuple will be transformed when the best permutation is found
            ClusterTuple bestPerm;
            
            // compare each new tuple with every tuple from previous column
            for (std::pair<ClusterTuple, ClusterEntry> predEntry : m[pos-1-start]) {
                
                // retrieve precomputed sorted vector over global ids
                std::vector<GlobalClusterId> prevTupleGlobal = sortedGlobalTuples[predEntry.first];
                
                // compute optimal switch cost
                Score s = predEntry.second.score + getSwitchCostAllPerms(prevTupleGlobal, rowTupleGlobal);
                
                if (s < minimum) {
                    minExists = true;
                    minimum = s;
                    // minDissim = d;
                    minimumPred = predEntry.first;
                }
            }
            
            if (minExists) {
                // in addition to best score over all predecessors, we need the best permutation of rowTuple to achieve this
                std::vector<GlobalClusterId> prevTuple = sortedGlobalTuples[minimumPred];
                std::vector<uint32_t> residualPosPrev; // positions in previous tuple, which could not be matched to position in current tuple
                std::vector<uint32_t> residualPosCur; // positions in current tuple, which could not be matched to position in previous tuple
                getSwitchCostAllPerms(prevTuple, rowTupleGlobal, residualPosPrev, residualPosCur);
                
                if (residualPosPrev.size() != residualPosCur.size()) {
                    std::cout<<"Residual sizes unequal"<<std::endl;
                    for (auto i = prevTuple.begin(); i != prevTuple.end(); ++i)
                        std::cout << *i << ' ';
                    std::cout<<std::endl;
                    for (auto i = rowTupleGlobal.begin(); i != rowTupleGlobal.end(); ++i)
                        std::cout << *i << ' ';
                    std::cout<<std::endl;
                }
                
                std::vector<GlobalClusterId> bestPermGlobal(minimumPred.asVector(ploidy, covMap[pos-1]));
                for (uint32_t i = 0; i < residualPosCur.size(); i++) {
                    GlobalClusterId residueCur = rowTupleGlobal[residualPosCur[i]];
                    GlobalClusterId residuePrev = prevTuple[residualPosPrev[i]];
                    for (uint32_t j = 0; j < ploidy; j++) {
                        if (bestPermGlobal[j] == residuePrev) {
                            bestPermGlobal[j] = residueCur;
                            break;
                        }
                    }
                }
            
                std::unordered_map<GlobalClusterId, LocalClusterId> globalToLocal;
                for (uint32_t i = 0; i < covMap[pos].size(); i++)
                    globalToLocal[covMap[pos][i]] = i;
                for (uint32_t i = 0; i < ploidy; i++)
                    bestPerm.set(globalToLocal[bestPermGlobal[i]], i);
            } else {
                bestPerm = rowTuple;
            }
            
            Score coverageCost = getCoverageCost(rowTuple, coverage[pos]);
            if (coverageCost != getCoverageCost(bestPerm, coverage[pos])) {
                std::cout<<"Row tuples have unequal coverage cost"<<std::endl;
                std::cout<<rowTuple.asString(ploidy, covMap[pos])<<std::endl;
                std::cout<<bestPerm.asString(ploidy, covMap[pos])<<std::endl;
            }
            
            // report best recursion
            if (minExists) {
                column[bestPerm] = ClusterEntry(minimum + coverageCost, minimumPred);
            } else {
                column[bestPerm] = ClusterEntry(coverageCost, ClusterTuple::INVALID_TUPLE);
            }
            firstUnthreadedPosition = pos+1;
            if (column[bestPerm].score < minimumInColumn) {
                minimumInColumn = column[bestPerm].score;
                minimumTupleInColumn = bestPerm;
                minimumPredTupleInColumn = minimumPred;
            }
            permedTuples.push_back(bestPerm);
        }
        
        // precompute the sorted vectors with global cluster ids for this column (will be reused in next column)
        sortedGlobalTuples.clear();
        for (ClusterTuple t : permedTuples) {
            std::vector<GlobalClusterId> tupleGlobal = t.asVector(ploidy, covMap[pos]);
            std::sort(tupleGlobal.begin(), tupleGlobal.end());
            sortedGlobalTuples[t] = tupleGlobal;
        }
        
        uint32_t allRows = column.size();
        if (symmetryOptimization) {
            // remove non-profitable entries
            std::vector<ClusterTuple> profitableTuples;
            std::vector<ClusterTuple> pivotTuples;
            profitableTuples.push_back(minimumTupleInColumn);
            pivotTuples.push_back(minimumTupleInColumn);
            
            uint32_t rounds = 2;
            for (uint32_t i = 0; i < rounds; i++) {
                for (ClusterTuple t : permedTuples) {
                    bool profitable = true;
                    bool pivot = true;
                    
                    for (ClusterTuple p : pivotTuples) {
                        if (p == t)
                            continue;
                        Score s = getSwitchCostAllPerms(sortedGlobalTuples[p], sortedGlobalTuples[t]);
                        if (column[t].score >= column[p].score + s) {
                            profitable = false;
                            pivot = false;
                            break;
                        } else if (s < (double)(rounds-i) * switchCost) {
                            pivot = false;
                        }
                    }
                    if (profitable) {
                        profitableTuples.push_back(t);
                        if (pivot && pivotTuples.size() < ploidy*ploidy) {
                            pivotTuples.push_back(t);
                        }
                    } else {
                        column.erase(t);
                    }
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
        
        //std::cout<<"Column "<<pos<<": "<<minimumInColumn<<"\t ("<<column.size()<<" rows, before "<<allRows<<")"<<std::endl;
        
        // write column into dp table(s)
        m.push_back(std::unordered_map<ClusterTuple, ClusterEntry>(column.begin(), column.end()));
    }
    
    // discard auxiliary data
    sortedGlobalTuples.clear();
    
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

Score HaploThreader::getSwitchCostAllPerms(const std::vector<GlobalClusterId>& prevTuple, const std::vector<GlobalClusterId>& curTuple) const {
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

Score HaploThreader::getSwitchCostAllPerms(const std::vector<GlobalClusterId>& prevTuple, const std::vector<GlobalClusterId>& curTuple,
                                           std::vector<uint32_t>& residualPosPrev, std::vector<uint32_t>& residualPosCur) const {
    uint32_t pIdx = 0;
    uint32_t cIdx = 0;
    // compare zig-zag-wise over sorted tuples
    while ((pIdx < ploidy) & (cIdx < ploidy)) {
        if (prevTuple[pIdx] == curTuple[cIdx]) {
            pIdx++; cIdx++;
        } else if (prevTuple[pIdx] < curTuple[cIdx]) {
            residualPosPrev.push_back(pIdx);
            pIdx++;
        } else {
            residualPosCur.push_back(cIdx);
            cIdx++;
        }
    }
    for (uint32_t i = pIdx; i < ploidy; i++)
        residualPosPrev.push_back(i);
    for (uint32_t j = cIdx; j < ploidy; j++)
        residualPosCur.push_back(j);
    
    return switchCost*residualPosPrev.size() + affineSwitchCost*(residualPosPrev.size() > 0);
}

std::vector<ClusterTuple> HaploThreader::computeGenotypeConformTuples (const std::vector<GlobalClusterId>& clusters,
                                                                       const std::vector<uint32_t>& consensus,
                                                                       const std::unordered_map<uint32_t, uint32_t>& genotype) const {
    std::vector<ClusterTuple> perfect = getGenotypeConformTuples (clusters, consensus, genotype);
    if (perfect.size() > 0) {
        return perfect;
    } else {
        const std::vector<uint32_t> consensusDummy(clusters.size(), 0);
        const std::unordered_map<uint32_t, uint32_t> genotypeDummy({ {0, ploidy} });
        std::vector<ClusterTuple> noMatch = getGenotypeConformTuples (clusters, consensusDummy, genotypeDummy);
        return noMatch;
    }
}

std::vector<ClusterTuple> HaploThreader::getGenotypeConformTuples (const std::vector<GlobalClusterId>& clusters, const std::vector<uint32_t>& consensus, 
                                                                const std::unordered_map<uint32_t, uint32_t>& genotype) const {
    std::vector<ClusterTuple> conformTuples;
    
    // convert genotype map to vector
    uint32_t maxAllele = 0;
    for (std::pair<uint32_t, uint32_t> entry : genotype) {
        maxAllele = std::max(maxAllele, entry.first+1);
    }
    std::vector<uint32_t> genotypeVec(maxAllele, 0);
    for (std::pair<uint32_t, uint32_t> entry : genotype) {
        genotypeVec[entry.first] = entry.second;
    }
    // split clusters by consensus-allele
    std::vector<std::vector<LocalClusterId>> clusterGroups(maxAllele, std::vector<LocalClusterId>(0));
    for (LocalClusterId i = 0; i < clusters.size(); i++) {
        clusterGroups[consensus[i]].push_back(i);
    }
    
    // if genotype not reachable: return empty vector
    for (uint32_t allele = 0; allele < maxAllele; allele++) {
        if (genotypeVec[allele] > 0 && clusterGroups[allele].size() == 0)
            return conformTuples;
    }

    /*
     * For each allele, store a vector, which contains all combinations (as vectors) from respective local cluster ids
     */
    std::vector<std::vector<std::vector<LocalClusterId>>> alleleWiseCombs;
    
    for (uint32_t allele = 0; allele < maxAllele; allele++) {
        // create vector of combinations for current allele
        std::vector<std::vector<LocalClusterId>> combsOfAllele;
        
        uint32_t numElem = genotypeVec[allele];
        uint32_t maxElem = clusterGroups[allele].size();
        if (numElem > 0) {
            std::vector<uint32_t> v(numElem, 0);
            
            /*
             * store vector, until maxElem-1 is in the last field, because then we must just
             * have stored the vector [maxElem-1, maxElem-1, ..., maxElem-1] and we are finished.
             */
            while (v[numElem-1] < maxElem) {
                // translate to local cluster ids and store in combsOfAllele
                std::vector<uint32_t> c(numElem, 0);
                for (uint32_t i = 0; i < numElem; i++)
                    c[i] = clusterGroups[allele][v[i]];
                combsOfAllele.push_back(c);
                
                // increment like a counter
                v[0]++;
                
                // if element i-1 overflowed, increase element i by 1 (which then might also overflow and so on)
                for (uint32_t i = 1; i < numElem; i++)
                    if (v[i-1] >= maxElem)
                        v[i]++;
                
                // any element i-1 which overflowed will be set to its minimum, i.e. the value of element i
                for (uint32_t i = numElem-1; i > 0; i--)
                    if (v[i-1] >= maxElem)
                        v[i-1] = v[i];
            }

        }
        alleleWiseCombs.push_back(combsOfAllele); // added vector may be empty
    }
    
    /*
     * Next step is to generate all combinations of size ploidy, so all combinations from
     * allele 0 times all combinations from allele 1 times ...
     */
    
    /*
     * This vector contains one entry per allele. indices[i] = j means, that we choose the j-th combination
     * from allele i to construct the next element of our result
     */
    std::vector<uint32_t> indices(maxAllele, 0);
    while (indices[maxAllele-1] < alleleWiseCombs[maxAllele-1].size()) {
        /*
         * As long as our last index is still valid, we can copy the current combination to result
         */
        std::vector<uint32_t> x;
        for (uint32_t allele = 0; allele < maxAllele; allele++) {
            if (alleleWiseCombs[allele].size() == 0)
                continue;
            // append the indices[allele]-th combination from alleleWiseCombs[allele] to x
            for (uint32_t c : alleleWiseCombs[allele][indices[allele]])
                x.push_back(c);
        }
        
        // append to solution
        conformTuples.push_back(ClusterTuple(x));
        
        // increment like a counter
        indices[0]++;
        
        /*
         * If element i-1 overflowed, increase element i by 1 (which then might also overflow and so on)
         * and set element i-1 to zero
         */
        for (uint32_t i = 1; i < maxAllele; i++) {
            // second condition indices[i-1] > 0 is necessary, because there might be alleles with no combinations
            if (indices[i-1] >= alleleWiseCombs[i-1].size() && indices[i-1] > 0) {
                indices[i]++;
                indices[i-1] = 0;
            }
        }
    }
    
    return conformTuples;
}
