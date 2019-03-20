#include "HaploThreader.h"
#include <limits>
#include <algorithm>
#include <unordered_set>

// using Score = HaploThreader::Score;
// using GlobalClusterId = HaploThreader::GlobalClusterId;
// using LocalClusterId = HaploThreader::LocalClusterId;
// using ClusterTuple = HaploThreader::ClusterTuple;
// using Position = HaploThreader::Position;

constexpr uint64_t ClusterTuple::TUPLE_MASKS[];
const ClusterTuple ClusterTuple::INVALID_TUPLE = ClusterTuple((TupleCode)-1);

HaploThreader::HaploThreader (uint32_t ploidy, double switchCost, double affineSwitchCost) :
    ploidy(ploidy),
    switchCost(switchCost),
    affineSwitchCost(affineSwitchCost)
{
    ;
}

std::vector<std::vector<GlobalClusterId>> HaploThreader::computePaths (Position num_vars, 
                    std::vector<std::vector<GlobalClusterId>> &covMap,
                    std::vector<std::vector<double>> &coverage, 
                    std::vector<std::vector<uint32_t>> consensus,
                    std::vector<uint32_t> genotypes
                   ) const {
    
    // setup data structures
    std::unordered_map<ClusterTuple, Score> pred;
    std::vector<std::unordered_map<ClusterTuple, Score>> d;
    std::vector<std::unordered_map<ClusterTuple, ClusterTuple>> b;
    Position firstUnthreadedPosition = 0;
    std::vector<std::vector<GlobalClusterId>> path;
    
    // compute genotype conform sets (first one without, second one with permutations)
    std::vector<std::vector<ClusterTuple>> confSets;
    std::vector<std::vector<ClusterTuple>> confSetsExt;
    for (uint32_t pos = 0; pos < num_vars; pos++) {
        std::vector<ClusterTuple> perfect = getGenotypeConformTuples (covMap[pos], consensus[pos], genotypes[pos], 0, true);
        if (perfect.size() > 0) {
//             std::cout<<"Position "<<pos<<": "<<perfect.size()<<" perfect matches"<<std::endl;
            confSetsExt.push_back(perfect);
            confSets.push_back(getGenotypeConformTuples(covMap[pos], consensus[pos], genotypes[pos], 0, false));
        } else {
            std::vector<ClusterTuple> oneOff = getGenotypeConformTuples (covMap[pos], consensus[pos], genotypes[pos], 1, true);
            if (oneOff.size() > 0) {
//                 std::cout<<"Position "<<pos<<": "<<oneOff.size()<<" approximate matches"<<std::endl;
                confSetsExt.push_back(oneOff);
                confSets.push_back(getGenotypeConformTuples(covMap[pos], consensus[pos], genotypes[pos], 1, false));
            } else {
                std::vector<ClusterTuple> noMatch = getCombinations(covMap[pos].size(), true);
//                 std::cout<<"Position "<<pos<<": "<<noMatch.size()<<" nonconform matches"<<std::endl;
                confSetsExt.push_back(noMatch);
                confSets.push_back(getCombinations(covMap[pos].size(), false));
            }
        }
    }
    
    // initialize first column
    std::vector<ClusterTuple> rowTuples = confSetsExt[0];
    std::unordered_map<ClusterTuple, Score> scoreColumn;
    std::unordered_map<ClusterTuple, ClusterTuple> btColumn;
    if (rowTuples.size() == 0) {
        std::cout<<"First variant has no clusters!"<<std::endl;
        return path;
    }
    Score minimumInColumn = std::numeric_limits<Score>::infinity();
    ClusterTuple minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
    ClusterTuple minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;
    for (ClusterTuple t : rowTuples) {
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
        // reset variables
        rowTuples.clear();
        scoreColumn.clear();
        btColumn.clear();
        
        Score minimum = std::numeric_limits<Score>::infinity();
        ClusterTuple minimumPred = ClusterTuple::INVALID_TUPLE;
        bool minExists = false;
        minimumInColumn = std::numeric_limits<Score>::infinity();
        minimumTupleInColumn = ClusterTuple::INVALID_TUPLE;
        minimumPredTupleInColumn = ClusterTuple::INVALID_TUPLE;
        
        // iterate over rows
        rowTuples = confSetsExt[pos];
        std::cout<<"Threading haplotypes through clusters .. ("<<pos<<"/"<<num_vars<<")\r"<<std::flush;
        for (ClusterTuple rowTuple : rowTuples) {
            minimum = std::numeric_limits<Score>::infinity();
            minimumPred = ClusterTuple::INVALID_TUPLE;
            
            // get relevant predecessors
            for (ClusterTuple predTup : confSetsExt[pos-1]) {
                Score s = d[pos-1][predTup] + getSwitchCost(rowTuple, predTup, covMap[pos], covMap[pos-1]);
                if (s < minimum) {
                    minExists = true;
                    minimum = s;
                    minimumPred = predTup;
                }
            }
            
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
    std::cout<<"Start backtracking in column "<<(firstUnthreadedPosition-1)<<std::endl;
    
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
        if (currentRow.asVector(ploidy, covMap[firstUnthreadedPosition-1]).size() == 0) {
            std::cout<<"Problem occured at position "<<(firstUnthreadedPosition-1)<<" in row "<<currentRow.asString(ploidy, covMap[firstUnthreadedPosition-1])<<std::endl;
        }
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
    }
    
    std::vector<ClusterTuple> uniqueTuples(tuples.begin(), tuples.end());
    return uniqueTuples;
}

std::vector<std::vector<uint32_t>> HaploThreader::getPermutationsVector(uint32_t maxElem) const {
    std::vector<std::vector<uint32_t>> tuples;
    std::vector<LocalClusterId> curPerm(ploidy, 0);
    
    // enumerate all permutations of clusters
    while (curPerm[ploidy-1] < maxElem) {
        tuples.push_back(curPerm);
        
        // increment permutation iterator
        curPerm[0]++;
        for (uint32_t i = 1; i < ploidy; i++) {
            if (curPerm[i-1] >= maxElem) {
                curPerm[i]++;
                curPerm[i-1] = 0;
            }
        }
    }
    return tuples;
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
