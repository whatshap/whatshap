#include "switchflipcalculator.h"
#include "staticsparsegraph.h"
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <random>

constexpr uint64_t Permutation::TUPLE_MASKS[];
const Permutation Permutation::INVALID = Permutation(0xf000000000000000);

SwitchFlipCalculator::SwitchFlipCalculator (uint32_t ploidy, double switchCost, double flipCost) :
    ploidy(ploidy),
    switchCost(switchCost),
    flipCost(flipCost)
{
}

std::pair<Score, Score> SwitchFlipCalculator::compare (const std::vector<std::vector<uint32_t>>& phasing0,
                    const std::vector<std::vector<uint32_t>>& phasing1,
                    std::vector<uint32_t>& switchesInColumn,
                    std::vector<std::vector<uint32_t>>& flippedHapsInColumn,
                    std::vector<std::vector<uint32_t>>& permInColumn
                    ) const {
    
    // setup data structures
    std::vector<std::unordered_map<Permutation, PermutationEntry>> m;
    std::vector<Permutation> permutations = getPermutations();
    Position numVars = phasing0.size();

    // initialize first column
    std::unordered_map<Permutation, PermutationEntry> column;
    Score minimumInColumn = std::numeric_limits<Score>::infinity();
    Permutation minimumPermInColumn = Permutation::INVALID;
    Permutation minimumPredPermInColumn = Permutation::INVALID;
    uint32_t allEntries = 0;
    uint32_t keptEntries = 0;
    for (Permutation p : permutations) {
        column[p] = PermutationEntry(this->flipCost * getNumFlips(p, phasing0[0], phasing1[0]), Permutation::INVALID);
        if (column[p].score < minimumInColumn) {
            minimumInColumn = column[p].score;
            minimumPermInColumn = p;
        }
    }
    
    m.push_back(std::unordered_map<Permutation, PermutationEntry>(column.begin(), column.end()));
    allEntries += permutations.size();
    keptEntries += column.size();
    
    // iterate over positions
    for (Position pos = 1; pos < numVars; pos++) {
        
        // reset variables
        column.clear();
        Score minimum = std::numeric_limits<Score>::infinity();
        Permutation minimumPred = Permutation::INVALID;
        bool minExists = false;
        minimumInColumn = std::numeric_limits<Score>::infinity();
        minimumPermInColumn = Permutation::INVALID;
        minimumPredPermInColumn = Permutation::INVALID;
        
        // iterate over rows
        for (Permutation rowPerm : permutations) {
            minimum = std::numeric_limits<Score>::infinity();
            minimumPred = Permutation::INVALID;
            
            // iterate over previous rows
            for (std::pair<Permutation, PermutationEntry> predEntry : m[pos-1]) {
                Score s = predEntry.second.score + switchCost * getNumSwitches(rowPerm, predEntry.first);
                if (s < minimum) {
                    minExists = true;
                    minimum = s;
                    minimumPred = predEntry.first;
                }
            }
            
            Score flipCost = this->flipCost * getNumFlips(rowPerm, phasing0[pos], phasing1[pos]);
            
            // report best recursion
            if (minExists) {
                column[rowPerm] = PermutationEntry(minimum + flipCost, minimumPred);
            } else {
                column[rowPerm] = PermutationEntry(flipCost, Permutation::INVALID);
            }
            if (column[rowPerm].score < minimumInColumn) {
                minimumInColumn = column[rowPerm].score;
                minimumPermInColumn = rowPerm;
                minimumPredPermInColumn = minimumPred;
            }
        }
        
        // remove non-profitable entries
        std::vector<Permutation> profitableTuples;
        std::vector<Permutation> openTuples;
        for (Permutation t : permutations) {
            if (column[t].score <= minimumInColumn) {
                profitableTuples.push_back(t);
            } else {
                openTuples.push_back(t);
            }
        }
        
        for (Permutation t : openTuples) {
            bool profitable = true;
            for (Permutation p : profitableTuples) {
                if (column[t].score >= column[p].score + switchCost * getNumSwitches(t, p)) {
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
        
        // write column into dp table(s)
        m.push_back(std::unordered_map<Permutation, PermutationEntry>(column.begin(), column.end()));
        allEntries += permutations.size();
        keptEntries += column.size();
    }
    
    // backtracking start
    Score flips = 0.0;
    Score switches = 0.0;
    Permutation currentRow = Permutation::INVALID;
    Score minimum = std::numeric_limits<Score>::infinity();
    for (std::pair<Permutation, PermutationEntry> entry : m[numVars-1]) {
        if (entry.second.score < minimum) {
            minimum = entry.second.score;
            currentRow = entry.first;
        }
    }
    if (currentRow == Permutation::INVALID) {
        std::cout<<"No minimum at end of compared block!"<<std::endl;
        return std::pair<Score, Score>(std::numeric_limits<Score>::infinity(), std::numeric_limits<Score>::infinity());
    } else {
        if (currentRow.asVector(ploidy).size() == 0)
            std::cout<<"Problem occured at variant "<<(numVars-1)<<" in row "<<currentRow.asString(ploidy)<<std::endl;
        permInColumn.push_back(currentRow.asVector(ploidy));
        Score localFlips = getNumFlips(currentRow, phasing0[numVars-1], phasing1[numVars-1]);
        Score localSwitches = getNumSwitches(currentRow, m[numVars-1][currentRow].pred);
        std::vector<uint32_t> localFlippedHaps = getFlippedHaps(currentRow, phasing0[numVars-1], phasing1[numVars-1]);
        flippedHapsInColumn.push_back(localFlippedHaps);
        switchesInColumn.push_back(localSwitches);
        flips += localFlips;
        switches += localSwitches;
    }
    
    // backtracking iteration
    for (Position pos = numVars-2; pos < numVars; pos--) {
        currentRow = m[pos+1][currentRow].pred;
        if (currentRow.asVector(ploidy).size() == 0) {
            std::cout<<"Problem occured at variant "<<(pos)<<" in row "<<currentRow.asString(ploidy)<<std::endl;
            return std::pair<Score, Score>(std::numeric_limits<Score>::infinity(), std::numeric_limits<Score>::infinity());
        } else {
            permInColumn.push_back(currentRow.asVector(ploidy));
            Score localFlips = getNumFlips(currentRow, phasing0[pos], phasing1[pos]);
            Score localSwitches = pos == 0 ? 0.0 : getNumSwitches(currentRow, m[pos][currentRow].pred);
            std::vector<uint32_t> localFlippedHaps = getFlippedHaps(currentRow, phasing0[pos], phasing1[pos]);
            flippedHapsInColumn.push_back(localFlippedHaps);
            switchesInColumn.push_back(localSwitches);
            flips += localFlips;
            switches += localSwitches;
        }
    }
    
    // reverse as we constructed it back to front
    std::reverse(permInColumn.begin(), permInColumn.end());
    std::reverse(flippedHapsInColumn.begin(), flippedHapsInColumn.end());
    std::reverse(switchesInColumn.begin(), switchesInColumn.end());
    
    return std::pair<Score, Score>(switches, flips);
}

Score SwitchFlipCalculator::getNumFlips(Permutation permutation,
                              const std::vector<uint32_t>& phase0, 
                              const std::vector<uint32_t>& phase1) const {
    uint32_t diffCount = 0;
    for (uint32_t i = 0; i < ploidy; i++) {
        diffCount += phase0[permutation.get(i)] != phase1[i];
    }
    return (Score)diffCount;
}

std::vector<uint32_t> SwitchFlipCalculator::getFlippedHaps(Permutation permutation,
                              const std::vector<uint32_t>& phase0, 
                              const std::vector<uint32_t>& phase1) const {
    std::vector<uint32_t> diffs;
    for (uint32_t i = 0; i < ploidy; i++) {
        if (phase0[permutation.get(i)] != phase1[i])
            diffs.push_back(i);
    }
    return diffs;
}

Score SwitchFlipCalculator::getNumSwitches(const Permutation p1, const Permutation p2) const {
    // count differences between permutations by some bit shift magic
    PermutationCode x = p1.code ^ p2.code;
    PermutationCode y = (x & 0x0111111111111111) | ((x >> 3) & 0x0111111111111111) | ((x >> 2) & 0x0111111111111111) | ((x >> 1) & 0x0111111111111111);
    
    // every 4-bit block contains now either a zero (permutation position equal) or a 1 (0001). Do popcount.
    return StaticSparseGraph::popcount(y);
}

std::vector<Permutation> SwitchFlipCalculator::getPermutations() const {        
    std::vector<Permutation> perms;
    Haplotype curPerm[ploidy];
    
    for (uint32_t i = 0; i < ploidy; i++) {
        curPerm[i] = i;
    }
    
    do {
        perms.push_back(Permutation(curPerm, ploidy));
    } while (std::next_permutation(curPerm, curPerm + ploidy));
    
    return perms;
}
