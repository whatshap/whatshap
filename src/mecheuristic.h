#ifndef MECHEURISTIC_H
#define MECHEURISTIC_H

#include "polyphase/allelematrix.h"
#include "readset.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <iostream>
#include <string> 
#include <sstream>
#include <algorithm>
#include <limits>
    
typedef uint32_t Position;
typedef double MecScore;
typedef std::vector<bool> Bipartition;
typedef uint32_t ReadId;
typedef int8_t Allele;
typedef uint16_t RowIndex;

struct BipartitionItem {
    Bipartition bp;
    Bipartition bpNew;
    MecScore score;
    RowIndex bt;
    
    BipartitionItem(const Bipartition bp, const Bipartition bpNew, const MecScore score, RowIndex bt) :
    bp(bp),
    bpNew(bpNew),
    score(score),
    bt(bt){}

    bool operator<(const BipartitionItem& other) const {
        return score < other.score;
    }
};

class MecHeuristic {

public:

    MecHeuristic(uint32_t rowLimit, bool weighted, bool allHet);
    std::vector<std::vector<Allele>> computeHaplotypes(ReadSet* rs) const;

private:
    uint32_t rowLimit;
    bool weighted;
    bool allHet;
    

    std::vector<std::pair<Bipartition, MecScore>> generateExtensions(std::vector<Allele>& alleleList, bool symmetric) const;
    std::vector<std::pair<Bipartition, MecScore>> generateExtensions(std::vector<Allele>& alleleList, bool reverse, uint32_t limit) const;
    std::vector<std::vector<uint32_t>> generateCombinations(const uint32_t n, const uint32_t k) const;
    bool bpEqual(Bipartition a, Bipartition b) const;
    
    void printColumnInfo(Position p, std::vector<ReadId> startIndex, std::vector<BipartitionItem> col) const;
};

#endif
