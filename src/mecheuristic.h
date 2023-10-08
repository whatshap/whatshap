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
typedef float MecScore;
typedef std::vector<bool> Bipartition;
typedef uint32_t ReadId;
typedef int8_t Allele;
typedef uint16_t RowIndex;
    
struct BipartitionItem {
    Bipartition bp;
    Bipartition bpNew;
    std::vector<std::vector<MecScore>> balances;
    MecScore score;
    RowIndex bt;
    
    BipartitionItem(const Bipartition bip, const Bipartition bpNew, std::vector<std::vector<MecScore>> balances, const MecScore score, RowIndex bt) :
    bp(bip),
    bpNew(bpNew),
    balances(balances),
    score(score),
    bt(bt) {
        bp.insert(bp.end(), bpNew.begin(), bpNew.end());
    }
    
    bool operator<(const BipartitionItem& other) const {
        return score < other.score;
    }
};

class MecHeuristic {

public:

    MecHeuristic(uint32_t rowLimit, bool allHet);
    std::vector<std::vector<Allele>> computeHaplotypes(ReadSet* rs, ReadSet* output) const;

private:
    uint32_t rowLimit;
    bool allHet;
    

    std::vector<std::pair<Bipartition, MecScore>> generateExtensions(std::vector<Allele>& alleleList, std::vector<MecScore>& qualityList, bool symmetric) const;
    std::vector<std::pair<Bipartition, MecScore>> generateExtensions(std::vector<Allele>& alleleList, std::vector<MecScore>& qualityList, bool reverse, uint32_t limit) const;
    std::vector<std::vector<uint32_t>> generateCombinations(const uint32_t n, const uint32_t k) const;
    bool bpEqual(Bipartition a, Bipartition b) const;
    Allele getAllele(ReadSet* rs, ReadId rid, uint32_t genPos) const;
    MecScore addBalance(std::vector<MecScore>& basis, std::vector<MecScore>& coBasis, std::vector<MecScore>& add) const;
    void printColumnInfo(Position p, std::vector<ReadId>& startIndex, std::vector<BipartitionItem>& col) const;
};

#endif
