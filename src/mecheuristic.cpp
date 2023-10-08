#include "mecheuristic.h"
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <cassert>
#include "binomial.h"
#include "assert.h"

MecHeuristic::MecHeuristic(uint32_t rowLimit, bool allHet) :
    rowLimit(rowLimit),
    allHet(allHet)
{
}

std::vector<std::vector<Allele>> MecHeuristic::computeHaplotypes(ReadSet* rs, ReadSet* output) const {
    ReadId m = rs->size();
    std::vector<uint32_t>* posptr = rs->get_positions();
    std::vector<uint32_t> positions(posptr->begin(), posptr->end());
    delete posptr;
    Position n = positions.size();
    std::unordered_map<uint32_t, Position> posMap;
    for (uint32_t i = 0; i < n; i++)
        posMap[positions[i]] = i;
    std::vector<ReadId> startIndex;
    std::vector<BipartitionItem> lastCol;
    std::vector<ReadId> active;
    std::vector<std::vector<RowIndex>> mBt (n, std::vector<RowIndex>());
    std::vector<std::pair<uint32_t, Bipartition>> mBp;
        
    // compute index of first read starting at position p (for p = 0, 1, ..., n - 1, n)
    startIndex.push_back(0);
    ReadId r = 0;
    for (Position p = 0; p < n; p++) {
        while (r < m && posMap[rs->get(r)->firstPosition()] <= p)
            r++;
        startIndex.push_back(r);
    }
    
    // create empty biparition as base for first column
    lastCol.emplace_back(Bipartition(), Bipartition(), std::vector<std::vector<MecScore>>(2, std::vector<MecScore>(1, 0)), 0, 0);
    
    // DP recurrence
    Position right = 0;
    for (Position p = 0; p < n; p++) {
        // compute intersection of active reads with previous column
        std::vector<ReadId> activeLast(active.begin(), active.end());
        std::vector<uint8_t> kept;
        active.clear();
        for (uint32_t i = 0; i < activeLast.size(); i++) {
            ReadId r = activeLast[i];
            if (posMap[rs->get(r)->lastPosition()] >= p) {
                active.push_back(r);
                kept.push_back(i);
            }
        }
        
        // copy bipartitions, but without the lost reads
        std::vector<MecScore> scoreTemp;
        std::vector<Bipartition> bpTemp;
        std::vector<Bipartition> bpNewTemp;
        std::vector<std::vector<std::vector<MecScore>>> balancesTemp;
        std::vector<RowIndex> btTemp;
        std::unordered_map<uint64_t, std::vector<uint32_t>> buckets;
        // i := index of BP in last column
        // j := index of BP in current column (at least those that exist yet)
        for (uint32_t i = 0; i < lastCol.size(); i++) {     
            Bipartition lastBp = lastCol[i].bp;
            Bipartition b(kept.size(), false);
            uint32_t bucket = 0;
            for (uint32_t a = 0; a < kept.size(); a++) {
                bool bo = lastBp[kept[a]];
                b[a] = bo;
                bucket = (bucket << 1) + bo;
            }
            bool duplicate = false;
            
            // use buckets (last 64 bits of bipartition bitvector) to avoid full pairwise check
            if (buckets.find(bucket) == buckets.end())
                buckets[bucket] = std::vector<uint32_t>();
            for (uint32_t j : buckets[bucket]) {
                if (bpEqual(b, bpTemp[j])) {
                    duplicate = true;
                    if (lastCol[i].score < lastCol[btTemp[j]].score) {
                        balancesTemp[j].clear();
                        for (auto& a: lastCol[btTemp[j]].balances)
                            balancesTemp[j].emplace_back(a.begin() + 1, a.end());
                        scoreTemp[j] = lastCol[i].score;
                        btTemp[j] = i;
                    }
                    break;
                }
            }
            if (!duplicate) {
                buckets[bucket].push_back(bpTemp.size());
                bpTemp.push_back(b);
                bpNewTemp.emplace_back(0);
                balancesTemp.emplace_back(0);
                for (auto& a: lastCol[i].balances)
                    balancesTemp[balancesTemp.size() - 1].emplace_back(a.begin() + 1, a.end());
                scoreTemp.push_back(lastCol[i].score);
                btTemp.push_back(i);
            }
        }
        buckets.clear();
        
        // from here on, we have bpTemp, balanceTemp, scoreTemp and btTemp available as basis for extensions
        
        // get last position of new reads and extend balance vectors
        for (ReadId r = startIndex[p]; r < startIndex[p + 1];  r++)
            right = std::max(right, posMap[rs->get(r)->lastPosition()]);
        for (auto& a: balancesTemp)
            for (auto& b: a)
                b.resize(right + 1 - p, 0);
        
        // for every new read: generate two extensions for every solution and keep only rowLimit best
        for (ReadId r = startIndex[p]; r < startIndex[p + 1];  r++) {
            active.push_back(r);
            // generate balance vector of read
            std::vector<MecScore> balance(right + 1 - p, 0);
            for (int32_t i = 0; i < rs->get(r)->getVariantCount(); i++) {
                Position o = posMap[rs->get(r)->getPosition(i)] - p;
                Allele a = rs->get(r)->getAllele(i);
                MecScore q = rs->get(r)->getVariantQuality(i);
                balance[o] += q * a - q * (1 - a);
            }
            std::vector<MecScore> scores;
            uint32_t solEnd = bpTemp.size();
            for (uint32_t sol = 0; sol < solEnd; sol++) {
                // create new solution with read in partiton 1
                bpTemp.emplace_back(bpTemp[sol].begin(), bpTemp[sol].end());
                bpNewTemp.emplace_back(bpNewTemp[sol].begin(), bpNewTemp[sol].end());
                balancesTemp.emplace_back(balancesTemp[sol].begin(), balancesTemp[sol].end());
                scoreTemp.push_back(scoreTemp[sol]);
                btTemp.push_back(btTemp[sol]);
                
                scoreTemp[scoreTemp.size() - 1] += addBalance(balancesTemp[balancesTemp.size() - 1][1], balancesTemp[balancesTemp.size() - 1][0], balance);
                scores.push_back(scoreTemp[scoreTemp.size() - 1]);
                bpNewTemp[bpTemp.size() - 1].push_back(true);

                // extend existing solution with read in partition 0
                
                scoreTemp[sol] += addBalance(balancesTemp[sol][0], balancesTemp[sol][1], balance);
                scores.push_back(scoreTemp[sol]);
                bpNewTemp[sol].push_back(false);
            }
            
            // prune solution set if too large
            if (scores.size() > rowLimit) {
                std::sort(scores.begin(), scores.end());
                MecScore tooHigh = scores.size() > rowLimit ? scores[rowLimit] : std::numeric_limits<MecScore>::infinity();
                std::vector<uint32_t> kept;
                for (uint32_t sol = 0; sol < solEnd * 2; sol++)
                    if (scoreTemp[sol] < tooHigh || scoreTemp[sol] == scores[0])
                        kept.push_back(sol);
                for (uint32_t sol = 0; sol < kept.size(); sol++) {
                    bpTemp[sol] = bpTemp[kept[sol]];
                    bpNewTemp[sol] = bpNewTemp[kept[sol]];
                    balancesTemp[sol] = balancesTemp[kept[sol]];
                    scoreTemp[sol] = scoreTemp[kept[sol]];
                    btTemp[sol] = btTemp[kept[sol]];
                }
                bpTemp.resize(kept.size());
                bpNewTemp.resize(kept.size());
                balancesTemp.resize(kept.size());
                scoreTemp.resize(kept.size());
                btTemp.resize(kept.size());
            }
        }
        
        // .. and then construct solutions only if score is sufficient (performance reasons)
        lastCol.clear();
        Bipartition btVector;
        for (uint32_t i = 0; i < bpTemp.size(); i++) {
            lastCol.emplace_back(bpTemp[i], bpNewTemp[i], balancesTemp[i], scoreTemp[i], btTemp[i]);
            mBt[p].push_back(btTemp[i]);
            btVector.insert(btVector.end(), bpNewTemp[i].begin(), bpNewTemp[i].end());
        }
        mBp.emplace_back(startIndex[p + 1] - startIndex[p], btVector);
        btVector.clear();
        printColumnInfo(p, startIndex, lastCol);
    }
    
    // find best score in last column
    MecScore s = std::numeric_limits<MecScore>::infinity();
    RowIndex ri = 0;
    for (uint32_t i = 0; i < lastCol.size(); i++) {
        if (lastCol[i].score < s) {
            s = lastCol[i].score;
            ri = i;
        }
    }
    
    // backtracking
    Bipartition full(m, false);
    for (Position p = n - 1;p < p + 1; p--) {
        ReadId offset = startIndex[p];
        ReadId newCount = mBp[p].first;
        Bipartition current(mBp[p].second.begin() + newCount * ri, mBp[p].second.begin() + newCount * (ri + 1));
        for (uint32_t i = 0; i < current.size(); i++)
            full[offset + i] = current[i];
        ri = mBt[p][ri];
    }
    
    // get allele votes
    MecScore actualScore = 0;
    std::vector<std::vector<MecScore>> vote(n, std::vector<MecScore>(4, 0));
    for (ReadId ri = 0; ri < m; ri++) {
        uint32_t part = full[ri];
        Read* r = rs->get(ri);
        for (int32_t i = 0; i < r->getVariantCount(); i++) {
            Allele a = r->getAllele(i);
            MecScore q = r->getVariantQuality(i);
            if (a >= 0)
                vote[posMap[r->getPosition(i)]][2 * part + a] += q;
        }
    }
    
    // compute consensus
    std::vector<std::vector<Allele>> haps(2, std::vector<Allele>(n, 0));
    for (Position p = 0; p < n; p++) {
        if (allHet) {
            haps[0][p] = (vote[p][0] + vote[p][3]) < (vote[p][1] + vote[p][2]);
            haps[1][p] = 1 - haps[0][p];
            actualScore += std::min(vote[p][0] + vote[p][3], vote[p][1] + vote[p][2]);
        } else {
            haps[0][p] = vote[p][0] < vote[p][1];
            haps[1][p] = vote[p][2] < vote[p][3];
            actualScore += std::min(vote[p][0], vote[p][1]);
            actualScore += std::min(vote[p][2], vote[p][3]);
        }
    }
    std::cout<<"Actual MEC score is "<<actualScore<<std::endl;
    
    Read* read0 = new Read("superread_0_0", -1, -1, 0);
    Read* read1 = new Read("superread_1_0", -1, -1, 0);
    
    for (uint32_t i = 0; i < n; i++)  {
        read0->addVariant(positions[i], haps[0][i], 30);
        read1->addVariant(positions[i], haps[1][i], 30);
    }
    
    output->add(read0);
    output->add(read1);
    
    return haps;
}


bool MecHeuristic::bpEqual(Bipartition a, Bipartition b) const {
    uint32_t m = a.size();
    if (m != b.size())
        return false;
    for (uint32_t i = 0; i < m; i++)
        if (a[i] != b[i])
            return false;
    return true;
}


Allele MecHeuristic::getAllele(ReadSet* rs, ReadId rid, uint32_t genPos) const {
    assert(rid < rs->size());
    Read* r = rs->get(rid);
    for (int32_t i = 0; i < r->getVariantCount(); i++)
        if (r->getPosition(i) == (int32_t)genPos)
            return r->getAllele(i);
    return -1;
}


MecScore MecHeuristic::addBalance(std::vector<MecScore>& basis, std::vector<MecScore>& coBasis, std::vector<MecScore>& add) const {
    MecScore penalty = 0;
    for (uint32_t i = 0; i < add.size(); i++) {
        if (allHet) {
            if (add[i] <= 0)
                penalty += std::min(-add[i], std::max(basis[i] - coBasis[i], (MecScore)0));
            else
                penalty += std::min(add[i], std::max(coBasis[i] - basis[i], (MecScore)0));
        } else {
            if (basis[i] * add[i] < 0)
                // penalty = how much did the added balance move the old one towards zero
                penalty += std::min(std::abs(basis[i]), std::abs(add[i]));
            // update balance vector
        }
        basis[i] += add[i];
    }
    return penalty;
}


void MecHeuristic::printColumnInfo(Position p, std::vector< ReadId >& startIndex, std::vector< BipartitionItem >& col) const
{
    MecScore s = col[0].score;
    for (uint32_t i = 1; i < col.size(); i++)
        if (col[i].score < s)
            s = col[i].score;
    std::cout<<"Column "<<p<<": ["<<startIndex[p]<<": "<<startIndex[p + 1] - 1<<"] with "<<col.size()<<" bipartitions and score "<<s<<std::endl;
}
