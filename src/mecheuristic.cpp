#include "mecheuristic.h"
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <cassert>
#include "binomial.h"
#include "assert.h"

MecHeuristic::MecHeuristic(uint32_t rowLimit, bool weighted, bool allHet) :
    rowLimit(rowLimit),
    weighted(weighted),
    allHet(allHet)
{
}

std::vector<std::vector<Allele>> MecHeuristic::computeHaplotypes(ReadSet* rs, ReadSet* output) {
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
    
    // fill first column of DP
    std::vector<Allele> newAlleles;
    for (ReadId r = 0; r < startIndex[1]; r++) {
        newAlleles.push_back(getAllele(rs, r, positions[0]));
        active.push_back(r);
    }
    std::vector<std::pair<Bipartition, MecScore>> extensions = generateExtensions(newAlleles, false);
    Bipartition btVector;
    for (uint32_t i = 0; i < extensions.size(); i++) {
        auto a = extensions[i];
        lastCol.emplace_back(a.first, a.first, a.second, 0);
        mBt[0].push_back(0);
        btVector.insert(btVector.end(), a.first.begin(), a.first.end());
    }
    mBp.emplace_back(newAlleles.size(), btVector);
    btVector.clear();
    printColumnInfo(0, startIndex, lastCol);
    
    // DP recurrence
    for (Position p = 1; p < n; p++) {
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
        std::unordered_map<uint64_t, std::vector<uint32_t>> buckets;
        std::vector<RowIndex> btTemp;
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
                        scoreTemp[j] = lastCol[i].score;
                        btTemp[j] = i;
                    }
                    break;
                }
            }
            if (!duplicate) {
                buckets[bucket].push_back(bpTemp.size());
                bpTemp.push_back(b);
                scoreTemp.push_back(lastCol[i].score);
                btTemp.push_back(i);
            }
        }
        buckets.clear();
        
        // precompute possible extensions based on new reads
        std::vector<Allele> newAlleles;
        std::vector<ReadId> ingoing;
        for (ReadId r = startIndex[p]; r < startIndex[p + 1]; r++) {
            newAlleles.push_back(getAllele(rs, r, positions[p]));
            ingoing.push_back(r);
        }
        std::vector<std::pair<Bipartition, MecScore>> extensions = generateExtensions(newAlleles, true);

        // precompute allele distribution (0 = #0 in part. 0, 1 = #1 in part 0, 2 = #0 in part 1, 3 = #1 in part 1)
        std::vector<std::vector<MecScore>> distBp;
        for (uint32_t i = 0; i < bpTemp.size(); i++) {
            Bipartition b = bpTemp[i];
            distBp.push_back(std::vector<MecScore>(4, 0.0));
            for (uint32_t j = 0; j < b.size(); j++) {
                // determine case 0, 1, 2 or 3
                Allele a = getAllele(rs, active[j], positions[p]);
                if (a >= 0) {
                    uint32_t c = 2 * b[j] + a;
                    distBp[i][c] += 1.0;
                }
            }
        }
        std::vector<std::vector<MecScore>> distExt;
        for (uint32_t i = 0; i < extensions.size(); i++) {
            Bipartition b = extensions[i].first;
            distExt.push_back(std::vector<MecScore>(4, 0.0));
            for (uint32_t j = 0; j < b.size(); j++) {
                // determine case 0, 1, 2 or 3
                Allele a = getAllele(rs, ingoing[j], positions[p]);
                if (a >= 0) {
                    uint32_t c = 2 * b[j] + a;
                    distExt[i][c] += 1.0;
                }
            }
        }
        
        // combine all old bipartitions with all extensions
        active.insert(active.end(), ingoing.begin(), ingoing.end());
        std::vector<MecScore> candScores;
        // only compute scores first ...
        for (uint32_t i = 0; i < bpTemp.size(); i++) {
            for (uint32_t j = 0; j < extensions.size(); j++) {
                // construct score (bt is taken from btTemp)
                MecScore dist00 = distBp[i][0] + distExt[j][0];
                MecScore dist01 = distBp[i][1] + distExt[j][1];
                MecScore dist10 = distBp[i][2] + distExt[j][2];
                MecScore dist11 = distBp[i][3] + distExt[j][3];
                MecScore sc = scoreTemp[i];
                if (allHet)
                    sc += std::min(dist00 + dist11, dist01 + dist10);
                else
                    sc += std::min(dist00, dist01) + std::min(dist10, dist11);
                candScores.push_back(sc);
            }
        }
        std::sort(candScores.begin(), candScores.end());
        MecScore tooHigh = candScores.size() > rowLimit ? candScores[rowLimit] : std::numeric_limits<MecScore>::infinity();
        
        // .. and then construct solutions only if score is sufficient (performance reasons)
        lastCol.clear();
        for (uint32_t i = 0; i < bpTemp.size(); i++) {
            for (uint32_t j = 0; j < extensions.size(); j++) {
                // TODO: Avoid duplicate score calculation without using complex structs as candidates
                MecScore dist00 = distBp[i][0] + distExt[j][0];
                MecScore dist01 = distBp[i][1] + distExt[j][1];
                MecScore dist10 = distBp[i][2] + distExt[j][2];
                MecScore dist11 = distBp[i][3] + distExt[j][3];
                MecScore sc = scoreTemp[i];
                if (allHet)
                    sc += std::min(dist00 + dist11, dist01 + dist10);
                else
                    sc += std::min(dist00, dist01) + std::min(dist10, dist11);
                if (sc < tooHigh || sc == candScores[0]) {
                    lastCol.emplace_back(bpTemp[i], extensions[j].first, sc, btTemp[i]);
                    mBt[p].push_back(lastCol[lastCol.size() - 1].bt);
                    btVector.insert(btVector.end(), extensions[j].first.begin(), extensions[j].first.end());
                }
            }
        }
        mBp.emplace_back(ingoing.size(), btVector);
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
    std::vector<std::vector<uint32_t>> vote(n, std::vector<uint32_t>(4, 0));
    for (ReadId ri = 0; ri < m; ri++) {
        uint32_t part = full[ri];
        Read* r = rs->get(ri);
        for (int32_t i = 0; i < r->getVariantCount(); i++) {
            Allele a = r->getAllele(i);
            if (a >= 0)
                vote[posMap[r->getPosition(i)]][2 * part + a]++;
        }
    }
    
    // compute consensus
    std::vector<std::vector<Allele>> haps(2, std::vector<Allele>(n, 0));
    for (Position p = 0; p < n; p++) {
        if (allHet) {
            haps[0][p] = (vote[p][0] + vote[p][3]) < (vote[p][1] + vote[p][2]);
            haps[1][p] = 1 - haps[0][p];
        } else {
            haps[0][p] = vote[p][0] < vote[p][1];
            haps[1][p] = vote[p][2] < vote[p][3];
        }
    }
    
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


std::vector<std::pair<Bipartition, MecScore>> MecHeuristic::generateExtensions(std::vector<Allele>& alleleList, bool symmetric) {
    if (symmetric && alleleList.size() > 0) {
        std::vector<std::pair<Bipartition, MecScore>> e1 = generateExtensions(alleleList, false, (rowLimit + 1) / 2);
        std::vector<std::pair<Bipartition, MecScore>> e2 = generateExtensions(alleleList, true, rowLimit / 2);
        e1.insert(e1.end(), e2.begin(), e2.end());
        return e1;
    } else {
        return generateExtensions(alleleList, false, rowLimit);
    }
}


std::vector<std::pair<Bipartition, MecScore>> MecHeuristic::generateExtensions(std::vector<Allele>& alleleList, bool reverse, uint32_t limit) {
    std::vector<std::pair<Bipartition, MecScore>> results;
    Bipartition opt(alleleList.size());
    uint32_t n = alleleList.size();
    
    // create conflict-free bipartition as pivot
    for (uint32_t i = 0; i < n; i++) {
        assert(i < alleleList.size());
        assert(alleleList[i] < 2);//*
        assert(i < opt.size());
        if (alleleList[i] == 0) {
            opt[i] = reverse;
        } else {
            opt[i] = !reverse;
        }
    }
    results.emplace_back(opt, 0);
    
    // generate all bipartitions, where exactly i elements are placed differently than in pivot
    for (uint32_t i = 1; i <= n; i++) {
        std::vector<std::vector<uint32_t>> candidates = generateCombinations(n, i);
        // if all combinations with i differences fit into row limit, add them all ...
        if (candidates.size() + results.size() <= limit) {
            for (std::vector<uint32_t>& changed : candidates) {
                Bipartition b(opt);
                for (uint32_t j : changed) {
                    assert(j < b.size());
                    b[j] = !b[j];
                }
                results.emplace_back(b, i);
            }
        } else { // ... otherwise stop
            break;
        }
    }    
    return results;
}


std::vector<std::vector<uint32_t>> MecHeuristic::generateCombinations(const uint32_t n, const uint32_t k) {
    assert(n >= k);
    std::vector<std::vector<uint32_t>> results;
    std::vector<uint32_t> v(k + 1, 0);
    for (uint32_t i = 0; i < k; i++)
        v[i] = k - i - 1;
    // Iterate like a counter. Components are not set to zero after overflow but to <next value> - 1. 
    // Exit on overflow on last component
    assert(k < v.size());
    while (v[k] == 0) {
        // report current combination
        results.emplace_back(v.begin(), v.end() - 1);
        // increment first component
        uint32_t i = 0;
        v[0]++;
        // propagate overflows
        while (i < k && v[i] >= n - i) {
            i++;
            v[i]++;
        }
        // i := first position without overflow. All preceeding indices are reset
        for (uint32_t j = i - 1; j < j + 1; j--) {
            assert(j + 1 < v.size());
            v[j] = v[j + 1] + 1;
        }
        assert(k < v.size());
    }
    return results;
}


bool MecHeuristic::bpEqual(Bipartition a, Bipartition b) {
    uint32_t m = a.size();
    if (m != b.size())
        return false;
    for (uint32_t i = 0; i < m; i++)
        if (a[i] != b[i])
            return false;
    return true;
}


Allele MecHeuristic::getAllele(ReadSet* rs, ReadId rid, uint32_t genPos) {
    assert(rid < rs->size());
    Read* r = rs->get(rid);
    for (int32_t i = 0; i < r->getVariantCount(); i++)
        if (r->getPosition(i) == (int32_t)genPos)
            return r->getAllele(i);
    return -1;
}


void MecHeuristic::printColumnInfo(Position p, std::vector<ReadId>& startIndex, std::vector<BipartitionItem>& col) {
    MecScore s = col[0].score;
    for (uint32_t i = 1; i < col.size(); i++)
        if (col[i].score < s)
            s = col[i].score;
    assert(p + 1 < startIndex.size());
    std::cout<<"Column "<<p<<": ["<<startIndex[p]<<": "<<startIndex[p + 1] - 1<<"] with "<<col.size()<<" bipartitions and score "<<s<<std::endl;
}


