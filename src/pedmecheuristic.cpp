#include "pedmecheuristic.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>
#include <string> 
#include <set>

PedMecHeuristic::PedMecHeuristic(ReadSet* rs, const std::vector<unsigned int>& recombcost, const Pedigree* pedigree, bool distrust_genotypes, const std::vector<unsigned int>* positions, uint32_t rowLimit, bool allowMutations, uint32_t verbosity) :
    rs(rs),
    recombCost(recombcost.size(), 0.0),
    pedigree(pedigree),
    distrustGenotypes(distrust_genotypes),
    positions(0),
    rowLimit(std::min(rowLimit, MAX_ROW_LIMIT)),
    allowMutations(allowMutations),
    verbosity(verbosity),
    solved(false),
    tmBits(2 * pedigree->triple_count()),
    numSamples(0),
    globalSampleIds(0),
    sampleMap(0),
    trios(pedigree->get_triples()),
    genotypes(0),
    alleles(0),
    posMap(0),
    mutationCost(recombcost.size(), std::numeric_limits<MecScore>::infinity()),
    optScore(0),
    optBipart(0),
    optHaps(0),
    mutations(0),
    optTrans(0) {
        // costs
        size_t n = recombcost.size();
        for (uint32_t i = 1; i < n; i++) {
            recombCost[i] = (MecScore)recombcost[i];
            if (allowMutations)
                mutationCost[i - 1] = 0.75 * (this->recombCost[i - 1] + this->recombCost[i]);
        }
        if (allowMutations)
            mutationCost[n - 1] = recombCost[n - 1] * 1.5;
        
        // positions
        if (positions == nullptr)
            this->positions = rs->get_positions();
        else
            this->positions = new std::vector<Position>(positions->begin(), positions->end());
        n = this->positions->size();
        for (uint32_t i = 0; i < n; i++)
            posMap[(*(this->positions))[i]] = i;
        
        // sample ids
        std::set<uint32_t> sampleSet;
        for (uint32_t i = 0; i < rs->size(); i++)
            sampleSet.insert(rs->get(i)->getSampleID());
        for (auto& trio: pedigree->get_triples())
            for (size_t s: trio)
                sampleSet.insert((uint32_t)s);
            
        globalSampleIds.insert(globalSampleIds.end(), sampleSet.begin(), sampleSet.end());
        std::sort(globalSampleIds.begin(), globalSampleIds.end());
        numSamples = globalSampleIds.size();
        for (uint32_t i = 0; i < numSamples; i++)
            sampleMap[globalSampleIds[i]] = i;
        
        for (auto& trio: trios) {
            trio[0] = sampleMap[trio[0]];
            trio[1] = sampleMap[trio[1]];
            trio[2] = sampleMap[trio[2]];
        }
        
        // genotypes
        for (uint32_t s = 0; s < numSamples; s++) {
            genotypes.emplace_back(n, 0);
            alleles.emplace_back(n, std::pair<Allele, Allele>(0, 1));
            for (uint32_t p = 0; p < n; p++) {
                std::vector<uint32_t> a = pedigree->get_genotype(s, p)->as_vector();
                //TODO Support different alleles
                genotypes[s][p] = a[0] + a[1];
            }
        }
}

PedMecHeuristic::~PedMecHeuristic() {
    delete positions;
}

MecScore PedMecHeuristic::getOptScore() const {
    return optScore;
}

Bipartition* PedMecHeuristic::getOptBipartition() const {
    return new Bipartition(optBipart.begin(), optBipart.end());
}

std::vector<Transmission>* PedMecHeuristic::getOptTransmission() const {
    return new std::vector<Transmission>(optTrans.begin(), optTrans.end());
}

std::vector<std::vector<std::vector<Allele>>> PedMecHeuristic::getOptHaplotypes() const {
    return optHaps;
}

void PedMecHeuristic::getSuperReads(std::vector<ReadSet*>* superReads) const {
    for (uint32_t sid = 0; sid < numSamples; sid++) {
        Read* read0 = new Read("superread_0", -1, -1, globalSampleIds[sid]);
        Read* read1 = new Read("superread_1", -1, -1, globalSampleIds[sid]);
        for (uint32_t p = 0; p < positions->size(); p++)  {
            read0->addVariant((*positions)[p], optHaps[sid][0][p], 30);
            read1->addVariant((*positions)[p], optHaps[sid][1][p], 30);
        }
        ReadSet* phasedSample = new ReadSet();
        phasedSample->add(read0);
        phasedSample->add(read1);
        superReads->push_back(phasedSample);
    }
}

std::vector<std::vector<std::pair<uint32_t, uint32_t>>>* PedMecHeuristic::getMutations() const {
    return new std::vector<std::vector<std::pair<uint32_t, uint32_t>>>(mutations.begin(), mutations.end());
}

void PedMecHeuristic::solve() {
    if (solved)
        return;
    ReadId m = rs->size();
    Position n = positions->size();
    
    // compute index of first read starting at position p (for p = 0, 1, ..., n - 1, n)
    std::vector<ReadId> startIndex;
    startIndex.push_back(0);
    ReadId q = 0;
    for (Position p = 0; p < n; p++) {
        while (q < m && posMap[rs->get(q)->firstPosition()] <= p)
            q++;
        startIndex.push_back(q);
    }

    // for yet unseen non-child samples, we want to put the first read always into partition 0
    std::vector<bool> seen(numSamples, false);
    for (auto& trio: trios)
        seen[trio[2]] = true;
    
    // allocate data structures for backtracing and computations
    std::vector<PedSolution> lastCol;
    std::vector<ReadId> active;
    std::vector<std::vector<RowIndex>> mBt (n, std::vector<RowIndex>());
    std::vector<std::pair<uint32_t, Bipartition>> mBp;
    std::vector<Bipartition> mTm;
    
    // create empty biparition with any transmission as base for first column
    lastCol.emplace_back(Bipartition(), 0, 0, numSamples);
    // since recombinations are free in first column, we don't care yet about transmission vectors
    
    // DP recurrence
    Position right = 0;
    for (Position p = 0; p < n; p++) {
        // compute intersection of active reads with previous column
        std::vector<ReadId> activeLast(active.begin(), active.end());
        std::vector<uint32_t> kept;
        active.clear();
        for (uint32_t i = 0; i < activeLast.size(); i++) {
            ReadId r = activeLast[i];
            if (posMap[rs->get(r)->lastPosition()] >= p) {
                active.push_back(r);
                kept.push_back(i);
            }
        }
        activeLast.clear();

        // copy bipartitions, but without the lost reads
        std::vector<PedSolution> sols;
        std::unordered_map<uint64_t, std::vector<RowIndex>> buckets;
        // i := index of BP in last column
        // j := index of BP in current column (at least those that exist yet)
        for (RowIndex i = 0; i < lastCol.size(); i++) {
            Bipartition lastBp = lastCol[i].bp;
            Bipartition b(kept.size(), false);
            uint64_t bucket = lastCol[i].trans;
            for (size_t a = 0; a < kept.size(); a++) {
                bool bo = lastBp[kept[a]];
                b[a] = bo;
                bucket = (bucket << 1) + bo;
            }
            bool duplicate = false;
            
            // use buckets (last 64 bits of bipartition bitvector) to avoid full pairwise check
            if (buckets.find(bucket) == buckets.end())
                buckets[bucket] = std::vector<RowIndex>();
            for (RowIndex j : buckets[bucket]) {
                if (lastCol[i].trans == sols[j].trans && bpEqual(b, sols[j].bp)) {
                    // found another previous solution with same remaining bipartition and transmission
                    duplicate = true;
                    updateSolution(sols[j], lastCol[i], i);
                    break;
                }
            }
            if (!duplicate) {
                buckets[bucket].push_back(sols.size());
                sols.emplace_back(b, lastCol[i].trans, std::numeric_limits<MecScore>::infinity(), numSamples);
                updateSolution(sols[sols.size() - 1], lastCol[i], i);
            }
        }
        buckets.clear();
        // from here on, we have solsTemp available as basis for extensions

        // get last position of new reads and extend balance vectors
        right = std::max(right, p);
        for (ReadId r = startIndex[p]; r < startIndex[p + 1];  r++)
            right = std::max(right, posMap[rs->get(r)->lastPosition()]);
        for (PedSolution& a: sols)
            for (auto& b: a.balances)
                b.resize(right + 1 - p, 0);
        
        // summarize equal reads
        size_t numNew = startIndex[p + 1] - startIndex[p];
        std::vector<int> equalTo(numNew, -1);
        std::vector<Balance> balances;
        std::vector<uint32_t> sampleIds;
        for (size_t i = 0; i < numNew;  i++) {
            ReadId r = startIndex[p] + i;
            active.push_back(r);
            // generate balance vector of read
            Balance b(right + 1 - p, 0);
            sampleIds.push_back(sampleMap[rs->get(r)->getSampleID()]);
            for (int32_t i = 0; i < rs->get(r)->getVariantCount(); i++) {
                Position o = posMap[rs->get(r)->getPosition(i)] - p;
                Allele a = rs->get(r)->getAllele(i);
                MecScore q = rs->get(r)->getVariantQuality(i);
                b[o] += q * a - q * (1 - a);
            }
            // check all previous balances and merge if alleles identical
            for (size_t j = 0; j < i; j++) {
                if (equalTo[j] != -1 || sampleIds[j] != sampleIds[i])
                    continue;
                bool equal = true;
                for (size_t k = 0; k < right + 1 - p; k++)
                    if (balances[j][k] * b[k] < 0 || (balances[j][k] != 0.0) != (b[k] != 0.0)) {
                        equal = false; break;
                    }
                if (equal) {
                    equalTo[i] = j;
                    for (size_t k = 0; k < right + 1 - p; k++)
                        balances[j][k] += b[k];
                    break;
                }
            }
            balances.push_back(b);
        }
        
        for (size_t i = 0; i < numNew;  i++) {
            if (equalTo[i] >= 0)
                continue;
            if (verbosity >= 2) {
                std::cout<<"  Read "<<i<<" ("<<sampleIds[i]<<"): ";
                for (size_t k = 0; k < right + 1 - p; k++)
                    std::cout<<balances[i][k]<<" ";
                std::cout<<std::endl;
            }
        }

        // for every new read: generate two extensions for every solution and keep only rowLimit best
        for (size_t i = 0; i < numNew;  i++) {
            Balance& balance = balances[i];
            uint32_t sampleId = sampleIds[i];
            std::vector<int> target(genotypes[sampleId].begin() + p, genotypes[sampleId].begin() + right + 1);

            uint32_t solEnd = sols.size();
            for (uint32_t sol = 0; sol < solEnd; sol++) {
                /* if read is identical to previous one, just put it in same partition and continue*/
                if (equalTo[i] >= 0) {
                    sols[sol].bpNew.push_back(sols[sol].bpNew[equalTo[i]]);
                    continue;
                }
                
                /* if read only contributes to positions with equal consensus and does not change it:
                 * do not branch, just put read where it fits best */
                bool useful = false;
                if (distrustGenotypes)
                    for (uint32_t j = 0; j < balance.size() && !useful; j++) {
                        MecScore s0 = sols[sol].balances[sampleId * 2][j], s1 = sols[sol].balances[sampleId * 2 + 1][j];
                        useful |= (balance[j] != 0 && s0 * s1 < 0) || ((balance[j] + s0) * s0 <= 0 && (balance[j] + s1) * s1 <= 0);
                    }
                else
                    for (uint32_t j = 0; j < balance.size() && !useful; j++)
                        useful |= (genotypes[sampleId][p + j] == 1 && balance[j] != 0);
                
                /* if sample seen before and read not identical to some previous one, 
                 * create new solution with read in partiton 1*/
                uint32_t sol1 = 0;
                if (seen[sampleId]) {
                    sols.emplace_back(sols[sol]);
                    sol1 = sols.size() - 1;
                    sols[sol1].score += addBalance(sols[sol1].balances[2 * sampleId + 1], sols[sol1].balances[2 * sampleId], balance, target);
                    sols[sol1].mutationScore = getMutationCost(sols[sol1].balances, sols[sol1].trans, p, true, 5);
                    sols[sol1].bpNew.push_back(true);
                }
                // update existing solution with read in partition 0
                sols[sol].score += addBalance(sols[sol].balances[2 * sampleId], sols[sol].balances[2 * sampleId + 1], balance, target);
                sols[sol].mutationScore = getMutationCost(sols[sol].balances, sols[sol].trans, p, true, 5);
                sols[sol].bpNew.push_back(false);

                if (sol1 && !useful) {
                    // store both solutions if read was useful, otherwise only keep the better
                    if (sols[sol].score + sols[sol].mutationScore > sols[sol1].score + sols[sol1].mutationScore)
                        sols[sol] = sols[sol1];
                    sols.pop_back();
                }
            }
            seen[sampleId] = true;
            
            // prune solution set if too large
            if (sols.size() > rowLimit)
                filterSolutions(sols);
        }
        // extend transmissions
        size_t solEnd = sols.size();
        for (uint32_t i = 0; i < solEnd; i++)
            extendSolutions(sols, i, p);
        if (sols.size() > rowLimit)
            filterSolutions(sols);
        
        // add costs for mutations
        for (uint32_t i = 0; i < sols.size(); i++) {
            // sols[i].score += getMutationCost(sols[i].balances, sols[i].trans, p, false, 0);
            std::vector<MecScore> firsts;
            for (Balance& b: sols[i].balances)
                firsts.push_back(b[0]);
            sols[i].score += getOptPhasing(firsts, sols[i].trans, p, nullptr);
        }
        
        // .. and then construct solutions only if score is sufficient (performance reasons)
        lastCol.clear();
        Bipartition btVector;
        Bipartition tmVector;
        for (uint32_t i = 0; i < sols.size(); i++) {
            sols[i].finalize();
            lastCol.emplace_back(sols[i]);
            mBt[p].push_back(sols[i].btRow);
            btVector.insert(btVector.end(), sols[i].bpNew.begin(), sols[i].bpNew.end());
            for (uint32_t j = 0; j < tmBits; j++)
                tmVector.push_back((sols[i].trans >> j) & 1);
        }
        mBp.emplace_back(startIndex[p + 1] - startIndex[p], btVector);
        mTm.push_back(tmVector);
        if (verbosity >= 1)
            printColumnInfo(p, startIndex, lastCol);
    }
    
    // find best score in last column
    optTrans.resize(n);
    MecScore s = std::numeric_limits<MecScore>::infinity();
    RowIndex ri = 0;
    for (uint32_t sol = 0; sol < lastCol.size(); sol++) {
        if (lastCol[sol].score < s) {
            s = lastCol[sol].score;
            ri = sol;
        }
    }

    // backtracking
    optBipart.resize(m);
    for (Position p = n - 1; p < p + 1; p--) {
        ReadId offset = startIndex[p];
        ReadId newCount = mBp[p].first;
        Bipartition current(mBp[p].second.begin() + newCount * ri, mBp[p].second.begin() + newCount * (ri + 1));
        for (uint32_t i = 0; i < current.size(); i++)
            optBipart[offset + i] = current[i];
        Bipartition tm(mTm[p].begin() + tmBits * ri, mTm[p].begin() + tmBits * (ri + 1));
        for (int t = tmBits - 1; t >= 0; t--)
            optTrans[p] = (optTrans[p] << 1) + tm[t];
        ri = mBt[p][ri];
    }

    // get allele votes
    std::vector<std::vector<MecScore>> balances(n, std::vector<MecScore>(2 * numSamples, 0));
    for (ReadId ri = 0; ri < m; ri++) {
        Read* r = rs->get(ri);
        for (int32_t i = 0; i < r->getVariantCount(); i++) {
            Allele a = r->getAllele(i);
            MecScore q = r->getVariantQuality(i);
            uint32_t sid = sampleMap[r->getSampleID()];
            if (a >= 0)
                balances[posMap[r->getPosition(i)]][2 * sid + optBipart[ri]] += (2 * a - 1) * q;
        }
    }
    
    // compute optimal phasing per position
    for (uint32_t sid = 0; sid < numSamples; sid++) {
        optHaps.emplace_back(2, std::vector<Allele>(n, -1));
    }
    mutations.resize(numSamples);
    for (uint32_t p = 0; p < n; p++) {
        std::vector<Allele>* posPhasing = new std::vector<Allele>(2 * numSamples, 0);
        std::vector<bool>* mut = new std::vector<bool>(2 * numSamples, false);
        MecScore s = getOptPhasing(balances[p], optTrans[p], p, posPhasing, mut);
        for (uint32_t sid = 0; sid < numSamples; sid++) {
            optHaps[sid][0][p] = (*posPhasing)[2 * sid];
            optHaps[sid][1][p] = (*posPhasing)[2 * sid + 1];
            if ((*mut)[2 * sid])
                mutations[sid].emplace_back(0, p);
            if ((*mut)[2 * sid + 1])
                mutations[sid].emplace_back(1, p);
        }
        delete posPhasing;
        delete mut;
    }
    solved = true;
}

bool PedMecHeuristic::bpEqual(const Bipartition a, const Bipartition b) const {
    uint32_t m = a.size();
    if (m != b.size())
        return false;
    for (uint32_t i = 0; i < m; i++)
        if (a[i] != b[i])
            return false;
    return true;
}

void PedMecHeuristic::updateSolution(PedSolution& newSol, const PedSolution& oldSol, RowIndex oldIdx) const {
    if (newSol.score > oldSol.score) {
        newSol.score = oldSol.score;
        newSol.btRow = oldIdx;
        newSol.balances.clear();
        for (auto& a: oldSol.balances)
            if (a.size() > 0)
                newSol.balances.emplace_back(a.begin() + 1, a.end());
            else
                newSol.balances.emplace_back();
    }
}

MecScore PedMecHeuristic::getRecombinationCost(const Transmission t1, const Transmission t2, Position p) const {
    return recombCost[p] * (MecScore)popcount((uint64_t)(t1 ^ t2));
}

MecScore PedMecHeuristic::getMutationCost(const std::vector<Balance>& balances,
                                          const Transmission& t,
                                          Position p,
                                          bool allowFlips,
                                          size_t ahead) const {
    MecScore cost = 0.0;
    size_t last = std::min(ahead, balances[0].size() - 1);
    for (size_t i = 0; i <= last; i++) {
        for (uint32_t k = 0; k < trios.size(); k++) {
            auto& trio = trios[k];
            uint32_t m2c = (t >> (2 * k)) & 1;
            uint32_t f2c = (t >> (2 * k + 1)) & 1;
            MecScore cm = balances[2 * trio[2]][i];
            MecScore cf = balances[2 * trio[2] + 1][i];
            MecScore m = balances[2 * trio[0] + m2c][i];
            MecScore f = balances[2 * trio[1] + f2c][i];
            if (allowFlips) {
                if (cm * m < 0)
                    cost += std::min(mutationCost[p], std::min(std::abs(cm), std::abs(m)));
                if (cf * f < 0)
                    cost += std::min(mutationCost[p], std::min(std::abs(cf), std::abs(f)));
            } else {
                cost += (cm * m < 0) * mutationCost[p];
                cost += (cf * f < 0) * mutationCost[p];
            }
        }
    }
    return cost;
}


MecScore PedMecHeuristic::getOptPhasing(const std::vector<MecScore>& balances,
                                                 const Transmission& t,
                                                 Position p,
                                                 std::vector<Allele>* optPhasing,
                                                 std::vector<bool>* mutated) const {
                                                     
    // precompute corresponding cost based on balance vector
    std::vector<std::vector<MecScore>> phaseCost(numSamples, std::vector<MecScore>(5));
    for (size_t s = 0; s < numSamples; s++) {
        MecScore a0 = balances[2 * s];
        MecScore a1 = balances[2 * s + 1];
        phaseCost[s][0] = (a0 * (a0 > 0) + a1 * (a1 > 0));
        phaseCost[s][1] = (-a0 * (a0 < 0) + a1 * (a1 > 0));
        phaseCost[s][2] = (a0 * (a0 > 0) - a1 * (a1 < 0));
        phaseCost[s][3] = (-a0 * (a0 < 0) - a1 * (a1 < 0));
        phaseCost[s][4] = *std::max_element(phaseCost[s].begin(), phaseCost[s].begin() + 4);
    }
    
    /* enumerate all allowed phasings per sample
     * 0 = 0|0, 1 = 0|1, 2 = 1|0, 3 = 1|1 */
    std::vector<std::vector<int>> phases(numSamples);
    if (distrustGenotypes) {
        for (size_t s = 0; s < numSamples; s++)
            for (size_t i = 0; i < 4; i++)
                if (phaseCost[s][i] < phaseCost[s][4] + 2 * mutationCost[p])
                    phases[s].push_back(i);
    } else {
        for (size_t s = 0; s < numSamples; s++) {
            if (genotypes[s][p] == 0)
                phases[s].push_back(0);
            else if (genotypes[s][p] == 2)
                phases[s].push_back(3);
            else {
                phases[s].push_back(1);
                phases[s].push_back(2);
            }
        }
    }
    
    // iterate over all allowed genotype combinations
    MecScore minCost = std::numeric_limits<MecScore>::infinity();
    std::vector<size_t> v(numSamples, 0);        
    while (v[numSamples - 1] < phases[numSamples - 1].size()) {
        MecScore cost = 0.0;
        std::vector<bool> mut(numSamples * 2, false);
        // iterate over trios to detect mutations
        for (uint32_t k = 0; k < trios.size(); k++) {
            auto& trio = trios[k];
            uint32_t m2c = (t >> (2 * k)) & 1;
            uint32_t f2c = (t >> (2 * k + 1)) & 1;
            // extract alleles according to current genotype combination
            Allele acm = phases[trio[2]][v[trio[2]]] & 1;
            Allele acf = (phases[trio[2]][v[trio[2]]] & 2) >> 1;
            Allele am = (phases[trio[0]][v[trio[0]]] & (1 + m2c)) >> m2c;
            Allele af = (phases[trio[1]][v[trio[1]]] & (1 + f2c)) >> f2c;
            // add costs for every violated inheritance rule
            cost += (am != acm) * mutationCost[p];
            cost += (af != acf) * mutationCost[p];
            mut[2 * trio[2]] = (am != acm);
            mut[2 * trio[2] + 1] = (af != acf);
        }
        // iterate over samples to penalize deviating genotype and allele balances
        for (uint32_t s = 0; s < numSamples; s++)
            cost += phaseCost[s][phases[s][v[s]]];
        
        if (cost < minCost) {
            minCost = cost;
            if (optPhasing != nullptr)
                for (uint32_t s = 0; s < numSamples; s++) {
                    (*optPhasing)[2 * s] = phases[s][v[s]] & 1;
                    (*optPhasing)[2 * s + 1] = (phases[s][v[s]] & 2) >> 1;
                }
            if (mutated != nullptr)
                for (uint32_t s = 0; s < 2 * numSamples; s++)
                    (*mutated)[s] = mut[s];
        }
        
        // increment counters
        v[0]++;
        for (size_t j = 0; j < numSamples - 1; j++)
            if (v[j] >= phases[j].size()) {
                v[j] = 0;
                v[j + 1]++;
            }
        
    }
    // report only cost for best genotype combination
    return minCost;
}

MecScore PedMecHeuristic::addBalance(Balance& basis,
                                     const Balance& coBasis,
                                     const Balance& add,
                                     const std::vector<int>& target) const {
    MecScore penalty = 0;
    for (uint32_t i = 0; i < add.size(); i++) {
        if (distrustGenotypes) {
            if (basis[i] * add[i] < 0)
                // penalty = how much did the added balance move the old one towards zero
                penalty += std::min(std::abs(basis[i]), std::abs(add[i]));
        } else if (target[i] == 1) {
            if (add[i] <= 0)
                penalty += std::min(-add[i], std::max(basis[i] - coBasis[i], (MecScore)0));
            else
                penalty += std::min(add[i], std::max(coBasis[i] - basis[i], (MecScore)0));
        } else {
            penalty += std::abs(add[i]) * (add[i] * (target[i] - 1) < 0);
        }
        // update balance vector
        basis[i] += add[i];
    }
    return penalty;
}

void PedMecHeuristic::extendSolutions(std::vector<PedSolution>& sols,
                                      uint32_t toExt,
                                      Position p) const {
    sols[toExt].mutationScore = getMutationCost(sols[toExt].balances, sols[toExt].trans, p, false, 0);
    if (sols[toExt].mutationScore > 0) {
        // if mutation score could be improved, try all other transmissions as well
        for (Transmission t = 0; t < std::pow(2, tmBits); t++) {
            if (t == sols[toExt].trans)
                continue;
            MecScore rc = getRecombinationCost(sols[toExt].trans, t, p);
            if (rc >= sols[toExt].mutationScore)
                continue;
            MecScore mut = getMutationCost(sols[toExt].balances, t, p, false, 0);
            if (mut + rc >= sols[toExt].mutationScore)
                continue;
            // if t good enough, add new solution
            sols.emplace_back(sols[toExt].bp, sols[toExt].bpNew, t, sols[toExt].score + rc, mut, sols[toExt].btRow, sols[toExt].balances);
        }
    }
}

void PedMecHeuristic::filterSolutions(std::vector<PedSolution>& sols) const {
    // collect and sort scores
    std::vector<MecScore> scores;
    for (uint32_t sol = 0; sol < sols.size(); sol++)
        scores.push_back(sols[sol].score + sols[sol].mutationScore);
    std::sort(scores.begin(), scores.end());
    MecScore tooHigh = scores.size() > rowLimit ? scores[rowLimit] : std::numeric_limits<MecScore>::infinity();
    
    // only keep what is better than rowLimit+1-th solution or optimal
    std::vector<uint32_t> kept;
    for (uint32_t sol = 0; sol < sols.size(); sol++) {
        MecScore score = sols[sol].score + sols[sol].mutationScore;
        if ((score < tooHigh || score == scores[0]) && kept.size() < MAX_ROW_LIMIT)
            kept.push_back(sol);
    }
    for (uint32_t sol = 0; sol < kept.size(); sol++)
        sols[sol] = sols[kept[sol]];
    sols.resize(kept.size());
}

void PedMecHeuristic::printColumnInfo(Position p, std::vector<ReadId>& startIndex, std::vector<PedSolution>& col) const {
    MecScore s = col[0].score;
    for (uint32_t i = 0; i < col.size(); i++)
        if (col[i].score < s)
            s = col[i].score;
    std::cout<<"Column "<<p<<": ["<<startIndex[p]<<": "<<startIndex[p + 1] - 1<<"] with "<<col.size()<<" bipartitions and score "<<s<<" (costs "<<recombCost[p]<<", "<<mutationCost[p]<<")"<<std::endl;
}
