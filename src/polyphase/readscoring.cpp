#include "readscoring.h"
#include "../binomial.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <assert.h>

void ReadScoring::scoreReadset(TriangleSparseMatrix* result, ReadSet* readset, const uint32_t minOverlap, const uint32_t ploidy, double err) const {
    
    if (ploidy < 2) {
        std::cout<<"Error: Ploidy < 2!"<<std::endl;
        return;
    }
            
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint8_t>> alleles;
    std::vector<uint32_t> posList;
    std::unordered_map<uint32_t, uint32_t> posMap;
    uint32_t longestReadSpan = 0;
    computeStartEnd(readset, begins, ends, positions, alleles, posList, posMap, longestReadSpan);
    
    std::unordered_set<Genotype> occGenotypesSet;
    uint8_t maxAllele = 0;
    
    // compute allele depths for every position
    std::vector<std::unordered_map<uint8_t, uint32_t>> alleleDepths(posMap.size(), std::unordered_map<uint8_t, uint32_t>());
    for (uint32_t read = 0; read < begins.size(); read++) {
        for (uint32_t pos = 0; pos < positions[read].size(); pos++) {
            uint32_t varPos = posMap[positions[read][pos]];
            uint8_t allele = alleles[read][pos];
            if (alleleDepths[varPos].find(allele) == alleleDepths[varPos].end())
                alleleDepths[varPos][allele] = 0;
            alleleDepths[varPos][alleles[read][pos]] += 1;
            maxAllele = std::max(maxAllele, (uint8_t)(allele+1));
        }
    }
    
    if (err == 0.0)
        err = estimateAlleleErrorRate(alleleDepths, ploidy);
    
    // compute genotype likelihoods
    std::vector<std::unordered_map<Genotype, double>> gl(posMap.size());
    for (uint32_t i = 0; i < posMap.size(); i++) {
        std::unordered_map<Genotype, double> l = computeGenotypeLikelihoods(alleleDepths[i], ploidy, err, true);
        gl[i] = l;
        for (auto& g: l) {
            occGenotypesSet.insert(g.first);
        }
    }
    
    // precompute allele-pair likelihood ratios for existing genotypes
    // aplp[g][a1][a2] = P(a1, a2 | g and alleles from same haplotype)
    // apld[g][a1][a2] = P(a1, a2 | g and alleles from different haplotypes)
    std::vector<Genotype> occGenotypes(occGenotypesSet.begin(), occGenotypesSet.end());
    std::unordered_map<Genotype, uint32_t> gMap;
    for (uint32_t i = 0; i < occGenotypes.size(); i++) {
        gMap[occGenotypes[i]] = i;
    }
    std::vector<double> apls(maxAllele*maxAllele*occGenotypes.size());
    std::vector<double> apld(maxAllele*maxAllele*occGenotypes.size());
    computeAllelePairLikelihoods(occGenotypes, apls, apld, maxAllele, err);
        
    // for each read pair, iterate over their overlapping positions and determine combined score
    for (uint32_t i = 0; i < begins.size(); i++) {
        // iterate until start position of read is behind required start
        for (uint32_t j = i+1; j < begins.size() && begins[j] <= ends[i]; j++) {
            float score = computeLogScore(positions[i], positions[j], alleles[i], alleles[j], maxAllele, posMap, gl, gMap, apls, apld, minOverlap);
            if (score != 0.0) {
                // multiply probability with (1/k) / ((1-1/k)) as prior probability for same/diff haplotype
                result->set(i, j, score - std::log(ploidy*(1.0-1.0/ploidy)));
            }
        }
    }
}

double ReadScoring::estimateAlleleErrorRate(std::vector<std::unordered_map<uint8_t, uint32_t>>& alleleDepths, uint32_t ploidy) const {  
    
    std::vector<std::unordered_map<Genotype, double>> gl(alleleDepths.size());
    double bestErr = 0.0;
    double bestSum = -std::numeric_limits<double>::infinity();
    
    for (double err = 0.01; err < 0.2; err += 0.01) {
        for (uint32_t i = 0; i < alleleDepths.size(); i++) {
            gl[i].clear();
            gl[i] = computeGenotypeLikelihoods(alleleDepths[i], ploidy, err, true);
        }
        double sum = evaluateGenotypeLikelihoods(gl);
        std::cout<<"Err="<<err<<" -> Sum="<<sum<<std::endl;
        if (sum > bestSum) {
            bestSum = sum;
            bestErr = err;
        }
    }
    std::cout<<"BestErr="<<bestErr<<std::endl;
    return bestErr;
}

double ReadScoring::evaluateGenotypeLikelihoods(std::vector<std::unordered_map<Genotype, double>>& gl) const {
    double logLikelihoodSum = 0.0;
    for (uint32_t i = 0; i < gl.size(); i++) {
        double max = 0.0;
        for (auto& l: gl[i]) {
            if (l.second > max)
                max = l.second;
        }
        logLikelihoodSum += std::log(max);
    }
    return logLikelihoodSum;
}

std::unordered_map<Genotype, double> ReadScoring::computeGenotypeLikelihoods (std::unordered_map<uint8_t, uint32_t> alleleDepth,
                                                                              const uint32_t ploidy,
                                                                              const double err,
                                                                              const bool normalize) const {
    std::unordered_map<Genotype, double> gl;
    uint32_t numAlleles = alleleDepth.size();
    assert(numAlleles <= 2);
    uint32_t numGenotypes = binomial_coefficient(ploidy+numAlleles-1, numAlleles-1);
    double weight = 0.0;
    
    // generate all possible genotypes
    for (uint32_t index = 0; index < numGenotypes; index++) {
        Genotype g(index, ploidy);
        // TODO: Make multi-allelic
        if (numAlleles == 1) {
            weight += 1;
            gl[g] = 1;
        } else {
            double fracAlt = (double)index / (double)ploidy;
            double l = binom_pmf(alleleDepth[0]+alleleDepth[1], alleleDepth[1], (1-fracAlt)*err + fracAlt*(1-err));
            weight += l;
            gl[g] = l;
        }
    }
    
    // normalize
    if (normalize) {
        for (uint32_t index = 0; index < numGenotypes; index++) {
            Genotype g(index, ploidy);
            gl[g] /= weight;
        }
    }
    return gl;
}

void ReadScoring::computeAllelePairLikelihoods(std::vector<Genotype>& genos,
                                               std::vector<double>& apls,
                                               std::vector<double>& apld,
                                               uint8_t maxAllele,
                                               double err) const {
    uint32_t numGenos = genos.size();
    for (uint32_t a1 = 0; a1 < maxAllele; a1++) {
        for (uint32_t a2 = a1; a2 < maxAllele; a2++) {
            for (uint32_t gi = 0; gi < numGenos; gi++) {
                uint32_t i1 = numGenos * (maxAllele*a1 + a2) + gi;
                uint32_t i2 = numGenos * (maxAllele*a2 + a1) + gi;
                
                double lEqual = 0.0;
                double lDiff = 0.0;
                std::vector<uint32_t> gv = genos[gi].as_vector();
                for (uint32_t gv1 = 0; gv1 < gv.size(); gv1++) {
                    for (uint32_t gv2 = 0; gv2 < gv.size(); gv2++) {
                        // for matching alleles, multiply with (1-err), for non-matching with err
                        double l = (1-err)*(gv[gv1] == a1) + err*(gv[gv1] != a1);
                        l *= (1-err)*(gv[gv2] == a2) + err*(gv[gv2] != a2);
                        // add to according sum and normalize to number of cases (prior probabilities are handles elsewhere)
                        if (gv1 == gv2)
                            lEqual += l / gv.size();
                        else
                            lDiff += l / (gv.size()*(gv.size()-1));
                    }
                }
                apls[i1] = apls[i2] = lEqual;
                apld[i1] = apld[i2] = lDiff;
            }
        }
    }
}

float ReadScoring::computeLogScore (std::vector<uint32_t>& posRead1,
                                    std::vector<uint32_t>& posRead2,
                                    std::vector<uint8_t>& alleleRead1,
                                    std::vector<uint8_t>& alleleRead2,
                                    const uint32_t maxAllele,
                                    std::unordered_map<uint32_t, uint32_t>& posMap,
                                    std::vector<std::unordered_map<Genotype, double>>& gl,
                                    std::unordered_map<Genotype, uint32_t>& gMap,
                                    std::vector<double>& apls,
                                    std::vector<double>& apld,
                                    const uint32_t minOverlap) const {
    ;
    uint32_t ov = 0;
    double logScore = 0.0;
    uint32_t k = 0;
    uint32_t l = 0;
    while (k < posRead1.size() && l < posRead2.size()) {
        if (posRead1[k] == posRead2[l]) {
            logScore += computeLogScoreSinglePos(alleleRead1[k], alleleRead2[l], maxAllele, gl[posMap[posRead1[k]]], gMap, apls, apld);
            ov += 1;
            k++; l++;
        } else if (posRead1[k] < posRead2[l]) {
            k++;
        } else {
            l++;
        }
    }
    if (ov >= minOverlap)
        return logScore;
    else
        return 0.0;
}

float ReadScoring::computeLogScoreSinglePos (uint8_t allele1,
                                             uint8_t allele2,
                                             const uint32_t maxAllele,
                                             std::unordered_map<Genotype, double>& gl,
                                             std::unordered_map<Genotype, uint32_t>& gMap,
                                             std::vector<double>& apls,
                                             std::vector<double>& apld) const {
    double same = 0.0;
    double diff = 0.0;
    for (std::pair<Genotype, double> p : gl) {
        uint32_t i = gMap.size()*(allele1*maxAllele + allele2) + gMap[p.first];
        same += p.second * apls[i];
        diff += p.second * apld[i];
    }
    if (same == 0.0 | diff == 0.0) {
        std::cout<<"Zero-score on "<<gl.size()<<" genotypes "<<"same="<<same<<" diff="<<diff<<std::endl;
        return 0.0;
    } else {
        return (float)std::log(same/diff);
    }
}

void ReadScoring::computeStartEnd (const ReadSet* readset,
                                   std::vector<uint32_t>& begins,
                                   std::vector<uint32_t>& ends,
                                   std::vector<std::vector<uint32_t>>& positions,
                                   std::vector<std::vector<uint8_t>>& alleles,
                                   std::vector<uint32_t>& posList,
                                   std::unordered_map<uint32_t, uint32_t>& posMap,
                                   uint32_t& longestReadSpan) const {
    // copy all relevant information from the readset into vectors for efficient access
    uint32_t numReads = readset->size();
    std::unordered_set<uint32_t> allPos;
    for (uint32_t i = 0; i < numReads; i++) {
        begins.push_back(readset->get(i)->firstPosition());
        ends.push_back(readset->get(i)->lastPosition());
        std::vector<uint32_t> pos;
        std::vector<uint8_t> all;
        for (int k = 0; k < readset->get(i)->getVariantCount(); k++) {
            pos.push_back(readset->get(i)->getPosition(k));
            all.push_back((uint8_t)(readset->get(i)->getAllele(k)));
            allPos.insert(readset->get(i)->getPosition(k));
        }
        positions.push_back(pos);
        alleles.push_back(all);
    }
    
    // create position map
    posList.clear();
    for (uint32_t pos : allPos) {
        posList.push_back(pos);
    }
    std::sort(posList.begin(), posList.end());
    for (uint32_t i = 0; i < posList.size(); i++) {
        posMap[posList[i]] = i;
    }
    
    // determine length of the longest read
    longestReadSpan = 0;
    for (uint32_t i = 0; i < numReads; i++) {
        longestReadSpan = std::max(longestReadSpan, ends[i] - begins[i]);
    }
}
