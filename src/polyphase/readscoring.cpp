#include "readscoring.h"
#include "../binomial.h"
#include "../multinomial.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <assert.h>

using Allele = AlleleMatrix::Allele;
using Position = AlleleMatrix::Position;
using AlleleRow = AlleleMatrix::AlleleRow;
using AlleleItem = AlleleMatrix::AlleleItem;

void ReadScoring::scoreReadset(TriangleSparseMatrix* result, AlleleMatrix* am, const uint32_t minOverlap, const uint32_t ploidy, double err) const {
    
    if (ploidy < 2) {
        std::cout<<"Error: Ploidy < 2!"<<std::endl;
        return;
    }

    uint32_t longestReadSpan = 0;
    for (Position i = 0; i < am->size(); i++)
        longestReadSpan = std::max(longestReadSpan, am->getLastPos(i) - am->getFirstPos(i));
    
    std::unordered_set<Genotype> occGenotypesSet;
    
    if (err == 0.0)
        err = estimateAlleleErrorRate(am, ploidy);
    
    // compute genotype likelihoods
    std::vector<std::unordered_map<Genotype, double>> gl(am->getNumPositions());
    for (uint32_t i = 0; i < gl.size(); i++) {
        std::unordered_map<Genotype, double> l = computeGenotypeLikelihoods(am->getAlleleDepths(i), ploidy, err, true);
        gl[i] = l;
        for (auto& g: l) {
            occGenotypesSet.insert(g.first);
        }
    }
    
    // precompute allele-pair likelihood ratios for existing genotypes
    // apls[g][a1][a2] = P(a1, a2 | g and alleles from same haplotype)
    // apld[g][a1][a2] = P(a1, a2 | g and alleles from different haplotypes)
    std::vector<Genotype> occGenotypes(occGenotypesSet.begin(), occGenotypesSet.end());
    std::unordered_map<Genotype, uint32_t> gMap;
    for (uint32_t i = 0; i < occGenotypes.size(); i++) {
        gMap[occGenotypes[i]] = i;
    }
    Allele numAlleles = am->getMaxNumAllele();
    std::vector<double> apls(numAlleles * numAlleles * occGenotypes.size());
    std::vector<double> apld(numAlleles * numAlleles * occGenotypes.size());
    computeAllelePairLikelihoods(occGenotypes, apls, apld, numAlleles, err);
        
    // for each read pair, iterate over their overlapping positions and determine combined score
    std::vector<uint32_t> sortedReads(am->size());
    for (uint32_t i = 0; i < sortedReads.size(); i++)
        sortedReads[i] = i;
    std::sort(sortedReads.begin(), sortedReads.end(), [am] (const uint32_t a, const uint32_t b) { return am->getFirstPos(a) < am->getFirstPos(b); });
    float offset = -std::log(ploidy * (1.0 - 1.0 / ploidy));
    for (uint32_t i = 0; i < am->size(); i++) {
        // iterate until start position of read is behind required start
        uint32_t terminal = am->getLastPos(sortedReads[i]) - minOverlap + 1;
        for (uint32_t j = i + 1; j < sortedReads.size() && am->getFirstPos(sortedReads[j]) <= terminal; j++) {
            float score = computeLogScore(am, sortedReads[i], sortedReads[j], gl, gMap, apls, apld, minOverlap);
            if (score != 0.0) {
                // multiply probability with (1/k) / ((1-1/k)) as prior probability for same/diff haplotype
                result->set(sortedReads[i], sortedReads[j], score + offset);
            }
        }
    }
}

double ReadScoring::estimateAlleleErrorRate(AlleleMatrix *am, uint32_t ploidy) const {
    
    std::vector<std::unordered_map<Genotype, double>> gl(am->getNumPositions());
    double bestErr = 0.0;
    double bestSum = -std::numeric_limits<double>::infinity();
    
    for (double err = 0.01; err < 0.2; err += 0.01) {
        for (Position i = 0; i < am->getNumPositions(); i++) {
            gl[i].clear();
            gl[i] = computeGenotypeLikelihoods(am->getAlleleDepths(i), ploidy, err, true);
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

std::unordered_map<Genotype, double> ReadScoring::computeGenotypeLikelihoods (std::vector<uint32_t> alleleDepth,
                                                                              const uint32_t ploidy,
                                                                              const double err,
                                                                              const bool normalize) const {
    std::unordered_map<Genotype, double> gl;
    uint32_t numAlleles = alleleDepth.size();
    uint32_t numGenotypes = binomial_coefficient(ploidy + numAlleles - 1, numAlleles - 1);
    double weight = 0.0;
    std::vector<Allele> alleles;
    uint32_t numExAlleles = 0;
    for (uint32_t i = 0; i < numAlleles; i++)
        if (alleleDepth[i] > 0) {
            alleles.push_back(i);
            numExAlleles++;
        }
    
    // generate all possible genotypes
    for (uint32_t index = 0; index < numGenotypes; index++) {
        Genotype g(index, ploidy);
        bool reachable = true;
        // skip genotypes containing alleles with zero support
        for (uint32_t a : g.as_vector())
            reachable &= (alleleDepth[a] > 0);
        if (!reachable) {
            continue;
        } else if (numExAlleles == 1) {
            weight += 1;
            gl[g] = 1;
        } else if (numExAlleles == 2) {
            double fracAlt = (double)index / (double)ploidy;
            double l = binom_pmf(alleleDepth[alleles[0]] + alleleDepth[alleles[1]], alleleDepth[alleles[1]], (1 - fracAlt) * err + fracAlt * (1 - err));
            weight += l;
            gl[g] = l;
        } else {
            std::vector<uint32_t> gv = g.as_vector();
            std::vector<double> p(numExAlleles);
            std::vector<uint32_t> n(numExAlleles);
            for (uint32_t a = 0; a < numExAlleles; a++) {
                double num = 0;
                for (uint32_t i = 0; i < gv.size(); i++)
                    if (gv[i] == (uint32_t)alleles[a])
                        num += 1.0;
                double freq = num / ploidy;
                p[a] = freq * (1 - err * (numExAlleles - 1)) + (1 - freq) * err;
                n[a] = alleleDepth[alleles[a]];
            }
            double l = multinom_pmf(n, p);
            weight += l;
            gl[g] = l;
        }
    }
    
    if (weight == 0.0) {
        uint32_t g = gl.size();
        for (std::pair<Genotype, double> p : gl)
            gl[p.first] = 1 / g;
    } else if (normalize) {
        uint32_t g = gl.size();
        for (std::pair<Genotype, double> p : gl)
            gl[p.first] = p.second / weight;
    }
    return gl;
}

void ReadScoring::computeAllelePairLikelihoods(std::vector<Genotype>& genos,
                                               std::vector<double>& apls,
                                               std::vector<double>& apld,
                                               uint8_t numAlleles,
                                               double err) const {
    uint32_t numGenos = genos.size();
    for (uint32_t a1 = 0; a1 < numAlleles; a1++) {
        for (uint32_t a2 = a1; a2 < numAlleles; a2++) {
            for (uint32_t gi = 0; gi < numGenos; gi++) {
                uint32_t i1 = numGenos * (numAlleles * a1 + a2) + gi;
                uint32_t i2 = numGenos * (numAlleles * a2 + a1) + gi;
                
                double lEqual = 0.0;
                double lDiff = 0.0;
                std::vector<uint32_t> gv = genos[gi].as_vector();
                for (uint32_t gv1 = 0; gv1 < gv.size(); gv1++) {
                    for (uint32_t gv2 = 0; gv2 < gv.size(); gv2++) {
                        // for matching alleles, multiply with (1 - err), for non-matching with err
                        double l = (1 - err) * (gv[gv1] == a1) + err * (gv[gv1] != a1);
                        l *= (1 - err) * (gv[gv2] == a2) + err * (gv[gv2] != a2);
                        // add to according sum and normalize to number of cases (prior probabilities are handles elsewhere)
                        if (gv1 == gv2)
                            lEqual += l / gv.size();
                        else
                            lDiff += l / (gv.size() * (gv.size() - 1));
                    }
                }
                apls[i1] = apls[i2] = lEqual;
                apld[i1] = apld[i2] = lDiff;
            }
        }
    }
}

float ReadScoring::computeLogScore (AlleleMatrix* am,
                                    uint32_t readId1,
                                    uint32_t readId2,
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
    uint32_t numAlleles = am->getMaxNumAllele();
    std::vector<AlleleItem> read1 = am->getRead(readId1);
    std::vector<AlleleItem> read2 = am->getRead(readId2);
    while (k < read1.size() && l < read2.size()) {
        if (read1[k].first == read2[l].first) {
            logScore += computeLogScoreSinglePos(read1[k].second, read2[l].second, numAlleles, gl[read1[k].first], gMap, apls, apld);
            ov += 1;
            k++; l++;
        } else if (read1[k].first < read2[l].first) {
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
                                             const uint32_t numAlleles,
                                             std::unordered_map<Genotype, double>& gl,
                                             std::unordered_map<Genotype, uint32_t>& gMap,
                                             std::vector<double>& apls,
                                             std::vector<double>& apld) const {
    double same = 0.0;
    double diff = 0.0;
    for (std::pair<Genotype, double> p : gl) {
        uint32_t i = gMap.size() * (allele1 * numAlleles + allele2) + gMap[p.first];
        same += p.second * apls[i];
        diff += p.second * apld[i];
    }
    if (same == 0.0 | diff == 0.0) {
        return 0.0;
    } else {
        return (float)std::log(same / diff);
    }
}
