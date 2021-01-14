#include "readscoring.h"
#include "../binomial.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <assert.h>

void ReadScoring::scoreReadsetGlobal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const {
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint8_t>> alleles;
    std::vector<uint32_t> posList;
    std::unordered_map<uint32_t, uint32_t> posMap;
    uint32_t longestReadSpan = 0;
    computeStartEnd(readset, begins, ends, positions, alleles, posList, posMap, longestReadSpan);
    
    // compute length of overlap and difference for all read pairs
    double hammingDistSame = 0;
    double hammingDistDiff = 0;
    computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, result, hammingDistSame, hammingDistDiff, minOverlap, ploidy, longestReadSpan);
    hammingDistSame = 0.10;
    hammingDistDiff = 0.40;
    
    // compute pair wise scores
    std::vector<std::pair<uint32_t, uint32_t>> entries = result->getEntries();
    std::unordered_map<uint64_t, float> cache;
    for (std::pair<uint32_t, uint32_t> p : entries) {
        uint32_t i = p.first;
        uint32_t j = p.second;
        TriangleSparseMatrix::DoubleInt d = result->getDoubleInt(i, j);
        uint32_t ov = d.u1;
        uint32_t di = d.u2;
        uint64_t ovdi = (ov*(ov+1))/2+di;
        std::unordered_map<uint64_t, float>::const_iterator it = cache.find(ovdi);
        if (it == cache.end()) {
            cache[ovdi] = logratioSim(ov, di, hammingDistSame, hammingDistDiff);
        }
        result->set(i, j, cache[ovdi]);
    }
    
}

void ReadScoring::scoreReadsetLocal(TriangleSparseMatrix* result, ReadSet* readset, const uint32_t minOverlap, const uint32_t ploidy) const {
    std::vector<std::vector<uint32_t>> emptyRef;
    scoreReadsetLocal(result, readset, emptyRef, minOverlap, ploidy);
}

void ReadScoring::scoreReadsetLocal(TriangleSparseMatrix* result, ReadSet* readset, std::vector<std::vector<uint32_t>>& refHaplotypes, const uint32_t minOverlap, const uint32_t ploidy) const {
    
    if (ploidy < 2) {
        std::cout<<"Error: Ploidy < 2!"<<std::endl;
        return;
    }
    
    uint32_t numReads = readset->size();
    
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint8_t>> alleles;
    std::vector<uint32_t> posList;
    std::unordered_map<uint32_t, uint32_t> posMap;
    uint32_t longestReadSpan = 0;
    computeStartEnd(readset, begins, ends, positions, alleles, posList, posMap, longestReadSpan);
    
    // check ref haplotypes
    if (refHaplotypes.size() > 0) {
        if (refHaplotypes.size() != ploidy) {
            std::cout<<"Error: Inconsistent ploidy in reference haplotypes! Was "
            <<refHaplotypes.size()<<" but expected "<<ploidy<<std::endl;
            return;
        } else if (refHaplotypes[0].size() != posList.size()) {
            std::cout<<"Error: Number of positions in reference haplotypes does not match number of positions in read set! Was "
            <<refHaplotypes[0].size()<<" but expected "<<posList.size()<<std::endl;
            return;
        }
    }
    
    // compute length of overlap and difference for all read pairs
    // reuse the result matrix to store overlaps and diffs. they will be overwritten later on
    double defaultSameDist = 0;
    double defaultDiffDist = 0;
    computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, result, defaultSameDist, defaultDiffDist, minOverlap, ploidy, longestReadSpan);
    
    // compute longest read length and average read length (in base pairs) and divide by 2
    uint32_t windowSize = 0;
    for (uint32_t i = 0; i < numReads; i++) {
        uint32_t len = readset->get(i)->lastPosition() - readset->get(i)->firstPosition();
        windowSize += len;
    }
    windowSize /= (4*numReads);
    
    //divide snp positions by window size
    std::vector<uint32_t> windowStarts;
    std::unordered_map<uint32_t, double> posToSameDist;
    std::unordered_map<uint32_t, double> posToDiffDist;
    
    uint32_t windowStartPosition = 0;
    for (uint32_t current = 0; current < posList.size(); current++) {
        if (posList[current] - windowStartPosition > windowSize || current == 0) {
            windowStarts.push_back(current);
            windowStartPosition = posList[current];
        }
    }
    windowStarts.push_back(posList.size());
    
    // determine relative hamming distance for same and different haplotypes for each window
    for (uint32_t windowIdx = 0; windowIdx < windowStarts.size()-1; windowIdx++) {
        // window bounds
        uint32_t startVariant = windowStarts[windowIdx];
        uint32_t endVariant = windowStarts[windowIdx+1];
        uint32_t start = posList[startVariant];
        uint32_t end = posList[endVariant-1];

        // compute length of overlap and difference for all read pairs
        TriangleSparseMatrix overlapsDiffsLocal;
        double localSameDist = 0;
        double localDiffDist = 0;
        computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, &overlapsDiffsLocal, 
                           localSameDist, localDiffDist, minOverlap, ploidy, longestReadSpan, start, end);
        
        if (overlapsDiffsLocal.getEntries().size() < ploidy) {
            // too few read pairs, use average over all reads instead
            localSameDist = defaultSameDist;
            localDiffDist = defaultDiffDist;
        } else if (refHaplotypes.size() == ploidy) {
            // use reference haplotypes to determine localDiffDist
            std::vector<double> pairDiffs;
            for (uint32_t h1 = 0; h1 < ploidy-1; h1++) {
                for (uint32_t h2 = h1+1; h2 < ploidy; h2++) {
                    double diffs = 0;
                    for (uint32_t pos = startVariant; pos < endVariant; pos++) {
                        if (refHaplotypes[h1][pos] != refHaplotypes[h2][pos])
                            diffs += 1.0;
                    }
                    pairDiffs.push_back(diffs/(double)(endVariant-startVariant));
                }
            }
            std::sort(pairDiffs.begin(), pairDiffs.end());
            bool found = false;
            double bestDiffDist = localDiffDist;
            for (uint32_t i = 0; !found && i < pairDiffs.size(); i++) {
                if (pairDiffs[i] > localSameDist/2) {
                    bestDiffDist = pairDiffs[i];
                    found = true;
                }
            }
            if (!found) {
                bestDiffDist = pairDiffs.back();
            }
            localSameDist = std::max(0.001, localSameDist);
            localDiffDist = std::min(localDiffDist, bestDiffDist*(1-localSameDist)+(1-bestDiffDist)*localSameDist);
        }
        
        // store values in a map for every snp position
        for (uint32_t j = windowStarts[windowIdx]; j < windowStarts[windowIdx+1]; j++) {
            posToSameDist[posList[j]] = localSameDist;
            posToDiffDist[posList[j]] = localDiffDist;
        }
    }
            
    // now, iterate over all overlapping read pairs and compute their score
    std::vector<std::pair<uint32_t, uint32_t>> entries = result->getEntries();
    for (std::pair<uint32_t, uint32_t> p : entries) {
        uint32_t i = p.first;
        uint32_t j = p.second;
        TriangleSparseMatrix::DoubleInt ovdi = result->getDoubleInt(i, j);
        uint32_t ov = ovdi.u1;
        uint32_t di = ovdi.u2;
        
        // zig zag over read positions to determine their individual same/diff dists
        double same = 0.0;
        double diff = 0.0;
        uint32_t k = 0;
        uint32_t l = 0;
        while (k < positions[i].size() && l < positions[j].size()) {
            if (positions[i][k] == positions[j][l]) {
                same += posToSameDist[positions[i][k]];
                diff += posToDiffDist[positions[i][k]];
                k++; l++;
            } else if (positions[i][k] < positions[j][l]) {
                k++;
            } else {
                l++;
            }
        }
        same /= ov;
        diff /= ov;
        same = std::max(same, 0.001);
		diff = std::min(0.999, std::max(diff, same + 0.001));
        result->set(i, j, logratioSim(ov, di, same, diff));
    }
}

void ReadScoring::scoreReadsetBayesian(TriangleSparseMatrix* result, ReadSet* readset, const uint32_t minOverlap, const uint32_t ploidy, double err) const {
    
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
    
    /*
    for (uint32_t pos = 0; pos < 100; pos++) {
        std::cout<<"AD("<<pos<<"): [";
        for (auto p: alleleDepths[pos])
            std::cout<<(uint32_t)(p.first)<<":"<<p.second<<" ";
        std::cout<<"]"<<std::endl;
    }*/
    
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
    
    /*
    for (uint32_t i = 0; i < 100; i++) {
        std::cout<<"GL("<<i<<")=[";
        for (auto p : gl[i])
            std::cout<<p.first.toString()<<":"<<p.second<<" ";
        std::cout<<"]"<<std::endl;
    }*/
    
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
    
    /*
    for (uint32_t a1 = 0; a1 < maxAllele; a1++) {
        for (uint32_t a2 = 0; a2 < maxAllele; a2++) {
            std::cout<<"ALPs("<<a1<<","<<a2<<")=[";
            uint32_t offset = occGenotypes.size() * (maxAllele*a1 + a2);
            for (auto p : gMap)
                std::cout<<p.first.toString()<<":"<<apls[p.second+offset]<<" ";
            std::cout<<"]"<<std::endl;
        }
    }
    for (uint32_t a1 = 0; a1 < maxAllele; a1++) {
        for (uint32_t a2 = 0; a2 < maxAllele; a2++) {
            std::cout<<"ALPd("<<a1<<","<<a2<<")=[";
            uint32_t offset = occGenotypes.size() * (maxAllele*a1 + a2);
            for (auto p : gMap)
                std::cout<<p.first.toString()<<":"<<apld[p.second+offset]<<" ";
            std::cout<<"]"<<std::endl;
        }
    }*/
        
    // for each read pair, iterate over their overlapping positions and determine combined score
    for (uint32_t i = 0; i < begins.size(); i++) {
        // iterate until start position of read is behind required start
        for (uint32_t j = i+1; j < begins.size() && begins[j] <= ends[i]; j++) {
//             if (ends[i] < begins[j] || ends[j] < begins[i])
//                 continue;
            float score = computeLogScore(positions[i], positions[j], alleles[i], alleles[j], maxAllele, posMap, gl, gMap, apls, apld, minOverlap, err);
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
    double LogLikelihoodSum = 0.0;
    for (uint32_t i = 0; i < gl.size(); i++) {
        double max = 0.0;
        for (auto& l: gl[i]) {
            if (l.second > max)
                max = l.second;
        }
        LogLikelihoodSum += std::log(max);
    }
    return LogLikelihoodSum;
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
            double l = binomPmf(alleleDepth[0]+alleleDepth[1], alleleDepth[1], (1-fracAlt)*err + fracAlt*(1-err));
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
                                    const uint32_t minOverlap,
                                    const double err) const {
    ;
    uint32_t ov = 0;
    double logScore = 0.0;
    uint32_t k = 0;
    uint32_t l = 0;
    while (k < posRead1.size() && l < posRead2.size()) {
        if (posRead1[k] == posRead2[l]) {
            logScore += computeLogScoreSinglePos(alleleRead1[k], alleleRead2[l], maxAllele, gl[posMap[posRead1[k]]], gMap, apls, apld, err);
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
                                             std::vector<double>& apld,
                                             const double err) const {
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

void ReadScoring::computeOverlapDiff (const ReadSet* readset,
                                      const std::vector<uint32_t>& begins,
                                      const std::vector<uint32_t>& ends,
                                      const std::vector<std::vector<uint32_t>>& positions,
                                      const std::vector<std::vector<uint8_t>>& alleles,
                                      const std::vector<uint32_t>& posList,
                                      std::unordered_map<uint32_t, uint32_t>& posMap,
                                      TriangleSparseMatrix* overlapDiffs,
                                      double& distSame,
                                      double& distDiff,
                                      const uint32_t minOverlap,
                                      const uint32_t ploidy,
                                      const uint32_t longestReadSpan,
                                      const uint32_t begin,
                                      const uint32_t end) const {
                                          
    std::vector<uint32_t> coveredReads;
    
    // determine coveredReads if begin and end position are specified (otherwise use all reads)
    if (begin == 0 && end == begins.size()) {
        for (uint32_t i = 0; i < begins.size(); i++) {
            coveredReads.push_back(i);
        }
    } else {
        // binary search to find first read, which can theoretically cover this window
        uint32_t firstIndex = std::lower_bound(begins.begin(), begins.end(), begin-longestReadSpan) - begins.begin();
        // iterate until start position of read is behind required start
        for (uint32_t j = firstIndex; begins[j] <= begin && j < begins.size(); j++) {
            if (ends[j] >= end) {
                coveredReads.push_back(j);
            }
        }
    }
    
    // iterate over all read pairs (efficiently omitting those who can certainly not overlap)
    std::vector<double> relativeDiffs;
    for (uint32_t i = 0; i < coveredReads.size(); i++) {
        // iterate until start position of read is behind required start
        uint32_t ci = coveredReads[i];
        for (uint32_t j = i+1; j < coveredReads.size() && begins[coveredReads[j]] <= ends[ci]; j++) {
            uint32_t cj = coveredReads[j];
            if (ends[ci] < begins[cj] || ends[cj] < begins[ci])
                continue;
            uint32_t ov = 0;
            uint32_t di = 0;
            uint32_t k = 0;
            uint32_t l = 0;
            while (k < positions[ci].size() && l < positions[cj].size()) {
                if (positions[ci][k] == positions[cj][l]) {
                    if (alleles[ci][k] != alleles[cj][l]) {
                        di++;
                    }
                    ov++; k++; l++; 
                } else if (positions[ci][k] < positions[cj][l]) {
                    k++;
                } else {
                    l++;
                }
            }
            if (ov >= minOverlap) {
                overlapDiffs->setDoubleInt(ci, cj, (uint16_t)ov, (uint16_t)di);
                relativeDiffs.push_back((double)di / (double)ov);
            }
        }
    }
    
    // estimate error rate of reads within equal haplotype
    computeCutoff(coveredReads, ploidy, relativeDiffs, distSame, distDiff);
}

void ReadScoring::computeOverlapDiff (const ReadSet* readset,
                                      const std::vector<uint32_t>& begins,
                                      const std::vector<uint32_t>& ends,
                                      const std::vector<std::vector<uint32_t>>& positions,
                                      const std::vector<std::vector<uint8_t>>& alleles,
                                      const std::vector<uint32_t>& posList,
                                      std::unordered_map<uint32_t, uint32_t>& posMap,
                                      TriangleSparseMatrix* overlapDiffs,
                                      double& distSame,
                                      double& distDiff,
                                      const uint32_t minOverlap,
                                      const uint32_t ploidy,
                                      const uint32_t longestReadSpan) const {
    computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, overlapDiffs, distSame, distDiff, minOverlap, ploidy, longestReadSpan, 0, begins.size());
}

void ReadScoring::computeCutoff(const std::vector<uint32_t>& coveredReads, const uint32_t ploidy, std::vector<double> relDiffs, double& distSame, double& distDiff) const {
    uint32_t numReads = coveredReads.size();
    std::sort(relDiffs.begin(), relDiffs.end());
    double sameSum = 0.0;
    int sameNum = 0;
    double diffSum = 0.0;
    int diffNum = 0;
    
    double p = (double)ploidy;
    double n = (double)numReads;
    uint32_t cutoff = std::min((uint32_t)1, (uint32_t)(relDiffs.size() - 1));
    if (ploidy < numReads) {
        cutoff = (uint32_t) std::ceil((p * (n/p) * (n/p-1) / 2) / 
                    ( (p * (n/p) * (n/p-1) / 2) + (p*(p-1)/2)*(n/p)*(n/p) ) * relDiffs.size());
        for (unsigned int i = 0; i < relDiffs.size(); i++) {
            if (i < cutoff) {
                sameSum += relDiffs[i];
                sameNum++;
            } else {
                diffSum += relDiffs[i];
                diffNum++;
            }
        }
        
        if (cutoff == 0) {
            distSame = 0.1;
        } else {
            distSame = sameSum/sameNum;
        }
        distDiff = diffSum/diffNum;
    }
}

float ReadScoring::logratioSim(const uint32_t overlap, const uint32_t diff, const double distSame, const double distDiff) const {
    double pSame = binomPmf(overlap, diff, distSame);
    double pDiff = binomPmf(overlap, diff, distDiff);
    if (pSame == 0)
        return -std::numeric_limits<float>::infinity();
    else if (pDiff == 0)
        return std::numeric_limits<float>::infinity();
    else
        return log(pSame / pDiff);
}

double ReadScoring::binomPmf(const uint32_t n, const uint32_t k, const double p) const {
    double coeff = 1.0;
    for (uint32_t i = 0; i < k; i++) {
        coeff *= ((n-i) / (k-i));
    }
    return coeff * pow(p, k) * pow(1-p, n-k);
}
