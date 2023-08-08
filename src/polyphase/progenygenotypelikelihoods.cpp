#include "progenygenotypelikelihoods.h"
#include <cmath>

ProgenyGenotypeLikelihoods::ProgenyGenotypeLikelihoods(uint32_t ploidy, uint32_t numSamples, uint32_t numPositions) :
    numSamples(numSamples),
    ploidy(ploidy),
    numPositions(numPositions),
    genotypePairs(6),
    likelihoodSameSN(genotypePairs.size()),
    likelihoodDiffSN(genotypePairs.size()),
    likelihoodSameS2(genotypePairs.size()),
    likelihoodDiffS2(genotypePairs.size()),
    likelihoodSameDN(genotypePairs.size()),
    likelihoodDiffDN(genotypePairs.size())
{
    setNumPositions(numPositions);

    // possible six signal counts for the three supported variant types
    // right number is observed signals for simplex-nulliplex variant
    // left number is observed signals for the other variant
    genotypePairs[0] = std::pair<uint32_t, uint32_t>(0, 0);
    genotypePairs[1] = std::pair<uint32_t, uint32_t>(0, 1);
    genotypePairs[2] = std::pair<uint32_t, uint32_t>(1, 0);
    genotypePairs[3] = std::pair<uint32_t, uint32_t>(1, 1);
    genotypePairs[4] = std::pair<uint32_t, uint32_t>(2, 0);
    genotypePairs[5] = std::pair<uint32_t, uint32_t>(2, 1);
    GenotypeLikelihood k = (GenotypeLikelihood)ploidy;

    // likelihoods: simplex-nulliplex vs simplex-nulliplex
    likelihoodSameSN[0] =  0.5;
    likelihoodSameSN[1] =  0;
    likelihoodSameSN[2] =  0;
    likelihoodSameSN[3] =  0.5;
    likelihoodSameSN[4] =  0;
    likelihoodSameSN[5] =  0;
    likelihoodDiffSN[0] = (k/2-1) / (2*(k-1));
    likelihoodDiffSN[1] =  k      / (4*(k-1));
    likelihoodDiffSN[2] =  k      / (4*(k-1));
    likelihoodDiffSN[3] = (k/2-1) / (2*(k-1));
    likelihoodDiffSN[4] =  0;
    likelihoodDiffSN[5] =  0;

    // likelihoods: simplex-simplex vs simplex-nulliplex
    likelihoodSameS2[0] =  likelihoodSameSN[0] / 2.0;
    likelihoodSameS2[1] =  likelihoodSameSN[1] / 2.0;
    likelihoodSameS2[2] = (likelihoodSameSN[2] + likelihoodSameSN[0]) / 2.0;
    likelihoodSameS2[3] = (likelihoodSameSN[3] + likelihoodSameSN[1]) / 2.0;
    likelihoodSameS2[4] = (likelihoodSameSN[4] + likelihoodSameSN[2]) / 2.0;
    likelihoodSameS2[5] = (likelihoodSameSN[5] + likelihoodSameSN[3]) / 2.0;
    likelihoodDiffS2[0] =  likelihoodDiffSN[0] / 2.0;
    likelihoodDiffS2[1] =  likelihoodDiffSN[1] / 2.0;
    likelihoodDiffS2[2] = (likelihoodDiffSN[2] + likelihoodDiffSN[0]) / 2.0;
    likelihoodDiffS2[3] = (likelihoodDiffSN[3] + likelihoodDiffSN[1]) / 2.0;
    likelihoodDiffS2[4] = (likelihoodDiffSN[4] + likelihoodDiffSN[2]) / 2.0;
    likelihoodDiffS2[5] = (likelihoodDiffSN[5] + likelihoodDiffSN[3]) / 2.0;

    // likelihoods: duplex-nulliplex vs simplex-nulliplex
    likelihoodSameDN[0] = (k/2-1) / (2*(k-1));
    likelihoodSameDN[1] =  0;
    likelihoodSameDN[2] =  k      / (4*(k-1));
    likelihoodSameDN[3] =  k      / (4*(k-1));
    likelihoodSameDN[4] =  0;
    likelihoodSameDN[5] = (k/2-1) / (2*(k-1));
    likelihoodDiffDN[0] = (k/2-2)*(k/2-1) / (2*(k-1)*(k-2));
    likelihoodDiffDN[1] = (k/2)  *(k/2-1) / (2*(k-1)*(k-2));
    likelihoodDiffDN[2] = (k/2)  *(k/2-1) / (k-1)   *(k-2);
    likelihoodDiffDN[3] = (k/2)  *(k/2-1) / (k-1)   *(k-2);
    likelihoodDiffDN[4] = (k/2)  *(k/2-1) / (2*(k-1)*(k-2));
    likelihoodDiffDN[5] = (k/2-2)*(k/2-1) / (2*(k-1)*(k-2));
}

GenotypeLikelihood ProgenyGenotypeLikelihoods::getGl(uint32_t pos, uint32_t sampleId, uint32_t genotype) const {
    if (pos >= numPositions)
        return 0.0;
    return gl[getIndex(pos, sampleId, genotype)];
}

std::vector<GenotypeLikelihood> ProgenyGenotypeLikelihoods::getGlv(uint32_t pos, uint32_t sampleId) const {
    if (pos >= numPositions)
        return std::vector<GenotypeLikelihood>(ploidy+1, 0.0);
    uint32_t index = getIndex(pos, sampleId, 0);
    return std::vector<GenotypeLikelihood>(gl.begin()+index, gl.begin()+(index+ploidy+1));
}

uint32_t ProgenyGenotypeLikelihoods::getPloidy() const {
    return ploidy;
}
uint32_t ProgenyGenotypeLikelihoods::getNumSamples() const {
    return numSamples;
}
uint32_t ProgenyGenotypeLikelihoods::getNumPositions() const {
    return numPositions;
}

void ProgenyGenotypeLikelihoods::setGl(uint32_t pos, uint32_t sampleId, uint32_t genotype, GenotypeLikelihood l) {
    if (pos >= numPositions) {
        setNumPositions(pos);
    }
    
    gl[getIndex(pos, sampleId, genotype)] = l;
}

void ProgenyGenotypeLikelihoods::setGlv(uint32_t pos, uint32_t sampleId, std::vector<GenotypeLikelihood> l) {
    if (pos >= numPositions) {
        setNumPositions(pos);
    }
    uint32_t start = getIndex(pos, sampleId, 0);
    for (uint32_t i = 0; i <= ploidy; i++)
        gl[start+i] = l[i];
}

void ProgenyGenotypeLikelihoods::setNumPositions(uint32_t pos) {
    gl.resize(getIndex(pos+1, 0, 0), -1.0);
}

double ProgenyGenotypeLikelihoods::getSimplexNulliplexScore(uint32_t pos1, uint32_t pos2) const {
    return getLogLikelihoodDifference(pos1, pos2, likelihoodSameSN, likelihoodDiffSN, 4);
}

double ProgenyGenotypeLikelihoods::getSimplexSimplexScore(uint32_t pos1, uint32_t pos2) const {
    return getLogLikelihoodDifference(pos1, pos2, likelihoodSameS2, likelihoodDiffS2, 6);
}

double ProgenyGenotypeLikelihoods::getDuplexNulliplexScore(uint32_t pos1, uint32_t pos2) const {
    return getLogLikelihoodDifference(pos1, pos2, likelihoodSameDN, likelihoodDiffDN, 6);
}

uint32_t ProgenyGenotypeLikelihoods::getIndex(uint32_t pos, u_int32_t sampleId, uint32_t genotype) const {
    return pos*numSamples*(ploidy+1) + sampleId*(ploidy+1) + genotype;
}

double ProgenyGenotypeLikelihoods::getLogLikelihoodDifference(uint32_t pos1, uint32_t pos2, const std::vector<GenotypeLikelihood>& likelihoodSame, const std::vector<GenotypeLikelihood>& likelihoodDiff, const uint32_t numCases) const {
    double log_likelihood_difference = std::log(1.0/(double)(ploidy-1));
    for (uint32_t i = 0; i < getNumSamples(); i++) {
        if (getGl(pos1, i, 0) < 0.0 or getGl(pos2, i, 0) < 0.0)
            continue;
        double likelihoodCooccur = 0.0;
        double likelihoodDisjoint = 0.0;
        for (uint32_t j = 0; j < numCases; j++) {
            double gl = getGl(pos1, i, genotypePairs[j].first) * getGl(pos2, i, genotypePairs[j].second);
            likelihoodCooccur += gl * likelihoodSame[j];
            likelihoodDisjoint += gl * likelihoodDiff[j];
        }
        if (likelihoodCooccur * likelihoodDisjoint > 0)
            log_likelihood_difference += std::log(likelihoodCooccur / likelihoodDisjoint);
    }
    
    return log_likelihood_difference;
}
