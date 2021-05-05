#include "progenygenotypelikelihoods.h"
#include <cmath>

ProgenyGenotypeLikelihoods::ProgenyGenotypeLikelihoods(uint32_t ploidy, uint32_t numSamples, uint32_t numPositions) :
    numSamples(numSamples),
    ploidy(ploidy),
    numPositions(numPositions)
{
    setNumPositions(numPositions);
}



GenotypeLikelihood ProgenyGenotypeLikelihoods::getGl(uint32_t pos, uint32_t sampleId, uint32_t genotype) const {
    if (pos >= numPositions)
        return 0.0;
    return gl[getIndex(pos, sampleId, genotype)];
}

std::vector<GenotypeLikelihood> ProgenyGenotypeLikelihoods::getGlv(uint32_t pos, uint32_t sampleId) const {
    if (pos >= numPositions)
        return std::vector<GenotypeLikelihood>(ploidy, 0.0);
    uint32_t index = getIndex(pos, sampleId, 0);
    return std::vector<GenotypeLikelihood>(gl.begin()+index, gl.begin()+(index+ploidy));
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
    for (uint32_t i = 0; i < ploidy; i++)
        gl[start+i] = l[i];
}

void ProgenyGenotypeLikelihoods::setNumPositions(uint32_t pos) {
    gl.resize(getIndex(pos+1, 0, 0), -1.0);
}

double ProgenyGenotypeLikelihoods::getLogLikelihoodDifference(uint32_t pos1, uint32_t pos2, std::vector<std::pair<uint32_t, uint32_t>> genotypePairs, std::vector<std::pair<GenotypeLikelihood, GenotypeLikelihood>> likelihoodPairs) const {
    double log_likelihood_difference = std::log(1.0/3.0);
    for (uint32_t i = 0; i < getNumSamples(); i++) {
        if (getGl(pos1, i, 0) < 0.0 or getGl(pos2, i, 0) < 0.0)
            continue;
        double likelihoodCooccur = 0.0;
        double likelihoodDisjoint = 0.0;
        for (uint32_t j = 0; j < genotypePairs.size(); j++) {
            double gl = getGl(pos1, i, genotypePairs[j].first) * getGl(pos2, i, genotypePairs[j].second);
            likelihoodCooccur += gl * likelihoodPairs[j].first;
            likelihoodDisjoint += gl * likelihoodPairs[j].second;
        }
        log_likelihood_difference += std::log(likelihoodCooccur / likelihoodDisjoint);
    }
    
    return log_likelihood_difference;
}



uint32_t ProgenyGenotypeLikelihoods::getIndex(uint32_t pos, u_int32_t sampleId, uint32_t genotype) const {
    return pos*numSamples*ploidy + sampleId*ploidy + genotype;
}
