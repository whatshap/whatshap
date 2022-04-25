#ifndef PROGENYGENOTYPELIKELIHOODS_H
#define PROGENYGENOTYPELIKELIHOODS_H

#include <cstdint>
#include <vector>
#include <unordered_map>

typedef double GenotypeLikelihood;

class ProgenyGenotypeLikelihoods {
    
public:

    ProgenyGenotypeLikelihoods(uint32_t ploidy, uint32_t numSamples, uint32_t numPositions);
    
    GenotypeLikelihood getGl(uint32_t pos, uint32_t sampleId, uint32_t genotype) const;
    std::vector<GenotypeLikelihood> getGlv(uint32_t pos, uint32_t sampleId) const;
    uint32_t getPloidy() const;
    uint32_t getNumSamples() const;
    uint32_t getNumPositions() const;
    
    void setGl(uint32_t pos, uint32_t sampleId, uint32_t genotype, GenotypeLikelihood gl);
    void setGlv(uint32_t pos, uint32_t sampleId, std::vector<GenotypeLikelihood> gl);
    void setNumPositions(uint32_t pos);
    
    double getSimplexNulliplexScore(uint32_t pos1, uint32_t pos2) const;
    double getSimplexSimplexScore(uint32_t pos1, uint32_t pos2) const;
    double getDuplexNulliplexScore(uint32_t pos1, uint32_t pos2) const;

private:
    
    std::vector<float> gl;
    uint32_t numSamples;
    uint32_t ploidy;
    uint32_t numPositions;
    std::vector<std::pair<uint32_t, uint32_t>> genotypePairs;
    std::vector<GenotypeLikelihood> likelihoodSameSN;
    std::vector<GenotypeLikelihood> likelihoodDiffSN;
    std::vector<GenotypeLikelihood> likelihoodSameS2;
    std::vector<GenotypeLikelihood> likelihoodDiffS2;
    std::vector<GenotypeLikelihood> likelihoodSameDN;
    std::vector<GenotypeLikelihood> likelihoodDiffDN;
    
    uint32_t getIndex(uint32_t pos, uint32_t sampleId, uint32_t genotype) const;
    double getLogLikelihoodDifference(uint32_t pos1, uint32_t pos2, const std::vector<GenotypeLikelihood>& likelihoodSame, const std::vector<GenotypeLikelihood>& likelihoodDiff, const uint32_t numCases) const;
};

#endif // PROGENYGENOTYPELIKELIHOODS_H
