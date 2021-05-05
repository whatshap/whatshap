#ifndef PROGENYGENOTYPELIKELIHOODS_H
#define PROGENYGENOTYPELIKELIHOODS_H

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
    
    double getLogLikelihoodDifference(uint32_t pos1, uint32_t pos2, std::vector<std::pair<uint32_t, uint32_t>> genotypePairs, std::vector<std::pair<GenotypeLikelihood, GenotypeLikelihood>> likelihoodPairs) const;

private:
    
     std::vector<float> gl;
     uint32_t numSamples;
     uint32_t ploidy;
     uint32_t numPositions;
     
     uint32_t getIndex(uint32_t pos, u_int32_t sampleId, uint32_t genotype) const;
};

#endif // PROGENYGENOTYPELIKELIHOODS_H
