#ifndef READSCORING_H
#define READSCORING_H

#include "../readset.h"
#include "../read.h"
#include "../genotype.h"
#include "trianglesparsematrix.h"

class ReadScoring {
 
public:
    /**
     * Computes pairwise scores for all reads in the readset and returns a sparse triangle matrix, where elements with a score of zero are not included.
     */
    void scoreReadsetGlobal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const;
    void scoreReadsetLocal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const;
    void scoreReadsetLocal(TriangleSparseMatrix *result, ReadSet *readset, std::vector<std::vector<uint32_t>>& refHaplotypes, const uint32_t minOverlap, const uint32_t ploidy) const;
    void scoreReadsetBayesian(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const;

private:
    double estimateAlleleErrorRate() const;
    
    std::unordered_map<Genotype, double> computeGenotypeLikelihoods (std::unordered_map<uint8_t, uint32_t> alleleDepths, uint32_t ploidy, const double err) const;
    
    void computeAllelePairLikelihoods(std::vector<Genotype>& genos,
                                      std::vector<double>& apls,
                                      std::vector<double>& apld,
                                      uint8_t maxAllele,
                                      double err) const;
    
    float computeLogScore (std::vector<uint32_t>& posRead1,
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
                           const double err) const;
                           
    float computeLogScoreSinglePos (uint8_t allele1,
                                    uint8_t allele2,
                                    const uint32_t maxAllele,
                                    std::unordered_map<Genotype, double>& gl,
                                    std::unordered_map<Genotype, uint32_t>& gMap,
                                    std::vector<double>& apls,
                                    std::vector<double>& apld,
                                    const double err) const;
    
    void computeStartEnd (const ReadSet* readset,
                          std::vector<uint32_t>& begins,
                          std::vector<uint32_t>& ends,
                          std::vector<std::vector<uint32_t>>& positions,
                          std::vector<std::vector<uint8_t>>& alleles,
                          std::vector<uint32_t>& posList,
                          std::unordered_map<uint32_t, uint32_t>& posMap,
                          uint32_t& longestReadSpan) const;

    void computeOverlapDiff (const ReadSet* readset,
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
                             const uint32_t end) const;
                             
    void computeOverlapDiff (const ReadSet* readset,
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
                             const uint32_t longestReadSpan) const;
    void computeCutoff(const std::vector<uint32_t>& coveredReads, const uint32_t ploidy, std::vector<double> relDiffs, double& distSame, double& distDiff) const;
    float logratioSim(const uint32_t overlap, const uint32_t diff, const double distSame, const double distDiff) const;
    double binomPmf(const uint32_t n, const uint32_t k, const double p) const;
};

#endif
