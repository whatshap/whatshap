#ifndef READSCORING_H
#define READSCORING_H

#include <cstdint>
#include "../readset.h"
#include "../read.h"
#include "../genotype.h"
#include "trianglesparsematrix.h"
#include "allelematrix.h"

class ReadScoring {
 
public:
    /**
     * Computes pairwise scores for all reads in the allele matrix and returns a sparse triangle matrix, where elements with a score of zero are not included.
     * 
     * @param result Pointer of the matrix, where the results are written
     * @param am Allele matrix to score
     * @param minOverlap Minimum number of positions two reads must share in order to receive a score (otherwise they are scored as zero)
     * @param err Allele error rate, which should be assumed for the scoring. If set to zero, the error rate is estimated from the reads themselves.
     */
    void scoreReadset(TriangleSparseMatrix *result, AlleleMatrix *am, const uint32_t minOverlap, const uint32_t ploidy, double err) const;
    
    /**
     * Computes an estimate for the allele error rate, based on the deviation of allele depths from valid genotypes.
     * 
     * @param am Allele matrix containing allele depths.
     * @param ploidy The underlying ploidy
     */
    double estimateAlleleErrorRate(AlleleMatrix *am, uint32_t ploidy) const;

private:
    /**
     * For a list of genotype likelihoods computes a score how plausible these likelihoods are. This is done by taking the logarithm of the highest
     * likelihood for each position and return the sum.
     *
     * @param gl Genotype likelihoods per position
     */
    double evaluateGenotypeLikelihoods(std::vector<std::unordered_map<Genotype, double>>& gl) const;

    /**
     * Computes a likelihood for every possible genotype based on the allele depths.
     * 
     * @param alleleDepths The allele depths to base the calculation on
     * @param ploidy Ploidy of the genotypes
     * @param err Allele error rate, which should be assumed
     * @param normalize If true, the likelihoods are normalized, such that they sum up to 1
     */
    std::unordered_map<Genotype, double> computeGenotypeLikelihoods (std::vector<uint32_t> alleleDepths,
                                                                     const uint32_t ploidy,
                                                                     const double err,
                                                                     const bool normalize) const;

    /**
     * For genotype and each pair of alleles, computes the likelihood to observe this exact allele pair under this exact genotype, under the
     * assumption that the alleles come from the same haplotype and that the alleles come from different haplotypes. The results are written
     * into the provided vectors respectively. The indices for the vectors are coded as
     * 
     * numGenos*(maxAllele*a1 + a2) + geno
     * 
     * with numGenos being the total number of genotypes, a1 and a2 being the two alleles and geno being the index of the genotype from the
     * input.
     * 
     * @param genos Vector of possible genotypes
     * @param apls Result vector for likelihoods under same haplotype assumption
     * @param apld Result vector for likelihoods under different haplotype assumption
     * @param numAlleles Number of different alleles
     * @param err Assumed allele error rate
     */
    void computeAllelePairLikelihoods(std::vector<Genotype>& genos,
                                      std::vector<double>& apls,
                                      std::vector<double>& apld,
                                      uint8_t numAlleles,
                                      double err) const;

    /**
     * For two reads, computes their score based on their alleles on overlapping positions.
     * 
     * @param am Allele matrix
     * @param readId1 id of first read
     * @param readId2 id of second read
     * @param gl Genotype likelihoods for each (indexed) position
     * @param gMap Mapping from genotype to its index in apls and apld
     * @param apls Allele pair likelihoods, assuming both alleles come from same haplotypes
     * @param apld Allele pair likelihoods, assuming both alleles come from different haplotypes
     * @param minOverlap Minimum number of positions in order to receive a score (otherwise they are scored as zero)
     */
    float computeLogScore(AlleleMatrix* am,
                          uint32_t readId1,
                          uint32_t readId2,
                          std::vector<std::unordered_map<Genotype, double>>& gl,
                          std::unordered_map<Genotype, uint32_t>& gMap,
                          std::vector<double>& apls,
                          std::vector<double>& apld,
                          const uint32_t minOverlap) const;

    /**
     * Computes the score contribution of a single position for two reads.
     * 
     * @param allele1 Allele of first read
     * @param allele2 Allele of second read
     * @param numAlleles Number of different alleles
     * @param gl Genotype likelihoods for given position
     * @param gMap Mapping from genotype to its index in apls and apld
     * @param apls Allele pair likelihoods, assuming both alleles come from same haplotypes
     * @param apld Allele pair likelihoods, assuming both alleles come from different haplotypes
     */
    float computeLogScoreSinglePos (uint8_t allele1,
                                    uint8_t allele2,
                                    const uint32_t numAlleles,
                                    std::unordered_map<Genotype, double>& gl,
                                    std::unordered_map<Genotype, uint32_t>& gMap,
                                    std::vector<double>& apls,
                                    std::vector<double>& apld) const;
};
#endif
