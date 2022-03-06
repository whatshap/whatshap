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
     * 
     * @param result Pointer of the matrix, where the results are written
     * @param readset Pointer to the scorable readset
     * @param minOverlap Minimum number of positions two reads must share in order to receive a score (otherwise they are scored as zero)
     * @param err Allele error rate, which should be assumed for the scoring. If set to zero, the error rate is estimated from the reads themselves.
     */
    void scoreReadset(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy, double err) const;
    
    /**
     * Computes an estimate for the allele error rate, based on the deviation of allele depths from valid genotypes.
     * 
     * @param alleleDepths A vector which contains the allele depths per position.
     * @param ploidy The underlying ploidy
     */
    double estimateAlleleErrorRate(std::vector<std::unordered_map<uint8_t, uint32_t>>& alleleDepths, uint32_t ploidy) const;

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
    std::unordered_map<Genotype, double> computeGenotypeLikelihoods (std::unordered_map<uint8_t, uint32_t> alleleDepths,
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
     * @param maxAllele Number of different alleles
     * @param err Assumed allele error rate
     */
    void computeAllelePairLikelihoods(std::vector<Genotype>& genos,
                                      std::vector<double>& apls,
                                      std::vector<double>& apld,
                                      uint8_t maxAllele,
                                      double err) const;

    /**
     * For two reads, computes their score based on their alleles on overlapping positions.
     * 
     * @param posRead1 Vector of (genome) positions for first read
     * @param posRead2 Vector of (genome) positions for second read
     * @param alleleRead1 Vector of alleles for the first read's positions
     * @param alleleRead2 Vector of alleles for the second read's positions
     * @param maxAllele Number of different alleles
     * @param posMap Mapping from genome positions to indexed positions, [0, l)
     * @param gl Genotype likelihoods for each (indexed) position
     * @param gMap Mapping from genotype to its index in apls and apld
     * @param apls Allele pair likelihoods, assuming both alleles come from same haplotypes
     * @param apld Allele pair likelihoods, assuming both alleles come from different haplotypes
     * @param minOverlap Minimum number of positions in order to receive a score (otherwise they are scored as zero)
     */
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
                           const uint32_t minOverlap) const;

    /**
     * Computes the score contribution of a single position for two reads.
     * 
     * @param allele1 Allele of first read
     * @param allele2 Allele of second read
     * @param maxAllele Number of different alleles
     * @param gl Genotype likelihoods for given position
     * @param gMap Mapping from genotype to its index in apls and apld
     * @param apls Allele pair likelihoods, assuming both alleles come from same haplotypes
     * @param apld Allele pair likelihoods, assuming both alleles come from different haplotypes
     */
    float computeLogScoreSinglePos (uint8_t allele1,
                                    uint8_t allele2,
                                    const uint32_t maxAllele,
                                    std::unordered_map<Genotype, double>& gl,
                                    std::unordered_map<Genotype, uint32_t>& gMap,
                                    std::vector<double>& apls,
                                    std::vector<double>& apld) const;

    /**
     * Writes the information of a readset into vectors for more efficient access.
     */
    void computeStartEnd (const ReadSet* readset,
                          std::vector<uint32_t>& begins,
                          std::vector<uint32_t>& ends,
                          std::vector<std::vector<uint32_t>>& positions,
                          std::vector<std::vector<uint8_t>>& alleles,
                          std::vector<uint32_t>& posList,
                          std::unordered_map<uint32_t, uint32_t>& posMap,
                          uint32_t& longestReadSpan) const;
};

#endif
