#ifndef READSCORING_H
#define READSCORING_H

#include "../readset.h"
#include "../read.h"
#include "trianglesparsematrix.h"

class ReadScoring {
 
public:
    /**
     * Computes pairwise scores for all reads in the readset and returns a sparse triangle matrix, where elements with a score of zero are not included.
     */
    void scoreReadsetGlobal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const;
    void scoreReadsetLocal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const;
    void scoreReadsetLocal(TriangleSparseMatrix *result, ReadSet *readset, std::vector<std::vector<uint32_t>>& refHaplotypes, const uint32_t minOverlap, const uint32_t ploidy) const;

private:
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
