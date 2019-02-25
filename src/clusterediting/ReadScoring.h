#ifndef READSCORING_H
#define READSCORING_H

#include "../readset.h"
#include "../read.h"
#include "TriangleSparseMatrix.h"

class ReadScoring {
 
public:
    /**
     * Computes pairwise scores for all reads in the readset and returns a sparse triangle matrix, where elements with a score of zero are not included.
     */
    void scoreReadsetGlobal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const;
    void scoreReadsetLocal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const;
    void scoreReadsetPatterns(TriangleSparseMatrix* result, ReadSet* readset, const uint32_t minOverlap, const uint32_t ploidy, const double errorrate, const uint32_t windowSize) const;

private:
    void computeStartEndOverlapDiff (const ReadSet* readset, std::vector< uint32_t >& begins, std::vector< uint32_t >& ends, std::vector< std::vector< uint32_t > >& positions, std::vector< std::vector< uint32_t > >& alleles, TriangleSparseMatrix& overlaps, TriangleSparseMatrix& diffs, double& distSame, double& distDiff, const uint32_t minOverlap, const uint32_t ploidy) const;
    void computeCutoff(const uint32_t numReads, const uint32_t ploidy, std::vector<double> relDiffs, double& distSame, double& distDiff) const;
    float logratioSim(const uint32_t overlap, const uint32_t diff, const double distSame, const double distDiff) const;
    double binomPmf(const uint32_t n, const uint32_t k, const double p) const;
    
    // used for inefficient popcount
    const uint64_t m1  = 0x5555555555555555;
    const uint64_t m2  = 0x3333333333333333;
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
    const uint64_t h01 = 0x0101010101010101;
    /**
     * Returns the number of set bits in a 64bit-word.
     */
    uint64_t popcount(uint64_t bitv) const;
};

#endif // READSCORING_H
