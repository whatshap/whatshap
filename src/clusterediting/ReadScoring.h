#ifndef READSCORING_H
#define READSCORING_H

#include "../readset.h"
#include "../read.h"
#include "TriangleSparseMatrix.h"
#include "Globals.h"

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
};

#endif // READSCORING_H
