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
    void scoreReadset(TriangleSparseMatrix *result, ReadSet *readset, const double errorrate, const uint32_t minOverlap, const uint32_t ploidy) const;

private:    
    float logratioSim(const uint32_t overlap, const uint32_t diff, const double distSame, const double distDiff) const;
    double binomPmf(const uint32_t n, const uint32_t k, const double p) const;
};

#endif // READSCORING_H
