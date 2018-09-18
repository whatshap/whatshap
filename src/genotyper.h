#ifndef GENOTYPER_H
#define GENOTYPER_H

#include <vector>

#include "readset.h"
#include "genotypelikelihoods.h"
#include "genotype.h"

/** Computes genotypes (and genotype likelihoods) for a given read set. 
 *  @param genotypes Vector is cleared and one value from {0,1,2} is added for each variant in the readset.
 *  @param genotype_likelihoods Vector is cleared and one set of likelihoods is added for each variant in the readset.
 */
void compute_genotypes(const ReadSet& readset, std::vector<Genotype>* genotypes, std::vector<GenotypeDistribution>* genotype_likelihoods, std::vector<unsigned int>* positions = nullptr);

#endif
