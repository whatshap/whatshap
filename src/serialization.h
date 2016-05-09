#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include "readset.h"
#include "pedigree.h"

void serialize(const ReadSet* readset, const std::vector<unsigned int>& recombcost, const Pedigree* pedigree, bool distrust_genotypes, const std::vector<unsigned int>& positions);

#endif
