#ifndef GENOTYPECOLUMNCOSTCOMPUTER_H
#define GENOTYPECOLUMNCOSTCOMPUTER_H

#include <array>
#include <vector>
#include <set>
#include <memory>
#include <map>
#include <utility>
#include <array>
#include "entry.h"
#include "pedigree.h"
#include "pedigreepartitions.h"
#include "columnindexingiterator.h"


class GenotypeColumnCostComputer
{
private:
  // the corresponding matrix column of reads
  const std::vector<const Entry*>& column;
  // the corresponding column index
  size_t column_index;
  const std::vector<unsigned int>& read_marks;
  // current partitioning of the reads
  unsigned int partitioning;
  // corresponding pedigree
  const Pedigree* pedigree;
  // stores the current costs of the partitions (for both alleles 0 and 1) = Z_j
  std::vector<std::array<long double, 2>> cost_partition;
  // the pedigree partitions
  const PedigreePartitions& pedigree_partitions;

public:
  GenotypeColumnCostComputer(const std::vector<const Entry*>& column, size_t column_index, const std::vector<unsigned int>& read_marks, const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions);
  // set partitioning to the given one
  void set_partitioning(unsigned int p);
  // update the partitioning by flipping read corresponding to the given bit
  void update_partitioning(int bit_to_flip);
  // returns the local cost for a given allele assignment prod_j Z_j
  long double get_cost(unsigned int allele_assignment);

};

#endif // GENOTYPECOLUMNCOSTCOMPUTER_H
