#ifndef PEDIGREE_COLUMN_COST_COMPUTER_H
#define PEDIGREE_COLUMN_COST_COMPUTER_H

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


  
class PedigreeColumnCostComputer {
private:
  const std::vector<const Entry*>& column;
  size_t column_index;
  const std::vector<unsigned int>& read_marks;  
  unsigned int partitioning;
  const Pedigree* pedigree;
  std::vector<std::array<unsigned int, 2>> cost_partition;
  const PedigreePartitions& pedigree_partitions;
  // all (bit-encoded) assignments of alleles to pedigree partititions that are
  // compatible with the genotypes
  std::vector<unsigned int> allele_assignments;
  
public:
  
  PedigreeColumnCostComputer(const std::vector<const Entry*>& column, size_t column_index, const std::vector<unsigned int>& read_marks, const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions);
  
  void set_partitioning(unsigned int partitioning);

  void update_partitioning(int bit_to_flip);
  
  std::vector<unsigned int> compute_roots(std::vector<Pedigree::triple_entry_t> triples);
  
  unsigned int get_cost();
  
  /** Returns a pair of alleles for all individuals in the ordered
   *  by their index in the pedigree given at construction time.
   */
  std::vector<std::pair<Entry::allele_t,Entry::allele_t>> get_alleles();

  // returns the weight (what will be the phred_score) at the current position for super-read s
  unsigned int get_weight(bool second_haplotype);
};

#endif
