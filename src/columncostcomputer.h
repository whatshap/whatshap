#ifndef COLUMN_COST_COMPUTER_H
#define COLUMN_COST_COMPUTER_H

#include <array>
#include <vector>
#include <set>
#include <memory>
#include "entry.h"
#include "columnindexingiterator.h"

class ColumnCostComputer {
private:
  const std::vector<const Entry*>& column;
  const std::vector<unsigned int>& read_marks;
  
  unsigned int inheritance_val;
  unsigned int partitioning;
  unsigned int num_of_roots;
  std::vector<unsigned int> roots:
  int num_of_triples;
  std::vector<int> triple;
  std::vector<std::vector<int>> triples;
  
  // cost_partitionX[Y] is the cost of flipping all entries in the partition X to Y for Y = 0,1.
  std::vector<std::pair<unsigned int,unsigned int>> cost_partition;
  std::map<unsigned, std::pair<unsigned int>> haps; 
  
public:
  struct ind_alleles_t {
    std::pair<Entry::allele_t,Entry::allele_t> ind;
    ind_alleles_t(Entry::allele_t ind1, Entry::allele_t ind2) :ind(std::make_pair(ind1, ind2)) {}
  };
  
  ColumnCostComputer(const std::vector<const Entry*>& column, const std::vector<unsigned int>& read_marks, unsigned int inheritance_val, std::vector<std::vector<int>>& triples);
  
  void set_partitioning(unsigned int partitioning);

  void update_partitioning(int bit_to_flip);
  
  std::vector<unsigned int> compute_roots(std::vector<triple>& triples);
  
  unsigned int get_cost(std::vector<std::vector<int>>& genotypes);
  
  /** Returns the six alleles for all three individuals at the current position.
   *  Alleles are ensured to agree with the given genotypes.
   */
  trio_alleles_t get_alleles(std::vector<std::vector<int>>& genotypes);

  // returns the weight (what will be the phred_score) at the current position for super-read s
  unsigned int get_weight(bool second_haplotype);

  // HALF_TABLE stuff ...

  // set partitioning to the dual of partitioning (for the backtracking)
  void set_dual_partitioning(unsigned int partitioning);
  
};

#endif
