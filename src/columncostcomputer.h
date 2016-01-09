#ifndef COLUMN_COST_COMPUTER_H
#define COLUMN_COST_COMPUTER_H

#include <array>
#include <vector>
#include <memory>
#include "entry.h"
#include "columnindexingiterator.h"

class ColumnCostComputer {
private:
  const std::vector<const Entry*>& column;
  const std::vector<unsigned int>& read_marks;
  
  // cost_partitionX[Y] is the cost of flipping all entries in the partition X to Y for Y = 0,1.
  // 1 and 2 -> mother
  unsigned int cost_partition_m1[2];
  unsigned int cost_partition_m2[2];
  // 3 and 4 -> father
  unsigned int cost_partition_f1[2];
  unsigned int cost_partition_f2[2];
  
  unsigned int inheritance_val;
  unsigned int partitioning;
  
public:
  typedef struct trio_alleles_t {
    std::pair<Entry::allele_t,Entry::allele_t> mother;
    std::pair<Entry::allele_t,Entry::allele_t> father;
    std::pair<Entry::allele_t,Entry::allele_t> child;
    trio_alleles_t(Entry::allele_t m1, Entry::allele_t m2, Entry::allele_t f1, Entry::allele_t f2, Entry::allele_t c1, Entry::allele_t c2) :
      mother(std::make_pair(m1, m2)), father(std::make_pair(f1, f2)), child(std::make_pair(c1, c2)) {}
  } trio_alleles_t;
  
  ColumnCostComputer(const std::vector<const Entry*>& column, const std::vector<unsigned int>& read_marks, unsigned int inheritance_val);
  
  void set_partitioning(unsigned int partitioning);

  void update_partitioning(int bit_to_flip);
  
  unsigned int get_cost(unsigned int genotypem, unsigned int genotypef, unsigned int genotypec);
  
  /** Returns the six alleles for all three individuals at the current position.
   *  Alleles are ensured to agree with the given genotypes.
   */
  trio_alleles_t get_alleles(unsigned int genotypem, unsigned int genotypef, unsigned int genotypec);

  // returns the weight (what will be the phred_score) at the current position for super-read s
  unsigned int get_weight(bool second_haplotype);

  // HALF_TABLE stuff ...

  // set partitioning to the dual of partitioning (for the backtracking)
  void set_dual_partitioning(unsigned int partitioning);
  
};

#endif
