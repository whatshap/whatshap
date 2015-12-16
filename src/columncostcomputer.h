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
  
    unsigned int cost_partition1[2];
  unsigned int cost_partition2[2];
  
  unsigned int inheritance_val;
  unsigned int partitioning;
  bool all_heterozygous;
  
public:
  ColumnCostComputer(const std::vector<const Entry*>& column, const std::vector<unsigned int>& read_marks, unsigned int inheritance_val, bool all_heterozygous = false);
  
  void set_partitioning(unsigned int partitioning);
  void set_partitioning_m(unsigned int partitioning);
  void set_partitioning_f(unsigned int partitioning);
  void set_partitioning_c(unsigned int partitioning);

  void update_partitioning(int bit_to_flip);
  
  unsigned int get_cost();
  
  /** Returns the allele at the current position for the given haplotype.
   *  @param second_haplotype If true the second haplotype is returned, otherwise the first.
   */
  Entry::allele_t get_allele(bool second_haplotype);

  // returns the weight (what will be the phred_score) at the current position for super-read s
  unsigned int get_weight(bool s);

  // HALF_TABLE stuff ...

  // set partitioning to the dual of partitioning (for the backtracking)
  void set_dual_partitioning(unsigned int partitioning);
  
};

#endif
