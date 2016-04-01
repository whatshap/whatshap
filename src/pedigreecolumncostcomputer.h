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
#include "columnindexingiterator.h"


  
class PedigreeColumnCostComputer {
private:
  const std::vector<const Entry*>& column;
  const std::vector<unsigned int>& read_marks;
  
  unsigned int inheritance_val;
  unsigned int partitioning;
 // unsigned int num_of_roots;
  std::vector<unsigned int> roots;
  unsigned int num_of_triples;
 // std::vector<unsigned int> triple;
  std::vector<Pedigree::triple_entry_t> triples;
  
  // cost_partitionX[Y] is the cost of flipping all entries in the partition X to Y for Y = 0,1.

  //  std::array<unsigned int,2> null_array{{0,0}};
  std::vector<std::array<unsigned int, 2>> cost_partition;
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> haps; 
  
public:
//    struct ind_alleles_t {
//     std::pair<Entry::allele_t,Entry::allele_t> ind;
//     ind_alleles_t(Entry::allele_t ind1, Entry::allele_t ind2) : ind(std::make_pair(ind1, ind2)) {}
//   };
//   
  PedigreeColumnCostComputer(const std::vector<const Entry*>& column, const std::vector<unsigned int>& read_marks, unsigned int inheritance_val, std::vector<Pedigree::triple_entry_t> triples);
  
  void set_partitioning(unsigned int partitioning);

  void update_partitioning(int bit_to_flip);
  
  std::vector<unsigned int> compute_roots(std::vector<Pedigree::triple_entry_t> triples);
  
  unsigned int get_cost(std::vector<Pedigree::triple_entry_t>& genotypes);
  
  /** Returns the six alleles for all three individuals at the current position.
   *  Alleles are ensured to agree with the given genotypes.
   */
  std::map<unsigned int, std::pair<Entry::allele_t,Entry::allele_t>> get_alleles(std::vector<Pedigree::triple_entry_t>& genotypes);

  // returns the weight (what will be the phred_score) at the current position for super-read s
  unsigned int get_weight(bool second_haplotype);

  // HALF_TABLE stuff ...

  // set partitioning to the dual of partitioning (for the backtracking)
  void set_dual_partitioning(unsigned int partitioning);
  
};

#endif
