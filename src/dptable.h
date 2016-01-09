#ifndef DP_TABLE_H
#define DP_TABLE_H

#include <array>
#include <vector>
#include <memory>

#include "columnindexingscheme.h"
#include "columniterator.h"
#include "entry.h"
#include "read.h"
#include "readset.h"

class DPTable {
private:
  ReadSet* read_set;
  
  const std::vector<unsigned int> read_marks;
  const std::vector<unsigned int> recombcost;
  const std::vector<unsigned int> genotypesm;
  const std::vector<unsigned int> genotypesf;
  const std::vector<unsigned int> genotypesc;
  // vector of indexingschemes
  std::vector<ColumnIndexingScheme*> indexers;
  // optimal score and its index in the rightmost DP table column
  unsigned int optimal_score;
  unsigned int optimal_score_index;
  unsigned int optimal_score_array_index;
  // backtrace_table[x][i] indicates the index in table x from which the
  // i-th entry in the forward projection of table x comes from
  std::vector<std::vector<std::array<unsigned int, 4>>* > backtrace_table;
  std::vector<std::vector<std::array<unsigned int, 4>>* > forrecomb;
  unsigned int read_count;
  // helper function to pull read ids out of read column
  std::unique_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);

  // helper function to compute the optimal path through the backtrace table
  std::unique_ptr<std::vector<unsigned int> > get_index_path();

  void compute_table();

public:
  /** Constructor.
   *  @param read_set DP table is constructed for the contained reads. Ownership is retained
   *                  by caller. Pointer must remain valid during the lifetime of this DPTable.
   *  @param all_heterozygous If true, then the "all heterozygous" assumption is made;
   *                          i.e., all positions are forced to be heterozygous even when
   *                          reads suggest a homozygous site. */
  DPTable(ReadSet* read_set, std::vector<unsigned int> read_marks, std::vector<unsigned int> recombcost, std::vector<unsigned int> genotypesm, std::vector<unsigned int> genotypesf, std::vector<unsigned int> genotypesc);
 
  ~DPTable();

  unsigned int get_optimal_score();
  
  /** Computes optimal haplotypes and adds them (in the form of "super reads") to 
   *  the given read_set.
   */
  void get_super_readsm(ReadSet* output_read_set);
  void get_super_readsf(ReadSet* output_read_set);
  void get_super_readsc(ReadSet* output_read_set);

  /** Performs a backtrace through the DP table and returns optimal partitioning of the reads.
   *  Pointer ownership is transferred to caller. */
  std::vector<bool>* get_optimal_partitioning();
};

#endif

