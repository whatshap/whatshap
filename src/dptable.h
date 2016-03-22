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

typedef struct index_and_inheritance_t {
  unsigned int index;
  unsigned int inheritance_value;
} index_and_inheritance_t;

typedef std::vector<unsigned int> num_of_recomb_uints_t;

class DPTable {
private:
  ReadSet* read_set;
  
  const std::vector<unsigned int> read_marks;
  const std::vector<unsigned int> recombcost;
  const vector<array<vector<unsigned int> , 3>> genotypes;
  const std::vector<std::vector<int>> triples;
  const unsigned int number_of_ind;
  // vector of indexingschemes
  std::vector<ColumnIndexingScheme*> indexers;
  // optimal score and its index in the rightmost DP table column
  unsigned int optimal_score;
  unsigned int optimal_score_index;
  unsigned int optimal_transmission_value;
  // transmission value preceeding the optimal one (in the column before / in the backtrace)
  unsigned int previous_transmission_value;
  // index_backtrace_table[c][i][t] indicates the index (=bipartition) in column c from which the
  // i-th entry in the FORWARD projection of column c comes from, assuming a transmission value of t
  std::vector<std::vector<num_of_recomb_uints_t>* > index_backtrace_table;
  // let x := index_backtrace_table[c][i][t] and dp[x][t] the corresponding DP entry
  // and j be the BACKWARD projection of x.
  // Then t' = transmission_backtrace_table[c][i][t] is the transmission index (from {0,1,2,3})
  // that gave rise to dp[x][t].
  std::vector<std::vector<num_of_recomb_uints_t>* > transmission_backtrace_table;
  unsigned int read_count;
  // helper function to pull read ids out of read column
  std::unique_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);
 //Read* r0m;
 //Read* r1m;
 //Read* r0f;
 //Read* r1f;
 //Read* r0c;
 //Read* r1c;
  std::vector<unsigned int> id_of_individuals;
  // helper function to compute the optimal path through the backtrace table
  std::unique_ptr<std::vector<index_and_inheritance_t> > get_index_path();

  void compute_table();

public:
  /** Constructor.
   *  @param read_set DP table is constructed for the contained reads. Ownership is retained
   *                  by caller. Pointer must remain valid during the lifetime of this DPTable.
   *  @param all_heterozygous If true, then the "all heterozygous" assumption is made;
   *                          i.e., all positions are forced to be heterozygous even when
   *                          reads suggest a homozygous site. */
  DPTable(ReadSet* read_set, vector<unsigned int> read_marks, vector<unsigned int> recombcost, vector<array<vector<unsigned int> , 3>> genotypes, std::vector<std::vector<int>> triples, unsigned int number_of_ind);
 
  ~DPTable();

  unsigned int get_optimal_score();
  
  /** Computes optimal haplotypes and adds them (in the form of "super reads") to 
   *  the given read_set.
   */
  void get_super_reads(std::vector<ReadSet*>* output_read_set, vector<unsigned int>* transmission_vector);

  /** Performs a backtrace through the DP table and returns optimal partitioning of the reads.
   *  Pointer ownership is transferred to caller. */
  std::vector<bool>* get_optimal_partitioning();
};

#endif

