#ifndef DP_TABLE_H
#define DP_TABLE_H

#include <vector>
#include <memory>

#include "columnindexingscheme.h"
#include "columniterator.h"
#include "entry.h"
#include "read.h"
#include "readset.h"

#ifdef FF_PARALLEL
#include <ff/parallel_for.hpp>

struct ff_task_t {
    ff_task_t(): 
        forward_index(-1),index(-1),cost(-1) {}
    ff_task_t(unsigned int forward_index, unsigned int index=-1, int cost=-1):
        forward_index(forward_index),index(index),cost(cost) {}

    unsigned int forward_index;
    unsigned int index;
    unsigned int cost;
};
#endif

class DPTable {
private:
  ReadSet* read_set;
  // vector of indexingschemes
  std::vector<ColumnIndexingScheme*> indexers;
  // optimal score and its index in the rightmost DP table column
  unsigned int optimal_score;
  unsigned int optimal_score_index;
  // backtrace_table[x][i] indicates the index in table x from which the
  // i-th entry in the forward projection of table x comes from
  std::vector<std::vector<unsigned int>* > backtrace_table;
  unsigned int read_count;
  // helper function to pull read ids out of read column
  std::shared_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);
  bool all_heterozygous;

  // helper function to compute the optimal path through the backtrace table
  std::shared_ptr<std::vector<unsigned int> > get_index_path();

  void compute_table();
    
#ifdef FF_PARALLEL
  ff::ParallelForPipeReduce<std::vector<ff_task_t>* > pf;
  int numthreads;
  long chunksize;
  unsigned int threshold;
#endif

public:
  /** Constructor.
   *  @param read_set DP table is constructed for the contained reads. Ownership is retained
   *                  by caller. Pointer must remain valid during the lifetime of this DPTable.
   *  @param all_heterozygous If true, then the "all heterozygous" assumption is made;
   *                          i.e., all positions are forced to be heterozygous even when
   *                          reads suggest a homozygous site. */
#ifdef FF_PARALLEL
  DPTable(ReadSet* read_set, bool all_heterozygous = false, int numthreads=-1, long chunksize=1024, unsigned int threshold=(1<<20));
#else
  DPTable(ReadSet* read_set, bool all_heterozygous = false);
#endif
    
  ~DPTable();

  unsigned int get_optimal_score();
  
  /** Computes optimal haplotypes and adds them (in the form of "super reads") to 
   *  the given read_set.
   */
  void get_super_reads(ReadSet* output_read_set);

  /** Performs a backtrace through the DP table and returns optimal partitioning of the reads.
   *  Pointer ownership is transferred to caller. */
  std::vector<bool>* get_optimal_partitioning();
};

#endif

