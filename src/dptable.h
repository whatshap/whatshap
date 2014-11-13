#ifndef DP_TABLE_H
#define DP_TABLE_H

#include <vector>
#include <memory>

#include "columnindexingscheme.h"
#include "columniterator.h"
#include "entry.h"

class DPTable {
private:
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
  std::auto_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);
  bool all_heterozygous;

public:
  /** Constructor. 
   *  @param all_heterozygous If true, then the "all heterozygous" assumption is made;
   *                          i.e., all positions are forced to be heterozygous even when
   *                          reads suggest a homozygous site. */
  DPTable(bool all_heterozygous = false);
  
  void compute_table(ColumnIterator * column_iterator);
  
  unsigned int get_optimal_score();
  
  // helper function to compute the optimal path through the backtrace table
  std::auto_ptr<std::vector<unsigned int> > get_index_path();

  // returns the haplotype for the reads (needs the input file again)
  typedef std::vector<std::vector<Entry::allele_t> > haplotype_t;
  std::auto_ptr<haplotype_t> get_haplotype(ColumnIterator * column_iterator);
  std::auto_ptr<haplotype_t> get_raw_haplotype(ColumnIterator * column_iterator);

  // returns the two super-reads for the reads (see above)
  typedef std::vector<std::vector<Entry *> > read_t;
  std::auto_ptr<read_t> get_super_reads(ColumnIterator * column_iterator);

  // returns optimal partitioning of the reads
  std::auto_ptr<std::vector<bool> > get_optimal_partitioning();
  
  // HALF_TABLE stuff ... when the half-table idea pans out

  // helper function to compute the optimal path through the backtrace table
  std::auto_ptr<std::vector<unsigned int> > get_index_path_HALF_TABLE();

  // returns the haplotype for the reads (needs the file again)
  std::auto_ptr<haplotype_t> get_haplotype_HALF_TABLE(ColumnIterator & column_iterator);

};

#endif

