/*
  reads and outputs, incrementally, the columns of a .wif file.
  ColumnReader will never output a column of height larger than
  'coverage_threshold', even if it reads in a column of height larger
  than this (in which case discards columns until height is within
  this threshold, according to some critieria) -- Murray Patterson,
  Sept 2013

  obacht : think of some sort of criteria.  For now, threshold is a
  cutoff (that is, it assumes that the preprocessing step took care of
  controlling the coverage) -- murray
*/

#ifndef COLUMN_READER_H
#define COLUMN_READER_H

#include "entry.h"
#include <string>
#include <fstream>
#include <vector>
#include <queue>
#include <list>
#include <memory>

class ColumnReader {

public:
  /** Constructor.
   *  @param remove_weights If true, all weights (i.e. phred scores) are removed, 
   *                        that is, set to one. */
  ColumnReader(std::string f, unsigned int c, bool remove_weights = false);
  ~ColumnReader();

  unsigned int num_cols(); // number of columns (snp positions)
  unsigned int num_rows(); // number of rows (so far)
  // true iff the reader has another column
  bool has_next();
  bool remove_weights;
  // get next column
  std::auto_ptr<std::vector<Entry *> > get_next();
  // return a const pointer to the positions
  const std::vector<unsigned int> * get_positions();

private:
  std::ifstream ifs;
  unsigned int coverage_threshold;

  std::vector<unsigned int> positions; // entry (snp) positions (indexed by column)
  typedef std::queue <Entry *> queue_t;
  typedef std::list<queue_t *> buffer_t;
  buffer_t buffer; // buffer on file

  unsigned int row; // current row (column) during the process
  unsigned int column;

  // auxiliary functions

  /*
    reads in .wif file, gets its entry (snp) positions and places
    them positions vector and returns true iff all of this was
    successful
  */
  bool compute_positions();

  /*
    returns true if buffer is not leftmost total (a buffer is
    'leftmost total' when all start positions of rows (row suffixes)
    in the buffer are the same (as current column), and does not
    read in this case, o.w. it reads line from .wif file into buffer
    and returns true if this causes the buffer to no longer be
    leftmost total, i.e., when the start position of this line is
    beyond the current column)

    precondition: the file stream has 'good' status
  */
  bool read_line();

};

std::ostream& operator<<(std::ostream& out, const std::vector<unsigned int>& v);

#endif
