#ifndef SLICER_H
#define SLICER_H

#include "activelistdelegator.h"
#include <vector>
#include <fstream>

class Slicer {

private:
  std::string filename;
  unsigned int coverage;

  std::ifstream ifs;
  std::vector<std::ofstream *> ofss; // set of active output streams,
  ActiveListDelegator delegator; // corresponding to this delegator

  std::string prefix; // filename prefix and extension
  std::string extension;

public:
  /** Constructor.
   * @ param filename The name of the file to slice
   * @ param coverage At what coverage to slice
   */
  Slicer(std::string filename, unsigned int coverage);

  /** Destructor. */
  ~Slicer();

  /** Slice up the files, generating files with names filenamecXXsYY
   *  where XX = coverage and YY \in {1,N} where N is the number of
   *  slices the file ends up in
   */
  void slice();

};

#endif
