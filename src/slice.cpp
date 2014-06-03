#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "slicer.h"

// effectively the max coverage we can handle (due to the size of a byte)
#define MAX_COVERAGE 32

using namespace std;

// MAIN
int main(int argc, char * const argv[]) {

  if(argc < 3) {
    cout << "usage : " << argv[0] << " input.wif coverage" << endl;
    exit(0);
  }

  string inputfilename = argv[1];
  unsigned int coverage = atoi(argv[2]);
  assert(coverage <= MAX_COVERAGE);

  Slicer slicer(inputfilename, coverage);
  slicer.slice();

  return 0;
}
