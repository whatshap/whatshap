#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "dicer.h"

using namespace std;

// MAIN
int main(int argc, char * const argv[]) {

  if(argc < 2) {
    cout << "usage : " << argv[0] << " input.wif" << endl;
    exit(0);
  }

  string inputfilename = argv[1];

  Dicer dicer(inputfilename);
  dicer.dice();

  return 0;
}
