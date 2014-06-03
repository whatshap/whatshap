#include <cassert>
#include <stdlib.h>
#include <sstream>
#include "dicer.h" // for its helper functions
#include "slicer.h"

using namespace std;

Slicer::Slicer(string filename, unsigned int coverage) : delegator(coverage) {
  this->filename = filename;
  this->coverage = coverage;

  assert(!filename.empty());
  try { ifs.open(filename.c_str(), ios::in); }
  catch(exception & e) { throw e; }

  Dicer::prefix_extension_tok(filename, prefix, extension); // get prefix and extension
}

Slicer::~Slicer() {

  ifs.close(); // close all file streams
  for(size_t i=0; i< ofss.size(); ++i) ofss[i]->close();
}

void Slicer::slice() {

  string line;
  while(getline(ifs,line)) {

    unsigned int first = Dicer::get_first(line);
    delegator.set_active_element(first);
    unsigned int index = delegator.delegate();
    assert(index <= ofss.size()); // sanity check

    // need to open a new output stream for the n+1th slice
    if(index == ofss.size()) {

      ostringstream oss;
      oss << prefix << "c" << coverage << "s" << (index+1) << extension;
      string slice_filename = oss.str();

      ofss.push_back(new ofstream());
      try { ofss[index]->open(slice_filename.c_str(), ios::out); }
      catch(exception & e) { throw e; }
    }
    
    *ofss[index] << line << endl; // write to "line delegated to" outputstream
    unsigned int last = Dicer::get_last(line);
    delegator.push(last); // and push "line" to active list to which it was delegated
  }
}
