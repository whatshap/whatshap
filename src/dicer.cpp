#include <cassert>
#include <stdlib.h>
#include <sstream>
#include <string>
#include "dicer.h"

using namespace std;

Dicer::Dicer(string filename) {
  this->filename = filename;
  this->index = -1;

  assert(!filename.empty());
  try { ifs.open(filename.c_str(), ios::in); }
  catch(exception & e) { throw e; }

  prefix_extension_tok(filename, prefix, extension); // get prefix (extension)
}

Dicer::~Dicer() {

  ifs.close(); // close all file streams
  ofs.close();
}

/*
  To "dice" a file is to cut it vertically (the slicer cut
  horizontally), i.e., into its connected components (in the graph
  sense).  We need this because the dp assumes the input is connected.

  This is implemented in a separate class because later we may add
  functionality for more refined dicing, such as into (a) bi-connected
  components, or (b) dices separated by intrisically heterozygous
  columns
*/
void Dicer::dice() {

  // the rightmost position in what we have read so far
  unsigned int end = 0; // detault -- NOTE : this assumes that positions start at 1

  string line;
  while(getline(ifs,line)) {

    unsigned int first = get_first(line);
    if(first > end) { // next read is disjoint from what we have read so far
      if(index >= 0) { // taper off current output stream
	ofs.close();
      }

      // now we open an output stream for the next (or first) dice
      ++index;
      ostringstream oss;
      oss << prefix << "d" << (index+1) << extension;
      string dice_filename = oss.str();

      try { ofs.open(dice_filename.c_str(), ios::out); }
      catch(exception & e) { throw e; }
    }

    ofs << line << endl; // write to the current output stream
    unsigned int last = get_last(line);    
    if(last > end) end = last; // update the end
  }
}

unsigned int Dicer::get_first(const string & line) {

  unsigned int first;
  istringstream iss(line);
  iss >> first;

  return first;
}

unsigned int Dicer::get_last(const string & line) {

  unsigned int last = 0; // default
  istringstream iss(line);
  string s; // temporary input holder

  while(1) {
    iss >> s; // either a snp or "#"
    if(s == "#") return last;
    last = atoi(s.c_str()); // else, cast it to an unsigned int
    iss >> s; // discard nucleotide,
    iss >> s; // allele
    iss >> s; // phred score, and
    iss >> s; // ":"
  }
}

void Dicer::prefix_extension_tok(const string & filename, string & prefix, string & extension) {

  unsigned int last = filename.find_last_of("/.");
  prefix = filename.substr(0,last);
  extension = filename.substr(last);
}
