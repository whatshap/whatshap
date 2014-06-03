#ifndef DICER_H
#define DICER_H

#include <vector>
#include <fstream>

class Dicer {

private:
  std::string filename;
  int index; // to keep track of what dice we are on

  std::ifstream ifs;
  std::ofstream ofs;

  std::string prefix; // filename prefix and extension
  std::string extension;

public:
  /** Constructor.
   * @ param filename The name of the file to dice up
   */
  Dicer(std::string filename);

  /** Destructor. */
  ~Dicer();

  /** Dice up the files, generating files with names filenamedXX where
   *  XX \in {1,N} where N is the number of dices the file ends up in
   */
  void dice();

  /* helper functions to extract first (last) snp position of a line (i.e., a line is a read) */
  static unsigned int get_first(const std::string & line);
  static unsigned int get_last(const std::string & line);

  /* helper functions to tokenize a filename into prefix and extension */
  static void prefix_extension_tok(const std::string & filename, std::string & prefix, std::string & extension);

};

#endif
