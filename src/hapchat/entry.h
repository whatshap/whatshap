#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>

typedef long int readid_t;

#define SENTINEL_READID (-1)

class Entry {

public:

  typedef enum { MAJOR_ALLELE = 0, MINOR_ALLELE = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;

  Entry(readid_t r, allele_t m, unsigned int p)
    :read_id(r), phred_score(p), allele_type(m), gap(false)
  {
    if(m == Entry::BLANK) {
      allele_type = Entry::MAJOR_ALLELE;
      gap = true;
      phred_score = 0;
    }
  }

  readid_t get_read_id() const { return read_id; }
  allele_t get_allele_type() const { return allele_type; }
  unsigned int get_phred_score() const { return phred_score; }
  bool is_gap() const {return gap;}

  void set_read_id(readid_t r) { read_id = r; }
  void set_allele_type(allele_t m) {
    if(m == Entry::BLANK) {
      allele_type = Entry::MAJOR_ALLELE;
      gap = true;
    } else {
      allele_type = m;
      gap = false;
    }
  }
  void set_phred_score(unsigned int p) { phred_score = p; }
  void set_gap(bool g) {gap = g;}

  friend std::ostream& operator<<(std::ostream& out, const Entry& e);

private:
  readid_t read_id;
  unsigned int phred_score;
  allele_t allele_type;
  bool gap;
};

#endif
