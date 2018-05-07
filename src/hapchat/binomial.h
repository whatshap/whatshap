#ifndef BINOMIAL_H
#define BINOMIAL_H

#include <vector>
#include <bitset>
#include <string.h>

#include "basic_types.h"

class binom_coeff {

 private:
  static std::vector<std::vector<unsigned int> > btable;
  static std::vector<std::vector<unsigned int> > ctable;

 public:
  static void
    initialize_binomial_coefficients(const unsigned int n, const unsigned int k);

//Note: if k > n the result is equal to 0
  static unsigned int
    binomial_coefficient(const unsigned int n, const unsigned int k) {
    return btable[n][k];
  }

//Note: if k > n then k is considered equal to n
  static unsigned int
    cumulative_binomial_coefficient(const unsigned int n, const unsigned int k) {
    return ctable[n][k];
  }

  // index of comb in the binomial coefficients
  static unsigned int indexof(std::bitset<MAX_COVERAGE> comb);

  // index of comb in the cumulative binomial coeffients
  static unsigned int cumulative_indexof(std::bitset<MAX_COVERAGE> comb, const unsigned int n_elements);

};


#endif
