/*

  Copyright (C) 2015-2018 Yuri Pirola, Simone Zaccaria

  Distributed under the MIT license.

  You should have received a copy of the MIT license along with this
  program.

*/

#include "binomialcoefficient.h"

std::vector<std::vector<unsigned int> > BinomialCoefficient::btable;
std::vector<std::vector<unsigned int> > BinomialCoefficient::ctable;

void
BinomialCoefficient::initialize_binomial_coefficients(const unsigned int n,
                                              const unsigned int k) {
  // binomial coefficients
  btable.clear();
  btable.resize(n+1, std::vector<unsigned int>(n + 1, 0));
  for (unsigned int i = 0; i <= n; ++i) {
    for (unsigned int j = 0; j <= i; j++) {
      if (i == 0 || j == 0 || j == i) {
        btable[i][j] = 1;
      } else {
        btable[i][j] = btable[i - 1][j - 1] + btable[i - 1][j];
      }
    }
  }
  // cumulative binomial coefficients
  ctable.clear();
  ctable.resize(n+1, std::vector<unsigned int>(n + 1, 0));
  for (unsigned int i = 0; i <= n; i++) {
    for (unsigned int j = 0; j <= k; j++) {
      for(unsigned int x = 0; x <= j; x++) {
        ctable[i][j] += btable[i][x];
      }
    }
  }
}


unsigned int BinomialCoefficient::indexof(std::bitset<MAX_COVERAGE> comb) {

  int k = 0;
  int c_k = 0;
  int temp = 0;
  int result = 0;

  while(comb.any()) {
    temp = ffsl(comb.to_ulong());
    c_k += temp;
    k++;
    result += BinomialCoefficient::binomial_coefficient(c_k - 1, k);
    comb >>= (temp);
  }

  return result;
}


unsigned int BinomialCoefficient::cumulative_indexof(std::bitset<MAX_COVERAGE> comb, const unsigned int n_elements) {

  unsigned int k = comb.count();
  unsigned int result = indexof(comb);

  for(unsigned int i = 0; i < k; i++) {
    result += BinomialCoefficient::binomial_coefficient(n_elements, i);
  }

  return result;
}
