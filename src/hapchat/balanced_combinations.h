#ifndef BALANCED_COMBINATIONS_H
#define BALANCED_COMBINATIONS_H

#include "basic_types.h"
#include "combinations.h"

class BalancedCombinations {

 public :

  typedef std::vector<std::vector<BitColumn> > Array;
  typedef std::vector<Counter> Mapping;

  // constructor
  BalancedCombinations();

  // initialize a generator
  void initialize(const Counter n, const Counter k,
		  const BitColumn & col, const double r);

  bool has_next();
  void next();

  void get_combination(BitColumn & result);

 private :

  // initial arguments
  Counter n_, k_;
  BitColumn col_; // column we are correcting
  double r_; // the ratio for computing c (below)

  Counter c_; // min support on a side: c = ceil(n*r)
  std::vector<Counter> p; // pi_0 and pi_1
  std::vector<Mapping> map; // map corrections to the right places in col_
  std::vector<Array> c; // C_B0 and C_B1

  // global configuration of the generator
  Counter t_; // 0 <= t <= k
  Counter i_, j_; // C_B0^i and C_B1^j
  Counter ii_, jj_; // which element of C_B0^i (resp. C_B1^j) we are on
  bool has_next_; // whether generator has a next

  // auxiliary (private) functions (and associated variables)
 private :

  // build mapping to col_
  void build_mapping();

  // initialize arrays C_B0, C_B1
  Array a;
  void initialize_arrays();

  // retrieve a C_B0, C_B1 array (building it should it be empty)
  Combinations generator;
  BitColumn comb;
  void retrieve_c0();
  void retrieve_c1();

  // compose the current combination
  Counter i,j;
  BitColumn current_;
  void make_current();

  // try to get the next one (may not exist)
  bool s_;
  void try_next();

};

#endif // BALANCED_COMBINATIONS_H
