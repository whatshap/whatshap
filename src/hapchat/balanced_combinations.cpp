#include "balanced_combinations.h"
#include <cmath>

using namespace std;


BalancedCombinations::BalancedCombinations() : generator() {}


void BalancedCombinations::initialize(const Counter n, const Counter k,
				      const BitColumn & col, const double r) {

  n_ = n;
  k_ = k;
  r_ = r;

  col_.reset();
  col_ |= col;

  c_ = (Counter)floor(n_ * r_);

  // pi_0 and pi_1
  p.clear();
  p.push_back(n_ - col_.count());
  p.push_back(col_.count());

  // build mapping for composing combinations and initialize arrays
  build_mapping();
  initialize_arrays();

  // initialize the counters
  t_ = 0;
  i_ = 0;
  j_ = 0;
  ii_ = 0;
  jj_ = 0;

  has_next_ = true;
  s_ = true; // prime the try loop
  try_next();
}


bool BalancedCombinations::has_next() {

  return has_next_;
}


void BalancedCombinations::next() {

  make_current();
  s_ = false;
  try_next();
}


void BalancedCombinations::get_combination(BitColumn & result) {

  result.reset();
  result |= current_;
}


// auxiliary (private) functions
/**********************************************************************/


void BalancedCombinations::build_mapping() {

  map.clear();
  map.resize(2);
  for(i_ = 0; i_ < n_; ++i_) {

    if(col_.test(i_))
      map[1].push_back(i_);
    else
      map[0].push_back(i_);
  }
}


void BalancedCombinations::initialize_arrays() {

  c.clear();

  // c[0][.]
  a.clear();
  a.resize(p[0]+1);
  c.push_back(a);

  // c[1][.]
  a.clear();
  a.resize(p[1]+1);
  c.push_back(a);
}


void BalancedCombinations::retrieve_c0() {

  if(c[0][i_].empty()) {

    generator.initialize(p[0], i_);
    while(generator.has_next()) {

      generator.next(); // should always be at least the empty comb
      generator.get_combination(comb);
      c[0][i_].push_back(comb);

    }
  }
}


void BalancedCombinations::retrieve_c1() {

  if(c[1][j_].empty()) {

    generator.initialize(p[1], j_);
    while(generator.has_next()) {

      generator.next(); // should always be at least the empty comb
      generator.get_combination(comb);
      c[1][j_].push_back(comb);

    }
  }
}


void BalancedCombinations::make_current() {

  current_.reset();

  // fill c0
  for(i = 0; i < p[0]; ++i)
    if(c[0][i_][ii_].test(i))
      current_.set(map[0][i]);

  // fill c1
  for(j = 0; j < p[1]; ++j)
    if(c[1][j_][jj_].test(j))
      current_.set(map[1][j]);
}


void BalancedCombinations::try_next() {

  // loop with switch, advancing exactly once each call to function
  while(t_ <= k_) {
    while(i_ <= min(p[0], t_)) {
      j_ = t_ - i_;

      // check if j_ is feasible
      if(j_ <= p[1]) {

	// check balance threshold
	if((p[0]-i_ + min(p[1],t_-i_) >= c_) and (p[1]-j_ + min(p[0],t_-j_) >= c_)) {

	  retrieve_c0(); // c[0][i_]
	  while(ii_ < c[0][i_].size()) {

	    retrieve_c1(); // c[1][j_]
	    while(jj_ < c[1][j_].size()) {

	      // at this point, jj_,ii_,j_,i_,t_ is a valid configuration
	      if(s_)
		return;

	      s_ = true;
	      ++jj_;
	    }
	    jj_ = 0;
	    ++ii_;
	  }
	  ii_ = 0;
	}
      }
      ++i_;
    }
    i_ = 0;
    ++t_;
  }

  has_next_ = false; // the end
}
