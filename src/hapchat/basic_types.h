#ifndef _BASIC_TYPES_H_
#define _BASIC_TYPES_H_

#include <bitset>
#include <vector>
#include <iostream>
#include <limits>

#include "entry.h"

#define MAX_COVERAGE 64
#define MAX_CORRECTIONS 63

typedef unsigned int Counter;

#define MAX_COUNTER std::numeric_limits<Counter>::max()

typedef int Pointer;
typedef std::bitset<MAX_COVERAGE> BitColumn;
typedef std::vector<Entry> Column;
typedef std::vector<Column> Block;


struct EntryRead {
  
  int position;
  bool allele;
  unsigned int phred_score;

  EntryRead(int pos, bool a, unsigned int p) 
  {
    position = pos;
    allele = a;
    phred_score = p;
  }

};

typedef std::vector<EntryRead> Fragment;


struct options_t {
  bool options_initialized;
  std::string input_filename;
  std::string haplotype_filename;
  bool unweighted;
  bool no_xs;
  bool all_heterozygous;
  bool unique;
  double error_rate;
  double alpha;
  bool balancing;
  unsigned int balance_cov;
  double balance_ratio;

  options_t()
  : options_initialized(false),
    input_filename(""),
    haplotype_filename(""),
    unweighted(false),
    no_xs(false),
    all_heterozygous(false),
    unique(false),
    error_rate(0.05),
    alpha(0.01),
    balancing(false),
    balance_cov(20),
    balance_ratio(0.5)
  {}

};

std::ostream& operator<<(std::ostream& out, const options_t& options);


struct constants_t
{

  BitColumn zeroes;
  BitColumn ones;

  constants_t() {
    ones.flip();
  };

};


struct Backtrace1
{
  Pointer jump;
  Counter index;

  Backtrace1()
    : jump(-1), index(0)
  {};
};



// A type for representing costs
class Cost {
public:
  typedef unsigned int cost_t;

private:
  static const cost_t infinity_;

  cost_t cost_;

  static bool is_addition_unsafe(const Cost& c1_, const Cost& c2_) {
    return (c1_.cost_ > (infinity_ - c2_.cost_));
  }

public:

  static const Cost INFTY;

  Cost(const cost_t cost= 0)
    :cost_(cost)
  {}

  Cost(const Cost& c)
    :cost_(c.cost_)
  {}

  Cost& operator=(const Cost& c) {
    cost_= c.cost_;
    return *this;
  }

  Cost& operator+=(const Cost& c) {
    if (is_addition_unsafe(*this, c))
      cost_ = infinity_;
    else
      cost_ += c.cost_;
    return *this;
  }

  Cost operator+(const Cost& c) const {
    if (is_addition_unsafe(*this, c))
      return INFTY;
    return Cost(cost_ + c.cost_);
  }

  bool operator<(const Cost& c) const {
    return cost_ < c.cost_;
  }

  bool operator<=(const Cost& c) const {
    return cost_ <= c.cost_;
  }

  bool operator>(const Cost& c) const {
    return cost_ > c.cost_;
  }

  bool operator>=(const Cost& c) const {
    return cost_ >= c.cost_;
  }

  bool operator==(const Cost& c) const {
    return cost_ == c.cost_;
  }

  friend std::ostream& operator<<(std::ostream& out, const Cost& c);
};

// Pretty-print costs
std::ostream& operator<<(std::ostream& out, const Cost& c);

// Pretty-print binary vectors
std::ostream& operator<<(std::ostream& out, const std::vector<bool>& v);

// Pretty-print char vectors
std::ostream& operator<<(std::ostream& out, const std::vector<char>& v);


options_t parse_arguments(int argc, char** argv);


#endif // _BASIC_TYPES_H_
