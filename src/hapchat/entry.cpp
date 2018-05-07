#include <cassert>
#include "entry.h"

std::ostream& operator<<(std::ostream& out, const Entry& e) {
  out << "Entry(" << e.read_id << ',';
  switch (e.allele_type) {
    case Entry::MAJOR_ALLELE:
      out << "MAJOR";
      break;
    case Entry::MINOR_ALLELE:
      out << "MINOR";
      break;
    case Entry::BLANK:
      out << "BLANK";
      break;
    case Entry::EQUAL_SCORES:
      out << "EQUAL_SCORES";
      break;
    default:
      assert(false);
  }
  out << ',' << ((int)e.phred_score) << ')';
  return out;
}
