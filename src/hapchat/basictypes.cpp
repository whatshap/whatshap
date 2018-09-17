/*

  Copyright (C) 2017-2018 Marco Dell'Acqua, Yuri Pirola, Simone
  Zaccaria

  Distributed under the MIT license.

  You should have received a copy of the MIT license along with this
  program.

*/

#include "basictypes.h"
#include <limits>

const Cost::cost_t Cost::infinity_(std::numeric_limits<Cost::cost_t>::max());
const Cost Cost::INFTY(Cost::infinity_);


std::ostream& operator<<(std::ostream& out, const Cost& c) {
  if (c == Cost::INFTY)
    out << "INFINITY";
  else
    out << c.cost_;
  return out;
}


std::ostream& operator<<(std::ostream& out, const std::vector<bool>& v) {
  for (std::vector<bool>::const_iterator it= v.begin(); it != v.end(); ++it) {
    out << (*it? '1' : '0');
  }
  return out;
}


std::ostream& operator<<(std::ostream& out, const std::vector<char>& v) {
  for (std::vector<char>::const_iterator it= v.begin(); it != v.end(); ++it) {
    out << *it;
  }
  return out;
}
