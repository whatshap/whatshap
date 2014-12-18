#include <cassert>
#include <stdexcept>
#include <sstream>

#include "indexset.h"

using namespace std;

IndexSet::IndexSet() {}

IndexSet::~IndexSet() {}

bool IndexSet::contains(size_t index) const {
	return set.find(index) != set.end();
}

void IndexSet::add(size_t index) {
	set.insert(index);
}

size_t IndexSet::size() const {
	return set.size();
}

string IndexSet::toString() {
	ostringstream oss;
	oss << '{';
	bool is_first = true;
	std::set<int>::const_iterator it = set.begin();
	for (; it != set.end(); ++it) {
		if (is_first) {
			oss << (*it);
			is_first = false;
		} else {
			oss << ',' << (*it);
		}
	}
	oss << '}';
	return oss.str();
}
