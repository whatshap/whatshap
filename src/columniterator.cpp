#include <cassert>
#include <unordered_set>

#include "columniterator.h"

using namespace std;

ColumnIterator::ColumnIterator(const ReadSet& set) : set(set) {
	assert(set.finalized);
	this->n = 0;
	this->next_read_index = 0;
}

ColumnIterator::~ColumnIterator() {}

unsigned int ColumnIterator::get_column_count() {
	return set.get_positions()->size();
}

unsigned int ColumnIterator::get_read_count() {
	return set.get_read_count();
}

const vector<unsigned int>* ColumnIterator::get_positions() {
	return set.get_positions();
}

bool ColumnIterator::has_next() {
	n < set.positions->size();
}

auto_ptr<vector<const Entry*> > ColumnIterator::get_next() {
	// check which of the current reads remain active
	// TODO
	// check which new reads might become active
	// TODO
	n += 1;
	//return result;
}
