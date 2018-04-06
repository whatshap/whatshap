#include <cassert>
#include <limits>
#include <unordered_set>
#include <stdexcept>

#include "backwardcolumniterator.h"

using namespace std;

BackwardColumnIterator::BackwardColumnIterator(const ReadSet& set, const std::vector<unsigned int>* positions) : set(set) {
	if (positions == nullptr) {
		this->positions = set.get_positions();
	} else {
		this->positions = new vector<unsigned int>(positions->begin(), positions->end());
	}

    this->n = (int)this->positions->size()-1;

    if(this->n < 0) return;

	// create a mapping of genomic positions to column indices
	std::unordered_map<unsigned int, size_t> position_map;
	for (size_t i=0; i<this->positions->size(); ++i) {
		position_map[this->positions->at(i)] = i;
	}
	// precompute first_reads
	first_reads.assign(this->positions->size(),  numeric_limits<size_t>::max());
	int pos = 0;
	for (size_t i=0; i<set.size(); ++i) {
		const Read* read = set.get(i);
		if (read->firstPosition() < pos) {
			throw std::runtime_error("ColumnIterator: reads in ReadSet are not sorted.");
		}
		if (!read->isSorted()) {
			throw std::runtime_error("ColumnIterator: encountered read with unsorted variants.");
		}
		auto first_column_it = position_map.find(read->firstPosition());
		auto last_column_it = position_map.find(read->lastPosition());
		assert(first_column_it != position_map.end());
		assert(last_column_it != position_map.end());
		assert(first_column_it->second < last_column_it->second);
		assert(last_column_it->second < this->positions->size());
		for (size_t j=first_column_it->second; j<=last_column_it->second; ++j) {
			if (first_reads[j] == numeric_limits<size_t>::max()) {
				first_reads[j] = i;
			}
		}
		pos = read->firstPosition();
	}
	// For positions not covered by any read, fill in the index of next read that will
	// become active in subsequent columns
	if (first_reads.size() >= 2) {
		size_t next_index = first_reads[first_reads.size()-1];
		for (int i=first_reads.size()-2; i>=0; --i) {
			if (first_reads[i] == numeric_limits<size_t>::max()) {
				first_reads[i] = next_index;
			} else {
				next_index = first_reads[i];
			}
		}
	}
    jump_to_column(this->n);
}


BackwardColumnIterator::~BackwardColumnIterator() {
	for (size_t i=0; i<blank_entries.size(); ++i) {
		delete blank_entries[i];
	}
	blank_entries.clear();
	delete positions;
}


unsigned int BackwardColumnIterator::get_column_count() {
	return positions->size();
}


unsigned int BackwardColumnIterator::get_read_count() {
	return set.size();
}


const vector<unsigned int>* BackwardColumnIterator::get_positions() {
	return positions;
}


bool BackwardColumnIterator::has_next() {
    return n >= 0;
}

unique_ptr<vector<const Entry*> > BackwardColumnIterator::get_next() {
	// genomic position of the column to be returned
	int next_pos = positions->at(n);
    jump_to_column(n);

	// gather entries from active reads
    list<active_read_t>::iterator list_it = active_reads.begin();
	unique_ptr<vector<const Entry*> > result(new vector<const Entry*>());
	for (list_it = active_reads.begin(); list_it != active_reads.end(); ++list_it) {
		const Read* read = set.get(list_it->read_index);
		// Does read cover the current position?
		if (read->getPosition(list_it->active_entry) == next_pos) {
			result->push_back(read->getEntry(list_it->active_entry));
		} else {
			// if not, generate a blank entry
			Entry* e = new Entry(read->getID(), Entry::BLANK, 0);
			blank_entries.push_back(e);
			result->push_back(e);
		}
	}

    n -= 1;
	return result;
}


void BackwardColumnIterator::jump_to_column(int k) {
	assert(k < positions->size());
    //if(n == k) return;

	active_reads.clear();
	n = k;
    size_t next_read_index = first_reads[k];
	int pos = positions->at(k);

	// determine set of active reads
	while (next_read_index < set.size()) {
		const Read* read = set.get(next_read_index);
		if (read->lastPosition() < pos) {
			next_read_index += 1;
			continue;
		}
		if (read->firstPosition() <= pos) {
			size_t active_entry = 0;
			while (read->getPosition(active_entry) < pos) {
				active_entry += 1;
				assert(active_entry < read->getVariantCount());
			}
			active_reads.push_back(active_read_t(next_read_index, active_entry));
			next_read_index += 1;
		} else {
			break;
		}
	}
}
