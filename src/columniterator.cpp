#include <cassert>
#include <unordered_set>
#include <stdexcept>

#include "columniterator.h"

using namespace std;

ColumnIterator::ColumnIterator(const ReadSet& set) : set(set) {
	this->n = 0;
	this->next_read_index = 0;
	this->positions = set.get_positions();
}

ColumnIterator::~ColumnIterator() {
	for (size_t i=0; i<blank_entries.size(); ++i) {
		delete blank_entries[i];
	}
	blank_entries.clear();
	delete positions;
}

unsigned int ColumnIterator::get_column_count() {
	return positions->size();
}

unsigned int ColumnIterator::get_read_count() {
	return set.size();
}

const vector<unsigned int>* ColumnIterator::get_positions() {
	return positions;
}

bool ColumnIterator::has_next() {
	return n < positions->size();
}

unique_ptr<vector<const Entry*> > ColumnIterator::get_next() {
	// genomic position of the column to be returned
	int next_pos = positions->at(n);
	
	// check which of the current reads remain active
	list<active_read_t>::iterator list_it = active_reads.begin();
	while (list_it != active_reads.end()) {
		const Read* read = set.get(list_it->read_index);
		if (read->lastPosition() < next_pos) {
			list_it = active_reads.erase(list_it);
			continue;
		}
		while (read->getPosition(list_it->active_entry) < next_pos) {
			list_it->active_entry += 1;
			assert(list_it->active_entry < read->getVariantCount());
		}
		++list_it;
	}

	// check which new reads might become active
	while (next_read_index < set.size()) {
		int read_start = set.get(next_read_index)->firstPosition();
		if (read_start < next_pos) {
			throw std::runtime_error("ColumnIterator: reads in ReadSet are not sorted.");
		}
		if (read_start == next_pos) {
			if (!set.get(next_read_index)->isSorted()) {
				throw std::runtime_error("ColumnIterator: encountered read with unsorted variants.");
			}
			active_reads.push_back(active_read_t(next_read_index));
			next_read_index += 1;
		} else {
			break;
		}
	}
	
	// gather entries from active reads
	unique_ptr<vector<const Entry*> > result = unique_ptr<vector<const Entry*> >(new vector<const Entry*>());
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

// 	cerr << "ColumnIterator::get_next():";
// 	for (size_t i=0; i<result->size(); ++i){
// 		cerr << " " << (*(result->at(i)));
// 	}
// 	cerr << endl;
	n += 1;
	return result;
}
