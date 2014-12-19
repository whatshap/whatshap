#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <iomanip>

#include "readset.h"

using namespace std;

ReadSet::ReadSet() {
	this->finalized = false;
	this->positions = 0;
}

ReadSet::~ReadSet() {
	for (size_t i=0; i<reads.size(); ++i) {
		delete reads[i];
	}
	if (positions != 0) delete positions;
}

void ReadSet::add(Read* read) {
	if (finalized) throw std::runtime_error("Cannot add to finalized ReadSet");
	reads.push_back(read);
}

string ReadSet::toString() {
	ostringstream oss;
	oss << "ReadSet:" << endl;
	for (size_t i=0; i<reads.size(); ++i) {
		oss << "  " << setw(5) << i << ' ' << reads[i]->toString() << endl;
	}
	return oss.str();
}

void ReadSet::finalize() {
	if (finalized) throw std::runtime_error("Cannot finalize a finalized ReadSet");
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->sortVariants();
	}
	sort(reads.begin(), reads.end(), read_comparator_t());
	unordered_set<unsigned int> position_set;
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->setID(i);
		reads[i]->addPositionsToSet(&position_set);
	}
	positions = new vector<unsigned int>(position_set.begin(), position_set.end());
	sort(positions->begin(), positions->end());
	finalized = true;
}

bool ReadSet::isFinalized() const {
	return finalized;
}

const vector<unsigned int>* ReadSet::get_positions() const {
	if (!finalized) throw std::runtime_error("ReadSet::get_positions: can only be called after finalization");
	return positions;
}

unsigned int ReadSet::size() const {
	return reads.size();
}

Read* ReadSet::get(int i) const {
	return reads[i];
}

ReadSet* ReadSet::subset(const IndexSet* indices) const {
	ReadSet* result = new ReadSet();
	IndexSet::const_iterator it = indices->begin();
	for (; it != indices->end(); ++it) {
		result->add(new Read(*(reads[*it])));
	}
	result->finalize();
	return result;
}
