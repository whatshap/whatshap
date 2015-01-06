#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <iomanip>

#include "readset.h"

using namespace std;

ReadSet::ReadSet() {
	this->finalized = false;
}

ReadSet::~ReadSet() {
	for (size_t i=0; i<reads.size(); ++i) {
		delete reads[i];
	}
}

void ReadSet::add(Read* read) {
	if (finalized) throw std::runtime_error("Cannot add to finalized ReadSet");
	if (read_name_map.find(read->getName()) != read_name_map.end()) {
		throw std::runtime_error("ReadSet::add: duplicate read name.");
	}
	reads.push_back(read);
	read_name_map[read->getName()] = reads.size() - 1;
}

string ReadSet::toString() {
	ostringstream oss;
	oss << "ReadSet:" << endl;
	for (size_t i=0; i<reads.size(); ++i) {
		oss << "  " << setw(5) << i << ' ' << reads[i]->toString() << endl;
	}
	return oss.str();
}

void ReadSet::sort() {
	// Sort the reads by position
	std::sort(reads.begin(), reads.end(), read_comparator_t());
	
	// Update read_name_map
	read_name_map.clear();
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->setID(i);
		read_name_map[reads[i]->getName()] = i;
	}
}

void ReadSet::finalize() {
	if (finalized) throw std::runtime_error("Cannot finalize a finalized ReadSet");

	// Sort the variants in the remaining reads
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->sortVariants();
	}

	this->sort();

	finalized = true;
}

bool ReadSet::isFinalized() const {
	return finalized;
}

vector<unsigned int>* ReadSet::get_positions() const {
	unordered_set<unsigned int> position_set;
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->addPositionsToSet(&position_set);
	}
	vector<unsigned int>* positions = new vector<unsigned int>(position_set.begin(), position_set.end());
	std::sort(positions->begin(), positions->end());
	return positions;
}

unsigned int ReadSet::size() const {
	return reads.size();
}

Read* ReadSet::get(int i) const {
	return reads[i];
}

Read* ReadSet::getByName(std::string name) const {
	read_name_map_t::const_iterator it = read_name_map.find(name);
	if (it == read_name_map.end()) {
		return 0;
	} else {
		return reads[it->second];
	}
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
