#include <sstream>
#include <stdexcept>
#include <algorithm>

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
	reads.push_back(read);
}

string ReadSet::toString() {
	ostringstream oss;
	oss << "ReadSet:" << endl;
	for (size_t i=0; i<reads.size(); ++i) {
		oss << "  " << reads[i]->toString() << endl;
	}
	return oss.str();
}

void ReadSet::finalize() {
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->sortVariants();
	}
	sort(reads.begin(), reads.end(), read_comparator_t());
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->setID(i);
	}
	finalized = true;
}
