#include <sstream>

#include "readset.h"

using namespace std;

ReadSet::ReadSet() {
}

ReadSet::~ReadSet() {
	for (size_t i=0; i<reads.size(); ++i) {
		delete reads[i];
	}
}

void ReadSet::add(Read* read) {
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