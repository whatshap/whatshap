#include <cassert>

#include "entry.h"

using namespace std;

Entry::Entry(unsigned int r, vector<allele_t> m, vector<int> p) {
	read_id = r;
	allele_type = m;
	phred_score = p;
}


unsigned int Entry::get_read_id() const {
	return read_id;
}


vector<Entry::allele_t> Entry::get_allele_type() const {
	return allele_type;
}


vector<int> Entry::get_phred_score() const {
	return phred_score;
}

void Entry::set_read_id(unsigned int r) {
	read_id = r;
}


void Entry::set_allele_type(vector<allele_t> m) {
	allele_type = m;
}


void Entry::set_phred_score(vector<int> p) {
	phred_score = p;
}


std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id << ", [";
	for (unsigned int i = 0; i < e.allele_type.size(); i++){
		if (i > 0) out << ",";
		switch (e.allele_type[i]) {
			case Entry::REF_ALLELE:
				out << "REF";
				break;
			case Entry::ALT_ALLELE:
				out << "ALT";
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
	}
	out << "],[";
	for (unsigned int i = 0; i < e.phred_score.size(); i++){
		if (i > 0) out << ",";
		out << e.phred_score[i];
	}
	out << "])";
}
