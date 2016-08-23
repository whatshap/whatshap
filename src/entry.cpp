#include <cassert>

#include "entry.h"

Entry::Entry(unsigned int r, allele_t m, unsigned int p) {
	read_id = r;
	allele_type = m;
	phred_score = p;
}


unsigned int Entry::get_read_id() const {
	return read_id;
}


Entry::allele_t Entry::get_allele_type() const {
	return allele_type;
}


unsigned int Entry::get_phred_score() const {
	return phred_score;
}

void Entry::set_read_id(unsigned int r) {
	read_id = r;
}


void Entry::set_allele_type(allele_t m) {
	allele_type = m;
}


void Entry::set_phred_score(unsigned int p) {
	phred_score = p;
}


std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id << ',';
	switch (e.allele_type) {
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
	out << ',' << ((int)e.phred_score) << ')';
}
