#include <cassert>


#include "entry.h"

Entry::Entry(unsigned int r, int m, std::vector<unsigned int> p) {
	read_id = r;
	allele = m;
	phred_score = p;
}


unsigned int Entry::get_read_id() const {
	return read_id;
}


int Entry::get_allele_type() const {
	return allele;
}


std::vector<unsigned int> Entry::get_phred_score() const {
	return phred_score;
}

void Entry::set_read_id(unsigned int r) {
	read_id = r;
}


void Entry::set_allele_type(int m) {
	allele = m;
}


void Entry::set_phred_score(std::vector<unsigned int> p) {
	phred_score = p;
}


std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id ;
	out << ","<< e.allele << ",";
	for (int i=0; i< e.phred_score.size(); i++){
		out << e.phred_score.at(i);
	}
	out << ")";
}
