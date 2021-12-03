#include <cassert>


#include "entry.h"

Entry::Entry(unsigned int r, int m, std::vector<float> e, int q) {
	read_id = r;
	allele = m;
	emission_score = e;
	quality = q;
}

Entry::Entry(unsigned int r, int m) {
	read_id = r;
	allele = m;
}

unsigned int Entry::get_read_id() const {
	return read_id;
}


int Entry::get_allele_type() const {
	return allele;
}


std::vector<float> Entry::get_emission_score() const {
	return emission_score;
}

int Entry::get_quality() const {
	return quality;
}

void Entry::set_read_id(unsigned int r) {
	read_id = r;
}


void Entry::set_allele_type(int m) {
	allele = m;
}


void Entry::set_emission_score(std::vector<float> e) {
	emission_score = e;
}

void Entry::set_quality(int q) {
	quality = q;
}

std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id ;
	out << ","<< e.allele << ",";
	return out;
}
