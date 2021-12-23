#include <cassert>
#include <bits/stdc++.h>
using namespace std;

#include "entry.h"

Entry::Entry(unsigned int r, int m, std::vector<unsigned int> e, int q) {
	read_id = r;
	allele = m;
	set_emission_score(e);
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


void Entry::set_emission_score(std::vector<unsigned int> e) {
	emission_score.resize(e.size());
	unsigned int d_max = *max_element(e.begin(), e.end());
	unsigned int d_min = *min_element(e.begin(), e.end());
	float normalization = 0.0;
	int i = 0;
	for (auto it = e.begin(); it != e.end(); it++, i++) {
		emission_score[i] = pow(0.1, (*it - d_min)) * pow (0.9, (d_max - *it));
		normalization += emission_score[i];
	}
	transform((emission_score).begin(), (emission_score).end(), (emission_score).begin(), std::bind2nd(std::divides<float>(), normalization));
}

void Entry::set_quality(int q) {
	quality = q;
}

std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id ;
	out << ","<< e.allele << ",(";
	for (auto i : e.emission_score) {
		out << i << ",";
	}
	out << ")," << e.quality << ")" << std::endl;
	return out;
}
