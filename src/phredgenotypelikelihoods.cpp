#include <sstream>
#include <cassert>
#include <algorithm>
#include "phredgenotypelikelihoods.h"

using namespace std;

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods(unsigned int ploidy, unsigned int n_alleles, const vector<double>& gl) :ploidy(ploidy), n_alleles(n_alleles), gl(gl) 
{}

double PhredGenotypeLikelihoods::get(Genotype genotype) const {
	unsigned int index = genotype.get_index(ploidy, n_alleles);
	assert(index < gl.size());
	return this->gl[index];
}

const vector<double>& PhredGenotypeLikelihoods::as_vector() const {
	return gl;
}

size_t PhredGenotypeLikelihoods::genotype_count() const {
	return gl.size();
}

std::string PhredGenotypeLikelihoods::toString() const {
	ostringstream oss;
	oss << "PhredGenotypeLikelihoods( ";
	for (size_t i = 0; i < gl.size(); i++){
		if (i > 0) oss << ",";
		oss << gl[i];
	}

	oss << ")" << endl;
	return oss.str();
}

unsigned int PhredGenotypeLikelihoods::get_ploidy() const {
	return ploidy;
}

unsigned int PhredGenotypeLikelihoods::get_n_alleles() const {
	return n_alleles;
}
