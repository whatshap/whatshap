#include <sstream>
#include <cassert>

#include "phredgenotypelikelihoods.h"
#include "binomial.h"

using namespace std;

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods(const vector<double>& gl, unsigned int ploidy, unsigned int nr_alleles) : gl(gl), ploidy(ploidy), nr_alleles(nr_alleles) {
	unsigned int expected_size = binomial_coefficient(ploidy + nr_alleles - 1, nr_alleles - 1);
	if (expected_size != this->gl.size()) {
		throw runtime_error("Error: wrong number of given genotype likelihoods given.");
	}
}


double PhredGenotypeLikelihoods::get(Genotype genotype) const {
	assert(this->ploidy == genotype.get_ploidy());
	unsigned int index = genotype.get_index();
	assert(index < this->gl.size());
	return this->gl[index];
}


std::string PhredGenotypeLikelihoods::toString() const {
	ostringstream oss;
	oss << "PhredGenotypeLikelihoods(";
	for (size_t i = 0; i < this->gl.size(); ++i) {
		if (i > 0) oss << ",";
		oss << gl[i];
	}
	return oss.str();
}

unsigned int PhredGenotypeLikelihoods::get_ploidy() const {
	return this->ploidy;
}

unsigned int PhredGenotypeLikelihoods::get_nr_alleles() const {
	return this->nr_alleles;
}

unsigned int PhredGenotypeLikelihoods::size() const {
	return this->gl.size();
}

const vector<double>& PhredGenotypeLikelihoods::as_vector() const {
	return this->gl;
}

void PhredGenotypeLikelihoods::get_genotypes(vector<Genotype>& genotypes) const {
	for (unsigned int i = 0; i < this->size(); ++i) {
		genotypes.push_back(Genotype(i, this->ploidy));
	}
} 
