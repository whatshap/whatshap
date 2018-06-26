#include <sstream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods() : gl({0,0,0}) {}

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods(const std::vector<double>& gl) : gl(gl) {}


double PhredGenotypeLikelihoods::get(size_t genotype) const {
	assert(genotype < gl.size());
	return this->gl[genotype];
}

const std::vector<double>& PhredGenotypeLikelihoods::as_vector() const {
	return gl;
}


size_t PhredGenotypeLikelihoods::genotype_count() const {
	return gl.size();
}


std::string PhredGenotypeLikelihoods::toString() const {
	ostringstream oss;
	oss << "PhredGenotypeLikelihoods(" << this->gl[0] << "," << this->gl[1] << "," << this->gl[2] << ")" << endl;
	return oss.str();
}
