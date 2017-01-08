#include <sstream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods(std::vector<unsigned int> gl) : gl({gl}) {}


unsigned int PhredGenotypeLikelihoods::get(size_t genotype) const {
	assert(genotype < gl.size());
	return this->gl[genotype];
}

std::vector<unsigned int> PhredGenotypeLikelihoods::get_gl() const{
	return gl;
}


std::string PhredGenotypeLikelihoods::toString() const {
	ostringstream oss;
	oss << "PhredGenotypeLikelihoods(" << this->gl[0] << "," << this->gl[1] << "," << this->gl[2] << ")" << endl;
	return oss.str();
}
