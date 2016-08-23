#include <sstream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods(unsigned int gl0, unsigned int gl1, unsigned int gl2) : gl({{gl0,gl1,gl2}}) {}


unsigned int PhredGenotypeLikelihoods::get(size_t genotype) const {
	assert(genotype < 3);
	return this->gl[genotype];
}


std::string PhredGenotypeLikelihoods::toString() const {
	ostringstream oss;
	oss << "PhredGenotypeLikelihoods(" << this->gl[0] << "," << this->gl[1] << "," << this->gl[2] << ")" << endl;
	return oss.str();
}
