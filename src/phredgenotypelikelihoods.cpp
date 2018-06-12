#include <sstream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods(const vector<double>& gl) : gl(gl) {}

double PhredGenotypeLikelihoods::get(size_t genotype) const {
	assert(genotype < gl.size());
	return this->gl[genotype];
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

	for(auto g_likelihood : this->gl){
		oss  << g_likelihood << " ";
	}
	oss << ")" << endl;
	return oss.str();
}
