#include <sstream>
#include <cassert>
#include <algorithm>
#include "phredgenotypelikelihoods.h"

using namespace std;

PhredGenotypeLikelihoods::PhredGenotypeLikelihoods(size_t ploidy, size_t, n_alleles, const vector<double>& gl) :ploidy(ploidy), n_alleles(n_alleles), gl(gl) 
{
	// make sure that the right number of likelihoods is given
	// use formula given here: https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
	size_t nr_of_genotypes = binomial_coeff(ploidy + n_alleles - 1, n_alleles - 1);
	assert(gl.size() == nr_of_genotypes);
}

double PhredGenotypeLikelihoods::get(PhredGenotypeLikelihoods::genotype alleles) const {
	assert(alleles.size() == ploidy);
	// sort the alleles numerically
	sort(alleles.begin(), alleles.end());
	// compute the index
	size_t index = genotype_to_index(alleles);
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

size_t genotype_to_index(PhredGenotypeLikelihoods::genotype){
	// use formula given here: https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
	size_t index = 0;
	for (size_t k = 1; k <= ploidy; k++){
		unsigned int allele = genotype[k];
		index += binomial_coeff(k + allele - 1, allele - 1);
	}
	return index;
}

// use implementation from here: https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient
int binomial_coeff(int n, int k){
	int result = 1;
	if (k > n-k) k = n-k;

	for (int i = 0; i < k; i++){
		result *= (n-i);
		result /= (i+1);
	}

	return result;
}
