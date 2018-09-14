#include <sstream>
#include <cassert>
#include <algorithm>
#include "phredgenotypelikelihoods.h"
#include "binomial.h"
#include<iostream>

using namespace std;

bool compare(const pair<double,unsigned int>& a, const pair<double, unsigned int>& b){
	return a.first < b.first;
}

Genotype index_to_genotype(unsigned int genotype_index, unsigned int ploidy){
	// use implementation here: https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
	vector<unsigned int> alleles(ploidy, 0);
	int pth = ploidy;
	int max_allele_index = genotype_index;
	int leftover_genotype_index = genotype_index;

	while(pth > 0){
		for (int allele_index = 0; allele_index <= max_allele_index; allele_index++){
			unsigned int i = binomial_coeff(pth + allele_index - 1, pth);
			if (i >= leftover_genotype_index || allele_index == max_allele_index) {
				if (i > leftover_genotype_index) --allele_index;
				leftover_genotype_index -= binomial_coeff(pth + allele_index - 1, pth);
				pth -= 1;
				max_allele_index = allele_index;
				alleles[pth] = allele_index;
				break;
			}
		}
	}
	return Genotype(alleles);
}

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

Genotype PhredGenotypeLikelihoods::get_likeliest_genotype(double threshold_prob) const {
	// vector of (likelihood,index) pairs
	vector< pair<double,unsigned int> > likelihoods;
	for ( size_t i = 0; i < gl.size(); i++){
		likelihoods.push_back(make_pair(gl[i],i));
	}
	// sort according to likelihoods
	sort(likelihoods.begin(),likelihoods.end(),compare);

	// check if there is a unique maximum and return best genotype
	size_t last_index = likelihoods.size() - 1;
	pair<double,unsigned int> best = likelihoods[last_index];
	pair<double,unsigned int> second_best = likelihoods[last_index-1];

	if ( (best.first > second_best.first) && (best.first > threshold_prob) ){
		return index_to_genotype(best.second, ploidy);
	} else {
		return Genotype();
	}
}
	
