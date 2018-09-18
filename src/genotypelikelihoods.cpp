#include <sstream>
#include <cassert>
#include <algorithm>
#include "genotypelikelihoods.h"
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

GenotypeLikelihoods::GenotypeLikelihoods(unsigned int ploidy, unsigned int n_alleles, const vector<double>& gl, bool is_phred_scaled) :ploidy(ploidy), n_alleles(n_alleles), gl(gl), is_phred_scaled(is_phred_scaled)
{}

double GenotypeLikelihoods::get(Genotype genotype) const {
	unsigned int index = genotype.get_index(ploidy, n_alleles);
	assert(index < gl.size());
	return this->gl[index];
}

const vector<double>& GenotypeLikelihoods::as_vector() const {
	return gl;
}

size_t GenotypeLikelihoods::genotype_count() const {
	return gl.size();
}

std::string GenotypeLikelihoods::toString() const {
	ostringstream oss;
	oss << "GenotypeLikelihoods( ";
	for (size_t i = 0; i < gl.size(); i++){
		if (i > 0) oss << ",";
		oss << gl[i];
	}

	oss << ")" << endl;
	return oss.str();
}

unsigned int GenotypeLikelihoods::get_ploidy() const {
	return ploidy;
}

unsigned int GenotypeLikelihoods::get_n_alleles() const {
	return n_alleles;
}

bool GenotypeLikelihoods::is_phred() const {
	return is_phred_scaled;
}

// TODO if is_phred, the likeliest genotype is the one with minimum score
Genotype GenotypeLikelihoods::get_likeliest_genotype(double threshold_prob) const {
	// vector of (likelihood,index) pairs
	vector< pair<double,unsigned int> > likelihoods;
	for ( size_t i = 0; i < gl.size(); i++){
		likelihoods.push_back(make_pair(gl[i],i));
	}
	// sort according to likelihoods
	sort(likelihoods.begin(),likelihoods.end(),compare);

	// if empty list of likelihoods, return None-genotype
	if (gl.size() == 0) return Genotype();

	// if only one likelihood, return corresponding gt
	if (gl.size()  == 1) return index_to_genotype(0,ploidy);

	// check if there is a unique maximum and return best genotype
	size_t last_index = likelihoods.size() - 1;
	pair<double,unsigned int> best;
	pair<double,unsigned int> second_best;

	// find the likeliest genotype
	if (is_phred()){
		// likeliest gt is the one with lowest (phred-scaled) probability
		best = likelihoods[0];
		second_best = likelihoods[1];

		if ((best.first < second_best.first) && (best.first < threshold_prob)){
			return index_to_genotype(best.second, ploidy);
		} else {
			return Genotype();
		}

	} else {
		// likeliest gt is the one with highest probability
		best = likelihoods[last_index];
		second_best = likelihoods[last_index-1];

		if ((best.first > second_best.first) && (best.first > threshold_prob)){
			return index_to_genotype(best.second, ploidy);
		} else {
			return Genotype();
		}
	}
}

