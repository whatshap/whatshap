#include <sstream>
#include "genotype.h"
#include "binomial.h"

using namespace std;

Genotype::Genotype(){}

Genotype::Genotype(vector<unsigned int> alleles) :alleles(alleles.begin(), alleles.end()){}

void Genotype::add_allele(unsigned int allele){
	alleles.insert(allele);
}

vector<unsigned int> Genotype::as_vector() const{
	return vector<unsigned int>(alleles.begin(),alleles.end());
}

bool Genotype::is_none() const{
	return alleles.size() == 0;
}

unsigned int Genotype::get_index(unsigned int ploidy, unsigned int n_alleles) const {
	// use formula given here: https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
	unsigned int index = 0;
	unsigned int k = 1;
	for (auto it = alleles.begin(); it != alleles.end(); ++it){
		unsigned int allele = (*it);
		index += binomial_coeff(k + allele - 1, allele - 1);
		k += 1;
	}
	return index;
}

string Genotype::toString() const{
	ostringstream oss;
	if (is_none()) {
		oss << ".";
		return oss.str();
	}

	for (auto it = alleles.begin(); it != alleles.end(); ++it){
		if (it != alleles.begin()) oss << '/';
		oss << *it;
	}
	return oss.str();
}

bool Genotype::is_homozygous() const{
	// homozygous <=> all alleles identical
	if (is_none()) return false;
	if (std::adjacent_find(alleles.begin(), alleles.end(), std::not_equal_to<unsigned int>()) == alleles.end()){
		return true;
	} else {
		return false;
	}
}

bool operator== (const Genotype &g1, const Genotype &g2){
	return g1.alleles == g2.alleles;
}

bool operator!= (const Genotype &g1, const Genotype &g2){
	return !(g1.alleles == g2.alleles);
}

bool operator< (const Genotype &g1, const Genotype &g2) {
	return g1.alleles < g2.alleles;
}
