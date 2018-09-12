#include <sstream>
#include "genotype.h"

using namespace std;

// use implementation from here: https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient
int binomial_coeff(int n, int k){
	if (k < 0 || n < 0) return 0;
	int result = 1;
	if (k > n-k) k = n-k;

	for (int i = 0; i < k; i++){
		result *= (n-i);
		result /= (i+1);
	}
	return result;
}

Genotype::Genotype(){}

Genotype::Genotype(vector<unsigned int> alleles) :alleles(alleles.begin(), alleles.end()){}

void Genotype::add_allele(unsigned int allele){
	alleles.insert(allele);
}

vector<unsigned int> Genotype::as_vector(){
	return vector<unsigned int>(alleles.begin(),alleles.end());
}

bool Genotype::is_none(){
	return alleles.size() == 0;
}

unsigned int Genotype::get_index(unsigned int ploidy, unsigned int n_alleles){
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

string Genotype::toString(){
	ostringstream oss;
	for (auto it = alleles.begin(); it != alleles.end(); ++it){
		if (it != alleles.begin()) oss << '/';
		oss << *it;
	}
	return oss.str();
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
