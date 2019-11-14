#include <sstream>
#include "genotype.h"
#include "binomial.h"

using namespace std;

Genotype::Genotype(){}

Genotype::Genotype(uint32_t index, uint32_t ploidy) {
	std::vector<uint32_t> genotype = convert_index_to_genotype(index, ploidy);
	this->alleles = std::multiset<uint32_t>(genotype.begin(), genotype.end());
}

Genotype::Genotype(vector<uint32_t> alleles) :alleles(alleles.begin(), alleles.end()){}

vector<uint32_t> Genotype::as_vector() const{
	return vector<uint32_t>(alleles.begin(),alleles.end());
}

bool Genotype::is_none() const{
	return alleles.size() == 0;
}

uint32_t Genotype::get_index() const {
	// use formula given here: https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
	uint32_t index = 0;
	uint32_t k = 1;
	for (auto it = alleles.begin(); it != alleles.end(); ++it){
		uint32_t allele = (*it);
		index += binomial_coefficient(k + allele - 1, allele - 1);
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
	if (is_none())
		return false;
	return std::adjacent_find(alleles.begin(), alleles.end(), std::not_equal_to<uint32_t>()) == alleles.end();
}

bool Genotype::is_diploid_and_biallelic() const{
	if (get_ploidy() != 2)
		return false;
	for (uint32_t allele : alleles) {
		if (allele > 1) {
			return false;
		}
	}
	return true;
}

uint32_t Genotype::get_ploidy() const {
	return alleles.size();
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

std::vector<uint32_t> Genotype::convert_index_to_genotype(uint32_t index, uint32_t ploidy) {
	/* The conversion code was taken from here: 
	 * https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
	 */
	
	std::vector<uint32_t> genotype(ploidy, 0);
	uint32_t pth = ploidy;
	uint32_t max_allele_index = index;
	uint32_t leftover_genotype_index = index;
	while (pth > 0)	{
	   for (uint32_t allele_index=0; allele_index <= max_allele_index; ++allele_index) {
		   uint32_t i = binomial_coefficient(pth+allele_index-1, pth);
		   if (i>=leftover_genotype_index || allele_index==max_allele_index) {
			   if (i>leftover_genotype_index)
				   --allele_index;
			   leftover_genotype_index -= binomial_coefficient(pth+allele_index-1, pth);
			   --pth;
			   max_allele_index = allele_index;
			   genotype[pth] = allele_index;
			   break;                
		   }
	   }
	}
	return genotype;
}