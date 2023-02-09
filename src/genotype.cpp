#include <sstream>
#include <stdexcept>
#include "genotype.h"
#include "binomial.h"
#include <iostream>

using namespace std;

Genotype::Genotype() : gt(0) {}

Genotype::Genotype(uint64_t index, uint32_t ploidy) {
	std::vector<uint32_t> genotype = convert_index_to_alleles(index, ploidy);
	
	std::sort(genotype.begin(), genotype.end());
	gt = 0UL;
	// copy to our representation
	for (uint32_t i = 0; i < ploidy; i++) {
		if (genotype[i] >= MAX_ALLELES) {
			throw std::runtime_error("Error: Maximum alleles for genotype exceeded!");
		}
		set_position(ploidy - i - 1, genotype[i]);
	}
	set_ploidy(ploidy);
	for (uint32_t i = 0; i < ploidy-1; i++) {
		if (get_position(i)<get_position(i+1)) {
			throw std::runtime_error("Error: Genotype not sorted! 0 ");
		}
	}
	//std::cout<<"Constructed genotype with index "<<index<<": "<<toString()<<std::endl;
}

Genotype::Genotype(vector<uint32_t> alleles) {
	// parameter check
	gt = 0;
	uint32_t ploidy = alleles.size();
	if (ploidy >= MAX_PLOIDY) {
		throw std::runtime_error("Error: Maximum ploidy for genotype exceeded!");
	}
	std::sort(alleles.begin(), alleles.end());
	for (uint32_t i = 0; i < ploidy; i++) {
		if (alleles[i] >= MAX_ALLELES) {
			throw std::runtime_error("Error: Maximum alleles for genotype exceeded!");
		}
		set_position(ploidy - i - 1, alleles[i]);
	}
	set_ploidy(ploidy);
	if (ploidy > 0) {
		for (uint32_t i = 0; i < ploidy-1; i++) {
			uint32_t first = get_position(i);
			uint32_t second = get_position(i+1);
			if (first < second) {
				std::cout<<"Not sorted at positions "<<i<<" and "<<(i+1)<<" with "<<first<<" < "<<second<<std::endl;
				std::cout<<"Genotype (vector): ";
				for (uint32_t i = 0; i < ploidy; i++) {
					std::cout<<alleles[i]<<" ";
				}
				std::cout<<std::endl;
				std::cout<<"Genotype (bits): ";
				for (uint32_t i = 0; i < ploidy; i++) {
					std::cout<<get_position(i)<<" ";
				}
				std::cout<<std::endl;
				throw std::runtime_error("Error: Genotype not sorted! 1 ");
			}
		}
	}
}

vector<uint32_t> Genotype::as_vector() const {
	vector<uint32_t> alleles;
	uint32_t ploidy = get_ploidy();
	for (uint32_t i = 0; i < ploidy; i++) {
		alleles.push_back(get_position(i));
	}
	return alleles;
}

bool Genotype::is_none() const {
	return get_ploidy() == 0;
}

uint64_t Genotype::get_index() const {
	// use formula given here: https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
	uint32_t ploidy = get_ploidy();
	uint32_t index = 0;
	uint32_t k = 1;
	for (uint32_t i = 0; i < ploidy; i++) {
		uint32_t allele = get_position(i);
		index += binomial_coefficient(k + allele - 1, allele - 1);
		k += 1;
	}
	return index;
}

string Genotype::toString() const {
	ostringstream oss;
	if (is_none()) {
		oss << ".";
		return oss.str();
	}

	uint32_t ploidy = get_ploidy();
	oss << get_position(ploidy-1);
	for (uint32_t i = 1; i < ploidy; i++) {
		oss << '/' << get_position(ploidy-i-1);
	}
	return oss.str();
}

bool Genotype::is_homozygous() const {
	// homozygous <=> all alleles identical
	if (is_none()) return false;
	
	uint32_t ploidy = get_ploidy();
	uint32_t allele = get_position(0);
	for (uint32_t i = 1; i < ploidy; i++) {
		if (get_position(i) != allele)
			return false;
	}
	
	return true;
}

bool Genotype::is_diploid_and_biallelic() const{
	uint32_t ploidy = get_ploidy();
	if (ploidy != 2)
		return false;
	for (uint32_t i = 0; i < ploidy; i++) {
		if (get_position(i) > 1) {
			return false;
		}
	}
	return true;
}

uint32_t Genotype::get_ploidy() const {
	return get_position(MAX_PLOIDY);
}

uint64_t Genotype::get_code() const {
	return gt;
}

bool operator== (const Genotype &g1, const Genotype &g2) {
	return !(g1.gt ^ g2.gt);
}

bool operator!= (const Genotype &g1, const Genotype &g2) {
	return g1.gt ^ g2.gt;
}

bool operator< (const Genotype &g1, const Genotype &g2) {
	return g1.get_index() < g2.get_index();
}

void Genotype::set_ploidy(const uint32_t ploidy) {
	set_position(MAX_PLOIDY, ploidy);
}
	
uint32_t Genotype::get_position(const uint32_t pos) const {
	if (pos < 0 || pos > MAX_PLOIDY)
		throw std::runtime_error("Error: Invalid get position");
	return (gt >> (pos*4)) & (uint64_t)15UL;
}

void Genotype::set_position(const uint32_t pos, const uint32_t allele) {
	if (pos < 0 || pos > MAX_PLOIDY)
		throw std::runtime_error("Error: Invalid set position");
	if (allele >= MAX_ALLELES)
		throw std::runtime_error("Error: Invalid set allele");
	
	uint64_t code = (uint64_t)(allele);
	
	uint64_t set_mask = (code << (pos*4));
	uint64_t delete_mask = ((15UL << (pos*4)) ^ (-1UL));
	
	gt &= delete_mask;
	gt |= set_mask;
}

std::vector<uint32_t> convert_index_to_alleles(uint64_t index, uint32_t ploidy) {
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

uint32_t get_max_genotype_ploidy() {
	return Genotype::MAX_PLOIDY;
}

uint32_t get_max_genotype_alleles() {
	return Genotype::MAX_ALLELES;
}
