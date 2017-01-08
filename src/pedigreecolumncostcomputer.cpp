#include <cassert>
#include <limits>
#include <utility>
#include <algorithm>
#include <array>
#include <map>
#include <unordered_set>
#include "vector2d.h"
#include "pedigreecolumncostcomputer.h"

using namespace std;
int power(unsigned int base, unsigned int exp)
{
	int result=1;
	while(exp)
	{
		if(exp & 1)
			result*=base;
		exp>>=1;
		base*=base;
	}
	return result;
}
PedigreeColumnCostComputer::PedigreeColumnCostComputer(const std::vector <const Entry *>&column, size_t column_index, const std::vector <unsigned int>& read_marks, const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions, bool distrust_genotypes, unsigned int n_alleles):
	column(column),
	column_index(column_index),
	read_marks(read_marks),
	partitioning(0),
	pedigree(pedigree),
	cost_partition(pedigree_partitions.count(), n_alleles, 0),
	pedigree_partitions(pedigree_partitions),
	n_alleles(n_alleles)
{
	// Enumerate all possible assignments of alleles to haplotypes and 
	// store those that are compatible with genotypes.
	for (unsigned int i = 0; i < power(n_alleles,pedigree_partitions.count()); ++i) {
		bool genotypes_compatible = true;
		unsigned int cost = 0;
		for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
			unsigned int partition0 = pedigree_partitions.haplotype_to_partition(individuals_index,0);
			unsigned int partition1 = pedigree_partitions.haplotype_to_partition(individuals_index,1);
			int allele0 = (i/power(n_alleles,partition0)) % n_alleles;
			int allele1 = (i/power(n_alleles,partition1)) % n_alleles;
			if (distrust_genotypes) {
				int genotype = allele0 + allele1;
				const PhredGenotypeLikelihoods* gls = pedigree->get_genotype_likelihoods(individuals_index, column_index);
				assert(gls != nullptr);
				cost += gls->get(genotype);
			} else {
				int genotype = pedigree->get_genotype(individuals_index, column_index);
				if (allele0 + allele1 != genotype) {
					genotypes_compatible = false;
					break;
				}
			}
		}
		if (genotypes_compatible) {
			allele_assignments.push_back(allele_assignment_t(i,cost));
		}
	}
}


void PedigreeColumnCostComputer::set_partitioning(unsigned int partitioning) {
	cost_partition = Vector2D<unsigned int>(pedigree_partitions.count(), n_alleles, 0);

	partitioning = partitioning;
	for (vector < const Entry * >::const_iterator it = column.begin(); it != column.end(); ++it) {
		auto & entry = **it;
		unsigned int    ind_id = read_marks[entry.get_read_id()];
		bool  entry_in_partition1 = (partitioning & ((unsigned int) 1)) == 0;
		auto atype = entry.get_allele_type();
		if(atype != -1){
			for(int i = 0; i < n_alleles; ++i){
				if(i != atype) {
				      (entry_in_partition1 ? cost_partition.at(pedigree_partitions.haplotype_to_partition(ind_id,0),i):cost_partition.at(pedigree_partitions.haplotype_to_partition(ind_id,1),i)) += entry.get_phred_score()[i];
			  
			      }
			}
		}
	partitioning = partitioning >> 1;
	}
}


void PedigreeColumnCostComputer::update_partitioning(int bit_to_flip) {
	const Entry & entry = *column[bit_to_flip];
	partitioning = partitioning ^ (((unsigned int) 1) << bit_to_flip);
	bool entry_in_partition1 = (partitioning & (((unsigned int) 1) << bit_to_flip)) == 0;
	unsigned int ind_id = read_marks[entry.get_read_id()];
	auto atype = entry.get_allele_type();
	if(atype == -1) return;
	for(int i = 0; i < n_alleles; ++i){
		if(i != atype) {
		(entry_in_partition1 ? cost_partition.at(pedigree_partitions.haplotype_to_partition(ind_id,1),i) : cost_partition.at(pedigree_partitions.haplotype_to_partition(ind_id,0), i)) -= entry.get_phred_score()[i];
		(entry_in_partition1 ? cost_partition.at(pedigree_partitions.haplotype_to_partition(ind_id,0),i) : cost_partition.at(pedigree_partitions.haplotype_to_partition(ind_id,1), i)) += entry.get_phred_score()[i];
		}
	}
}


unsigned int PedigreeColumnCostComputer::get_cost() {
	unsigned int best_cost = numeric_limits < unsigned int >::max();
	for (const allele_assignment_t& a : allele_assignments) {
		unsigned int cost = a.cost;
		for (size_t p = 0; p < pedigree_partitions.count(); ++p) {
			int allele = (a.assignment/power(n_alleles,p)) % n_alleles;
			cost += cost_partition.at(p, allele);
		}
		if (cost < best_cost) {
			best_cost = cost;
		}
	}
	return best_cost;
}


vector <PedigreeColumnCostComputer::phased_variant_t> PedigreeColumnCostComputer::get_alleles() {
	unsigned int best_cost = numeric_limits < unsigned int >::max();
	unsigned int second_best_cost = numeric_limits < unsigned int >::max();
	vector<phased_variant_t> pop_haps(pedigree->size(), phased_variant_t());
	// best_cost_for_allele[individual][haplotype][to_allele] is the best cost for flipping
	// "haplotype" in "individual" to "to_allele"
	std::vector<unsigned int> vec(n_alleles,  numeric_limits < unsigned int >::max());
	Vector2D<std::vector<unsigned int>> best_cost_for_allele(pedigree->size(), 2, vec);
	for (const allele_assignment_t& a : allele_assignments) {
		unsigned int cost = a.cost;
		for (size_t p = 0; p < pedigree_partitions.count(); ++p) {
			int allele = (a.assignment/power(n_alleles,p)) % n_alleles;
			cost += cost_partition.at(p, allele);
		}
		bool new_best = false;
		if (cost <= best_cost) {
			best_cost = cost;
			new_best = true;
		}
		for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
			unsigned int partition0 = pedigree_partitions.haplotype_to_partition(individuals_index,0);
			unsigned int partition1 = pedigree_partitions.haplotype_to_partition(individuals_index,1);
			int allele0 = (a.assignment/power(n_alleles,partition0)) % n_alleles;
			int allele1 = (a.assignment/power(n_alleles,partition1)) % n_alleles;
			if (new_best) {
				pop_haps[individuals_index] = phased_variant_t(allele0,allele1);
				}
			if (cost < best_cost_for_allele.at(individuals_index,0)[allele0]) {
				best_cost_for_allele.at(individuals_index,0)[allele0] = cost;
			}
			if (cost < best_cost_for_allele.at(individuals_index,1)[allele1]) {
				best_cost_for_allele.at(individuals_index,1)[allele1] = cost;
			}
		}
	}

	if (best_cost == numeric_limits < unsigned int >::max()) {
		throw std::runtime_error("Error: Mendelian conflict");
	}

	// Test whether some of the allele assignments (in particular first_best and second_best cost) are ambiguous
	for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
		for (size_t haplotype = 0; haplotype < 2; ++haplotype) {
			unsigned int first_best = numeric_limits < unsigned int >::max();
			unsigned int second_best =  numeric_limits < unsigned int >::max();
			for (int i=0; i< n_alleles;i++){
				int quality = (int)(best_cost_for_allele.at(individuals_index,haplotype)[i]);
				if (quality < first_best)
				{
					second_best = first_best;
					first_best = quality;
				}else if (quality < second_best)
					 second_best = quality;
			}
			//TODO: compute the phasing qualities
			pop_haps[individuals_index].quality = 0;
			if (first_best == second_best && first_best != numeric_limits < unsigned int >::max()){
				if (haplotype == 0) {
					pop_haps[individuals_index].allele0 = -2;
				} else {
					pop_haps[individuals_index].allele1 = -2;
				}
			}
		}
	}

	return pop_haps;
}


unsigned int PedigreeColumnCostComputer::get_weight(bool second_haplotype) {
	throw std::runtime_error("Not yet implemented for pedigrees.");
}
