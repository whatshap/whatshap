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

PedigreeColumnCostComputer::PedigreeColumnCostComputer(const std::vector <const Entry *>&column, size_t column_index, const std::vector <unsigned int>& read_marks, const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions, bool distrust_genotypes):
	column(column),
	column_index(column_index),
	read_marks(read_marks),
	partitioning(0),
	pedigree(pedigree),
	cost_partition(pedigree_partitions.count(), {0,0}),
	pedigree_partitions(pedigree_partitions)
{
	// Enumerate all possible assignments of alleles to haplotypes and 
	// store those that are compatible with genotypes.
	for (unsigned int i = 0; i < (1<<pedigree_partitions.count()); ++i) {
		bool genotypes_compatible = true;
		unsigned int cost = 0;
		for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
			unsigned int partition0 = pedigree_partitions.haplotype_to_partition(individuals_index,0);
			unsigned int partition1 = pedigree_partitions.haplotype_to_partition(individuals_index,1);
			unsigned int allele0 = (i >> partition0) & 1;
			unsigned int allele1 = (i >> partition1) & 1;
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
	cost_partition.assign(pedigree_partitions.count(), {0,0});

	partitioning = partitioning;
	for (vector < const Entry * >::const_iterator it = column.begin(); it != column.end(); ++it) {
		auto & entry = **it;
		bool  entry_in_partition1 = (partitioning & ((unsigned int) 1)) == 0;
		unsigned int    ind_id = read_marks[entry.get_read_id()];
		switch (entry.get_allele_type()) {

		case Entry::REF_ALLELE:
			(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)] :cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)])[1] += entry.get_phred_score();
			break;
		case Entry::ALT_ALLELE:
			(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)] :cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)])[0] += entry.get_phred_score();
			break;
		case Entry::BLANK:
			break;
		default:
			assert(false);
		}
		partitioning = partitioning >> 1;
	}
}


void PedigreeColumnCostComputer::update_partitioning(int bit_to_flip) {
	const Entry & entry = *column[bit_to_flip];
	partitioning = partitioning ^ (((unsigned int) 1) << bit_to_flip);
	bool entry_in_partition1 = (partitioning & (((unsigned int) 1) << bit_to_flip)) == 0;
	unsigned int ind_id = read_marks[entry.get_read_id()];
	switch (entry.get_allele_type()) {
	case Entry::REF_ALLELE:
		(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)] : cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)])[1] -= entry.get_phred_score();
		(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)] :  cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)])[1] += entry.get_phred_score();
		break;
	case Entry::ALT_ALLELE:
		(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)] : cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)])[0] -= entry.get_phred_score();
		(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)] :  cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)])[0] += entry.get_phred_score();
		break;
    case Entry::BLANK:
		break;
	default:
		assert(false);
	}
}


unsigned int PedigreeColumnCostComputer::get_cost() {
	unsigned int best_cost = numeric_limits < unsigned int >::max();
	for (const allele_assignment_t& a : allele_assignments) {
		unsigned int cost = a.cost;
		for (size_t p = 0; p < pedigree_partitions.count(); ++p) {
			unsigned int allele = (a.assignment >> p) & 1;
			cost += cost_partition[p][allele];
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
	Vector2D<array<unsigned int,2>> best_cost_for_allele(pedigree->size(), 2, {numeric_limits<unsigned int>::max(),numeric_limits<unsigned int>::max()});
	for (const allele_assignment_t& a : allele_assignments) {
		unsigned int cost = a.cost;
		for (size_t p = 0; p < pedigree_partitions.count(); ++p) {
			unsigned int allele = (a.assignment >> p) & 1;
			cost += cost_partition[p][allele];
		}
		bool new_best = false;
		if (cost <= best_cost) {
			best_cost = cost;
			new_best = true;
		}
		for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
			unsigned int partition0 = pedigree_partitions.haplotype_to_partition(individuals_index,0);
			unsigned int partition1 = pedigree_partitions.haplotype_to_partition(individuals_index,1);
			unsigned int allele0 = (a.assignment >> partition0) & 1;
			unsigned int allele1 = (a.assignment >> partition1) & 1;
			if (new_best) {
				pop_haps[individuals_index] = phased_variant_t(
					(allele0 == 0)?Entry::REF_ALLELE:Entry::ALT_ALLELE,
					(allele1 == 0)?Entry::REF_ALLELE:Entry::ALT_ALLELE
				);
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

	// Test whether some of the allele assignments are ambiguous
	for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
		for (size_t haplotype = 0; haplotype < 2; ++haplotype) {
			int quality = abs(((int)(best_cost_for_allele.at(individuals_index,haplotype)[0])) - ((int)(best_cost_for_allele.at(individuals_index,haplotype)[1])));
			pop_haps[individuals_index].quality = (unsigned int)quality;
			if (quality == 0) {
				if (haplotype == 0) {
					pop_haps[individuals_index].allele0 = Entry::EQUAL_SCORES;
				} else {
					pop_haps[individuals_index].allele1 = Entry::EQUAL_SCORES;
				}
			}
		}
	}

	return pop_haps;
}


unsigned int PedigreeColumnCostComputer::get_weight(bool second_haplotype) {
	throw std::runtime_error("Not yet implemented for pedigrees.");
}
