#include <cassert>
#include <limits>
#include <utility>
#include <algorithm>
#include <array>
#include <map>
#include <unordered_set>
#include "pedigreecolumncostcomputer.h"

using namespace std;

PedigreeColumnCostComputer::PedigreeColumnCostComputer(const std::vector <const Entry *>&column, size_t column_index, const std::vector <unsigned int>& read_marks,
						const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions):
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
		for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
			unsigned int partition0 = pedigree_partitions.haplotype_to_partition(individuals_index,0);
			unsigned int partition1 = pedigree_partitions.haplotype_to_partition(individuals_index,1);
			unsigned int allele0 = (i >> partition0) & 1;
			unsigned int allele1 = (i >> partition1) & 1;
			int genotype = pedigree->get_genotype(individuals_index, column_index);
			if (allele0 + allele1 != genotype) {
				genotypes_compatible = false;
				break;
			}
		}
		if (genotypes_compatible) {
			allele_assignments.push_back(i);
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

		case Entry::MAJOR_ALLELE:
			(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)] :cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)])[1] += entry.get_phred_score();
			break;
		case Entry::MINOR_ALLELE:
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

void PedigreeColumnCostComputer::update_partitioning(int bit_to_flip)
{
    // update cost based on the changed bit

	const Entry & entry = *column[bit_to_flip];
	partitioning = partitioning ^ (((unsigned int) 1) << bit_to_flip);
	bool entry_in_partition1 = (partitioning & (((unsigned int) 1) << bit_to_flip)) == 0;
	unsigned int ind_id = read_marks[entry.get_read_id()];
	switch (entry.get_allele_type()) {
	case Entry::MAJOR_ALLELE:
		(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)] : cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)])[1] -= entry.get_phred_score();
		(entry_in_partition1 ? cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)] :  cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)])[1] += entry.get_phred_score();
		break;
	case Entry::MINOR_ALLELE:
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
	for (unsigned int& i : allele_assignments) {
		unsigned int cost = 0;
		for (size_t p = 0; p < pedigree_partitions.count(); ++p) {
			unsigned int allele = (i >> p) & 1;
			cost += cost_partition[p][allele];
		}
		if (cost < best_cost) {
			best_cost = cost;
		}
	}
	return best_cost;
}

std::vector <std::pair <Entry::allele_t,Entry::allele_t >> PedigreeColumnCostComputer::get_alleles() {
    // TODO: avoid code duplication
	unsigned int best_cost = numeric_limits < unsigned int >::max();
	unsigned int second_best_cost = numeric_limits < unsigned int >::max();
	std::vector <std::pair < Entry::allele_t, Entry::allele_t >> pop_haps(pedigree->size(), std::make_pair(Entry::EQUAL_SCORES, Entry::EQUAL_SCORES));
	for (unsigned int& i : allele_assignments) {
		unsigned int cost = 0;
		for (size_t p = 0; p < pedigree_partitions.count(); ++p) {
			unsigned int allele = (i >> p) & 1;
			cost += cost_partition[p][allele];
		}
		if (cost <= best_cost) {
			second_best_cost = best_cost;
			best_cost = cost;
			for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
				unsigned int partition0 = pedigree_partitions.haplotype_to_partition(individuals_index,0);
				unsigned int partition1 = pedigree_partitions.haplotype_to_partition(individuals_index,1);
				unsigned int allele0 = (i >> partition0) & 1;
				unsigned int allele1 = (i >> partition1) & 1;
				pop_haps[individuals_index] = std::make_pair(
					(allele0 == 0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE,
					(allele1 == 0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE
				);
			}
		}
	}

	if (best_cost == numeric_limits < unsigned int >::max()) {
		throw std::runtime_error("Error: Mendelian conflict");
	}

	if (second_best_cost == best_cost) {
		// TODO: Having too equal scores does not necissarily imply that
		// the case is undecidable for a pedigree individuals. Be smarter here.
		pop_haps.assign(pedigree->size(), std::make_pair(Entry::EQUAL_SCORES, Entry::EQUAL_SCORES));
	}

	return pop_haps;
}

unsigned int PedigreeColumnCostComputer::get_weight(bool second_haplotype)
{
    throw std::runtime_error("Not yet implemented for pedigrees.");
}
