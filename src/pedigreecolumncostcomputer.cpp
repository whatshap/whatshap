#include <cassert>
#include <limits>
#include <utility>
#include <algorithm>
#include <array>
#include <map>
#include <unordered_set>
#include "vector2d.h"
#include "pedigreecolumncostcomputer.h"
#include <math.h>

using namespace std;

PedigreeColumnCostComputer::PedigreeColumnCostComputer(const std::vector <const Entry *>&column, size_t column_index, const std::vector <unsigned int>& read_marks, const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions, bool distrust_genotypes):
	column(column),
	column_index(column_index),
	read_marks(read_marks),
	partitioning(0),
	pedigree(pedigree),
	ploidy(pedigree_partitions.get_ploidy()),
	cost_partition(pedigree_partitions.count(), {0,0}),
	pedigree_partitions(pedigree_partitions)
{
	// Enumerate all possible assignments of alleles to haplotypes and
	// store those that are compatible with genotypes.
	for (unsigned int i = 0; i < (1<<pedigree_partitions.count()); ++i) {
		bool genotypes_compatible = true;
		unsigned int cost = 0;
		for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
			// determine the individuals genotype
			int genotype = 0;
			for (unsigned int haplotype = 0; haplotype < ploidy; ++haplotype) {
				unsigned int partition = pedigree_partitions.haplotype_to_partition(individuals_index, haplotype);
				unsigned int allele = (i >> partition) & 1;
				genotype += allele;
			}
			if (distrust_genotypes) {
				const PhredGenotypeLikelihoods* gls = pedigree->get_genotype_likelihoods(individuals_index, column_index);
				assert(gls != nullptr);
				cost += gls->get(genotype);
			} else {
				int true_genotype = pedigree->get_genotype(individuals_index, column_index);
				if (genotype != true_genotype) {
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


void PedigreeColumnCostComputer::set_partitioning(unsigned int p) {
	cost_partition.assign(pedigree_partitions.count(), {0,0});
	partitioning = p;
	for (vector < const Entry * >::const_iterator it = column.begin(); it != column.end(); ++it) {
		auto & entry = **it;
		
		// determine which parition the read is in
		unsigned int partition = p % ploidy;
		unsigned int ind_id = read_marks[entry.get_read_id()];

		for (unsigned int a = 0; a < entry.get_allele_type().size(); a++){
			switch (entry.get_allele_type()[a]) {
			case Entry::REF_ALLELE: 
				cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,partition)][1] += entry.get_phred_score()[a];
				break;
			case Entry::ALT_ALLELE:
				cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,partition)][0] += entry.get_phred_score()[a];
				break;
			case Entry::BLANK: break;
			default:
				assert(false);
			}
		}
		p /= ploidy;
	}
}


void PedigreeColumnCostComputer::update_partitioning(int bit_to_flip, int new_partition) {
	// determine the new partitioning
	unsigned int factor = pow(ploidy, bit_to_flip);
	unsigned int old_partition = (partitioning / factor) % ploidy;
	unsigned int tmp = partitioning - (old_partition * factor);
	partitioning = tmp + new_partition * factor;

	// read entry to consider
	const Entry & entry = *column[bit_to_flip];
	unsigned int ind_id = read_marks[entry.get_read_id()];
	for (unsigned int a = 0; a < entry.get_allele_type().size(); a++){
		switch (entry.get_allele_type()[a]) {
		case Entry::REF_ALLELE:
			cost_partition[pedigree_partitions.haplotype_to_partition(ind_id, old_partition)][1] -= entry.get_phred_score()[a];
			cost_partition[pedigree_partitions.haplotype_to_partition(ind_id, new_partition)][1] += entry.get_phred_score()[a];
			break;
		case Entry::ALT_ALLELE:
			cost_partition[pedigree_partitions.haplotype_to_partition(ind_id, old_partition)][0] -= entry.get_phred_score()[a];
                	cost_partition[pedigree_partitions.haplotype_to_partition(ind_id, new_partition)][0] += entry.get_phred_score()[a];
			break;
	    	case Entry::BLANK:
			break;
		default:
			assert(false);
		}
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
	vector<phased_variant_t> pop_haps(pedigree->size(), phased_variant_t(ploidy));
	// best_cost_for_allele[individual][haplotype][to_allele] is the best cost for flipping
	// "haplotype" in "individual" to "to_allele"
	Vector2D<array<unsigned int,2>> best_cost_for_allele(pedigree->size(), ploidy, {numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max()});
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
			vector<unsigned int> alleles;
			for (unsigned int haplotype = 0; haplotype < ploidy; ++haplotype) {
				unsigned int partition = pedigree_partitions.haplotype_to_partition(individuals_index, haplotype);
				alleles.push_back((a.assignment >> partition) & 1);
			}
			if (new_best) {
				vector<Entry::allele_t> entries;
				for (auto allele: alleles){
					(allele == 0) ? entries.push_back(Entry::REF_ALLELE) : entries.push_back(Entry::ALT_ALLELE);
				}
				pop_haps[individuals_index] = phased_variant_t(entries);
			}
			for (unsigned int haplotype = 0; haplotype < ploidy; ++haplotype) {
				if (cost < best_cost_for_allele.at(individuals_index, haplotype)[alleles[haplotype]]) {
					best_cost_for_allele.at(individuals_index, haplotype)[alleles[haplotype]] = cost;
				}
			}
		}
	}

	if (best_cost == numeric_limits < unsigned int >::max()) {
		throw std::runtime_error("Error: Mendelian conflict");
	}

	// Test whether some of the allele assignments are ambiguous
	for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
		for (size_t haplotype = 0; haplotype < ploidy; ++haplotype) {
			int quality = abs(((int)(best_cost_for_allele.at(individuals_index,haplotype)[0])) - ((int)(best_cost_for_allele.at(individuals_index,haplotype)[1])));
			pop_haps[individuals_index].quality = (unsigned int)quality;
			if (quality == 0) {
				pop_haps[individuals_index].alleles[haplotype] = Entry::EQUAL_SCORES;
			}
		}
	}

	return pop_haps;
}


unsigned int PedigreeColumnCostComputer::get_weight(bool second_haplotype) {
	throw std::runtime_error("Not yet implemented for pedigrees.");
}
