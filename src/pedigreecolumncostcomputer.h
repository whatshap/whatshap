#ifndef PEDIGREE_COLUMN_COST_COMPUTER_H
#define PEDIGREE_COLUMN_COST_COMPUTER_H

#include <array>
#include <vector>
#include <set>
#include <memory>
#include <map>
#include <utility>
#include <array>
#include "entry.h"
#include "pedigree.h"
#include "pedigreepartitions.h"
#include "columnindexingiterator.h"


  
class PedigreeColumnCostComputer {
private:
	const std::vector<const Entry*>& column;
	size_t column_index;
	const std::vector<unsigned int>& read_marks;  
	unsigned int partitioning;
	const Pedigree* pedigree;
	std::vector<std::array<unsigned int, 2>> cost_partition;
	const PedigreePartitions& pedigree_partitions;
	typedef struct allele_assignment_t {
		/** The i-th bit in assignment gives the allele assigned to pedigree partition i. */
		unsigned int assignment;
		/** Cost of this assignment incurred by genotype changes. */
		unsigned int cost;
		allele_assignment_t() : assignment(0), cost(0) {}
		allele_assignment_t(unsigned int assignment, unsigned int cost) : assignment(assignment), cost(cost) {}
	} allele_assignment_t;
	/** All allowed assignments and their costs. */
	std::vector<allele_assignment_t> allele_assignments;
  
public:
  
	PedigreeColumnCostComputer(const std::vector<const Entry*>& column, size_t column_index, const std::vector<unsigned int>& read_marks, const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions, bool distrust_genotypes);

	void set_partitioning(unsigned int partitioning);

	void update_partitioning(int bit_to_flip);

	std::vector<unsigned int> compute_roots(std::vector<Pedigree::triple_entry_t> triples);

	unsigned int get_cost();

	typedef struct phased_variant_t {
		Entry::allele_t allele0;
		Entry::allele_t allele1;
		unsigned int quality;
		phased_variant_t() : allele0(Entry::BLANK), allele1(Entry::BLANK), quality(0) {}
		phased_variant_t(Entry::allele_t allele0, Entry::allele_t allele1) : allele0(allele0), allele1(allele1), quality(0) {}
		phased_variant_t(Entry::allele_t allele0, Entry::allele_t allele1, unsigned int quality) : allele0(allele0), allele1(allele1), quality(quality) {}
	} phased_variant_t;

	/** Returns a phased variants ordered by their index in the pedigree given at construction time. */
	std::vector<phased_variant_t> get_alleles();

	/* Returns the weight (what will be the phred_score) at the current position for super-reads. */
	unsigned int get_weight(bool second_haplotype);
};

#endif
