#ifndef PEDIGREE_PARTITIONS_H
#define PEDIGREE_PARTITIONS_H

#include <array>
#include <vector>

#include "pedigree.h"

/** Class to keep track partitions of reads that are IBD in a pedigree based on 
 *  a specific transmission pattern (given at construction time).
 */
class PedigreePartitions {
private:
	const Pedigree& pedigree;
	unsigned int transmission_vector;
	unsigned int partition_count;
	unsigned int ploidy;
	std::vector<std::vector<int>> haplotype_to_partition_map;
	/** (Recursively) compute entry haplotype_to_partition_map[i]. */
	void compute_haplotype_to_partition_rec(size_t i, const std::vector<int>& triple_indices);
public:
	PedigreePartitions(const Pedigree& pedigree, unsigned int transmission_vector, unsigned int ploidy);

	/** Returns the number of partitions. */
	unsigned int count() const;

	/** Returns an index of the partition for a given individual and haplotype (haplotype must be < ploidy). */
	size_t haplotype_to_partition(size_t individual_index, size_t haplotype) const;

	/** Returns the ploidy of the samples **/
	unsigned int get_ploidy() const;

	friend std::ostream& operator<<(std::ostream& out, const PedigreePartitions& pp);

};

std::ostream& operator<<(std::ostream& out, const PedigreePartitions& pp);

#endif
