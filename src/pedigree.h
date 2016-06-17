#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <iostream>
#include <unordered_map>
#include <vector>

#include "phredgenotypelikelihoods.h"

/*
 * how to use:
 * - build up a pedigree by adding individuals and their genotypes with addIndividual(),
 * - then add the relationship between them with addRelationship().
 */
class Pedigree {
public:
	Pedigree();
	virtual ~Pedigree();
	//TODO: maybe named: mother, father, child
	typedef std::array<size_t, 3> triple_entry_t;

	typedef struct parent_indices_t {
		size_t mother_index;
		size_t father_index;
		parent_indices_t(size_t mother_index, size_t father_index) : mother_index(mother_index), father_index(father_index) {}
	} parent_indices_t;
	//typedef std::array<std::vector<unsigned int>, 3> genotype_entry_t;

	/** Add an individual with associated genotypes and genotype_likelihoods.
	 *  Ownership of pointer given in genotype_likelihoods are transferred to Pedigree object. */
	void addIndividual(unsigned int individual_id, std::vector<unsigned int> genotypes, std::vector<PhredGenotypeLikelihoods*> genotype_likelihoods);

	// add a relationship (a mother/father/child triple)
	void addRelationship(unsigned int mother_id, unsigned int father_id, unsigned int child_id);

	/** Returns the genotype of individual with given index for the given variant_index.
	 *  Note that index of an individual is not its id (in general), see id_to_index().
	 */
	unsigned int get_genotype(size_t individual_index, size_t variant_index) const;

	/** Returns a genotype based on an individuals id. */
	unsigned int get_genotype_by_id(unsigned int individual_id, unsigned int variant_index) const;

	/** Returns the genotype likelihoods of individual with given index for the given variant_index.
	 *  Note that index of an individual is not its id (in general), see id_to_index().
	 */
	const PhredGenotypeLikelihoods* get_genotype_likelihoods(size_t individual_index, size_t variant_index) const;

	/** Returns a genotype based on an individuals id. */
	const PhredGenotypeLikelihoods* get_genotype_likelihoods_by_id(unsigned int individual_id, unsigned int variant_index) const;

	/** Turns the id of an individual into its index. */
	size_t id_to_index(unsigned int individual_id) const;

	/** Returns the id of individual with given index. */
	unsigned int index_to_id(size_t individual_index) const;

	size_t size() const {
		return individual_ids.size();
	}

	size_t get_variant_count() const {
		return variant_count;
	}

	size_t triple_count() const {
		return triples.size();
	}

	/** Get triples of indices in a trio relationship. Note that each individual
	 *  in the triple is represented by its index (and not its id), see id_to_index().
	 */
	const std::vector<triple_entry_t>& get_triples() const;

	std::string toString() const;

private:
	int variant_count;
	std::vector<triple_entry_t> triples;
	std::vector<unsigned int> individual_ids;
	std::unordered_map<unsigned int, size_t> id_to_index_map;
	// genotypes[i][j] is the genotype of individual with index i at locus j.
	std::vector<std::vector<unsigned int>> genotypes;
	std::vector<std::vector<PhredGenotypeLikelihoods*>> genotype_likelihoods;
};

#endif
