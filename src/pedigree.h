#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <iostream>
#include <unordered_map>
#include <vector>

/*
 * how to use:
 * - build up a pedigree by adding individuals and their genotypes with addIndividual(),
 * - then add the relationship between them with addRelationship().
 */
class Pedigree {
public:
	Pedigree();
	//TODO: maybe named: mother, father, child
	typedef std::array<size_t, 3> triple_entry_t;

	typedef struct parent_indices_t {
		size_t mother_index;
		size_t father_index;
		parent_indices_t(size_t mother_index, size_t father_index) : mother_index(mother_index), father_index(father_index) {}
	} parent_indices_t;
	//typedef std::array<std::vector<unsigned int>, 3> genotype_entry_t;

	// add genotypes of a single individual
	void addIndividual(unsigned int individual_id, std::vector<unsigned int> genotypes);

	// add a relationship (a mother/father/child triple)
	void addRelationship(unsigned int mother_id, unsigned int father_id, unsigned int child_id);

	/** Returns the genotype of individual with given index in the given column.
	 *  Note that index of an individual is not its id (in general), see id_to_index().
	 */
	unsigned int get_genotype(size_t individual_index, size_t column) const;
	
	/** Turns the id of an individual into its index. */
	size_t id_to_index(unsigned int individual_id) const;

	size_t size() const {
		return individual_ids.size();
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
	std::vector<triple_entry_t> triples;
	std::vector<unsigned int> individual_ids;
	std::unordered_map<unsigned int, size_t> id_to_index_map;
	// genotypes[i][j] is the genotype of individual with index i at locus j.
	std::vector<std::vector<unsigned int>> genotypes;
	
	// maps individual ids to genotypes

};

#endif
