#ifndef GENOTYPE_H
#define GENOTYPE_H

#include<vector>
#include<set>
#include<string>
#include<algorithm>

/**
 * Representation of a genotype of arbitrary ploidy with multi-allelic variants. Genotypes are 
 * unordered multisets of alleles and have a canonical index in the VCF-format:
 * 
 * https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
 * 
 * Given the ploidy, there is a bijective mapping between genotypes and non-negative integer
 * numbers. In the diploid, bi-allelic case, the index is equal to the number of alternative
 * allles (either 0, 1 or 2).
 */
class Genotype{
	public:
	
		const static uint32_t DIPLOID = 2;
	
		/**
		 * Creates an empty genotype with no alleles.
		 */
		Genotype();
	
		/**
		 * Creates a genotype of given ploidy using the canonical index (see class description).
		 */
		Genotype(uint32_t index, uint32_t ploidy);
	
		/**
		 * Creates a genotype from a list of given alleles.
		 */
		Genotype(std::vector<uint32_t> alleles);
	
		/**
		 * Returns the genotype's alleles as a vector.
		 */
		std::vector<uint32_t> as_vector() const;
	
		/**
		 * Returns whether the genotype is empty (i.e. invalid).
		 */
		bool is_none() const;
	
		/**
		 * Returns the canonical index of the genotype (see class description).
		 */
		uint32_t get_index() const;
	
		/**
		 * Returns the genotype as readable string.
		 */
		std::string toString() const;
	
		/**
		 * Returns whether the genotype is homozygous.
		 */
		bool is_homozygous() const;
	
		/**
		 * Returns whether the genotype has ploidy 2 and only alleles 0 and 1.
		 */
		bool is_diploid_and_biallelic() const;
	
		/**
		 * Returns the ploidy of the genotype.
		 */
		uint32_t get_ploidy() const;
	
		// operators
		friend bool operator== (const Genotype &g1, const Genotype &g2);
		friend bool operator!= (const Genotype &g1, const Genotype &g2);
		friend bool operator< (const Genotype &g1, const Genotype &g2);

	private:
		std::multiset<uint32_t> alleles;
	
		/**
		 * Creates a sorted vector of alleles from a given canonical index and ploidy.
		 */
		std::vector<uint32_t> convert_index_to_genotype(uint32_t index, uint32_t ploidy);
};

#endif // GENOTYPE_H
