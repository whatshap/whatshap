#ifndef GENOTYPE_H
#define GENOTYPE_H

#include<vector>
#include<set>
#include<string>
#include<algorithm>

/**
 * Representation of a genotype of arbitrary ploidy with multi-allelic variants. Genotypes are 
 * unordered multisets of alleles, which are stored in 64bit word, divied into 16 chunks of 4 bits
 * each. The highest 4 bits (left most ones) encode the ploidy, ranging from 0 (invalid genotype) to
 * a maximum of 15. The remaining 15 4-bit-blocks encode the alleles in ascending order, i.e. the
 * lowest 4 bits (right most ones) encode the highest allele. There is a maximum of 16 different
 * alleles, which can be stored by this method.
 *
 * Genotypes have a canonical index in the VCF-format:
 * 
 * https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
 * 
 * Given the ploidy, there is a bijective mapping between genotypes and non-negative integer
 * numbers. In the diploid, bi-allelic case, the index is equal to the number of alternative
 * allles (either 0, 1 or 2):
 * 
 * 0 -> 0/0, 1 -> 0/1, 2 -> 1/1
 * 
 * For the bi-allelic case, this can easily be generalized for higher ploidy, e.g. 4:
 *
 * 0 -> 0/0/0/0, 1 -> 0/0/0/1, 2 -> 0/0/1/1, 3 -> 0/1/1/1, 4 -> 1/1/1/1
 *
 * For multi-allelic variants, there are more possible indices, which are enumerated in such way, that
 * first we have all genotypes, which only use the first allele (i.e. 0/0 for diploid genotypes. Then,
 * all genotypes with the first two alleles follow, then all genotypes using the first three alleles,
 * and so on. This way, the index can be interpreted without knowing the highest allele, as this can
 * be derived from the number itself (the ploidy is mandatory, though!):
 *
 * 0 -> 0/0
 * 1 -> 0/1, 2 -> 1/1
 * 3 -> 0/2, 4 -> 1/2, 5 -> 2/2
 * 6 -> 0/3, 7 -> 1/3, 8 -> 2/3, 9 -> 3/3
 * 
 * 0 -> 0/0/0/0
 * 1 -> 0/0/0/1, 2 -> 0/0/1/1, 3 -> 0/1/1/1, 4 -> 1/1/1/1
 * 5 -> 0/0/0/2, 6 -> 0/0/1/2, 7 -> 0/1/1/2, 8 -> 1/1/1/2, 9 -> 0/0/2/2, 10 -> 0/1/2/2, 11 -> 1/1/2/2, 12 -> 0/2/2/2, 13 -> 1/2/2/2, 14 -> 2/2/2/2
 * etc.
 */
class Genotype{
	public:
		/**
		 * The maximum supported number of alleles
		 */
		const static uint32_t MAX_ALLELES = 16;
	
		/**
		 * The maximum supported ploidy
		 */
		const static uint32_t MAX_PLOIDY = 15;
	
		/**
		 * Constant for diploid ploidy.
		 */
		const static uint32_t DIPLOID = 2;
	
		/**
		 * Creates an empty genotype with no alleles.
		 */
		Genotype();
	
		/**
		 * Creates a genotype of given ploidy using the canonical index (see class description).
		 */
		Genotype(uint64_t index, uint32_t ploidy);
	
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
		uint64_t get_index() const;
	
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
        
        /**
         * Returns the internal representation.
         */
        uint64_t get_code() const;
	
		// operators
		friend bool operator== (const Genotype &g1, const Genotype &g2);
		friend bool operator!= (const Genotype &g1, const Genotype &g2);
		friend bool operator< (const Genotype &g1, const Genotype &g2);

	private:
		/**
		 * Bitstring for storage. Each allele is encoded in 4 bits with a total ploidy
		 * of 15 (=60 bits). The following 4 bits encode the ploidy, with 0 indicating
		 * an invalid genotype. The remaining bits are special flags.
		 */
		uint64_t gt;
	
		// general manipulation methods
		void set_ploidy(const uint32_t ploidy);
		uint32_t get_position(const uint32_t pos) const;
		void set_position(const uint32_t pos, const uint32_t allele);
};

/**
 * Creates a sorted vector of alleles from a given canonical index and ploidy.
 */
std::vector<uint32_t> convert_index_to_alleles(uint64_t index, uint32_t ploidy);

/**
 * Returns the maximum supported ploidy for genotypes
 */
uint32_t get_max_genotype_ploidy();

/**
 * Returns the maximum supported number of alleles per variant for genotypes
 */
uint32_t get_max_genotype_alleles();

/**
 * Provide hash function to use Genotypes as keys for maps.
 */
namespace std {
    template <>
    struct hash<Genotype> {
        size_t operator()(const Genotype& g) const {
            return hash<uint64_t>()(g.get_code());
        }
    };
}

#endif // GENOTYPE_H
