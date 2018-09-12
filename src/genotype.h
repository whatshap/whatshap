#ifndef GENOTYPE_H
#define GENOTYPE_H

#include<vector>
#include<set>
#include<string>

class Genotype{
	public:
		// construct empty genotype
		Genotype();
		// construct genotype from the given alleles
		Genotype(std::vector<unsigned int> alleles);
		// add allele to the genotype
		void add_allele(unsigned int allele);
		// return the alleles as (sorted) vector
		std::vector<unsigned int> as_vector();
		// check if genotype is empty
		bool is_none();
		// compute index of genotype in sorted list of all genotypes
		unsigned int get_index(unsigned int ploidy, unsigned int n_alleles);
		// return string representation of the genotype
		std::string toString();
		// operators
		friend bool operator== (const Genotype &g1, const Genotype &g2);
		friend bool operator!= (const Genotype &g1, const Genotype &g2);
		friend bool operator< (const Genotype &g1, const Genotype &g2);

	private:
		std::multiset<unsigned int> alleles;
};

#endif // GENOTYPE_H
