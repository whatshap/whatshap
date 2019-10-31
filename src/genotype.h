#ifndef GENOTYPE_H
#define GENOTYPE_H

#include<vector>
#include<set>
#include<string>
#include<algorithm>

class Genotype{
	public:
		// construct empty genotype
		Genotype();
	
		// construct genotype from given index and ploidy
		Genotype(unsigned int index, unsigned int ploidy);
	
		// construct genotype from the given alleles
		Genotype(std::vector<unsigned int> alleles);
	
		// add allele to the genotype
		void add_allele(unsigned int allele);
	
		// return the alleles as (sorted) vector
		std::vector<unsigned int> as_vector() const;
	
		// check if genotype is empty
		bool is_none() const;
	
		// compute index of genotype in sorted list of all genotypes
		unsigned int get_index() const;
	
		// return string representation of the genotype
		std::string toString() const;
	
		// check if genotype is homozygous
		bool is_homozygous() const;
	
		// return the ploidy of the genotype
		unsigned int get_ploidy() const;
	
		// operators
		friend bool operator== (const Genotype &g1, const Genotype &g2);
		friend bool operator!= (const Genotype &g1, const Genotype &g2);
		friend bool operator< (const Genotype &g1, const Genotype &g2);

	private:
		std::multiset<unsigned int> alleles;
};

#endif // GENOTYPE_H
