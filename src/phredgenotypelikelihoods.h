#ifndef PHRED_GENOTYPE_LIKELIHOODS_H
#define PHRED_GENOTYPE_LIKELIHOODS_H

#include <array>
#include "genotype.h"

class PhredGenotypeLikelihoods {
public:
	PhredGenotypeLikelihoods(const std::vector<double>& gl, unsigned int ploidy, unsigned int nr_alleles=2);

	double get(Genotype genotype) const;

	std::string toString() const;

	unsigned int get_ploidy() const;

	unsigned int get_nr_alleles() const;

	unsigned int size() const;

	const std::vector<double>& as_vector() const;

	void get_genotypes(std::vector<Genotype>& genotypes) const;

private:
	std::vector<double> gl;
	unsigned int ploidy;
	unsigned int nr_alleles;
};


#endif
