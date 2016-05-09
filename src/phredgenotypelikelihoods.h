#ifndef PHRED_GENOTYPE_LIKELIHOODS_H
#define PHRED_GENOTYPE_LIKELIHOODS_H

#include <array>

#include <cereal/types/array.hpp>

class PhredGenotypeLikelihoods {
public:
	PhredGenotypeLikelihoods(unsigned int gl0 = 0, unsigned int gl1 = 0, unsigned int gl2 = 0);

	unsigned int get(size_t genotype) const;

	std::string toString() const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(gl);
	}

private:
	std::array<unsigned int, 3> gl;
};


#endif
