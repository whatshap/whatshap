#include <iostream>
#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include "readset.h"

#include "serialization.h"

using namespace std;

void serialize(const ReadSet* readset, const std::vector<unsigned int>& recombcost, const Pedigree* pedigree, bool distrust_genotypes, const std::vector<unsigned int>& positions) {
	ofstream ofs;
	ofs.open("test.ser");
	cereal::BinaryOutputArchive oarchive(ofs);
	oarchive(*readset, recombcost, *pedigree, distrust_genotypes, positions);
	ofs.close();
}
