#include <iostream>
#include <fstream>

#include <cereal/archives/binary.hpp>

#include "readset.h"

#include "serialization.h"

using namespace std;

void serialize(const ReadSet* readset, const Pedigree* pedigree) {
	ofstream ofs;
	ofs.open("test.ser");
	cereal::BinaryOutputArchive oarchive(ofs);
	oarchive(*readset, *pedigree);
	ofs.close();
}
