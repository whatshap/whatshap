#include <iostream>
#include <fstream>

#include <cereal/archives/binary.hpp>

#include "readset.h"

#include "serialization.h"

using namespace std;

void serialize(const ReadSet* readset) {
	ofstream ofs;
	ofs.open("test.ser");
	cereal::BinaryOutputArchive oarchive(ofs);
	oarchive(*readset);
	ofs.close();
}
