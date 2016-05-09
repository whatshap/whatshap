#include <iostream>
#include <fstream>

#include <cereal/archives/binary.hpp>

#include "serialization.h"

using namespace std;

void serialize(const Read* read) {
	ofstream ofs;
	ofs.open("test.ser");
	cereal::BinaryOutputArchive oarchive(ofs);
	oarchive(*read);
	ofs.close();
}
