#include <iostream>
#include <fstream>
#include <vector>

#include <cereal/archives/binary.hpp>

#include "readset.h"
#include "pedigree.h"
#include "pedigreedptable.h"

using namespace std;

void usage(const char* name) {
	cerr << "Usage: " << name << " [options] <whatshap.instance>" << endl;
	cerr << endl;
	cerr << "De-serializes and runs a WhatsHap problem instance." << endl;
	exit(1);
}

int main(int argc, char* argv[]) {

	ifstream ifs;
	ifs.open("test.ser");
	cereal::BinaryInputArchive iarchive(ifs);
	ReadSet readset;
	Pedigree pedigree;
	std::vector<unsigned int> recombcost;
	bool distrust_genotypes;
	std::vector<unsigned int> positions;
	cout << "Deserializing data structures" << endl;
	iarchive(readset, recombcost, pedigree, distrust_genotypes, positions);
// 	cout << readset.toString() << endl;
// 	cout << pedigree.toString() << endl;
// 	cout << "recombcost: [";
// 	for (auto i : recombcost) {
// 		cout << ' '<< i;
// 	}
// 	cout << ']' << endl;
//	cout << pedigree.toString();
	cout << "ReadSet has " << readset.size() << " reads." << endl;
	ifs.close();

	cout << "Constructing DP table" << endl;
	PedigreeDPTable dptable(&readset, recombcost, &pedigree, distrust_genotypes, &positions);
	cout << "PedMEC cost:" << dptable.get_optimal_score() << endl;

	return 0;
}
