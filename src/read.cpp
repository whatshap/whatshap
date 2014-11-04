#include <sstream>

#include "read.h"

using namespace std;

Read::Read(const std::string& name, int mapq) : name(name), mapq(mapq) {
}

string Read::toString() {
	ostringstream oss;
	oss << name << " " << mapq << " (";
	for (size_t i=0; i<variants.size(); ++i) {
		if (i>0) oss << ";";
		oss << "[" << variants[i].position << "," << variants[i].base << "," << variants[i].entry << "]";
	}
	oss << ")";
	return oss.str();
}

void Read::addVariant(int position, char base, int allele, int quality) {
	variants.push_back(enriched_entry_t(position, base, allele, quality));
}
