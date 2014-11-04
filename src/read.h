#ifndef READ_H
#define READ_H

#include <string>
#include <vector>

#include "entry.h"
  
class Read {
public:
	Read(const std::string& name, int mapq);
	virtual ~Read() {}
	std::string toString();
	void addVariant(int position, char base, int allele, int quality);
private:
	typedef struct enriched_entry_t {
		Entry entry;
		int position;
		char base;
		enriched_entry_t(int position, char base, int allele, int quality) : entry(0,Entry::allele_t(allele),quality), position(position), base(base) {}
	} enriched_entry_t;
	
	std::string name;
	int mapq;
	std::vector<enriched_entry_t> variants;
};

#endif