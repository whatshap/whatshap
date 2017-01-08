#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>
#include<vector>

class Entry {
public:
	
	Entry(unsigned int r, int m, std::vector<unsigned int> p);

	unsigned int get_read_id() const;
	int get_allele_type() const;
	std::vector<unsigned int> get_phred_score() const;

	void set_read_id(unsigned int r);
	void set_allele_type(int m);
	void set_phred_score(std::vector<unsigned int> p);

	friend std::ostream& operator<<(std::ostream& out, const Entry& e);

private:
	unsigned int read_id;
	int allele;
	std::vector<unsigned int> phred_score;
};

#endif
