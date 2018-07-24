#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>
#include <vector>

class Entry {
public:
	typedef enum { REF_ALLELE = 0, ALT_ALLELE = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;
	Entry(unsigned int r, std::vector<allele_t> m, std::vector<int> p);

	unsigned int get_read_id() const;
	std::vector<Entry::allele_t> get_allele_type() const;
	std::vector<int> get_phred_score() const;

	void set_read_id(unsigned int r);
	void set_allele_type(std::vector<allele_t> m);
	void set_phred_score(std::vector<int> p);

	friend std::ostream& operator<<(std::ostream& out, const Entry& e);

private:
	unsigned int read_id;
	std::vector<allele_t> allele_type;
	std::vector<int> phred_score;
};

#endif
