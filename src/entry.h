#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>

class Entry {
public:
	typedef enum { REF_ALLELE = 0, ALT_ALLELE = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;
	Entry(unsigned int r, allele_t m, unsigned int p);

	unsigned int get_read_id() const;
	allele_t get_allele_type() const;
	unsigned int get_phred_score() const;

	void set_read_id(unsigned int r);
	void set_allele_type(allele_t m);
	void set_phred_score(unsigned int p);

	friend std::ostream& operator<<(std::ostream& out, const Entry& e);

private:
	unsigned int read_id;
	allele_t allele_type;
	unsigned int phred_score;
};

#endif
