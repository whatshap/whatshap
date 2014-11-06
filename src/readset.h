#ifndef READSET_H
#define READSET_H

#include <string>
#include <vector>

#include "read.h"
  
class ReadSet {
public:
	ReadSet();
	virtual ~ReadSet();
	/** Ownership of pointer is transferred from caller to the ReadSet. */
	void add(Read* read);
	/** Sorts reads and variants within reads and assigns unique read identifiers. 
	 *  After calling finalize(), the read set becomes "frozen" and cannot be modified. **/
	void finalize();
	std::string toString();
private:
	typedef struct read_comparator_t {
		read_comparator_t() {}
		bool operator()(const Read* r1, const Read* r2) {
			return r1->firstPosition() < r2->firstPosition();
		}
	} read_comparator_t;

	std::vector<Read*> reads;
	bool finalized;
};

#endif