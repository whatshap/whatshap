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
	std::string toString();
private:
	std::vector<Read*> reads;
};

#endif