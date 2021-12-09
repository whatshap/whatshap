#include <bitset>
#include <limits>
#include <cassert>

#include "graycodes.h"

using namespace std;

GrayCodes::GrayCodes(int l) {
	assert(l <= numeric_limits<unsigned int>::digits);
	this->length = l;
	this->s = ~((unsigned int)0);
	this->c = 0;
	this->i = -1;
	this->changed_bit = -1;
}


bool GrayCodes::has_next() {
	return i < length;
}


unsigned int GrayCodes::get_next(int* changed_bit) {
	unsigned int result = c;
	if (changed_bit != 0) {
		*changed_bit = this->changed_bit;
	}
	i = 0;
	while (i < this->length) {
		unsigned int mask = ((unsigned int)1) << i;
		if (((c&mask) ^ (s&mask)) != 0) {
			c = c ^ mask;
			this->changed_bit = i;
			break;
		}
		s = s ^ mask;
		i += 1;
	}
	return result;
}

vector<int> GrayCodes::toBinary(int n) {
    vector<int> binaryVector;
	binaryVector.resize(length);
	for (int index = 0; index < length; index++) {
        binaryVector[index] = n % 2;
        n = n / 2;
    }
	return binaryVector;
}