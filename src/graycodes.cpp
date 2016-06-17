#include <bitset>
#include <limits>
#include <cassert>

#include "graycodes.h"

using namespace std;

GrayCodes::GrayCodes(int length) {
	//cout << "length is : " << length << endl;
	//cout << "max is : " << numeric_limits<GrayCodes::int_t>::digits << endl;
	assert(length <= numeric_limits<GrayCodes::int_t>::digits);
	this->length = length;
	this->s = ~((int_t)0);
	this->c = 0;
	this->i = -1;
	this->changed_bit = -1;
}


bool GrayCodes::has_next() {
	return i < length;
}


GrayCodes::int_t GrayCodes::get_next(int* changed_bit) {
	GrayCodes::int_t result = c;
	if (changed_bit != 0) {
		*changed_bit = this->changed_bit;
	}
	i = 0;
	while (i<length) {
		int_t mask = ((GrayCodes::int_t)1) << i;
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
