#ifndef GRAYCODES_H
#define GRAYCODES_H

#include <iostream>

/** A class to generate Gray codes. 
  * Implementation is based on
  * "An Algorithm for Gray Codes", S. Mossige, Computing (18), pp. 89-92, 1977.
  */
class GrayCodes {
	public:
		typedef unsigned int int_t;

		GrayCodes(int length);

		bool has_next();

		/** Return the next Gray code.
		  * @param changed_bit If not null, the index of the changed bit is
		  *                    returned via this variable.
		  */
		int_t get_next(int* changed_bit = 0);
	private:
		int length;
		int i;
		int_t s;
		int_t c;
		int changed_bit;
};

#endif
