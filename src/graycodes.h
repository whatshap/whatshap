#ifndef GRAYCODES_H
#define GRAYCODES_H

#include <iostream>
#include <vector>

/** A class to generate Gray codes. 
  * Implementation is based on
  * "An Algorithm for Gray Codes", S. Mossige, Computing (18), pp. 89-92, 1977.
  */
class GrayCodes {
	public:

		GrayCodes(int length);

		bool has_next();

		/** Return the next Gray code.
		  * @param changed_bit If not null, the index of the changed bit is
		  *                    returned via this variable.
		  */
		unsigned int get_next(int* changed_bit = 0);

		// Returns the binary form of bipartiton based on number of active reads.
		std::vector<int> toBinary(int n);

	private:
		int length;
		int i;
		unsigned int s;
		unsigned int c;
		int changed_bit;
};

#endif
