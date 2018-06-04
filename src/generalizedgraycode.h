#ifndef  GENERALIZED_GRAYCODE
#define GENERALIZED_GRAYCODE

#include <vector>

/**
the algorithm of Guan is used to generate (n,k)-Gray Codes.
http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=DEDDD5F10CCC1F497AE62CF67A76E8D4?doi=10.1.1.119.1344&rep=rep1&type=pdf
**/

class GeneralizedGrayCodes {
	public:
		typedef unsigned int int_t;
		// generate gray codes of given length and base
		GeneralizedGrayCodes(unsigned int length, unsigned int base);
		// check if there is a next code
		bool has_next();
		// get the next gray code (represented as decimal number)
		// update position and partition to the position that changed and
		// the updated value there
		int_t get_next(int* position, int* partition );
	private:
		unsigned int length;
		unsigned int base;

		std::vector<int> n;
		std::vector<int> g;
		std::vector<int> u;
		int current_index;
		int changed_position;
		unsigned int changed_partition;
		
};

#endif // GENERALIZED_GRAYCODE
