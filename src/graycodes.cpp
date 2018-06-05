#include<iostream>
#include "graycodes.h"
#include <limits>
#include <cassert>

using namespace std;

GrayCodes::GrayCodes(unsigned int length, unsigned int base):
	length(length),
	base(base),
	n(length + 1,base),
	g(length + 1,0),
	u(length + 1,1),
	current_index(-1),
	changed_position(-1),
	changed_partition(-1)
{
	assert(length <= numeric_limits<GrayCodes::int_t>::digits);
}

bool GrayCodes::has_next()
{
	return g[length] == 0;
}

GrayCodes::int_t GrayCodes::get_next(int* position, int* partition)
{
	assert(this->has_next());
	unsigned int factor = 1;
	GrayCodes::int_t result = 0;
	for(unsigned int j = 0; j < length; j++){
		result += g[j] * factor;
		factor *= base;
	}

	if(position != 0){
		*position = changed_position;
	}

	if(partition != 0){
		*partition = changed_partition;
	}

	// print to check
//	for(int j = length-1; j >= 0; j--){
//		cout << g[j] <<  " ";
//	}

	unsigned int i = 0;
	int k = g[0] + u[0];
	while( (k >= n[i]) || (k < 0) ){
		u[i] = -u[i];
		i += 1;
		k = g[i] + u[i];
	}
	g[i] = k;
	changed_position = length - i - 1 ;
	//changed_position = i;
	changed_partition = k;

//	std::cout << "length: " << length << " base: " << base << std::endl;
//	std::cout << "GrayCodes: " << *position << " " << *partition << " " << result << std::endl;
	return result;
}
