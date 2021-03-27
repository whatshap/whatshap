#include "reverseComplement.h"

char rev_comp::complement[int('t') + 1];
rev_comp* rev_comp::_instance = nullptr;

std::string rev_comp::run(std::string seq) {
	rev_comp::build_instance();
	
	auto first = seq.begin(), last = seq.end();
        
	while(true) {
		if(first == last || first == --last) {
    	    if(seq.length() % 2) {
    		    *first = rev_comp::complement[(unsigned char) *first];
            }
    	    return seq;
        } else {
        	*first = rev_comp::complement[(unsigned char) *first];
    	    *last = rev_comp::complement[(unsigned char) *last];
    	    std::iter_swap(first, last);
    	    ++first;
        }
    }
}

void rev_comp::build_instance() {
    if(_instance == nullptr) {
        _instance = new rev_comp();
    }
}

rev_comp::rev_comp() {
    this->complement['A'] = 'T';
    this->complement['T'] = 'A';
    this->complement['C'] = 'G';
    this->complement['G'] = 'C';
    this->complement['a'] = 't';
    this->complement['t'] = 'a';
    this->complement['c'] = 'g';
    this->complement['g'] = 'c';
}

rev_comp::~rev_comp(){

}
