#include <cassert>
#include "columnindexingiterator.h"
#include "column.h"
#include <math.h>
#include <cmath>

using namespace std;

Column::Column(const unsigned int index, const unsigned int* n_ref, const std::vector<unsigned int>& read_ids, const std::vector<unsigned int>& next_read_ids) : 

read_ids(read_ids),
next_read_ids(next_read_ids) {

	bool read_found;
	n_references = *n_ref; // each reference sample has 2 haplotype paths.
	for (auto read = begin(read_ids); read != end(read_ids); read++) {
		read_found = false;
		for (auto next_read = begin(next_read_ids); next_read != end(next_read_ids); next_read++){
			if (*read == *next_read) {
				act_nonterminating_read_ids.push_back(*read);
				read_found = true;
				break;
			}
		}
		if (!read_found) {
			act_terminating_read_ids.push_back(*read);
		}
	}
}


unsigned int Column::get_index(vector<unsigned int>& b1, vector<unsigned int>& b2, unsigned int& r1, unsigned int& r2) {
	unsigned int index = (n_references*n_references*bipartition_to_index(b1,b2))+reference_allele_to_index(r1, r2);
	assert (index < this->get_column_size());
	return index;
}

unsigned int Column::get_index(unsigned int b_index, unsigned int r_index) {
	unsigned int index = (n_references*n_references*b_index)+r_index;
	assert (index < this->get_column_size());
	return index;
}

unique_ptr<ColumnIndexingIterator> Column::get_iterator() {
	return unique_ptr<ColumnIndexingIterator>(new ColumnIndexingIterator(this));
}

unsigned int Column::get_column_size() {
	return pow(2 ,read_ids.size()) * pow(n_references,2);
}

vector<unsigned int> * Column::get_read_ids() {
	return &(this->read_ids);
}

vector <unsigned int> * Column::get_active_nonterminating_read_ids() {
	return &(this->act_nonterminating_read_ids);
}

vector <unsigned int> * Column::get_active_terminating_read_ids() {
	return &(this->act_terminating_read_ids);
}

vector <unsigned int> * Column::get_next_read_ids() {
	return &(this->next_read_ids);
}

vector<vector<unsigned int>> Column::index_to_bipartition(unsigned int& index, int column_type) {
	int c_size;
	vector<unsigned int>* ri;
	if (column_type == 0) {
		ri = &read_ids;
		c_size = this->get_column_size();
	}
	else {
		ri = &act_nonterminating_read_ids;
		c_size = pow(2, ri->size())*pow(n_references, 2);
	}
	assert (index < c_size);
	unsigned int b_index = (unsigned int)floor((double)index/pow(n_references,2));
	vector<vector<unsigned int>> bipartition;
	bipartition.resize(2);
	for (int i = 0; i < ri->size(); i++){
		if (b_index == 0) {
			bipartition[0].push_back(ri->at(i));
			continue;
		}
		bipartition[(int)(b_index%2)].push_back(ri->at(i));
		b_index = b_index / 2;
	}
	return bipartition;
}

vector<unsigned int> Column::index_to_reference_allele(unsigned int& index, int column_type) {
	int c_size;
	vector<unsigned int>* ri;
	if (column_type == 0) {
		ri = &read_ids;
		c_size = this->get_column_size();
	}
	else {
		ri = &act_nonterminating_read_ids;
		c_size = pow(2, ri->size())*pow(n_references, 2);
	}
	assert (index < c_size);
	vector<unsigned int> ref;
	ref.resize(2);
	unsigned int r_index = (unsigned int)(index%(int)pow(n_references,2));
	for (int i = 0; i < 2; i++) {
		ref[1-i] = r_index%n_references;
		r_index = r_index / n_references;
	}
	return ref;
}

unsigned int Column::bipartition_to_index(vector<unsigned int>& b1, vector<unsigned int>& b2) {
	unsigned int index = 0;
	unsigned int count = 0;
	// TODO: Put code to check if all elements in b1 and b2 are in read_ids and b1 and b2 contain all elements of read_ids
	for (int i = 0; i < read_ids.size(); i++) {
		for (int j = 0; j < b2.size(); j++) {
			if (read_ids[i] == b2[j]) {
				index = index + pow(2,i);
				count++;
			}
		}
	}
	assert (read_ids.size() == (count + b1.size()));
	return index;
}

unsigned int Column::reference_allele_to_index(unsigned int& r1, unsigned int& r2) {
	return n_references*r1 + r2;
}

// Gives the compatible bipartitions in the left column (at position index) given the bipartition index of a bipartition of the right column (at position index + 1)

// TODO: You can precompute the free positions and use it somehow?
vector<unsigned int> Column::get_backward_compatible_bipartitions(int b_index) {
	vector<int> free_positions;
    vector<unsigned int> compatible_bipartition = {0};
	int base = 0;
    int count = 0;
    for (int i = 0; i < next_read_ids.size(); i++) {
        int ri = next_read_ids.at(i);
		// This break statement if the last count++ step happened outside the while loop
		if (count >= read_ids.size()) break;
        while (ri > read_ids.at(count)) {
            free_positions.push_back(count);
            compatible_bipartition.resize(2*compatible_bipartition.size());
            for (int j = 0; j < compatible_bipartition.size()/2; j++) {
                compatible_bipartition.at((compatible_bipartition.size()/2)+j) = compatible_bipartition.at(j) + pow(2, count);
            }
            count++;
			// Break out of while loop (otherwise read_ids.at(count) is not defined)
			if (count >= read_ids.size()) break;
        }
		// This break statement if the last count++ step happened inside the while loop and we don't want the base value to increase.
		if (count >= read_ids.size()) break;
        base = base + (pow(2,count)*(b_index%2));
        b_index = b_index/2;
        count++;
    }
    for (int i = 0; i < compatible_bipartition.size(); i++) {
        compatible_bipartition.at(i) += base;
    }
	return compatible_bipartition;
}

// // Gives the compatible bipartitions in the right column (at position index + 1) given the bipartition index of a bipartition of the left column (at position index)
// vector<unsigned int> Column::get_forward_compatible_bipartitions(int b_index) {
// 	vector<unsigned int> compatible_bipartition = {0};
//     int base = 0;
//     int count = 0;
// 	// Determining the bipartition of the non terminating reads
//     for (int i = 0; i < read_ids.size(); i++) {
//         int ri = read_ids.at(i);
//         if (ri == next_read_ids.at(count)) {
//             base = base + (pow(2,count)*(b_index%2));
//             count++;
//         }
//         b_index = b_index/2;
//     }
//     compatible_bipartition.at(0) = base;
// 	// Iterating over the new reads and accordingly getting the compatible bipartition indices
//     while (count < next_read_ids.size()) {
//         compatible_bipartition.resize(2*compatible_bipartition.size());
//         for (int j = 0; j < compatible_bipartition.size()/2; j++) {
//             compatible_bipartition.at((compatible_bipartition.size()/2)+j) = compatible_bipartition.at(j) + pow(2, count);
//         }
//         count++;
//     }
// 	return compatible_bipartition;
// }