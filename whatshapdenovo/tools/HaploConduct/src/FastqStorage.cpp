//============================================================================
// Name        : FastqStorage.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Read and store the reads from input fastq files
//============================================================================


#include <assert.h>
#include <fstream>
#include <sstream>
#include <utility>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "FastqStorage.h"

// write the (overlap)IDs of reads stored in this fastq storage to a file
void FastqStorage::writeIDsToFile(std::string filename) {
    std::string path = PATH + filename;
    std::ofstream id_file(path);
	std::vector< Read* >::const_iterator it;
	for (it = m_read_vec.begin(); it != m_read_vec.end(); it++) {
	    std::string line = read_id_to_str((*it)->get_read_id()) + "\n";
	    id_file << line;
    }
    id_file.close();
}

Read* FastqStorage::get_read(read_id_t ID) {
    unsigned int index = m_ID_to_index.at(ID);
    return m_read_vec[index];
}

unsigned int FastqStorage::get_readcount() {
    return m_read_vec.size();
}

// read fastq file into a (temporary) vector
void FastqStorage::fastq_to_stream(std::string filename, std::vector< std::string > &tmp_vector) {
    std::ifstream fastqfile (filename.c_str());
    unsigned int count = 0;
    std::string line;
	if (fastqfile.is_open()) {
		while (getline(fastqfile, line) && count < 4*MAX) {
			tmp_vector.push_back(line);
			count++;
		}
		fastqfile.close();
	}
	else {
	    std::cerr << "Unable to open fastq file " << filename << std::endl;
        exit(1);
    }
}

// create a dictionary for transforming fastq IDs to overlapfile IDs
void FastqStorage::read_new_ids() {
    std::string tupleline;
    std::stringstream ss;
    std::ifstream new_ids (m_new_id_file.c_str());
    unsigned int max_id = 0;
    if (new_ids.is_open()) {
        std::string new_read_ID;
        std::string old_read_ID;
        while (getline(new_ids, tupleline)) {
            ss << tupleline;
            getline(ss, new_read_ID, '\t');
            getline(ss, old_read_ID, '\t');
            if (old_read_ID.at(0) == '>') {
                std::string tmp = old_read_ID.substr(1);
                old_read_ID = tmp;
            }
            read_id_t newid = str_to_read_id(new_read_ID);
            if (newid > max_id) { max_id = newid; }
//            m_new_readIDs.insert(std::pair< std::string, std::string >(old_read_ID.substr(1), new_read_ID));
            m_new_readIDs.insert(std::pair< std::string, std::string >(old_read_ID, new_read_ID));
	        ss << "";
	        ss.clear();
        }
        new_ids.close();
    }
    else {
        std::cerr << "Unable to open read-to-overlapID file";
        exit(1);
    }
    m_largest_read_id = max_id;
}

void FastqStorage::read_singles() {
    std::vector< std::string > tmp_fastqvector;
	fastq_to_stream(m_single_file, tmp_fastqvector);

	std::vector<std::string>::const_iterator it;
	it = tmp_fastqvector.begin();
	int c = 0;
	read_id_t id = 0;
	std::string seq;
	std::string phred;
	while (it != tmp_fastqvector.end()) {
		switch (c%4) {
			case 0: { // ID-line
			    if (*(it->begin()) != '@') {
			        std::cerr << "Read ID does not start with @. Exiting read_singles.\n";
                    exit(1);
			    }
			    std::stringstream stream(it->substr(1));
			    std::string id1;
			    stream >> id1;
			    if (m_new_id_file.length() > 0) {
			        id = str_to_read_id(m_new_readIDs.at(id1)); // change read ID into corresponding overlaps read ID
			    }
			    else {
			        id = str_to_read_id(id1);
			    }
				it++;
				c++;
				break;
		    }
			case 1: { // seq-line
			    seq = boost::to_upper_copy(*it);
				it++;
				c++;
				break;
			}
			case 2: { // + line
				it++;
				c++;
				break;
			}
			case 3: { // qual-line
			    phred = *it;
				it++;
				c++;
				bool is_paired = false;
				bool is_super = false;
				if (seq.length() == 0) {
				    std::cerr << "single read with ID " << id << " has an empty sequence... exiting.\n";
				    exit(1);
				}
                else {
                    Read current_read(is_paired, is_super, id, seq, "", phred, "");
                    (m_singles_vec).push_back(current_read);
                }
				break;
			}
		}
	}
	tmp_fastqvector.clear();
}

void FastqStorage::read_pairs() {
    std::vector< std::string > tmp_fastqvector1;
    std::vector< std::string > tmp_fastqvector2;
	fastq_to_stream(m_paired1_file, tmp_fastqvector1);
	fastq_to_stream(m_paired2_file, tmp_fastqvector2);

	std::vector<std::string>::const_iterator it1, it2;
	it1 = tmp_fastqvector1.begin();
	it2 = tmp_fastqvector2.begin();
	int c = 0;
	read_id_t id = 0;
	std::string seq1, seq2;
	std::string phred1, phred2;
	while (it1 != tmp_fastqvector1.end() && it2 != tmp_fastqvector2.end()) {
		switch (c%4) {
			case 0: {
			    if (*(it1->begin()) != '@') {
			        std::cerr << "Read ID does not start with @. Exiting read_pairs.\n";
                    exit(1);
			    }
			    std::stringstream stream1(it1->substr(1));
			    std::string id1;
			    stream1 >> id1;
			    std::stringstream stream2(it2->substr(1));
			    std::string id2;
			    stream2 >> id2;
			    if (id1 != id2) {
			        std::cerr << "Fastq files /1 /2 are not ordered identically. Exiting read_pairs.\n";
                    exit(1);
			    }
			    if (m_new_id_file.length() > 0) {
			        id = str_to_read_id(m_new_readIDs.at(id1)); // change read ID into corresponding overlaps read ID
			    }
			    else {
			        id = str_to_read_id(id1);
			    }
				it1++;
				it2++;
				c++;
				break;
		    }
			case 1: {
			    seq1 = *it1;
			    seq2 = *it2;
//			    seq2 = Read::build_rev_comp(*it2); // TRANSFORM REVERSE (/2) READ TO FORWARD STRAND
				it1++;
				it2++;
				c++;
				break;
			}
			case 2: {
				it1++;
				it2++;
				c++;
				break;
			}
			case 3: {
			    phred1 = *it1;
			    phred2 = *it2;
//			    std::string rev_phred(it2->rbegin(), it2->rend()); // ALSO REVERSE THE QUALITY SEQUENCE
//			    phred2 = rev_phred;
				it1++;
				it2++;
				c++;
				bool is_paired = true;
				bool is_super = false;
				if (seq1.length() == 0 || seq2.length() == 0) {
				    std::cerr << "paired read with ID " << id << " has an empty sequence... exiting.\n";
				    exit(1);
				}
//                else if (seq1.length() > min_overlap_length && seq2.length() > min_overlap_length) {
                else {
                    Read current_read(is_paired, is_super, id, seq1, seq2, phred1, phred2);
	                (m_paired_vec).push_back(current_read);
	            }
				break;
			}
		}
	}
	tmp_fastqvector1.clear();
	tmp_fastqvector2.clear();
}
