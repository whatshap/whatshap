//============================================================================
// Name        : FastqStorage.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Read and store the reads from input fastq files
//============================================================================

#ifndef FASTQSTORAGE_H_
#define FASTQSTORAGE_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Read.h"
#include "Types.h"


// A class to store the initial reads from the input fastq files. Single- and paired-end reads
// are stored in separate vectors, but united in a vector of pointers for further usage.
// The read IDs by which the reads are stored are those used in the overlaps file.
class FastqStorage
{
private:
//    std::string PATH = "/ufs/baaijens/ViralQuasispecies/";
    std::string PATH;
    std::string m_single_file;
    std::string m_paired1_file;
    std::string m_paired2_file;
    std::string m_new_id_file;
    unsigned int MAX;
    std::map<std::string, std::string> m_new_readIDs; // dictionary for transforming fastq IDs to overlapfile IDs

//    std::vector< Read > m_singles_vec;
//    std::vector< Read > m_paired_vec;

    void read_singles();
    void read_pairs();
    void read_new_ids();
    void fastq_to_stream(std::string filename, std::vector< std::string > &tmp_vector);

//    unsigned int min_overlap_length = 0; // minimal absolute overlap length hence minimum read length

public:
    std::vector< Read > m_singles_vec;
    std::vector< Read > m_paired_vec;
    std::vector< Read* > m_read_vec;
    std::map< read_id_t, unsigned int > m_ID_to_index; // dictionary for transforming a read ID to the corresponding index of the (pointer)read vector

    unsigned int m_readcount_single;
    unsigned int m_readcount_paired;
    unsigned int m_largest_read_id;


    FastqStorage(const ProgramSettings program_settings)
    {
//        std::cout << "FastqStorage is being created.\n";
        PATH = program_settings.output_dir;
        m_single_file = program_settings.singles_file;
        m_paired1_file = program_settings.paired1_file;
        m_paired2_file = program_settings.paired2_file;
        m_new_id_file = program_settings.id_correspondence;
        MAX = program_settings.max_reads;

        if (m_new_id_file.length() > 0) {
            read_new_ids(); // create a dictionary for transforming fastq IDs to overlapfile IDs
        }
    	if (m_single_file.length() > 0 && m_single_file != "None") {
    		read_singles(); // store single-end reads in m_singles_vec
    	}
    	m_readcount_single = m_singles_vec.size();
    	if (m_paired1_file.length() > 0 && m_paired1_file != "None") {
        	read_pairs(); // store paired-end reads in m_paired_vec
    	}
    	m_readcount_paired = m_paired_vec.size();
        if (program_settings.verbose) {
        	std::cout << "Singles: " << m_readcount_single << std::endl;
        	std::cout << "Pairs: " << m_readcount_paired << std::endl;
        }
    	unsigned int count = 0;

    	// create a vector of pointers to all reads (both single and paired) together with an index m_ID_to_index
    	assert (m_read_vec.empty());
    	assert (m_ID_to_index.empty());
    	for (auto& read_it : m_singles_vec) {
    	    m_read_vec.push_back( &read_it );
    	    m_ID_to_index.insert( std::make_pair( read_it.get_read_id(), count) );
    	    count++;
    	}
    	for (auto& read_it : m_paired_vec) {
    	    m_read_vec.push_back( &read_it );
    	    m_ID_to_index.insert( std::make_pair( read_it.get_read_id(), count) );
    	    count++;
    	}
    }

    ~FastqStorage(void) {
//        std::cout << "FastqStorage is being deleted." << std::endl;
    }

    Read* get_read(read_id_t ID);
    unsigned int get_readcount();
    void writeIDsToFile(std::string filename); // write the (overlap)IDs of reads stored in this fastq storage to a file
};



#endif /* FASTQSTORAGE_H_ */
