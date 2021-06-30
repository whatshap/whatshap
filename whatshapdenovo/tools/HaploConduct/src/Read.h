//============================================================================
// Name        : Read.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Read class for overlap graph construction
//============================================================================

#ifndef READ_H_
#define READ_H_

#include <string>
#include <iostream>
#include <list>
#include <assert.h>
#include <set>
#include <unordered_map>

#include "Types.h"

class Read
{
private:
    bool m_is_paired; // indicates if paired-end read
    bool m_is_super; // indicates if superread
	read_id_t m_read_id; // read identifier
    node_id_t m_vertex_id_N; // vertex id in overlap graph, NORMAL ORIENTATION
    node_id_t m_vertex_id_R; // vertex id in overlap graph, REVERSE ORIENTATION
	std::string m_seq1, m_seq2; // nucleotide sequences (/1 and /2, respectively)
	std::string m_phred1, m_phred2; // Phred scores (/1 and /2, respectively)
	bool node_N_set;
	bool node_R_set;
    bool m_is_tip;

    // only for superreads:
    std::list< node_id_t > sorted_vertices1; // list of clique nodes ordered by index in superread for /1 seq
    std::list< node_id_t > sorted_vertices2; // similar for /2 seq
//    std::set< read_id_t > original_read_set; // set containing all read IDs of ORIGINAL reads that are part of the superread
    std::unordered_map< read_id_t, OriginalIndex > original_read_indexes; // dict mapping ORIGINAL read IDs that are part of the superread to the corresponding index(es)
    std::unordered_map< node_id_t, SubreadInfo > subreadMap; // dict mapping node to subread details

//	std::list<unsigned int> m_clique_sorted1; // sorted clique (in current graph) corresponding to read indexes, given by read IDs, for /1 seq
//	std::list<unsigned int> m_clique_sorted2; // similar for /2 seq
//	std::list<int> m_read_indexes1; // list containing the indexes of each of the reads in m_clique in the superread for /1 seq
//	std::list<int> m_read_indexes2; // similar for /2 seq
//	std::list<int> m_read_startpos1; // list containing the subsequence start position of each of the reads in m_clique for /1 seq
//	std::list<int> m_read_startpos2; // similar for /2 seq

//	bool clique0_set;
//	bool clique1_set;
//	bool clique2_set;

public:
	Read(bool is_paired, bool is_super, read_id_t id, std::string seq1, std::string seq2, std::string phred1, std::string phred2) :
        m_is_paired(is_paired),
        m_is_super(is_super),
        m_read_id(id),
        m_seq1(seq1),
        m_seq2(seq2),
        m_phred1(phred1),
        m_phred2(phred2)
        {
            node_N_set = false;
            node_R_set = false;
            m_is_tip = false;
//            clique0_set = false;
//	        clique1_set = false;
//	        clique2_set = false;
        }

    bool is_tip() {
        return m_is_tip;
    }

    void set_tip() {
        m_is_tip = true;
    }

    void set_vertex_id(bool normal, read_id_t id) {
        if (normal) {
            m_vertex_id_N = id;
            node_N_set = true;
        }
        else {
            m_vertex_id_R = id;
            node_R_set = true;
        }
    }

    void set_read_id(read_id_t id) {
        m_read_id = id;
    }

	read_id_t get_read_id() const {
		return m_read_id;
    }

    read_id_t get_vertex_id(bool normal) const {
        if (normal) {
            assert (node_N_set);
            return m_vertex_id_N;
        }
        else {
            assert (node_R_set);
            return m_vertex_id_R;
        }
    }

    void set_super(bool b) {
        m_is_super = b;
    }

    bool is_paired() const {
        return m_is_paired;
    }

    bool is_super() const {
        return m_is_super;
    }

	// static std::string build_rev_comp(std::string seq) {
	//     std::string::reverse_iterator it;
	//     std::string rev_comp = "";
    //     for (it = seq.rbegin(); it != seq.rend(); it++) {
    //         if (*it == 'A')
    //             rev_comp.append("T");
    //         else if (*it == 'T')
    //             rev_comp.append("A");
    //         else if (*it == 'C')
    //             rev_comp.append("G");
    //         else if (*it == 'G')
    //             rev_comp.append("C");
    //         else if (*it == 'N')
    //             rev_comp.append("N");
    //         else {
    //             std::cerr << "Invalid sequence character. Aborting.\n";
    //             exit(1);
    //         }
    //     }
    //     return rev_comp;
    // }

	std::string get_seq(int i) const {
	    if (!m_is_paired) {
	        assert (i == 0);
		    return m_seq1;
		}
	    else assert (i == 1 || i == 2);
	    if (i == 1) {
		    return m_seq1;
		}
		else {
		    return m_seq2;
		}
	}

	std::string get_phred(int i) const {
	    if (!m_is_paired) {
	        assert (i == 0);
		    return m_phred1;
		}
	    else assert (i == 1 || i == 2);
	    if (i == 1) {
		    return m_phred1;
		}
		else {
		    return m_phred2;
		}
    }

	std::string get_rev_comp(int i) const {
	    if (!m_is_paired) {
	        assert (i == 0);
		    return build_rev_comp(m_seq1);
		}
	    else assert (i == 1 || i == 2);
	    if (i == 1) {
		    return build_rev_comp(m_seq1);
		}
		else {
		    return build_rev_comp(m_seq2);
		}
	}

	std::string get_rev_phred(int i) const {
	    if (!m_is_paired) {
	        assert (i == 0);
		    std::string rev_phred(m_phred1.rbegin(), m_phred1.rend());
		    return rev_phred;
		}
	    else assert (i == 1 || i == 2);
	    if (i == 1) {
		    std::string rev_phred(m_phred1.rbegin(), m_phred1.rend());
		    return rev_phred;
		}
		else {
		    std::string rev_phred(m_phred2.rbegin(), m_phred2.rend());
		    return rev_phred;
		}
	}

    unsigned int get_len() {
        unsigned int len;
        if (m_is_paired) {
            len = m_seq1.size() + m_seq2.size();
        }
        else {
            len = m_seq1.size();
        }
        return len;
    }

    bool test_N_rate() {
        unsigned int N_count;
        unsigned int len;
        if (m_is_paired) {
            std::string total_seq = m_seq1 + m_seq2;
            N_count = std::count(total_seq.begin(), total_seq.end(), 'N');
            len = total_seq.size();
        }
        else {
            N_count = std::count(m_seq1.begin(), m_seq1.end(), 'N');
            len = m_seq1.size();
        }
        bool pass;
        if (N_count < 0.05*len) {
            pass = true;
        }
        else {
            pass = false;
        }
        return pass;
    }

    // only for superreads:
	std::list< node_id_t > get_sorted_clique(int i) const {
	    assert (m_is_super);
	    if (!m_is_paired) {
	        assert (i == 0);
	        assert (sorted_vertices1.size() > 0);
		    return sorted_vertices1;
		}
	    else {
	        assert (i == 1 || i == 2);
	        if (i == 1) {
	            assert (sorted_vertices1.size() > 0);
		        return sorted_vertices1;
		    }
		    else {
		        assert (sorted_vertices2.size() > 0);
		        return sorted_vertices2;
		    }
		}
    }

    void set_sorted_clique(int i, std::list< node_id_t > sorted_clique) {
	    assert (m_is_super);
	    if (!m_is_paired) {
	        assert (i == 0);
		    sorted_vertices1 = sorted_clique;
		}
	    else {
	        assert (i == 1 || i == 2);
	        if (i == 1) {
		        sorted_vertices1 = sorted_clique;
		    }
		    else {
		        sorted_vertices2 = sorted_clique;
		    }
		}
    }

	std::unordered_map< read_id_t, OriginalIndex > get_original_reads() const {
	    assert (m_is_super);
	    return original_read_indexes;
	}

	void set_original_reads(std::unordered_map< read_id_t, OriginalIndex > read_ids) {
	    assert (m_is_super);
	    original_read_indexes = read_ids;
    }

    void set_subread_map(std::unordered_map< node_id_t, SubreadInfo > map) {
	    assert (m_is_super);
        subreadMap = map;
    }

    std::unordered_map< node_id_t, SubreadInfo > get_subreadMap() {
        assert (m_is_super);
        return subreadMap;
    }

    SubreadInfo get_subread_info(node_id_t node) const {
	    assert (m_is_super);
        assert (!subreadMap.empty());
        if (subreadMap.find(node) == subreadMap.end()) {
            std::cout << "node " << node << " not found in subreadMap :( \n";
        }
        return subreadMap.at(node);
    }

//	std::list<int> get_read_indexes(int i) const {
//	    assert (m_is_super);
//	    if (!m_is_paired) {
//	        assert (i == 0);
//		    return m_read_indexes1;
//		}
//	    else {
//	        assert (i == 1 || i == 2);
//	        if (i == 1) {
//		        return m_read_indexes1;
//		    }
//		    else {
//		        return m_read_indexes2;
//		    }
//		}
//    }
//
//    void set_read_indexes(int i, std::list<int> indexes) {
//	    assert (m_is_super);
//	    if (!m_is_paired) {
//	        assert (i == 0);
//		    m_read_indexes1 = indexes;
//		}
//	    else {
//	        assert (i == 1 || i == 2);
//	        if (i == 1) {
//		        m_read_indexes1 = indexes;
//		    }
//		    else {
//		        m_read_indexes2 = indexes;
//		    }
//		}
//    }
};

#endif /* READ_H_ */
