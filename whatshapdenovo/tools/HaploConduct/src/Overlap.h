//============================================================================
// Name        : Overlap.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Overlap class for edge calculator
//============================================================================

#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <string>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>

#include "Types.h"

class Overlap
{
private:
	read_id_t m_id1;
	read_id_t m_id2;
	unsigned int m_pos1;
	unsigned int m_pos2;
	std::string m_ord;
	std::string m_ori1;
	std::string m_ori2;
	unsigned int m_perc1;
	unsigned int m_perc2;
	unsigned int m_len1;
	unsigned int m_len2;
	std::string m_type1;
	std::string m_type2;

public:

	Overlap(std::vector< std::string > ov_array) :
		m_id1(str_to_read_id(ov_array[0])),
		m_id2(str_to_read_id(ov_array[1])),
		m_pos1(atoi(ov_array[2].c_str())),
		m_pos2(atoi(ov_array[3].c_str())),
		m_ord(ov_array[4]),
		m_ori1(ov_array[5]),
		m_ori2(ov_array[6]),
		m_perc1(atoi(ov_array[7].c_str())),
		m_perc2(atoi(ov_array[8].c_str())),
		m_len1(atoi(ov_array[9].c_str())),
		m_len2(atoi(ov_array[10].c_str())),
		m_type1(ov_array[11]),
		m_type2(ov_array[12])
	{
	    // check for every entry if it has the right type and format:
	    if (ov_array[3] == "-") { // if pos2 irrelevant
	        m_pos2 = 0;
	        m_perc2 = 0;
	        m_len2 = 0;
	    }
	    check_id(m_id1);
	    check_id(m_id2);
	    check_pos(m_pos1);
	    check_pos(m_pos2);
	    check_ori(m_ori1);
	    check_ori(m_ori2);
	    check_perc(m_perc1);
	    check_perc(m_perc2);
	    check_len(m_len1);
	    check_len(m_len2);
	    check_type(m_type1);
	    check_type(m_type2);
	    check_ord(m_ord);
	}

	Overlap(read_id_t id1, read_id_t id2, unsigned int pos1, unsigned int pos2, std::string ord, std::string ori1, std::string ori2, unsigned int perc1, unsigned int perc2, unsigned int len1, unsigned int len2, std::string type1, std::string type2) :
		m_id1(id1),
		m_id2(id2),
		m_pos1(pos1),
		m_pos2(pos2),
		m_ord(ord),
		m_ori1(ori1),
		m_ori2(ori2),
		m_perc1(perc1),
		m_perc2(perc2),
		m_len1(len1),
		m_len2(len2),
		m_type1(type1),
		m_type2(type2)
	{
	    check_id(m_id1);
	    check_id(m_id2);
	    check_pos(m_pos1);
	    check_pos(m_pos2);
	    check_ori(m_ori1);
	    check_ori(m_ori2);
	    check_perc(m_perc1);
	    check_perc(m_perc2);
	    check_len(m_len1);
	    check_len(m_len2);
	    check_type(m_type1);
	    check_type(m_type2);
	    check_ord(m_ord);
	}

	~Overlap(void) {}

	void check_id(read_id_t id) {}

	void check_pos(int s) {
	    if (s < 0) {
	        std::cerr << "overlap.m_pos < 0; Exiting.\n";
	        exit(1);
	    }
	}

	void check_ord(std::string &s) {
	    if (s.length() != 1) {
            s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
//            std::cout << s.length() << " ";
            assert (s.length() == 1);
	    }
		assert (s == "1" || s == "2" || s == "-");
	    if (m_type1 == "s" || m_type2 == "s") { assert (s == "-"); }
	    else { assert (s == "1" || s == "2"); }
	}

	void check_ori(std::string &s) {
	    if (s.length() != 1) {
            s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
            assert (s.length() == 1);
	    }
	    if (s != "+" && s != "-") {
	        std::cerr << "overlap.m_ori not of the right format (+, -). Exiting.\n";
	        exit(1);
	    }
	}

	void check_perc(int s) {
	    if (s < 0 || s > 100) {
	        std::cerr << s << "\n";
            std::cerr << "overlap.m_perc not of the right format (0 <= perc <= 100). Exiting.\n";
	        exit(1);
	    }
	}

	void check_len(int s) {
		if (s < 0) {
            std::cerr << "overlap.m_len < 0. Exiting.\n";
	        exit(1);
	    }
	}

	void check_type(std::string &s) {
	    if (s.length() != 1) {
            s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());
            s.erase(std::remove(s.begin(), s.end(), '\t'), s.end());
            s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
	    }
	    assert (s.length() == 1);
	    if (s != "s" && s != "p") {
	        std::cerr << s << " not of the form 's' or 'p'. Exiting.\n";
	        exit(1);
	    }
	}

	read_id_t get_id(int i) const {
	    if (i == 1)
	        return m_id1;
	    else if (i == 2)
	        return m_id2;
	    else {
	        std::cerr <<  "id requested must be for 1 or 2. Exiting.\n";
	        exit(1);
	    }
	}

	int get_pos(int i) const {
	    if (i == 1)
	        return m_pos1;
	    else if (i == 2)
	        return m_pos2;
	    else {
	        std::cerr <<  "pos requested must be for 1 or 2. Exiting.\n";
	        exit(1);
	    }
	}

	std::string get_ord() const {
	    return m_ord;
	}

	std::string get_ori(int i) const {
	    if (i == 1)
	        return m_ori1;
	    else if (i == 2)
	        return m_ori2;
	    else {
	        std::cerr <<  "ori requested must be for 1 or 2. Exiting.\n";
	        exit(1);
	    }
	}

	unsigned int get_perc() const {
	    if (m_perc2 > 0) {
	        return 0.5*(m_perc1 + m_perc2);
	    }
	    else {
	        return m_perc1;
	    }
	}

	unsigned int get_len(int i) const {
	    if (i == 1)
	        return m_len1;
	    else if (i == 2)
	        return m_len2;
	    else {
	        std::cerr <<  "len requested must be for 1 or 2. Exiting.\n";
	        exit(1);
	    }
	}

	std::string get_type(int i) const {
	    if (i == 1)
	        return m_type1;
	    else if (i == 2)
	        return m_type2;
	    else {
	        std::cerr <<  "type requested must be for 1 or 2. Exiting.\n";
	        exit(1);
	    }
	}

	std::string get_overlap_line() const {
	    std::string overlap_line = std::to_string(m_id1) + "\t" + std::to_string(m_id2) + "\t" + std::to_string(m_pos1) + "\t" + std::to_string(m_pos2) + "\t" + m_ord + "\t" + m_ori1 + "\t" + m_ori2 + "\t" + std::to_string(m_perc1) + "\t" + std::to_string(m_perc2) + "\t" + std::to_string(m_len1) + "\t" + std::to_string(m_len2) + "\t" + m_type1 + "\t" + m_type2 + "\n";
	    return overlap_line;
	}
};


#endif /* OVERLAP_H_ */
