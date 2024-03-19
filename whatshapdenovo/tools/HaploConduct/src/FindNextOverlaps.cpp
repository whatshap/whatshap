//============================================================================
// Name        : FindNextOverlaps.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : An algorithm for finding the new overlaps:
//               it induces all superread edges from the existing edges and overlaps, and writes them to a new overlaps file.
//============================================================================

#include <set>
#include <list>
#include <utility>
#include <fstream>
#include <sstream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "SRBuilder.h"

// Compute overlaps to be considered at next iteration
void SRBuilder::updateOverlap(Edge edge_info, edge_count_t& copied_count, edge_count_t& u2SR_count, edge_count_t& v2SR_count, edge_count_t& SR2SR_count, std::set< std::string >& overlap_set) {
//    std::cout << "updateOverlap...\n";
    node_id_t u = edge_info.get_vertex(1);
    node_id_t v = edge_info.get_vertex(2);
    Read* read1 = edge_info.get_read(1);
    Read* read2 = edge_info.get_read(2);
    int pos1, pos2;
    char ord1, ord2;
    std::string ori1, ori2;
    if (program_settings.resolve_orientations && edge_info.get_score() == 0) { // nonedge overlap
        ori1 = (edge_info.get_ori(1) == overlap_graph->getOrientation(u)) ? "+" : "-"; // orientations may not agree with current labelling,
        ori2 = (edge_info.get_ori(2) == overlap_graph->getOrientation(v)) ? "+" : "-"; // so we allow for negative orientations
    }
    else { // existing edge or no labelling applied
        ori1 = "+"; // should agree with current labelling
        ori2 = "+";
    }
    int overlap_perc;
    int overlap_len1;
    int overlap_len2;
    std::string type1, type2;
    bool optimize = true;
    if (!visited[u] && !visited[v]) { // both vertices are not in any superread, so add edge to overlaps
//        std::cout << "u, v both not in superread\n";
        read_id_t id1 = nodes_to_new_IDs.at(u);
        read_id_t id2 = nodes_to_new_IDs.at(v);
        std::string overlap_line;
        pos1 = edge_info.get_pos(1);
        pos2 = edge_info.get_pos(2);
        overlap_line += read_id_to_str(id1) + "\t" + read_id_to_str(id2) + "\t"; // read IDs have changed!!
        overlap_line += std::to_string(pos1) + "\t" + std::to_string(pos2) + "\t";
        if (edge_info.get_ord() == '-') { overlap_line += "-\t"; }
        else if (edge_info.get_ord() == '1') { overlap_line += "1\t"; }
        else {
            assert (edge_info.get_ord() == '2');
            overlap_line += "2\t";
        }
        overlap_line += ori1 + "\t" + ori2 + "\t";
        overlap_line += std::to_string(edge_info.get_perc()) + "\t0\t";
        overlap_line += std::to_string(edge_info.get_len(1)) + "\t" + std::to_string(edge_info.get_len(2)) + "\t";
        (read1->is_paired()) ? type1="p" : type1="s";
        (read2->is_paired()) ? type2="p" : type2="s";
        overlap_line += type1 + "\t" + type2;
        if (!(program_settings.no_inclusions && edge_info.get_perc() == 100)) {
            overlap_set.insert(overlap_line);
            copied_count++;
        }
    }
    else if (!visited[u]) {
//        std::cout << "u not in superread\n";
        // consider all superreads S containing v
        std::vector< Read* > SR_list = nodes_to_SR.at(v);
        read_id_t id1 = nodes_to_new_IDs.at(u);
        for (auto it_SR : SR_list) {
            // for each S compute the overlap position of u and v in S
            read_id_t id2 = it_SR->get_read_id();
//            if ((it_SR->get_subread_info(v)).startpos1 > 0 || (it_SR->get_subread_info(v)).startpos2 > 0) {
//                continue;
//            }
            if (optimize) {
                // check if overlap was already found
                assert (id1 != id2);
                read_id_t smallest = std::min(id1, id2);
                read_id_t largest = std::max(id1, id2);
                #pragma omp flush(overlaps_found)
                if (overlaps_found.at(smallest).count(largest) != 0) { // overlap was found before
                    continue;
                }
                #pragma omp critical(insert_overlap)
                {
                    overlaps_found.at(smallest).insert(largest);
                }
            }
            bool leftside = true;
            bool second_occ = false;
            int idx2l, idx2r;
            if (it_SR->is_paired()) {
                idx2l = findCliqueIndex(v, it_SR, leftside, second_occ);
                leftside = false;
                idx2r = findCliqueIndex(v, it_SR, leftside, second_occ);
            }
            else if (read2->is_paired()) {
                idx2l = findCliqueIndex(v, it_SR, leftside, second_occ);
                second_occ = true;
                idx2r = findCliqueIndex(v, it_SR, leftside, second_occ); // find second occurrence of v!!
            }
            else {
                idx2l = findCliqueIndex(v, it_SR, leftside, second_occ);
                idx2r = idx2l;
            }
            bool success = computeOverlapData(read1, it_SR, 0, 0, idx2l, idx2r, edge_info, pos1, pos2, ord1, ord2, type1, type2, overlap_perc, overlap_len1, overlap_len2);
            if (!success) {
                continue;
            }
            // else if (pos1 == 0 && id1 > id2) {
            //     continue;
            // }
            std::string t1, t2;
            std::string overlap_line;
            if (ord1 == '1') {
                overlap_line += read_id_to_str(id1) + "\t" + read_id_to_str(id2) + "\t";
                t1 = type1;
                t2 = type2;
            }
            else {
                overlap_line += read_id_to_str(id2) + "\t" + read_id_to_str(id1) + "\t";
                t1 = type2;
                t2 = type1;
            }
            overlap_line += std::to_string(pos1) + "\t" + std::to_string(pos2) + "\t";
            if (ord2 == '-') { overlap_line += "-\t"; }
            else if (ord2 == '1') { overlap_line += "1\t"; }
            else {
                assert (ord2 == '2');
                overlap_line += "2\t";
            }
            overlap_line += ori1 + "\t" + ori2 + "\t";
            overlap_line += std::to_string(overlap_perc) + "\t0\t";
            overlap_line += std::to_string(overlap_len1) + "\t" + std::to_string(overlap_len2) + "\t";
            overlap_line += t1 + "\t" + t2;
            if (!(program_settings.no_inclusions && overlap_perc == 100)) {
                overlap_set.insert(overlap_line);
                u2SR_count++;
            }
        }
    }
    else if (!visited[v]) {
//        std::cout << "v not in superread\n";
        // consider all superreads S containing u
        std::vector< Read* > SR_list = nodes_to_SR.at(u);
        read_id_t id1 = nodes_to_new_IDs.at(v);
        for (auto it_SR : SR_list) {
            // for each S compute the overlap position of u and v in S
            read_id_t id2 = it_SR->get_read_id();
//            if ((it_SR->get_subread_info(u)).startpos1 > 0 || (it_SR->get_subread_info(u)).startpos2 > 0) {
//                continue;
//            }
            if (optimize) {
                // check if overlap was already found
                assert (id1 != id2);
                read_id_t smallest = std::min(id1, id2);
                read_id_t largest = std::max(id1, id2);
                #pragma omp flush(overlaps_found)
                if (overlaps_found.at(smallest).count(largest) != 0) { // overlap was found before
                    continue;
                }
                #pragma omp critical(insert_overlap)
                {
                    overlaps_found.at(smallest).insert(largest);
                }
            }
            int idx1l, idx1r;
            bool leftside = true;
            bool second_occ = false;
            if (it_SR->is_paired()) {
                idx1l = findCliqueIndex(u, it_SR, leftside, second_occ);
                leftside = false;
                idx1r = findCliqueIndex(u, it_SR, leftside, second_occ);
            }
            else if (read1->is_paired()) {
                idx1l = findCliqueIndex(u, it_SR, leftside, second_occ);
                second_occ = true;
                idx1r = findCliqueIndex(u, it_SR, leftside, second_occ);
            }
            else {
                idx1l = findCliqueIndex(u, it_SR, leftside, second_occ);
                idx1r = idx1l;
            }
            bool success = computeOverlapData(it_SR, read2, idx1l, idx1r, 0, 0, edge_info, pos1, pos2, ord1, ord2, type1, type2, overlap_perc, overlap_len1, overlap_len2);
            if (!success) {
                continue;
            }
            // else if (pos1 == 0 && id1 > id2) { // BAD IDEA
            //     continue;
            // }
            std::string t1, t2;
            std::string overlap_line;
            if (ord1 == '1') {
                overlap_line += read_id_to_str(id2) + "\t" + read_id_to_str(id1) + "\t";
                t1 = type1;
                t2 = type2;
            }
            else {
                overlap_line += read_id_to_str(id1) + "\t" + read_id_to_str(id2) + "\t";
                t1 = type2;
                t2 = type1;
            }
            overlap_line += std::to_string(pos1) + "\t" + std::to_string(pos2) + "\t";
            if (ord2 == '-') { overlap_line += "-\t"; }
            else if (ord2 == '1') { overlap_line += "1\t"; }
            else {
                assert (ord2 == '2');
                overlap_line += "2\t";
            }
            overlap_line += ori1 + "\t" + ori2 + "\t";
            overlap_line += std::to_string(overlap_perc) + "\t0\t";
            overlap_line += std::to_string(overlap_len1) + "\t" + std::to_string(overlap_len2) + "\t";
            overlap_line += t1 + "\t" + t2;
            if (!(program_settings.no_inclusions && overlap_perc == 100)) {
                overlap_set.insert(overlap_line);
                v2SR_count++;
            }
        }
    }
    else {
//        std::cout << "u, v both in superreads\n";
        // consider all superreads SR1 containing u
        std::vector< Read* > SR_list1 = nodes_to_SR.at(u);
        for (auto it_SR1 : SR_list1) {
            read_id_t id1 = it_SR1->get_read_id();
//            if ((it_SR1->get_subread_info(u)).startpos1 > 0 || (it_SR1->get_subread_info(u)).startpos2 > 0) {
//                continue;
//            }
            int idx1l, idx1r;
            if (it_SR1->is_paired()) {
                idx1l = findCliqueIndex(u, it_SR1, /*leftside*/ true, /*second_occ*/ false);
                idx1r = findCliqueIndex(u, it_SR1, /*leftside*/ false, /*second_occ*/ false);
            }
            else if (read1->is_paired()) {
                idx1l = findCliqueIndex(u, it_SR1, /*leftside*/ true, /*second_occ*/ false);
                idx1r = findCliqueIndex(u, it_SR1, /*leftside*/ true, /*second_occ*/ true);
            }
            else {
                idx1l = findCliqueIndex(u, it_SR1, /*leftside*/ true, /*second_occ*/ false);
                idx1r = idx1l;
            }
            // consider all superreads SR2 containing v
            std::vector< Read* > SR_list2 = nodes_to_SR.at(v);
            for (auto it_SR2 : SR_list2) {
                read_id_t id2 = it_SR2->get_read_id();
                if (id1 == id2) {
                    continue;
                }
//                if ((it_SR2->get_subread_info(v)).startpos1 > 0 || (it_SR2->get_subread_info(v)).startpos2 > 0) {
//                    continue;
//                }
                if (optimize) {
                    // check if overlap was already found
                    read_id_t smallest = std::min(id1, id2);
                    read_id_t largest = std::max(id1, id2);
                    #pragma omp flush(overlaps_found)
                    if (overlaps_found.at(smallest).count(largest) != 0) { // overlap was found before
                        continue;
                    }
                    #pragma omp critical(insert_overlap)
                    {
                        overlaps_found.at(smallest).insert(largest);
                    }
                }
                int idx2l, idx2r;
                if (it_SR2->is_paired()) {
                    idx2l = findCliqueIndex(v, it_SR2, /*leftside*/ true, /*second_occ*/ false);
                    idx2r = findCliqueIndex(v, it_SR2, /*leftside*/ false, /*second_occ*/ false);
                }
                else if (read2->is_paired()) {
                    idx2l = findCliqueIndex(v, it_SR2, /*leftside*/ true, /*second_occ*/ false);
                    idx2r = findCliqueIndex(v, it_SR2, /*leftside*/ true, /*second_occ*/ true);
                }
                else {
                    idx2l = findCliqueIndex(v, it_SR2, /*leftside*/ true, /*second_occ*/ false);
                    idx2r = idx2l;
                }
                bool success = computeOverlapData(it_SR1, it_SR2, idx1l, idx1r, idx2l, idx2r, edge_info, pos1, pos2, ord1, ord2, type1, type2, overlap_perc, overlap_len1, overlap_len2);
                if (!success) {
                    continue;
                }
                // else if (pos1 == 0 && id1 > id2) {
                //     continue;
                // }
                // else if (overlap_perc == 100) {
                //     continue;
                // }
                std::string t1, t2;
                std::string overlap_line;
                if (ord1 == '1') {
                    overlap_line += read_id_to_str(id1) + "\t" + read_id_to_str(id2) + "\t";
                    t1 = type1;
                    t2 = type2;
                }
                else {
                    overlap_line += read_id_to_str(id2) + "\t" + read_id_to_str(id1) + "\t";
                    t1 = type2;
                    t2 = type1;
                }
                overlap_line += std::to_string(pos1) + "\t" + std::to_string(pos2) + "\t";
                if (ord2 == '-') { overlap_line += "-\t"; }
                else if (ord2 == '1') { overlap_line += "1\t"; }
                else {
                    assert (ord2 == '2');
                    overlap_line += "2\t";
                }
                overlap_line += ori1 + "\t" + ori2 + "\t";
                overlap_line += std::to_string(overlap_perc) + "\t0\t";
                overlap_line += std::to_string(overlap_len1) + "\t" + std::to_string(overlap_len2) + "\t";
                overlap_line += t1 + "\t" + t2;
                if (!(program_settings.no_inclusions && overlap_perc == 100)) {
                    overlap_set.insert(overlap_line);
                    SR2SR_count++;
                }
            }
        }
    }
}


// Given a node, find the index (position) of the corresponding read in a superread
int SRBuilder::findCliqueIndex(node_id_t node, Read* superread, bool leftside, bool second_occ) {
//    std::cout << "findCliqueIndex... \n";
    SubreadInfo sub_info = superread->get_subread_info(node);
    if (leftside && !second_occ) {
        assert (sub_info.index1 >= 0 && sub_info.startpos1 >= 0);
        assert (!(sub_info.index1 > 0 && sub_info.startpos1 > 0));
        return sub_info.index1 - sub_info.startpos1;
    }
    else {
        assert (leftside || !second_occ);
        assert (sub_info.index2 >= 0 && sub_info.startpos2 >= 0);
        if (superread->is_paired()) {
            assert (!(sub_info.index2 > 0 && sub_info.startpos2 > 0)); // not true when self-overlap was merged
        }
        return sub_info.index2 - sub_info.startpos2;
    }
}


// computes overlap positions, vertex order and read type
bool SRBuilder::computeOverlapData(Read* superread1, Read* superread2, int idx1l, int idx1r, int idx2l, int idx2r, Edge edge, int &new_pos1, int &new_pos2, char &ord1, char &ord2, std::string &type1, std::string &type2, int& overlap_perc, int& overlap_len1, int& overlap_len2) {
//    std::cout << "in computeOverlapData... \n";
    int pos1 = edge.get_pos(1);
    int pos2 = edge.get_pos(2);
    int len;
    // S-S overlap
    if (!superread1->is_paired() && !superread2->is_paired()) {
        type1 = "s";
        type2 = "s";
        new_pos1 = (pos1 + idx1l) - idx2l;
        int len1 = (superread1->get_seq(0)).length();
        int len2 = (superread2->get_seq(0)).length();
        assert (len1 > 0 && len2 > 0);
        if (new_pos1 < 0) {
            ord1 = '2';
            new_pos1 = -new_pos1;
            len = len2;
        }
        else {
            ord1 = '1';
            len = len1;
        }
        overlap_len1 = std::min({len - new_pos1, len1, len2});
        overlap_len2 = 0;
        overlap_perc = (int)floor(std::max(overlap_len1/float(len1), overlap_len1/float(len2))*100);
        ord2 = '-';
        new_pos2 = 0;
        if (new_pos1 >= len) {
//            std::cout << new_pos1 << " " << len << "\n";
//            std::cout << pos1 << " " << idx1l << " " << idx2l << " " << len1 << " " << len2 << "\n";
//            std::cout << superread1->get_seq(0) << "\n";
//            assert ((superread1->get_seq(0)).size() > 0);
            return 0; // failure: too much of the read has been removed from the superread
        }
    }
    // P-S overlap
    else if (superread1->is_paired() && !superread2->is_paired()) {
        type1 = "p";
        type2 = "s";
        int len1 = (superread1->get_seq(1)).length() + (superread1->get_seq(2)).length();
        int len2 = (superread2->get_seq(0)).length();
        assert (len1 > 0 && len2 > 0);
        new_pos1 = (pos1 + idx1l) - idx2l;
        if (new_pos1 < 0) {
            ord1 = '2';
            new_pos1 = -new_pos1;
//            assert (new_pos1 < (int)(superread2->get_seq(0)).length());
            if (new_pos1 >= (int)(superread2->get_seq(0)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len1 = (superread1->get_seq(1)).length();
        }
        else {
            ord1 = '1';
//            assert (new_pos1 < (int)(superread1->get_seq(1)).length());
            if (new_pos1 >= (int)(superread1->get_seq(1)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len1 = (superread1->get_seq(1)).length() - new_pos1;
        }
        if (edge.get_ord() == '1') {
            new_pos2 = idx2r - (idx1r + pos2);
        }
        else {
            new_pos2 = (pos2 + idx2r) - idx1r;
        }
//        assert (new_pos2 < (int)(superread2->get_seq(0)).length());
        if (new_pos2 >= (int)(superread2->get_seq(0)).length()) {
            return 0; // failure: too much of the read has been removed from the superread
        }
        else if (new_pos2 < 0) {
//            std::cout << "P-S: new_pos2 < 0" << std::endl;
            return 0; // failure: paired-end read should be merged?
        }
        ord2 = '-';
        overlap_len2 = std::min((superread2->get_seq(0)).length() - new_pos2, (superread1->get_seq(2)).length());
        int total_overlap = overlap_len1 + overlap_len2;
//        overlap_len = total_overlap;
        overlap_perc = (int)floor(std::max(total_overlap/float(len1), total_overlap/float(len2))*100);
        overlap_perc = std::min(overlap_perc, 100);

        if (new_pos2 < 0) {
            std::cout << "new_pos2: " << new_pos2 << "\n";
            std::cout << "pos: " << pos1 << " " << pos2 << "\n";
            std::cout << "len1l, len1r, len2: " << (superread1->get_seq(1)).length() << " " << (superread1->get_seq(2)).length() << " " << len2 << "\n";
            std::cout << type1 << " " << type2 << "\n";
            std::cout << "indexes: " << idx1l << " " << idx1r << " " << idx2l << " " << idx2r << "\n";
        }
        assert (new_pos2 >= 0);
        if (overlap_perc < 0) {
            std::cout << type1 << type2 << "\n";
        }
    }
    // S-P overlap
    else if (!superread1->is_paired() && superread2->is_paired()) {
        type1 = "s";
        type2 = "p";
        int len1 = (superread1->get_seq(0)).length();
        int len2 = (superread2->get_seq(1)).length() + (superread2->get_seq(2)).length();
        assert (len1 > 0 && len2 > 0);
        new_pos1 = pos1 + idx1l - idx2l;
        if (new_pos1 < 0) {
            ord1 = '2';
            new_pos1 = -new_pos1;
//            assert (new_pos1 < (int)(superread2->get_seq(1)).length());
            if (new_pos1 >= (int)(superread2->get_seq(1)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len1 = (superread2->get_seq(1)).length() - new_pos1;
        }
        else {
            ord1 = '1';
//            assert (new_pos1 < (int)(superread1->get_seq(0)).length());
            if (new_pos1 >= (int)(superread1->get_seq(0)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len1 = (superread2->get_seq(1)).length();
        }
        if (edge.get_ord() == '2') {
            new_pos2 = idx1r - (pos2 + idx2r);
        }
        else {
            new_pos2 = idx1r + pos2 - idx2r;
        }
//        assert (new_pos2 < (int)(superread1->get_seq(0)).length());
        if (new_pos2 >= (int)(superread1->get_seq(0)).length()) {
            return 0; // failure: too much of the read has been removed from the superread
        }
        else if (new_pos2 < 0) {
//            std::cout << "S-P: new_pos2 < 0" << std::endl;
            return 0; // failure: paired-end read should be merged?
        }
        ord2 = '-';
        overlap_len2 = std::min((superread1->get_seq(0)).length() - new_pos2, (superread2->get_seq(2)).length());
        int total_overlap = overlap_len1 + overlap_len2;
//        overlap_len = total_overlap;
        overlap_perc = (int)floor(std::max(total_overlap/float(len1), total_overlap/float(len2))*100);
        overlap_perc = std::min(overlap_perc, 100);
    }
    // P-P overlap
    else if (superread1->is_paired() && superread2->is_paired()) {
        type1 = "p";
        type2 = "p";
        new_pos1 = (pos1 + idx1l) - idx2l;
        if (new_pos1 < 0) {
            ord1 = '2';
            new_pos1 = -new_pos1;
//            assert (new_pos1 < (int)(superread2->get_seq(1)).length());
            if (new_pos1 >= (int)(superread2->get_seq(1)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len1 = std::min((superread1->get_seq(1)).length(), (superread2->get_seq(1)).length()-new_pos1);
        }
        else {
            ord1 = '1';
//            assert (new_pos1 < (int)(superread1->get_seq(1)).length());
            if (new_pos1 >= (int)(superread1->get_seq(1)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len1 = std::min((superread1->get_seq(1)).length()-new_pos1, (superread2->get_seq(1)).length());
        }
        if (edge.get_ord() == '1') {
            new_pos2 = (pos2 + idx1r) - idx2r;
        }
        else {
            new_pos2 = idx1r - (pos2 + idx2r); // pos w.r.t. read1 -> if pos < 0 then ord = 2, UNLESS ord1=2
        }
        if (new_pos2 < 0) {
            if (ord1 == '1') {
                ord2 = '2';
            }
            else {
                ord2 = '1';
            }
            new_pos2 = -new_pos2;
//            assert (new_pos2 < (int)(superread2->get_seq(2)).length());
            if (new_pos2 >= (int)(superread2->get_seq(2)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len2 = std::min((superread1->get_seq(2)).length(), (superread2->get_seq(2)).length()-new_pos2);
        }
        else {
            if (ord1 == '1') {
                ord2 = '1';
            }
            else {
                ord2 = '2';
            }
//            assert (new_pos2 < (int)(superread1->get_seq(2)).length());
            if (new_pos2 >= (int)(superread1->get_seq(2)).length()) {
                return 0; // failure: too much of the read has been removed from the superread
            }
            overlap_len2 = std::min((superread1->get_seq(2)).length()-new_pos2, (superread2->get_seq(2)).length());
        }
        int total_overlap = overlap_len1 + overlap_len2;
//        overlap_len = total_overlap;
        int total_len1 = (superread1->get_seq(1)).length() + (superread1->get_seq(2)).length();
        int total_len2 = (superread2->get_seq(1)).length() + (superread2->get_seq(2)).length();
        overlap_perc = (int)floor(std::max(total_overlap/float(total_len1), total_overlap/float(total_len2))*100);
        overlap_perc = std::min(overlap_perc, 100);
    }
    assert (new_pos1 >= 0);
    if (new_pos2 < 0) {
        std::cout << "new_pos2: " << new_pos2 << "\n";
        std::cout << "pos: " << pos1 << " " << pos2 << "\n";
        std::cout << type1 << " " << type2 << "\n";
        std::cout << "indexes: " << idx1l << " " << idx1r << " " << idx2l << " " << idx2r << "\n";
    }
    assert (new_pos2 >= 0);
    if (overlap_perc < 0) {
        std::cout << type1 << type2 << "\n";
    }
    assert (overlap_perc >= 0 && overlap_perc <= 100);
    return 1; // success
}


void SRBuilder::processOverlaps(const std::vector<Edge>& edge_vec, edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set)
{
    if (program_settings.verbose) {
        std::cout << "processOverlaps\n";
        std::cout << "size: " << edge_vec.size() << std::endl;
    }
    edge_count_t size = edge_vec.size();
	#pragma omp parallel num_threads(N_THREADS) shared(total_copied_count, total_u2SR_count, total_v2SR_count, final_overlap_set)
	{
	    // for each thread, create a set of overlap lines not allowing duplicates
        std::set< std::string > overlap_set;
        edge_count_t copied_count = 0;
        edge_count_t u2SR_count = 0;
        edge_count_t v2SR_count = 0;
        edge_count_t SR2SR_count = 0;
	    #pragma omp for
        for (edge_count_t i=0; i < size; i++) {
            Edge edge = edge_vec.at(i);
            updateOverlap(edge, copied_count, u2SR_count, v2SR_count, SR2SR_count, overlap_set);
        }
        #pragma omp critical
        {
//            std::cout << "Merging overlap sets...\n";
            final_overlap_set.insert(overlap_set.begin(), overlap_set.end());
            total_copied_count += copied_count;
            total_u2SR_count += u2SR_count;
            total_v2SR_count += v2SR_count;
            total_SR2SR_count += SR2SR_count;
        }
    }
    if (program_settings.verbose) {
        std::cout << "overlaps: " << final_overlap_set.size() << std::endl;
    }
}


// Reconsider edge overlaps
void SRBuilder::reconsiderEdgeOverlaps(edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set) {
    if (program_settings.verbose) {
        std::cout << "Reconsider existing edges...\n";
    }
//    std::cout << "visited count: " << visited.count() << "\n";
    std::vector<Edge> edge_vec;
    const edge_count_t overlaps_per_vec = 1000000;
    for (auto it_adj_list : overlap_graph->adj_out) {
        if (edge_vec.size() > overlaps_per_vec) {
            // process the currently collected overlaps
            processOverlaps(edge_vec, total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
            edge_vec.clear(); // empty the vector
        }
        edge_vec.insert(edge_vec.end(), it_adj_list.begin(), it_adj_list.end());
    }
    if (edge_vec.size() > 0) {
        processOverlaps(edge_vec, total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
        edge_vec.clear();
    }
    // also reconsider the tips/branching edges that were removed without
    // checking for evidence (i.e. if minimum evidence = 0)
//    if (program_settings.branch_min_ev == 0) {
    if (true) {
        std::vector< Edge > additional_edges = overlap_graph->branching_edges;
        processOverlaps(additional_edges, total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
    }
}


// Reconsider non-edge overlaps
void SRBuilder::reconsiderNonedgeOverlaps(edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set) {
    if (program_settings.verbose) {
        std::cout << "Reconsider non-edge overlaps...\n";
    }
    std::string old_overlapsfile = PATH + "nonedge_overlaps.txt";
    std::ifstream overlapsfile (old_overlapsfile.c_str());
    edge_count_t i = 0;
    const edge_count_t overlaps_per_vec = 1000000;
    if (overlapsfile.is_open()) {
        std::stringstream ss;
        std::string tupleline;
        std::vector<Edge> edge_vec;
        while (getline(overlapsfile, tupleline)) {
            if (i%1000000 == 0 && i > 0 && program_settings.verbose) {
                std::cout << i << " overlaps done...\n";
            }
            i++;
            boost::trim_if(tupleline, boost::is_any_of("\t "));
            std::vector< std::string > tmp_vec;
            ss << tupleline;
            std::string tmp;
            while (getline(ss, tmp, '\t')) {
                tmp_vec.push_back(tmp);
            }
			ss << "";
			ss.clear();
            Overlap overlap(tmp_vec);
            // build a temporary edge such that we can re-use the procedure for updating edges
            read_id_t id1 = overlap.get_id(1);
	        read_id_t id2 = overlap.get_id(2);
            node_id_t index1 = (fastq_storage->m_ID_to_index).at(id1);
            node_id_t index2 = (fastq_storage->m_ID_to_index).at(id2);
            Read* read1 = (fastq_storage->m_read_vec).at(index1);
            Read* read2 = (fastq_storage->m_read_vec).at(index2);
            bool ori1 = (overlap.get_ori(1) == "+") ? true : false;
            bool ori2 = (overlap.get_ori(2) == "+") ? true : false;
            node_id_t vertex1, vertex2;
            if (program_settings.add_duplicates) {
                vertex1 = read1->get_vertex_id(ori1);
                vertex2 = read2->get_vertex_id(ori2);
//                ori1 = 1; // orientation in duplicates-graph always +
//                ori2 = 1;
            }
            else {
                vertex1 = read1->get_vertex_id(true);
                vertex2 = read2->get_vertex_id(true);
//                bool new_ori1 = (ori1 == overlap_graph->getOrientation(vertex1)); // orientation corresponding to labelled read
//                bool new_ori2 = (ori2 == overlap_graph->getOrientation(vertex2));
//                ori1 = new_ori1;
//                ori2 = new_ori2;
            }
            Edge edge(0, overlap.get_pos(1), overlap.get_pos(2), ori1, ori2, overlap.get_ord(), read1, read2);
            edge.set_perc(overlap.get_perc());
            int overlap_len1 = overlap.get_len(1);
            int overlap_len2 = overlap.get_len(2);
            edge.set_len(overlap_len1, overlap_len2);
            edge.set_vertices(vertex1, vertex2);

	        // check if there is an edge between this pair of vertices
	        if (overlap_graph->checkEdge(vertex1, vertex2, /*reverse_allowed*/ true) > 0) {
	            continue;
	        }
            edge_vec.push_back(edge);
            // when adding duplicates, also add opposite overlap
            if (program_settings.add_duplicates) {
                node_id_t v1 = read1->get_vertex_id(!ori1);
                node_id_t v2 = read2->get_vertex_id(!ori2);
                if (!(read1->is_paired()) && !(read2->is_paired())) { // S-S
                    int pos1 = (read1->get_seq(0)).size() - overlap.get_pos(1) - (read2->get_seq(0)).size();
                    if (pos1 < 0) {
                        Edge opposite_edge(0, -pos1, 0, !ori2, !ori1, overlap.get_ord(), read2, read1);
                        opposite_edge.set_vertices(v2, v1);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                    else {
                        Edge opposite_edge(0, pos1, 0, !ori1, !ori2, overlap.get_ord(), read1, read2);
                        opposite_edge.set_vertices(v1, v2);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                }
                else if ((read1->is_paired()) && !(read2->is_paired())) { // P-S
                    int pos1 = (read1->get_seq(2)).size() + overlap.get_pos(2) - (read2->get_seq(0)).size();
                    int pos2 = (read2->get_seq(0)).size() + overlap.get_pos(1) - (read1->get_seq(1)).size();
                    if (pos1 < 0) {
                        Edge opposite_edge(0, -pos1, pos2, !ori2, !ori1, overlap.get_ord(), read2, read1);
                        opposite_edge.set_vertices(v2, v1);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                    else {
                        Edge opposite_edge(0, pos1, pos2, !ori1, !ori2, overlap.get_ord(), read1, read2);
                        opposite_edge.set_vertices(v1, v2);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                }
                else if (!(read1->is_paired()) && (read2->is_paired())) { // S-P
                    int pos1 = (read1->get_seq(0)).size() - overlap.get_pos(2) - (read2->get_seq(2)).size();
                    int pos2 = (read1->get_seq(0)).size() - overlap.get_pos(1) - (read2->get_seq(1)).size();
                    if (pos1 < 0) {
                        Edge opposite_edge(0, -pos1, pos2, !ori2, !ori1, overlap.get_ord(), read2, read1);
                        opposite_edge.set_vertices(v2, v1);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                    else {
                        Edge opposite_edge(0, pos1, pos2, !ori1, !ori2, overlap.get_ord(), read1, read2);
                        opposite_edge.set_vertices(v1, v2);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                }
                else { // P-P
                    assert((read1->is_paired()) && (read2->is_paired()));
                    int pos1;
                    if (overlap.get_ord() == "1") {
                        pos1 = (read1->get_seq(2)).size() - overlap.get_pos(2) - (read2->get_seq(2)).size();
                    }
                    else {
                        assert (overlap.get_ord() == "2");
                        pos1 = (read1->get_seq(2)).size() + overlap.get_pos(2) - (read2->get_seq(2)).size();
                    }
                    int pos2 = (read1->get_seq(1)).size() - overlap.get_pos(1) - (read2->get_seq(1)).size();
                    std::string ord;
                    if (pos1 < 0) {
                        if (pos2 < 0) {
                            pos2 = -pos2;
                            ord = "1";
                        }
                        else {
                            ord = "2";
                        }
                        Edge opposite_edge(0, -pos1, pos2, !ori2, !ori1, ord, read2, read1);
                        opposite_edge.set_vertices(v2, v1);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                    else {
                        if (pos2 < 0) {
                            pos2 = -pos2;
                            ord = "2";
                        }
                        else {
                            ord = "1";
                        }
                        Edge opposite_edge(0, pos1, pos2, !ori1, !ori2, ord, read1, read2);
                        opposite_edge.set_vertices(v1, v2);
                        opposite_edge.set_perc(overlap.get_perc());
                        opposite_edge.set_len(overlap.get_len(1), overlap.get_len(2));
                        edge_vec.push_back(opposite_edge);
                    }
                }
	        }
            if (edge_vec.size() == overlaps_per_vec) {
                // process the currently collected overlaps
                processOverlaps(edge_vec, total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
                edge_vec.clear(); // empty the vector
            }
        }
        if (edge_vec.size() > 0) {
            processOverlaps(edge_vec, total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
            edge_vec.clear();
        }
        overlapsfile.close();
    }
    else {
        std::cerr << "Unable to open non-edge overlaps file\n";
        exit(1);
    }
}


void SRBuilder::findInclusionOverlaps(edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set) {
    // deduce overlaps that have been missed due to inclusions
    if (program_settings.verbose) {
        std::cout << "findInclusionOverlaps...\n";
    }
    std::vector<Edge> edge_vec;
    for (auto edge_list : overlap_graph->inclusion_edges) {
        unsigned int l = edge_list.size();
        for (unsigned int i=0; i<l; i++) {
            for (unsigned int j=i+1; j<l; j++) {
                Edge edge1 = edge_list.at(i);
                Edge edge2 = edge_list.at(j);
                double score = program_settings.edge_threshold; // overlap score
                int pos1;
                int pos2 = 0; // no PE-overlaps allowed
                bool ori1; // 1 if NORMAL, 0 if REVERSE (vertex 1)
                bool ori2; // similar for vertex 2
                Read* read1; // pointer to read corresponding to vertex1
                Read* read2; // pointer to read corresponding to vertex2
                std::string ord = "-";
                node_id_t node1;
                node_id_t node2;
                int len;
                int perc;
                // build induced edge
                if (edge1.get_vertex(1) == edge2.get_vertex(1)) {
                    continue; // irrelevant: does not induce edge u->v
                }
                else if (edge1.get_vertex(1) == edge2.get_vertex(2)) {
                    node1 = edge2.get_vertex(1);
                    node2 = edge1.get_vertex(2);
                    read1 = edge2.get_read(1);
                    read2 = edge1.get_read(2);
                    pos1 = edge2.get_pos(1);
                    ori1 = edge2.get_ori(1);
                    ori2 = edge1.get_ori(2);
                }
                else if (edge1.get_vertex(2) == edge2.get_vertex(1)) {
                    node1 = edge1.get_vertex(1);
                    node2 = edge2.get_vertex(2);
                    read1 = edge1.get_read(1);
                    read2 = edge2.get_read(2);
                    pos1 = edge1.get_pos(1);
                    ori1 = edge1.get_ori(1);
                    ori2 = edge2.get_ori(2);
                }
                else {
                    assert (edge1.get_vertex(2) == edge2.get_vertex(2));
                    continue; // irrelevant: does not induce edge u->v
                }
                // check if there are paired-end reads involved
                if (read1->is_paired() || read2->is_paired()) {
                    continue; // TODO: paired-end overlaps
                }
                len = std::min(read1->get_len()-pos1, read2->get_len());
                perc = (int)floor(100 * len / std::min(read1->get_len(), read2->get_len()));
                Edge new_edge(score, pos1, pos2, ori1, ori2, ord, read1, read2);
                new_edge.set_vertices(node1, node2);
                new_edge.set_perc(perc);
                new_edge.set_len(len, 0);
                if (overlap_graph->checkEdge(node1, node2, true) == -1) {
                    // edge doesn't exist yet
                    edge_vec.push_back(new_edge);
                }
            }
        }
    }
    if (edge_vec.size() > 0) {
        processOverlaps(edge_vec, total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
        edge_vec.clear();
    }
}


unsigned long SRBuilder::findNextOverlaps() {
    if (program_settings.verbose) {
        std::cout << "findNextOverlaps...\n";
    }
    // keep track of superread edges already found
//    std::set< read_id_t > empty_set = {};
    overlaps_found = std::vector< std::set< read_id_t >> (new_read_count, std::set< read_id_t >());
    // build an adjacency list mapping vertices superreads
    std::deque<Read>::iterator it;
    std::list<node_id_t>::const_iterator itv;
    for (it = single_SR_vec.begin(); it != single_SR_vec.end(); it++) {
        Read* read_ptr = & (*it);
        std::list< node_id_t > clique = it->get_sorted_clique(0);
        for (auto node_it : clique) {
            nodes_to_SR.at(node_it).push_back(read_ptr);
        }
    }
    for (it = paired_SR_vec.begin(); it != paired_SR_vec.end(); it++) {
        std::list< node_id_t > clique = it->get_sorted_clique(1);
        Read* read_ptr = & (*it);
        for (auto node_it : clique) {
            nodes_to_SR.at(node_it).push_back(read_ptr);
        }
    }
    edge_count_t total_copied_count = 0;
    edge_count_t total_u2SR_count = 0;
    edge_count_t total_v2SR_count = 0;
    edge_count_t total_SR2SR_count = 0;
    std::set< std::string > final_overlap_set;

    reconsiderEdgeOverlaps(total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
    if (program_settings.verbose) {
        std::cout << "Current edges have been updated...\n";
        std::cout << "Number of overlaps so far: " << final_overlap_set.size() << " of which " << total_copied_count << " are copied edges.\n";
        std::cout << "u2SR, v2SR, SR2SR: " << total_u2SR_count << " " << total_v2SR_count << " " << total_SR2SR_count << "\n";
    }
    if (!program_settings.optimize) {
        reconsiderNonedgeOverlaps(total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
        if (program_settings.verbose) {
            std::cout << "Old overlaps have been reconsidered...\n";
        }
    }
    findInclusionOverlaps(total_copied_count, total_u2SR_count, total_v2SR_count, total_SR2SR_count, final_overlap_set);
    if (program_settings.verbose) {
        std::cout << "Edges missed due to inclusions have been deduced...\n";
        std::cout << "Number of final overlaps: " << final_overlap_set.size() << std::endl;
    }
    std::set< std::string >::const_iterator it_ss;
    std::string filename = PATH + "overlaps.txt";
    std::ofstream outfile(filename);
    if (program_settings.verbose) {
        std::cout << "Number of overlap lines to write: " << final_overlap_set.size() << " of which " << total_copied_count << " are copied edges.\n";
        std::cout << "u2SR, v2SR, SR2SR: " << total_u2SR_count << " " << total_v2SR_count << " " << total_SR2SR_count << "\n";
    }
    clock_t t1, t2;
    t1 = clock();
    for (auto it_line : final_overlap_set) {
        outfile << it_line << "\n";
    }
    t2 = clock();
    if (program_settings.verbose) {
        std::cout << "Writing overlaps took " << ((float)(t2-t1))/CLOCKS_PER_SEC << " seconds.\n";
    }
    outfile.close();

    next_overlaps_count = final_overlap_set.size();

    return final_overlap_set.size();
}
