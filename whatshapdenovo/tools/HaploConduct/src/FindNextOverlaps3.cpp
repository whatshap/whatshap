//============================================================================
// Name        : FindNextOverlaps3.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : An alternative algorithm for finding the new overlaps:
//               instead of reconsidering all existing edges, we check all pairs of superreads that have an ORIGINAL subread in common.
//============================================================================

#include <list>
#include <boost/dynamic_bitset.hpp>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm> // std::max, std::min

#include "SRBuilder.h"

void SRBuilder::findNextOverlaps3() {
    if (program_settings.verbose) {
        std::cout << "FindNextOverlaps3...\n";
        std::cout << single_SR_vec.size() << " " << paired_SR_vec.size() << "\n";
    }
    // build an adjacency list mapping original reads to superreads
    std::deque<Read>::iterator it;
    std::unordered_map< read_id_t, node_id_t > original_to_index;
    node_id_t index = 0;
    for (it = single_SR_vec.begin(); it != single_SR_vec.end(); it++) {
        Read* read_ptr = & (*it);
        std::unordered_map< read_id_t, OriginalIndex > original_reads = it->get_original_reads();
        for (auto read_it : original_reads) {
            read_id_t original_id = read_it.first;
            auto existing_index = original_to_index.find(original_id);
            if (existing_index == original_to_index.end()) {
                original_to_index.insert(std::make_pair(original_id, index));
                nodes_to_SR.at(index).push_back(read_ptr);
                index++;
            }
            else {
                nodes_to_SR.at(existing_index->second).push_back(read_ptr);
            }
        }
    }
    for (it = paired_SR_vec.begin(); it != paired_SR_vec.end(); it++) {
        Read* read_ptr = & (*it);
        std::unordered_map< read_id_t, OriginalIndex > original_reads = it->get_original_reads();
        for (auto read_it : original_reads) {
            read_id_t original_id = read_it.first;
            auto existing_index = original_to_index.find(original_id);
            if (existing_index == original_to_index.end()) {
                original_to_index.insert(std::make_pair(original_id, index));
                nodes_to_SR.at(index).push_back(read_ptr);
                index++;
            }
            else {
                nodes_to_SR.at(existing_index->second).push_back(read_ptr);
            }
        }
    }
    for (it = trivial_SR_vec.begin(); it != trivial_SR_vec.end(); it++) {
        Read* read_ptr = & (*it);
        std::unordered_map< read_id_t, OriginalIndex > original_reads = it->get_original_reads();
        for (auto read_it : original_reads) {
            read_id_t original_id = read_it.first;
            auto existing_index = original_to_index.find(original_id);
            if (existing_index == original_to_index.end()) {
                original_to_index.insert(std::make_pair(original_id, index));
                nodes_to_SR.at(index).push_back(read_ptr);
                index++;
            }
            else {
                nodes_to_SR.at(existing_index->second).push_back(read_ptr);
            }
        }
    }
    node_id_t SR_count = single_SR_vec.size() + paired_SR_vec.size();
    unsigned int c = 0;
    for (auto it : nodes_to_SR) {
        if (it.size() > c) {
            c = it.size();
        }
    }
    if (program_settings.verbose) {
        std::cout << "SR_count = " << SR_count << ", max #superreads per original = " << c << "\n";
    }
    nodeDictApproach(original_to_index);
}

void SRBuilder::nodeDictApproach(std::unordered_map< read_id_t, node_id_t > original_to_index) {
//    std::cout << "nodeDictApproach...\n";
    // keep track of superread edges already found
    overlaps_found = std::vector< std::set< read_id_t >> (new_read_count, std::set< read_id_t >());
    // build a list containing every pair of superreads having a subread in common
    std::list< OverlapCandidate > overlaps_list;
    node_id_t node_count = nodes_to_SR.size();
    int status = 0;
    if (program_settings.verbose) {
        std::cout << "Building list of overlapping superreads: \n";
    }
    for (auto original_it : original_to_index) {
        read_id_t original_id = original_it.first;
        node_id_t node = original_it.second; // index of original read in nodes_to_SR
        std::vector< Read* > SR_list = nodes_to_SR.at(node);
        int perc = static_cast<int>((node*100)/node_count);
        if (perc % 10 == 0 && perc > status) { // TODO: this does not always print because perc is rounded
            status = perc;
            if (program_settings.verbose) {
                std::cout << status << "%.. \n";
            }
        }
        node_id_t count = SR_list.size();
        for (node_id_t i=0; i < count; i++) {
            Read* SR1 = SR_list.at(i);
            read_id_t id1 = SR1->get_read_id();
            for (node_id_t j=i+1; j < count; j++) {
                Read* SR2 = SR_list.at(j);
                read_id_t id2 = SR2->get_read_id();
                // check if overlap was already found
                node_id_t smallest = std::min(id1, id2);
                node_id_t largest = std::max(id1, id2);
                if (overlaps_found.at(smallest).count(largest) != 0) { // overlap was found before
                    continue;
                }
                overlaps_found.at(smallest).insert(largest);
                OverlapCandidate SR_overlap;
                SR_overlap.SR1 = SR1;
                SR_overlap.SR2 = SR2;
                SR_overlap.common_node = node;
                SR_overlap.original_id = original_id;
                overlaps_list.push_back(SR_overlap);
            }
        }
    }
    if (program_settings.verbose) {
        std::cout << "100%\n";
    }
    // for every entry in this list, find the corresponding superread overlap
    std::string filename = PATH + "overlaps.txt";
    std::ofstream outfile(filename);
    int overlaps_count = 0;
    if (program_settings.verbose) {
        std::cout << "Deducing and writing overlaps: \n";
    }
    node_id_t total_count = overlaps_list.size();
    status = 0;
    for (auto overlap_candidate : overlaps_list) {
        int perc = static_cast<int>((overlaps_count*100)/total_count);
        if (perc % 10 == 0 && perc > status) {
            status = perc;
            if (program_settings.verbose) {
                std::cout << status << "%.. \n";
            }
        }
        read_id_t common_id = overlap_candidate.original_id;
        Overlap overlap = deduceOverlap(overlap_candidate, common_id);
        if (program_settings.no_inclusions && overlap.get_perc() == 100) {
            // ignore inclusion overlap
            continue;
        }
        if (overlap.get_len(1) > 0) {
            std::string line = overlap.get_overlap_line();
            outfile << line;
            overlaps_count++;
        }
    }
    if (program_settings.verbose) {
        std::cout << "100% \n";
        std::cout << "Number of overlaps found: " << overlaps_count << "\n";
    }
    next_overlaps_count += overlaps_count;
    outfile.close();
}


Overlap SRBuilder::deduceOverlap(OverlapCandidate candidate, read_id_t original_id) {
//    std::cout << "deduceOverlap\n";
    // overlap parameters to determine:
    read_id_t id1;
    read_id_t id2;
    int pos1;
    int pos2;
    std::string ord;
    std::string ori1;
    std::string ori2;
    unsigned int perc1;
    unsigned int perc2;
    int len1;
    int len2;
    std::string type1;
    std::string type2;

    Read* SR1 = candidate.SR1;
    Read* SR2 = candidate.SR2;

//    unsigned int node = candidate.common_node;
    std::unordered_map< read_id_t, OriginalIndex > originals1 = SR1->get_original_reads();
    std::unordered_map< read_id_t, OriginalIndex > originals2 = SR2->get_original_reads();
    OriginalIndex original_index_SR1 = originals1.at(original_id);
    OriginalIndex original_index_SR2 = originals2.at(original_id);
    if (!(SR1->is_paired()) && !(SR2->is_paired())) { // S-S overlap
//        int idx1 = findCliqueIndex(node, SR1, true, false);
//        int idx2 = findCliqueIndex(node, SR2, true, false);
        int idx1 = original_index_SR1.index1;
        int idx2 = original_index_SR2.index1;
        int lenA = (SR1->get_seq(0)).length();
        int lenB = (SR2->get_seq(0)).length();
        if (idx1-idx2 >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1 - idx2;
            if (pos1 > lenA) { // no overlap (due to opposite subread ends being removed)
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "s"); // this overlap will be ignored
            }
            len1 = std::min(lenA - pos1, lenB);
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2 - idx1;
            if (pos1 > lenB) { // no overlap (due to opposite subread ends being removed)
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "s"); // this overlap will be ignored
            }
            len1 = std::min(lenA, lenB - pos1);
        }
        perc1 = (int)floor(std::max(len1/float(lenA), len1/float(lenB))*100);
        if (perc1 > 100) {
            std::cout << "len1=" << len1 << ", lenA=" << lenA << ", lenB=" << lenB << ", idx1=" << idx1 << ", idx2=" << idx2 << ", pos=" << pos1 << "\n";
        }
        // the remaining parameters are redundant here
        pos2 = 0;
        ord = "-";
        ori1 = "+";
        ori2 = "+";
        perc2 = 0;
        len2 = 0;
        type1 = "s";
        type2 = "s";
    }
    else if ((SR1->is_paired()) && !(SR2->is_paired())) { // P-S overlap
        // note: this case actually cannot occur because we insert singles into SR_list before pairs
        int idx1l = original_index_SR1.index1;
        int idx1r = original_index_SR1.index2;
        int idx2l = original_index_SR2.index1;
        int idx2r = original_index_SR2.index2;
        int lenA1 = (SR1->get_seq(1)).length();
        int lenA2 = (SR1->get_seq(2)).length();
        int lenB = (SR2->get_seq(0)).length();
        if (idx1l-idx2l >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1l - idx2l;
            len1 = lenA1 - pos1;
            if (len1 <= 0) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "p", "s"); // this overlap will be ignored
            }
            type1 = "p";
            type2 = "s";
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2l - idx1l;
            len1 = std::min(lenA1, lenB - pos1);
            if (len1 <= 0) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "p", "s"); // this overlap will be ignored
            }
            type1 = "s";
            type2 = "p";
        }
        perc1 = (int)floor(len1/float((SR1->get_seq(1)).length())*100);
        pos2 = idx2r - idx1r;
        len2 = std::min(lenA2, lenB-pos2);
        if (len2 <= 0 || pos2 < 0) {
            return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "p", "s"); // this overlap will be ignored
        }
        perc2 = (int)floor(len2/float(lenA2)*100);
        ord = "-";
        ori1 = "+";
        ori2 = "+";
    }
    else if (!(SR1->is_paired()) && (SR2->is_paired())) { // S-P overlap
        int idx1l = original_index_SR1.index1;
        int idx1r = original_index_SR1.index2;
        int idx2l = original_index_SR2.index1;
        int idx2r = original_index_SR2.index2;
        int lenA = (SR1->get_seq(0)).length();
        int lenB1 = (SR2->get_seq(1)).length();
        int lenB2 = (SR2->get_seq(2)).length();
        if (idx1l-idx2l >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1l - idx2l;
            len1 = std::min(lenB1, lenA-pos1);
            if (len1 <= 0) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "p"); // this overlap will be ignored
            }
            type1 = "s";
            type2 = "p";
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2l - idx1l;
            len1 = lenB1 - pos1;
            if (len1 <= 0) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "p"); // this overlap will be ignored
            }
            type1 = "p";
            type2 = "s";
        }
        perc1 = (int)floor(len1/float((SR2->get_seq(1)).length())*100);
        pos2 = idx1r - idx2r;
        len2 = std::min(lenB2, lenA-pos2);
        if (len2 <= 0 || pos2 < 0) {
            return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "p"); // this overlap will be ignored
        }
        perc2 = (int)floor(len2/float(lenB2)*100);
        ord = "-";
        ori1 = "+";
        ori2 = "+";
    }
    else { // P-P overlap
        int idx1l = original_index_SR1.index1;
        int idx1r = original_index_SR1.index2;
        int idx2l = original_index_SR2.index1;
        int idx2r = original_index_SR2.index2;
        int lenA = (SR1->get_seq(1)).length();
        int lenB = (SR2->get_seq(1)).length();
        int lenC = (SR1->get_seq(2)).length();
        int lenD = (SR2->get_seq(2)).length();
        bool front_ord, back_ord;
        if (idx1l-idx2l >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1l - idx2l;
            len1 = std::min(lenA - pos1, lenB);
            front_ord = 1;
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2l - idx1l;
            len1 = std::min(lenA, lenB - pos1);
            front_ord = 0;
        }
        if (idx1r-idx2r >= 0) {
            pos2 = idx1r - idx2r;
            len2 = std::min(lenC - pos2, lenD);
            back_ord = 1;
        }
        else {
            pos2 = idx2r - idx1r;
            len2 = std::min(lenC, lenD - pos2);
            back_ord = 0;
        }
        if (len1 <= 0 || len2 <= 0) {
            return Overlap(0, 0, 0, 0, "1", "-", "-", 0, 0, 0, 0, "p", "p"); // this overlap will be ignored
        }
        perc1 = (int)floor(std::max(len1/float(lenA), len1/float(lenB))*100);
        perc2 = (int)floor(std::max(len2/float(lenC), len2/float(lenD))*100);
        if (perc1 > 100) {
            std::cout << "perc1 > 100: " << perc1 << std::endl;
            std::cout << "idx1l, idx2l: " << idx1l << " " << idx2l << std::endl;
            std::cout << "pos1, pos2: " << pos1 << " " << pos2 << std::endl;
            std::cout << "len1, len2: " << len1 << " " << len2 << std::endl;
            std::cout << "lenA " << lenA << std::endl;
            std::cout << "lenB " << lenB << std::endl;
            std::cout << "lenC " << lenC << std::endl;
            std::cout << "lenD " << lenD << std::endl;
        }
        if (perc2 > 100) {
            std::cout << "perc2 > 100: " << perc2 << std::endl;
            std::cout << "pos1, pos2: " << pos1 << " " << pos2 << std::endl;
            std::cout << "len1, len2: " << len1 << " " << len2 << std::endl;
            std::cout << "lenA " << lenA << std::endl;
            std::cout << "lenB " << lenB << std::endl;
            std::cout << "lenC " << lenC << std::endl;
            std::cout << "lenD " << lenD << std::endl;
        }
        if (perc1 < 0 || perc2 < 0) {
            std::cout << "perc1: " << perc1 << std::endl;
            std::cout << "perc2: " << perc2 << std::endl;
            std::cout << "pos1, pos2: " << pos1 << " " << pos2 << std::endl;
            std::cout << "len1, len2: " << len1 << " " << len2 << std::endl;
            std::cout << "lenA " << lenA << std::endl;
            std::cout << "lenB " << lenB << std::endl;
            std::cout << "lenC " << lenC << std::endl;
            std::cout << "lenD " << lenD << std::endl;
        }
        assert (perc1 >= 0 && perc1 <= 100 && perc2 >= 0 && perc2 <= 100);
        if ((front_ord && back_ord) || (!front_ord && !back_ord)) {
            ord = "1";
        }
        else {
            ord = "2";
        }
        ori1 = "+";
        ori2 = "+";
        type1 = "p";
        type2 = "p";
    }
    // write data to overlap
    Overlap overlap(id1, id2, pos1, pos2, ord, ori1, ori2, perc1, perc2, len1, len2, type1, type2);
    return overlap;
}
