//============================================================================
// Name        : SRBuilder.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Cluster reads in the overlap graph and construct super-reads
//============================================================================

#include <unordered_map>
#include <list>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <utility>
#include <fstream>
#include <sstream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <iterator> // std::next
#include <algorithm> // std::count

#include "SRBuilder.h"
#include "Overlap.h"
#include "Edge.h"


/*  Sort vertices w.r.t. base_node: sort_vertices fills a list of sequences and their qualities to be inserted in the superread,
    while at the same time filling a list with the corresponding positions of these sequences in the superread.
    When finished, it returns the length of superread.
 */
int SRBuilder::sort_vertices(std::vector< node_id_t > vertices, char type, node_id_t base_node, std::list<int> &pos_list, std::list<std::string> &seq_list, std::list<std::string> &qual_list, std::list< node_id_t > &sorted_vertices, int thread_id)
{
//    std::cout << thread_id << " sort_vertices\n";
    assert (type == 'l'||type == 'r'||type == 's'); // type of the superread
    assert (pos_list.empty());
    assert (seq_list.empty());
    assert (qual_list.empty());
    assert (sorted_vertices.empty());

    read_id_t base_ID = (overlap_graph->vertex_to_read).at(base_node);
    Read* base_read = fastq_storage->get_read(base_ID);
    std::string base_seq;
    std::string base_qual;
//    if (base_node < fastq_storage->get_readcount()) { // forward read
    if (overlap_graph->getOrientation(base_node)) { // forward read
        if (type == 'l') {
            base_seq = base_read->get_seq(1);
            base_qual = base_read->get_phred(1);
        }
        else if (type == 'r') {
            base_seq = base_read->get_seq(2);
            base_qual = base_read->get_phred(2);
        }
        else {
            assert (type == 's');
            base_seq = base_read->get_seq(0);
            base_qual = base_read->get_phred(0);
        }
    }
    else { // reverse read
        if (type == 'l') {
            base_seq = base_read->get_rev_comp(2);
            base_qual = base_read->get_rev_phred(2);
        }
        else if (type == 'r') {
            base_seq = base_read->get_rev_comp(1);
            base_qual = base_read->get_rev_phred(1);
        }
        else {
            assert (type == 's');
            base_seq = base_read->get_rev_comp(0);
            base_qual = base_read->get_rev_phred(0);
        }
    }
    seq_list.push_back(base_seq);
    qual_list.push_back(base_qual);
    pos_list.push_back(0);
    sorted_vertices.push_back(base_node);
    int total_len = base_seq.length(); // to be updated iteratively
    int l_ext = 0;
    int r_ext = 0;

//    if (thread_id == 1) std::cout << thread_id << " starting for loop\n";
    std::vector< node_id_t >::const_iterator it1;
    for (it1 = vertices.begin(); it1 != vertices.end(); it1++) {
//        if (thread_id == 1) std::cout << thread_id << " node " << *it1 << "\n";
        node_id_t node = *it1;
        if (node == base_node)
            continue;
        Edge* edge = overlap_graph->getEdgeInfo(base_node, node); // Edge u->v if it exists; else v->u!!
        bool current_ori = overlap_graph->getOrientation(node);
//        if (node < fastq_storage->get_readcount()) {
//            current_ori = true;
//        }
//        else {
//            current_ori = false;
//        }
        read_id_t current_id;
        read_id_t id1 = (edge->get_read(1))->get_read_id();
        read_id_t id2 = (edge->get_read(2))->get_read_id();
        char ord = edge->get_ord();
        if (id1 == base_ID) {
            current_id = id2;
        }
        else {
            assert (id2 == base_ID);
            current_id = id1;
        }
        Read* current_read;
        char current_type;
        current_read = fastq_storage->get_read(current_id);
        if (type == 's') {
            if (current_read->is_paired())
                current_type = 'p';
            else
                current_type = 's';
        }
        else {
            current_type = type;
        }
        assert (current_type == 'l' || current_type == 'r' || current_type == 'p' || current_type == 's');

//        if (thread_id == 1) std::cout << thread_id << " type fixed\n";
        std::string current_seq;
        std::string current_qual;
        int new_pos;
        std::string current_seq1;
        std::string current_qual1;
        int new_pos1;
        if (current_type == 's') {
            int pos = edge->get_pos(1);
            if (current_ori) {
                current_seq = current_read->get_seq(0);
                current_qual = current_read->get_phred(0);
            }
            else {
                current_seq = current_read->get_rev_comp(0);
                current_qual = current_read->get_rev_phred(0);
            }
            if (base_ID == id1) {
                new_pos = pos;
            }
            else {
                new_pos = -pos;
            }
        }
        else if (current_type == 'l' || current_type == 'p') {
            int pos = edge->get_pos(1);
            if (current_ori) {
                current_seq = current_read->get_seq(1);
                current_qual = current_read->get_phred(1);
            }
            else {
                current_seq = current_read->get_rev_comp(2);
                current_qual = current_read->get_rev_phred(2);
            }
            if (base_ID == id1) {
                new_pos = pos;
            }
            else {
                new_pos = -pos;
            }
            if (current_type == 'p') {
                current_seq1 = current_seq;
                current_qual1 = current_qual;
                new_pos1 = new_pos;
            }
        }
        if (current_type == 'r' || current_type == 'p') {
            int pos = edge->get_pos(2);
            if (current_ori) {
                current_seq = current_read->get_seq(2);
                current_qual = current_read->get_phred(2);
            }
            else {
                current_seq = current_read->get_rev_comp(1);
                current_qual = current_read->get_rev_phred(1);
            }
            if ((current_type == 'p') || (base_ID == id1 && ord == '1') || (base_ID == id2 && ord == '2')) {
                new_pos = pos;
            }
            else {
                assert (current_type == 'r' && type == 'r');
                new_pos = -pos;
            }
        }
//        if (thread_id == 1) std::cout << thread_id << " positions fixed\n";
        std::list<int>::iterator it2;
        std::list<std::string>::iterator it3;
        std::list<std::string>::iterator it4;
        std::list<node_id_t>::iterator it5;
        if (current_type == 'p') {
            it3 = seq_list.begin();
            it4 = qual_list.begin();
            it5 = sorted_vertices.begin();
            for (it2 = pos_list.begin(); (*it2 < new_pos1) && (it2 != pos_list.end()); it2++) {
                it3++;
                it4++;
                it5++;
            }
            pos_list.insert(it2, new_pos1);
            seq_list.insert(it3, current_seq1);
            qual_list.insert(it4, current_qual1);
            sorted_vertices.insert(it5, node);
        }
//        if (thread_id == 1) std::cout << thread_id << " first insert done\n";
        it3 = seq_list.begin();
        it4 = qual_list.begin();
        it5 = sorted_vertices.begin();
        for (it2 = pos_list.begin(); (*it2 < new_pos) && (it2 != pos_list.end()); it2++) {
//            if (thread_id == 1) std::cout << thread_id << " second insert...\n";
//            if (thread_id == 1) std::cout << *it2 << " " << new_pos << "\n";
            it3++;
            it4++;
            it5++;
        }
        pos_list.insert(it2, new_pos);
        seq_list.insert(it3, current_seq);
        qual_list.insert(it4, current_qual);
        sorted_vertices.insert(it5, node);
//        if (thread_id == 1) std::cout << thread_id << " inserted in lists\n";
        // compute total length of consensus sequence
        int len1, len2;
        if (current_type == 'p') {
            assert (new_pos >= 0);
            len1 = -new_pos1;
            len2 = current_seq.length() + new_pos - base_seq.length();
            int seq1_len2 = current_seq1.length() + new_pos1 - base_seq.length();
            if (seq1_len2 > len2) {
                len2 = seq1_len2;
            }
        }
        else {
            len1 = -new_pos;
            len2 = current_seq.length() + new_pos - base_seq.length();
        }
        if (len1 > l_ext) { l_ext = len1; }
        if (len2 > r_ext) { r_ext = len2; }
//        if (thread_id == 1) std::cout << thread_id << " length checked\n";
    }
//    if (thread_id == 1) std::cout << thread_id << " finished for loop\n";
    total_len += (l_ext + r_ext);
    assert (pos_list.size() == seq_list.size());
    assert (total_len > pos_list.back());
    // make sure positions are >= 0 (these then correspond to the consensus positions)
    int min = *pos_list.begin();
    assert (min == *std::min_element(pos_list.begin(), pos_list.end()));
    if (min < 0) {
        transform(pos_list.begin(), pos_list.end(), pos_list.begin(), bind2nd(std::plus<int>(), -min));
    }
    std::list<int>::const_iterator posit = pos_list.begin();
    int c_pos = *posit;
    int n_pos;
    assert (c_pos == 0);
    posit++;
    std::list<std::string>::const_iterator seqit = seq_list.begin();
    assert ((int)(seqit->length()) <= total_len);
    seqit++;
    for (; posit != pos_list.end(); posit++) {
        n_pos = *posit;
        assert (n_pos >= 0);
        assert (c_pos <= n_pos);
        if (*posit + (int)(seqit->length()) > total_len) {
            std::cout << *posit << " " << seqit->length() << " " << total_len << "\n";
            std::cout << "pos_list" << std::endl;
            for (auto tmp1 : pos_list) {
                std::cout << tmp1 << std::endl;
            }
            std::cout << "seq_list" << std::endl;
            for (auto tmp2 : seq_list) {
                std::cout << tmp2.length() << std::endl;
            }
            std::cout << "vertices" << std::endl;
            for (auto tmp3 : vertices) {
                std::cout << tmp3 << std::endl;
            }
            std::cout << "base seq length: " << base_seq.length() << std::endl;
        }
        assert (*posit + (int)(seqit->length()) <= total_len);
        seqit++;
        c_pos = n_pos;
    }
    return total_len;
}


double SRBuilder::phred_to_prob(const int phred) {
    double P = pow(10, -phred/10.0);
    assert (P >= 0 && P <= 1);
    return P;
}


// compute a consensus nucleotide and its corresponding quality
bool SRBuilder::consensus_pos(std::string nucleotides, std::string qualities, std::string &cons_seq, std::string& cons_qual) {
//    std::cout << "consensus_pos\n";
    assert (nucleotides.length() == qualities.length());
    std::string nuc;
    std::string qual;
    double score_A = 0;
    double score_C = 0;
    double score_T = 0;
    double score_G = 0;
    unsigned int k = 0; // count overlap length excluding N's
    assert (nucleotides.length() > 0);
    assert (nucleotides.length() == qualities.length());
    for (unsigned int i = 0; i < nucleotides.length(); i++) {
        char n = nucleotides.at(i);
        char q = qualities.at(i);
        int Q = static_cast<int>(q) - 33;
        assert (Q >= 0);
        double p = phred_to_prob(Q);
        if (n == 'A') {
            score_A += log10(1-p);
            score_C += log10(p/3.0);
            score_T += log10(p/3.0);
            score_G += log10(p/3.0);
            k++;
        }
        else if (n == 'C') {
            score_C += log10(1-p);
            score_A += log10(p/3.0);
            score_T += log10(p/3.0);
            score_G += log10(p/3.0);
            k++;
        }
        else if (n == 'T') {
            score_T += log10(1-p);
            score_C += log10(p/3.0);
            score_A += log10(p/3.0);
            score_G += log10(p/3.0);
            k++;
        }
        else if (n == 'G') {
            score_G += log10(1-p);
            score_C += log10(p/3.0);
            score_T += log10(p/3.0);
            score_A += log10(p/3.0);
            k++;
        }
        else {
            assert (n == 'N');
//            cons_seq.push_back('N');
//            cons_qual.push_back('$');
//            return 1;
        }
    }
//    std::cout << score_A << " " << score_T << " " << score_C << " " << score_G << "\n";
    double max_score = std::max({score_A, score_T, score_C, score_G});
    double max_prob = std::pow(10.0, max_score);
    double total_prob = std::pow(10.0, score_A) + std::pow(10.0, score_T) + std::pow(10.0, score_C) + std::pow(10.0, score_G);
    if (max_score == 0 || total_prob == 0.0) {
        cons_seq.push_back('N');
        cons_qual.push_back('$');
        return 1;
//        return 0;
    }
    else {
        double p_incorrect = 1 - (max_prob / total_prob);
        if (nucleotides.length() > 1 && (1-p_incorrect) < minQual) { //default minQual = 0.99
//            std::cout << "p_incorrect = " << p_incorrect << "\n";
//            return 0;
            cons_seq.push_back('N');
            cons_qual.push_back('$');
            return 1;
        }
        if (p_incorrect != p_incorrect) {
            std::cout << "p_incorrect NaN\n";
            return 0;
        }
        int phred;
        if (p_incorrect < std::pow(10.0, -9.3)) {
            phred = 93;
        }
        else {
            assert (p_incorrect >= 0 && p_incorrect <= 1);
            phred = (int)round(-10 * log10(p_incorrect));
        }
        if (phred < 0) {
            std::cout << phred << ", " << max_prob/total_prob << ", " << total_prob << "\n";
            phred = 0;
        }
        else if (phred > 93) {
            std::cout << phred << ", " << max_prob/total_prob << ", " << total_prob << "\n";
            phred = 93;
        }
        char new_nuc;
        if (max_score == score_A) { new_nuc = 'A'; }
        else if (max_score == score_T) { new_nuc = 'T'; }
        else if (max_score == score_C) { new_nuc = 'C'; }
        else if (max_score == score_G) { new_nuc = 'G'; }
        else {
            std::cerr << "error computing consensus nucleotide, exiting.\n";
            exit(1);
        }
        cons_seq.push_back(new_nuc);
        cons_qual.push_back(static_cast<char>(phred + 33));
        return 1;
    }
}


// compute the whole consensus sequence and qualities by calling consensus_pos for every position
int SRBuilder::consensus(int total_len, std::list<int> &pos_list, std::list<std::string> &seq_list, std::list<std::string> &qual_list, std::string &cons_seq, std::string &cons_qual, bool subreads_needed, bool error_correction)
{
//    std::cout << "consensus\n";
    boost::dynamic_bitset<> active_nodes(pos_list.size()); // keep track of which sequences to use for computing the current consensus position, default entry is false
    std::vector<unsigned int> active_pos; // keep track of the position in each sequence
    unsigned int minimumSupport;
    if (subreads_needed) {
        minimumSupport = 2;
    }
    else {
        minimumSupport = program_settings.min_clique_size;
    }
    unsigned int currentSupport = 1;
    int trim_pos;
    if (error_correction) {
        std::list<int>::const_iterator it = pos_list.begin();
        while (currentSupport < minimumSupport and it != pos_list.end()) {
            assert (*it >= 0);
            currentSupport++;
            it++;
        }
        if (it == pos_list.end()) {
            std::cout << "Not enough support for super-read.\n";
            cons_seq = "";
            cons_qual = "";
            return -1;
        }
        trim_pos = *it;
    }
    else {
        trim_pos = 0;
    }
    assert (trim_pos >= 0);
    for (auto pos_it : pos_list) {
        if (pos_it < trim_pos) {
            active_pos.push_back(trim_pos - pos_it);
        }
        else {
            active_pos.push_back(0);
        }
    }
    int current_pos;
    int idx1 = 0;
    std::list<int>::const_iterator pos_it = pos_list.begin();
    assert (*pos_it == 0);
    bool prefix_removed = false;

    for (current_pos = 0; current_pos < total_len; current_pos++) {
//        std::cout << "current_pos, pos_it: " << current_pos << " " << *pos_it << "\n";
        while (current_pos == *pos_it && pos_it != pos_list.end()) {
            active_nodes[idx1] = true; // mark as active node
            idx1++;
            pos_it++;
        }
        std::vector<bool>::const_iterator it; // to iterate through active_nodes vector
        int idx2 = 0; // keep track of current index in seq- and qual_list;
        std::string nucleotides; // add all nucleotides to be used for consensus at this position
        std::string qualities; // add all quality scores corresponding to the nucleotides (in same order)
        std::list<std::string>::const_iterator seq_it = seq_list.begin();
        std::list<std::string>::const_iterator qual_it = qual_list.begin();
        if (error_correction && active_nodes.count() < minimumSupport) { // not enough support to compute accurate consensus
            if (pos_it == pos_list.end()) { // remove sequence suffix
                break;
            }
            else if (!prefix_removed) { // remove sequence prefix
                continue;
            }
        }
        prefix_removed = true;
        for (unsigned int node_idx = 0; node_idx < pos_list.size(); node_idx++) {
            if (active_nodes[node_idx] == true) {
                unsigned int pos = active_pos[idx2];
                if (pos >= seq_it->length() || pos >= qual_it->length()) {
                    cons_seq = "";
                    cons_qual = "";
                    return 0;
                }
                char nuc = seq_it->at(pos);
                char qual = qual_it->at(pos);
                nucleotides.push_back(nuc);
                qualities.push_back(qual);
                if (pos+1 < seq_it->length()) {
                    active_pos[idx2] = pos+1;
                }
                else {
                    active_nodes[idx2] = false;
                }
            }
            idx2++;
            seq_it++;
            qual_it++;
        }
        if (nucleotides.length() == 0) {
            cons_seq = "";
            cons_qual = "";
            return 0;
            // assert (trim_pos == 0);
            // std::cout << "current_pos, total_len, pos_it: " << current_pos << " " << total_len << " " << *pos_it << "\n";
            // std::cout << "sequence lengths: ";
            // for (std::list<std::string>::const_iterator seqit = seq_list.begin(); seqit != seq_list.end(); seqit++) {
            //     std::cout << seqit->length() << " ";
            // }
            // std::cout << "\n";
            // std::cout << "pos_list: ";
            // for (std::list<int>::const_iterator posit = pos_list.begin(); posit != pos_list.end(); posit++) {
            //     std::cout << *posit << " ";
            // }
            // std::cout << "\n";
        }
        bool result = consensus_pos(nucleotides, qualities, cons_seq, cons_qual);
        if (!result) {
            cons_seq = "";
            cons_qual = "";
            break;
        }
    }
    return trim_pos;
//    if (cons_seq != "") { // check consensus properties
//        std::list<std::string>::const_iterator seq_it = seq_list.begin();
//        for (unsigned int i = 0; i < pos_list.size(); i++) {
//            assert (seq_it != seq_list.end());
//            std::string seq = *seq_it;
//            assert (active_pos[i] == seq.length()-1);
//            seq_it++;
//        }
//        assert (cons_seq.length() == (unsigned int)total_len);
//    }
}


std::unordered_map< node_id_t, SubreadInfo > SRBuilder::calcSubreadInfo(int trim_pos1, int trim_pos2, std::list<int> pos_list1, std::list<int> pos_list2, std::list<node_id_t> sorted_vertices1, std::list<node_id_t> sorted_vertices2) {
    assert (pos_list1.size() == sorted_vertices1.size());
    assert (pos_list2.size() == sorted_vertices2.size());
    std::unordered_map< node_id_t, SubreadInfo > subread_map;
    std::list< node_id_t >::const_iterator node_it = sorted_vertices1.begin();
    for (auto pos : pos_list1) {
        std::unordered_map< node_id_t, SubreadInfo >::iterator it;
        it = subread_map.find(*node_it);
        if (it != subread_map.end()) { // left index of node already in dict, hence
            assert (trim_pos2 == -1);  // it must be a single-end superread
 //           assert (pos >= trim_pos1);
            SubreadInfo& sub_info = it->second;
            if (trim_pos1 > pos) {
                sub_info.startpos2 = trim_pos1 - pos;
                sub_info.index2 = 0;
            }
            else {
                sub_info.startpos2 = 0;
                sub_info.index2 = pos - trim_pos1;
            }
  //          sub_info.index2 = pos - trim_pos1;
  //          sub_info.startpos2 = 0;
        }
        else { // node not in dict yet so add new item
            SubreadInfo sub_info;
            if (trim_pos1 > pos) {
                sub_info.startpos1 = trim_pos1 - pos;
                sub_info.index1 = 0;
            }
            else {
                sub_info.startpos1 = 0;
                sub_info.index1 = pos - trim_pos1;
            }
            sub_info.index2 = -1; // dummy
            sub_info.startpos2 = -1; // dummy
            subread_map.insert(std::make_pair( *node_it, sub_info ));
        }
        node_it++;
    }
    if (trim_pos2 >= 0) { // this indicates a paired-end superread, hence also add the /2 subread info
        assert (pos_list1.size() == pos_list2.size());
        std::list< node_id_t >::const_iterator node_it2 = sorted_vertices2.begin();
        for (auto pos : pos_list2) {
            std::unordered_map< node_id_t, SubreadInfo >::iterator it;
            it = subread_map.find(*node_it2);
            SubreadInfo& sub_info = it->second;
            assert (it != subread_map.end()); // left index of node should already exist in dict
            if (trim_pos2 > pos) {
                sub_info.startpos2 = trim_pos2 - pos;
                sub_info.index2 = 0;
            }
            else {
                sub_info.startpos2 = 0;
                sub_info.index2 = pos - trim_pos2;
            }
            node_it2++;
        }
    }
    return subread_map;
}

void SRBuilder::filter_subreads(int num, node_id_t base_node, std::list< node_id_t > & sorted_vertices, std::list<int> & pos_list, std::list< std::string > & seq_list, std::list< std::string > & qual_list, std::list<int> & new_pos_list, std::list< std::string> & new_seq_list, std::list< std::string > & new_qual_list) {
    // filter subread set to a subset of reads
    auto pos_it = pos_list.begin();
    auto seq_it = seq_list.begin();
    auto qual_it = qual_list.begin();
    assert (sorted_vertices.size() == pos_list.size());
    std::set< node_id_t > selected_nodes;
    selected_nodes.insert(sorted_vertices.begin(), std::next(sorted_vertices.begin(), num/2)); // take leftmost set of reads
    selected_nodes.insert(base_node); // always include base node
    // sort vertices by end position
    std::vector< std::pair<node_id_t, int> > pairs;
    for (auto node_it : sorted_vertices) {
        int seq_len = static_cast<int>(seq_it->size());
        int endpos = *pos_it + seq_len;
        pairs.push_back(std::make_pair(node_it, endpos));
        pos_it++;
        seq_it++;
        qual_it++;
    }
    std::vector< node_id_t > sorted_by_endpos = sortVerticesByEndpos(pairs);
    auto last_node = sorted_by_endpos.rbegin();
    while ((int)selected_nodes.size() < num) {
        assert (last_node != sorted_by_endpos.rend());
        selected_nodes.insert(*last_node);
        last_node++;
    }
    pos_it = pos_list.begin();
    seq_it = seq_list.begin();
    qual_it = qual_list.begin();
    for (auto node_it : sorted_vertices) {
        if (selected_nodes.find(node_it) != selected_nodes.end()) { // TODO: smarter search: first sort selected nodes
            new_pos_list.push_back(*pos_it);
            new_seq_list.push_back(*seq_it);
            new_qual_list.push_back(*qual_it);
        }
        pos_it++;
        seq_it++;
        qual_it++;
    }
}


std::vector< node_id_t > SRBuilder::sortVerticesByEndpos(std::vector< std::pair<node_id_t, int> > pairs) {
    std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, int>& a, const std::pair<node_id_t, int>& b)
    {
        return a.second < b.second;
    }
    );
    std::vector< node_id_t > sorted_vertices;
    for (auto pair_it : pairs) {
        node_id_t node = pair_it.first;
        sorted_vertices.push_back(node);
    }
    return sorted_vertices;
}


Read SRBuilder::constructSuperread(std::vector< node_id_t > clique, read_id_t id, int thread_id) // construct superreads from maximal cliques
{
    assert (clique.size() > 1);
    // sort clique
    std::sort(clique.begin(), clique.end());
//    std::cout << "constructSuperread\n";
    // 1. Determine the superread type and find a base read:
    //     Superread is single-end iff any of its subreads is single-end; in this case
    //     we use a single-end subreads as base read. If the superread is paired-end,
    //     any of its subreads can be the base read.)
    //    Also build a set of original read IDs.
//    clock_t t1, t2;
//    t1 = clock();
    std::unordered_map< node_id_t, std::unordered_map< read_id_t, long> > original_reads_per_subread;
    std::vector< node_id_t >::const_iterator it;
    char superread_type = 'p';
    node_id_t base_node = clique.at(0);
    // Find a base node: any single-end read, or if there is none, then any paired-end read
    for (it = clique.begin(); it != clique.end(); it++) {
        read_id_t ID = (overlap_graph->vertex_to_read).at(*it);
        Read* read = fastq_storage->get_read(ID);
        if (superread_type == 'p' && !read->is_paired()) {
            base_node = *it;
            superread_type = 's';
        }
    }
//    t2 = clock();
//    std::cout << "Step 1 took " << ((float)(t2-t1))/CLOCKS_PER_SEC << " seconds.\n";

    // 2. Order the vertices from left to right
//    t1 = clock();
    std::list<int> pos_list1, pos_list2;
    int len1 = 0;
    int len2 = 0; // consensus sequence total length
    std::list<std::string> seq_list1, seq_list2;
    std::list<std::string> qual_list1, qual_list2;
    std::list<node_id_t> sorted_vertices1, sorted_vertices2;
    if (superread_type == 'p') {
        len1 = sort_vertices(clique, 'l', clique[0], pos_list1, seq_list1, qual_list1, sorted_vertices1, thread_id);
        len2 = sort_vertices(clique, 'r', clique[0], pos_list2, seq_list2, qual_list2, sorted_vertices2, thread_id);
        assert (pos_list1.size() == pos_list2.size());
    }
    else {
        len1 = sort_vertices(clique, 's', base_node, pos_list1, seq_list1, qual_list1, sorted_vertices1, thread_id);
    }
//    t2 = clock();
//    std::cout << "Step 2 took " << ((float)(t2-t1))/CLOCKS_PER_SEC << " seconds.\n";
//        std::cout << pos_list1.size() << " " << seq_list1.size() << " " << qual_list1.size() << "\n";
    assert (pos_list1.size() == seq_list1.size());
    assert (seq_list1.size() == qual_list1.size());
    assert (pos_list2.size() == seq_list2.size());
    assert (seq_list2.size() == qual_list2.size());

    // 3. For each position in consensus, compute its nucleotide and quality.
    //      At the first iteration, error correction is performed, i.e. we only compute a
    //      consensus nuc if the support is at least the minimum clique size. In other words,
    //      the read ends are removed because those have not been error corrected.
    //      Note that this also could cause a single-end superread to be split into a
    //      paired-end superread if there is not enough support for the internal segment.
    //      However, for now we ignore this and only check the dangling ends.
    //      TODO: check internal segment support
//    t1 = clock();
    std::string cons_seq1, cons_seq2;
    std::string cons_qual1, cons_qual2;
    int trim_pos1, trim_pos2;
    unsigned int min_clique_size = program_settings.min_clique_size;
    bool subreads_needed = false;
    if (clique.size() > 3 * min_clique_size) { // clique too large, filter clique down to size 2*min_clique_size
        std::list<int> new_pos_list1, new_pos_list2;
        std::list<std::string> new_seq_list1, new_seq_list2;
        std::list<std::string> new_qual_list1, new_qual_list2;
        filter_subreads(2*min_clique_size, base_node, sorted_vertices1, pos_list1, seq_list1, qual_list1, new_pos_list1, new_seq_list1, new_qual_list1);
        trim_pos1 = consensus(len1, new_pos_list1, new_seq_list1, new_qual_list1, cons_seq1, cons_qual1, subreads_needed, program_settings.error_correction);
        if (superread_type == 'p') {
            filter_subreads(2*min_clique_size, base_node, sorted_vertices2, pos_list2, seq_list2, qual_list2, new_pos_list2, new_seq_list2, new_qual_list2);
            trim_pos2 = consensus(len2, new_pos_list2, new_seq_list2, new_qual_list2, cons_seq2, cons_qual2, subreads_needed, program_settings.error_correction);
        }
        else {
            trim_pos2 = -1;
        }
    }
    else {
        // if (clique.size() < min_clique_size) { // min_clique_size only satisfied when considering subreads
        //     subreads_needed = true;
        // }
        trim_pos1 = consensus(len1, pos_list1, seq_list1, qual_list1, cons_seq1, cons_qual1, subreads_needed, program_settings.error_correction);
        if (superread_type == 'p') {
            trim_pos2 = consensus(len2, pos_list2, seq_list2, qual_list2, cons_seq2, cons_qual2, subreads_needed, program_settings.error_correction);
        }
        else {
            trim_pos2 = -1;
        }
    }
    // Add subread information corresponding to the consensus read trimming
    std::unordered_map< node_id_t, SubreadInfo > subreads_map = calcSubreadInfo(trim_pos1, trim_pos2, pos_list1, pos_list2, sorted_vertices1, sorted_vertices2);

    std::unordered_map< read_id_t, OriginalIndex > original_reads;
    // update existing subread (originals) indexes
    for (auto node_it : clique) {
        read_id_t subID = (overlap_graph->vertex_to_read).at(node_it);
        Read* subread_ptr = fastq_storage->get_read(subID);
        bool forward = (overlap_graph->getOrientation(node_it));
        std::unordered_map< read_id_t, OriginalIndex > subreads = (overlap_graph->original_ID_dict).at(subID);
        SubreadInfo sub_info = subreads_map.at(node_it);
        int idx1 = sub_info.index1 - sub_info.startpos1;
        int idx2 = sub_info.index2 - sub_info.startpos2;
        for (auto it : subreads) {
            read_id_t original_ID = it.first;
            if (original_reads.find(original_ID) != original_reads.end()) {
                continue; // subread already inserted by other node
            }
            OriginalIndex original_index = it.second;
            original_index.forward = (original_index.forward == forward);
            if (program_settings.first_it) {
                original_index.index1 = idx1;
                if (original_index.is_paired) {
                    original_index.index2 = idx2;
                }
            }
            else if (forward) {
                original_index.index1 += idx1;
                if (original_index.is_paired) {
                    if (sub_info.index2 >= 0) {
                        original_index.index2 += idx2;
                    }
                    else {
                        original_index.index2 += idx1;
                    }
                }
            }
            else {
                if (original_index.is_paired) {
                    if (subread_ptr->is_paired()) {
                        original_index.index1 = subread_ptr->get_seq(1).size() + idx1 - (original_index.len1 + original_index.index1);
                        if (len2 > 0 || sub_info.index2 >= 0) {
                            original_index.index2 = subread_ptr->get_seq(2).size() + idx2 - (original_index.len2 + original_index.index2);
                        }
                        else {
                            original_index.index2 = subread_ptr->get_seq(2).size() + idx1 - (original_index.len2 + original_index.index2);
                        }
                    }
                    else {
                        original_index.index1 = subread_ptr->get_seq(0).size() + idx1 - (original_index.len1 + original_index.index1);
                        original_index.index2 = subread_ptr->get_seq(0).size() + idx1 - (original_index.len2 + original_index.index2);
                    }
                }
                else {
                    original_index.index1 = subread_ptr->get_seq(0).size() + idx1 - (original_index.len1 + original_index.index1);

                }
            }
            original_reads.insert(std::make_pair(original_ID, original_index));
        }
        // if (!program_settings.first_it) { // update existing subread (originals) indexes
        //     std::unordered_map< read_id_t, OriginalIndex > subreads = (overlap_graph->original_ID_dict).at(ID);
        //     SubreadInfo sub_info = subreads_map.at(node_it);
        //     int idx1 = sub_info.index1 - sub_info.startpos1;
        //     int idx2 = sub_info.index2 - sub_info.startpos2;
        //     for (auto it : subreads) {
        //         read_id_t original_ID = it.first;
        //         OriginalIndex original_index = it.second;
        //         original_index.index1 += idx1;
        //         original_index.forward = (original_index.forward == forward);
        //         if (original_index.is_paired) {
        //             if (sub_info.index2 >= 0) {
        //                 original_index.index2 += idx2;
        //             }
        //             else {
        //                 original_index.index2 += idx1;
        //             }
        //         }
        //         original_reads.insert(std::make_pair(original_ID, original_index));
        //     }
        // }
        // else { // create original subread index information
        //     SubreadInfo sub_info = subreads_map.at(node_it);
        //     OriginalIndex original_index;
        //     if (sub_info.index2 == -1) { // single-end read
        //         original_index.is_paired = 0;
        //         original_index.index1 = sub_info.index1 - sub_info.startpos1;
        //     }
        //     else { // paired-end read
        //         assert (sub_info.index2 >= 0);
        //         original_index.is_paired = 1;
        //         original_index.index1 = sub_info.index1 - sub_info.startpos1;
        //         original_index.index2 = sub_info.index2 - sub_info.startpos2;
        //     }
        //     original_reads.insert(std::make_pair(ID, original_index)); // at first iteration all subreads are trivial
        // }
    }

//    t2 = clock();
//    std::cout << "Step 3 took " << ((float)(t2-t1))/CLOCKS_PER_SEC << " seconds.\n";

    // 4. Construct a superread and append to vector for storage
    bool is_super = true;
    if (superread_type == 'p') {
        bool is_paired = true;
        Read superread(is_paired, is_super, id, cons_seq1, cons_seq2, cons_qual1, cons_qual2);
        superread.set_sorted_clique(1, sorted_vertices1);
        superread.set_sorted_clique(2, sorted_vertices2);
//        superread.set_read_indexes(1, pos_list1);
//        superread.set_read_indexes(2, pos_list2);
        superread.set_subread_map(subreads_map);
        superread.set_original_reads(original_reads);
        return superread;
    }
    else {
        bool is_paired = false;
        Read superread(is_paired, is_super, id, cons_seq1, "", cons_qual1, "");
        superread.set_sorted_clique(0, sorted_vertices1);
//        superread.set_read_indexes(0, pos_list1);
        superread.set_subread_map(subreads_map);
        superread.set_original_reads(original_reads);
        return superread;
    }
}

Read SRBuilder::merge_self_overlap(Read superread, EdgeCalculator & edge_calculator) {
    int min_overlap = 15;
    double min_score = 0.99;
    std::string seq1 = superread.get_seq(1);
    std::string seq2 = superread.get_seq(2);
    std::string qual1 = superread.get_phred(1);
    std::string qual2 = superread.get_phred(2);
    int max_pos = seq1.length() - min_overlap;
    for (int pos = 0; pos < max_pos; pos++) {
        // start with smallest allowed overlap, since this is most likely
        int overlap_pos = seq1.length() - min_overlap - pos;
        assert (overlap_pos >= 0);
        assert (overlap_pos < static_cast<int>(seq1.length()));
        double mismatch_rate;
        double score = edge_calculator.overlap_score(seq1, seq2, qual1, qual2, overlap_pos, mismatch_rate);
        // once we find an overlap that is good enough, stop searching and merge
        if (score > min_score) {
            // prepare the input for computing the merged read
            int total_len = seq2.length() + overlap_pos;
            std::string cons_seq;
            std::string cons_qual;
            std::list< std::string > seq_list;
            seq_list.push_back(seq1);
            seq_list.push_back(seq2);
            std::list< std::string > qual_list;
            qual_list.push_back(qual1);
            qual_list.push_back(qual2);
            std::list< int > pos_list;
            pos_list.push_back(0);
            pos_list.push_back(overlap_pos);
            // find the consensus sequence and quality scores
            consensus(total_len, pos_list, seq_list, qual_list, cons_seq, cons_qual, /*subreads_needed*/ false, /*error_correction*/ false);
            if (cons_seq != "") {
                // now build the new superread
                bool is_paired = false;
                bool is_super = true;
                read_id_t id = superread.get_read_id();
                Read merged_superread(is_paired, is_super, id, cons_seq, "", cons_qual, "");
                // reformat the clique and subread info
                std::unordered_map< node_id_t, SubreadInfo > current_subreadMap = superread.get_subreadMap();
                std::unordered_map< node_id_t, SubreadInfo > new_subreadMap;
                std::list< node_id_t > new_sorted_clique;
                std::vector< std::pair< node_id_t, int > > pairs; // use (node_id, index) pairs for sorting clique nodes
                for (auto subread_it : current_subreadMap) {
                    SubreadInfo subread_info = subread_it.second;
                    pairs.push_back(std::make_pair(subread_it.first, subread_info.index1));
                    if (subread_info.index2 >= 0) {
                        int overlap_len = seq1.length() - overlap_pos;
                        int new_index2 = subread_info.index2 + seq1.length() - overlap_len;
                        subread_info.index2 = new_index2;
                        pairs.push_back(std::make_pair(subread_it.first, new_index2));
                    }
                    new_subreadMap.insert(std::make_pair(subread_it.first, subread_info));
                }

                std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, int>& a, const std::pair<node_id_t, int>& b)
                { // sort clique nodes (per read end) by index in super-read
                    return a.second < b.second;
                }
                );
                for (auto pair_it : pairs) {
                    new_sorted_clique.push_back(pair_it.first);
                }
                merged_superread.set_sorted_clique(0, new_sorted_clique);
                merged_superread.set_subread_map(new_subreadMap);
                // reformat the original read indexes in super-read
                std::unordered_map< read_id_t, OriginalIndex > current_originalsMap = superread.get_original_reads();
                std::unordered_map< read_id_t, OriginalIndex > new_originalsMap;
                for (auto read_it : current_originalsMap) {
                    OriginalIndex original_index = read_it.second;
                    if (original_index.is_paired) {
                        long idx2 = original_index.index2;
                        int overlap_len = seq1.length() - overlap_pos;
                        original_index.index2 = idx2 + seq1.length() - overlap_len;
                    }
                    new_originalsMap.insert(std::make_pair(read_it.first, original_index));
                }
                merged_superread.set_original_reads(new_originalsMap);
                return merged_superread;
            }
        }
    }
    return superread;
}


node_id_t SRBuilder::process_cliques(const std::vector< std::vector<node_id_t> >& clique_vec, read_id_t& count)
{
//    std::cout << "process_cliques\n";
    std::vector<Read> pairs;
    std::vector<Read> singles;
    node_id_t SR_count = clique_vec.size();
    if (program_settings.verbose) {
        std::cout << SR_count << " superread computations to do.\n";
    }
    clock_t t1, t2;
    t1 = clock();

    EdgeCalculator edge_calculator(fastq_storage, overlap_graph, program_settings); // for merging self-overlaps

	#pragma omp parallel num_threads(N_THREADS) shared(clique_vec)
	{
	    int tid = omp_get_thread_num();
        std::vector<Read> pairs_this_thread;
        std::vector<Read> singles_this_thread;
	    #pragma omp for
	    for (node_id_t i = 0; i < SR_count; i++)
	    {
            std::vector<node_id_t> clique = clique_vec.at(i);
            Read superread = constructSuperread(clique, 0, tid);
            if (superread.is_paired()) {
                if (superread.get_seq(1) != "" && superread.get_seq(2) != "") {
                    // test for self-overlap
                    Read new_superread = merge_self_overlap(superread, edge_calculator);
                    if (new_superread.test_N_rate() == false) {
                        // percentage of ambiguous base calls ('N's) too high
                        continue;
                    }
                    if (new_superread.is_paired()) {
                        pairs_this_thread.push_back(new_superread);
                    }
                    else {
                        singles_this_thread.push_back(new_superread);
                    }
                }
            }
            else {
                if (superread.get_seq(0) != "" && superread.test_N_rate() == true) {
                    singles_this_thread.push_back(superread);
                }
            }
        }
        #pragma omp critical(write_singles)
        {
            singles.insert(singles.end(), singles_this_thread.begin(), singles_this_thread.end());
        }
        #pragma omp critical(write_pairs)
        {
            pairs.insert(pairs.end(), pairs_this_thread.begin(), pairs_this_thread.end());
        }
    }
    t2 = clock();
    if (program_settings.verbose) {
        std::cout << "Superread computations took " << ((float)(t2-t1))/CLOCKS_PER_SEC << " seconds.\n"; // Note: this is N_THREADS*time
        std::cout << "Writing results per thread to single_SR_vec / paired_SR_vec...\n";
    }
//    std::cout << "m_largest_read_id " << fastq_storage->m_largest_read_id << "\n";
    node_id_t size_singles = singles.size();
    node_id_t size_pairs = pairs.size();
//    writeSinglesToFile(singles, count);
//    writePairsToFile(pairs, count);
    single_SR_vec.insert(single_SR_vec.end(), singles.begin(), singles.end());
    paired_SR_vec.insert(paired_SR_vec.end(), pairs.begin(), pairs.end());
    if (program_settings.verbose) {
        std::cout << "Current number of superreads (single, paired): " << size_singles << " " << size_pairs << "\n";
    }
    return size_singles + size_pairs;
}

void SRBuilder::cliquesToSuperreads() // construct superreads from maximal cliques
{
    // clear superread fastq files
    std::string filename0 = PATH + "singles.fastq";
    std::string filename1 = PATH + "paired1.fastq";
    std::string filename2 = PATH + "paired2.fastq";
    std::string originals = PATH + "subreads.txt";
    remove(filename0.c_str());
    remove(filename1.c_str());
    remove(filename2.c_str());
    remove(originals.c_str());

	read_id_t count = 0;
	int singleton_count = 0;

	std::ifstream cliquefile (PATH + "cliques.txt");
    const unsigned long int cliques_per_vec = 10000000;
    if (cliquefile.is_open()) {
        std::stringstream ss;
        std::string cliqueline;
        std::vector< std::vector< node_id_t> > SR_clique_vec;
        node_id_t totalcount = 0;
        node_id_t SRcount = 0;
        boost::dynamic_bitset<> bitvec(overlap_graph->getVertexCount());
        boost::dynamic_bitset<> used_nodes(overlap_graph->getVertexCount());
        while (getline(cliquefile, cliqueline)) {
            clique_count++;
            std::vector< node_id_t> clique;
            std::istringstream iss(cliqueline);
            node_id_t node;
            while (iss >> node) { // this automatically filters out comment lines because they can't be written to int
                clique.push_back(node);
            }

            if (program_settings.remove_multi_occ) { // filter out nodes that have occurred in previous cliques
                std::vector< node_id_t > filtered_clique;
                for (auto node_it : clique) {
                    if (used_nodes[node_it] == 0) {
                        filtered_clique.push_back(node_it);
                    }
                }
                clique = filtered_clique;
            }

            if (clique.size() == 1) {
                singleton_count++;
            }
            else if (clique.size() >= minCliqueSize) {
                totalcount++;
                SR_clique_vec.push_back(clique);
                for (auto node_it : clique) {
                    used_nodes[node_it] = 1;
                }
            }
//            else { // take original subreads into account for computing clique size
//                unsigned int clique_size_originals = 0;
//                for (auto node_it : clique) {
//                    read_id_t ID = (overlap_graph->vertex_to_read).at(node_it);
//                    std::unordered_map< read_id_t, OriginalIndex > originals = original_ID_dict.at(ID);
//                    assert (originals.size() > 0);
//                    clique_size_originals += originals.size();
//                }
//                if (clique_size_originals >= minCliqueSize) {
//                    totalcount++;
//                    SR_clique_vec.push_back(clique);;
//                    for (auto node_it : clique) {
//                        used_nodes[node_it] = 1;
//                    }
//                }
//            }

            if (SR_clique_vec.size() == cliques_per_vec) {
                SRcount += process_cliques(SR_clique_vec, count); // process the currently collected overlaps
                SR_clique_vec.clear(); // empty the vector
            }
        }
        if (SR_clique_vec.size() > 0) {
            SRcount += process_cliques(SR_clique_vec, count);
            SR_clique_vec.clear();
        }
        cliquefile.close();
        if (program_settings.verbose) {
            std::cout << "Total number of cliques considered: " << totalcount << "\n";
            std::cout << "Number of superreads constructed: " << SRcount << "\n";
            std::cout << "Number of cliques dismissed: " << totalcount - SRcount << "\n";
            std::cout << "Number of size one cliques (singletons): " << singleton_count << "\n";
        }
        SR_singles_count = single_SR_vec.size();
        SR_paired_count = paired_SR_vec.size();
        SR_trivials_count = trivial_SR_vec.size();

        // mark vertices that appear in superreads
        std::deque<Read>::const_iterator it;
        std::list<node_id_t>::const_iterator itv;
        for (it = single_SR_vec.begin(); it != single_SR_vec.end(); it++) {
            std::list< node_id_t > clique = it->get_sorted_clique(0);
            for (auto node : clique) {
                bitvec[node] = 1;
            }
        }
        for (it = paired_SR_vec.begin(); it != paired_SR_vec.end(); it++) {
            std::list< node_id_t > clique = it->get_sorted_clique(1);
            for (auto node : clique) {
                bitvec[node] = 1;
            }
        }
        // Save bitvec as visited of class
        visited = bitvec;

        // Write single-end super-reads to file
        count = 0;
        writeSinglesToFile(single_SR_vec, count);

        // Run through bitvector of nodes: if not visited, add to trivial_SR_vec
        for (node_id_t v = 0; v < (overlap_graph->getVertexCount()); v++) {
            if (visited[v] == 0) {
                read_id_t ID = (overlap_graph->vertex_to_read).at(v);
                Read* read = fastq_storage->get_read(ID);
                if (read->get_len() < program_settings.keep_singletons) { // too short to add as a trivial superread
                    visited[v] = 1; // by marking as visited we avoid further processing:
                                       // its superreads will be considered for new overlaps,
                                       // but since there are none also no overlaps will be added.
                    continue;
                }
                else if (read->test_N_rate() == false) { // too much N's in this read
                    visited[v] = 1; // by marking as visited we avoid further processing:
                                       // its superreads will be considered for new overlaps,
                                       // but since there are none also no overlaps will be added.
                    continue;
                }
                std::unordered_map< read_id_t, OriginalIndex > subreads;
                subreads = (overlap_graph->original_ID_dict).at(ID);
                // if (program_settings.first_it) {
                //     OriginalIndex original_index;
                //     original_index.index1 = 0;
                //     if (read->is_paired()) {
                //         original_index.is_paired = 1;
                //         original_index.index2 = 0;
                //     }
                //     else {
                //         original_index.is_paired = 0;
                //     }
                //     subreads.insert(std::make_pair(ID, original_index)); // trivial subread
                // }
                // else {
                //     subreads = (overlap_graph->original_ID_dict).at(ID);
                // }
                assert (subreads.size() > 0);
//                if (v == read->get_vertex_id(/*normal*/ true)) {
                if (overlap_graph->getOrientation(v)) { // forward read
                    read->set_read_id(count); // update read ID
                    read->set_super(true);
                    read->set_original_reads(subreads);
                    trivial_SR_vec.push_back(*read);
                }
                else {
                    // reverse reads are only allowed at the first iteration, so we copy the read into a new (forward) read
//                    assert (v == read->get_vertex_id(/*normal*/ false));
                    std::unordered_map< read_id_t, OriginalIndex > updated_subreads;
                    if (read->is_paired()) {
                        Read rev_read(true, true, count, read->get_rev_comp(2), read->get_rev_comp(1), read->get_rev_phred(2), read->get_rev_phred(1));
                        // update subread indexes because of reversal
                        for (auto sub_it : subreads) {
                            OriginalIndex original_index = sub_it.second;
                            original_index.forward = !original_index.forward;
                            original_index.index1 = read->get_seq(1).size() - (original_index.index1 + original_index.len1);
                            original_index.index2 = read->get_seq(2).size() - (original_index.index2 + original_index.len2);
                            updated_subreads.insert(std::make_pair(sub_it.first, original_index));
                        }
                        rev_read.set_original_reads(updated_subreads);
                        trivial_SR_vec.push_back(rev_read);
                    }
                    else {
                        Read rev_read(false, true, count, read->get_rev_comp(0), "", read->get_rev_phred(0), "");
                        // update subread indexes because of reversal
                        for (auto sub_it : subreads) {
                            OriginalIndex original_index = sub_it.second;
                            original_index.forward = !original_index.forward;
                            original_index.index1 = read->get_seq(0).size() - (original_index.index1 + original_index.len1);
                            if (original_index.is_paired) {
                                original_index.index2 = read->get_seq(0).size() - (original_index.index2 + original_index.len2);
                            }
                            updated_subreads.insert(std::make_pair(sub_it.first, original_index));
                        }
                        rev_read.set_original_reads(updated_subreads);
                        trivial_SR_vec.push_back(rev_read);
                    }
                }
                nodes_to_new_IDs.insert(std::pair<node_id_t, read_id_t>(v, count)); // create a map from old IDs to new IDs
                count++;
            }
        }
    }
    else {
        std::cerr << "Unable to open clique-file";
        exit(1);
    }
    if (program_settings.verbose) {
        std::cout << "Number of trivial superreads: " << trivial_SR_vec.size() << "\n";
    }
    writeTrivialsToFile();
    // finally write paired-end super-reads to file
    writePairsToFile(paired_SR_vec, count);
    new_read_count = count;
}


void SRBuilder::mergeAlongEdges() // construct superreads from high quality edges
{
    // clear superread fastq files
    std::string filename0 = PATH + "singles.fastq";
    std::string filename1 = PATH + "paired1.fastq";
    std::string filename2 = PATH + "paired2.fastq";
    std::string originals = PATH + "subreads.txt";
    remove(filename0.c_str());
    remove(filename1.c_str());
    remove(filename2.c_str());
    remove(originals.c_str());

	read_id_t count = 0;
	std::vector< std::vector< node_id_t > > SR_merge_vec = overlap_graph->getEdgesForMerging();
    node_id_t totalcount = SR_merge_vec.size();
    node_id_t SRcount = process_cliques(SR_merge_vec, count);
    SR_merge_vec.clear();
    if (program_settings.verbose) {
        std::cout << "Total number of edges considered: " << totalcount << "\n";
        std::cout << "Number of superreads constructed: " << SRcount << "\n";
        std::cout << "Number of superreads dismissed due to quality score: " << totalcount - SRcount << "\n";
    }
    // mark vertices that appear in superreads
    boost::dynamic_bitset<> bitvec(overlap_graph->getVertexCount());
    std::deque<Read>::const_iterator it;
    std::list<node_id_t>::const_iterator itv;
    for (it = single_SR_vec.begin(); it != single_SR_vec.end(); it++) {
        std::list< node_id_t > clique = it->get_sorted_clique(0);
        for (auto node : clique) {
            bitvec[node] = 1;
        }
    }
    for (it = paired_SR_vec.begin(); it != paired_SR_vec.end(); it++) {
        std::list< node_id_t > clique = it->get_sorted_clique(1);
        for (auto node : clique) {
            bitvec[node] = 1;
        }
    }
    // Save bitvec as visited of class
    visited = bitvec;
    // write single-end super-reads to file
    count = 0;
    writeSinglesToFile(single_SR_vec, count);
    // Run through bitvector of nodes: if not visited, add to trivial_SR_vec
    for (node_id_t v = 0; v < (overlap_graph->getVertexCount()); v++) {
        if (visited[v] == 0) {
            read_id_t ID = (overlap_graph->vertex_to_read).at(v);
            Read* read = fastq_storage->get_read(ID);
            if (read->get_len() < program_settings.keep_singletons) { // too short to add as a trivial superread
                visited[v] = 1; // by marking as visited we avoid further processing:
                                   // its superreads will be considered for new overlaps,
                                   // but since there are none also no overlaps will be added.
                continue;
            }
            else if (read->test_N_rate() == false) { // too much N's in this read
                visited[v] = 1; // by marking as visited we avoid further processing:
                                   // its superreads will be considered for new overlaps,
                                   // but since there are none also no overlaps will be added.
                continue;
            }
            else if (program_settings.ignore_inclusions && overlap_graph->inclusions[v] == 1) {
                // corresponds to tip node in overlap graph
                // store separately and eventually write sequences to fastq
                visited[v] = 1;
                tips_vec.push_back(*read);
                continue;
            }
            else if (read->is_tip() && program_settings.store_tips_separately) {
                // corresponds to tip node in overlap graph
                // store separately and eventually write sequences to fastq
                visited[v] = 1;
                tips_vec.push_back(*read);
                continue;
            }
            std::unordered_map< read_id_t, OriginalIndex > subreads;
            subreads = (overlap_graph->original_ID_dict).at(ID);
            // if (program_settings.first_it) {
            //     OriginalIndex original_index;
            //     original_index.index1 = 0;
            //     if (read->is_paired()) {
            //         original_index.is_paired = 1;
            //         original_index.index2 = 0;
            //     }
            //     else {
            //         original_index.is_paired = 0;
            //     }
            //     subreads.insert(std::make_pair(ID, original_index)); // trivial subread
            // }
            // else {
            //     subreads = (overlap_graph->original_ID_dict).at(ID);
            // }
            assert (subreads.size() > 0);
//                if (v == read->get_vertex_id(/*normal*/ true)) {
            if (overlap_graph->getOrientation(v)) { // forward read
                read->set_read_id(count); // update read ID
                read->set_super(true);
                read->set_original_reads(subreads);
                trivial_SR_vec.push_back(*read);
            }
            else {
                // reverse reads are only allowed at the first iteration, so we copy the read into a new (forward) read
//                    assert (v == read->get_vertex_id(/*normal*/ false));
                std::unordered_map< read_id_t, OriginalIndex > updated_subreads;
                if (read->is_paired()) {
                    Read rev_read(true, true, count, read->get_rev_comp(2), read->get_rev_comp(1), read->get_rev_phred(2), read->get_rev_phred(1));
                    // update subread indexes because of reversal
                    for (auto sub_it : subreads) {
                        OriginalIndex original_index = sub_it.second;
                        original_index.forward = !original_index.forward;
                        original_index.index1 = read->get_seq(1).size() - (original_index.index1 + original_index.len1);
                        original_index.index2 = read->get_seq(2).size() - (original_index.index2 + original_index.len2);
                        updated_subreads.insert(std::make_pair(sub_it.first, original_index));
                    }
                    rev_read.set_original_reads(updated_subreads);
                    trivial_SR_vec.push_back(rev_read);
                }
                else {
                    Read rev_read(false, true, count, read->get_rev_comp(0), "", read->get_rev_phred(0), "");
                    // update subread indexes because of reversal
                    for (auto sub_it : subreads) {
                        OriginalIndex original_index = sub_it.second;
                        original_index.forward = !original_index.forward;
                        original_index.index1 = read->get_seq(0).size() - (original_index.index1 + original_index.len1);
                        if (original_index.is_paired) {
                            original_index.index2 = read->get_seq(0).size() - (original_index.index2 + original_index.len2);
                        }
                        updated_subreads.insert(std::make_pair(sub_it.first, original_index));
                    }
                    rev_read.set_original_reads(updated_subreads);
                    trivial_SR_vec.push_back(rev_read);
                }
            }
            nodes_to_new_IDs.insert(std::pair<node_id_t, read_id_t>(v, count)); // create a map from old IDs to new IDs
            count++;
        }
    }
    if (program_settings.verbose) {
        std::cout << "Number of trivial superreads: " << trivial_SR_vec.size() << "\n";
    }
    writeTrivialsToFile();
    writePairsToFile(paired_SR_vec, count);
    writeTipsToFile();
    new_read_count = count;
    SR_singles_count = single_SR_vec.size();
    SR_paired_count = paired_SR_vec.size();
    SR_trivials_count = trivial_SR_vec.size();
}

void SRBuilder::writeTipsToFile() {
    if (tips_vec.empty()) { // nothing to be done
        return;
    }
    std::string tips_filename = PATH + "removed_tip_sequences.fastq";
    std::ofstream tips_file(tips_filename, std::fstream::app);
    read_id_t new_ID = 0;
    for (auto read : tips_vec) {
        if (read.is_paired()) {
            tips_file << "@" << new_ID << "_1\n";
            tips_file << read.get_seq(1) << "\n";
            tips_file << "+\n";
            tips_file << read.get_phred(1) << "\n";
            tips_file << "@" << new_ID << "_2\n";
            tips_file << read.get_seq(2) << "\n";
            tips_file << "+\n";
            tips_file << read.get_phred(2) << "\n";
            new_ID++;
        }
        else {
            tips_file << "@" << new_ID << "\n";
            tips_file << read.get_seq(0) << "\n";
            tips_file << "+\n";
            tips_file << read.get_phred(0) << "\n";
            new_ID++;
        }
    }
    tips_file.close();
}

void SRBuilder::writeTrivialsToFile() { // note that the trivial superreads have already been given new read IDs
//    std::cout << "writetrivials\n";
    std::string filename0 = PATH + "singles.fastq";
    std::string filename1 = PATH + "paired1.fastq";
    std::string filename2 = PATH + "paired2.fastq";
    std::string originals_file = PATH + "subreads.txt";
    std::ofstream outfile0(filename0, std::fstream::app);
    std::ofstream outfile1(filename1, std::fstream::app);
    std::ofstream outfile2(filename2, std::fstream::app);
    std::ofstream originals(originals_file, std::fstream::app);

    std::deque<Read>::const_iterator it;
    std::deque<Read> paired_trivials;
    for (it = trivial_SR_vec.begin(); it != trivial_SR_vec.end(); it++) {
        read_id_t ID = it->get_read_id();
        std::string ori;
        int len1;
        if (it->is_paired()) {
            outfile1 << "@" << ID << "\n";
            outfile1 << it->get_seq(1) << "\n";
            outfile1 << "+\n";
            outfile1 << it->get_phred(1) << "\n";
            outfile2 << "@" << ID << "\n";
            outfile2 << it->get_seq(2) << "\n";
            outfile2 << "+\n";
            outfile2 << it->get_phred(2) << "\n";
        }
        else {
            outfile0 << "@" << ID << "\n";
            outfile0 << it->get_seq(0) << "\n";
            outfile0 << "+\n";
            outfile0 << it->get_phred(0) << "\n";
        }
        std::unordered_map< read_id_t, OriginalIndex > subreads = it->get_original_reads();
        originals << ID;
        for (auto subread_it : subreads) {
            OriginalIndex original_index = subread_it.second;
            ori = original_index.forward ? "+" : "-";
            len1 = original_index.len1;
            originals << "\t" << subread_it.first << ":" << ori << ":" << original_index.index1;
            if (original_index.is_paired) {
                originals << "," << original_index.index2 << ":" << len1 << "," << original_index.len2;
            }
            else {
                originals << ":" << len1;
            }
        }
        originals << "\n";
    }
    outfile0.close();
    outfile1.close();
    outfile2.close();
    originals.close();
}

void SRBuilder::writeSinglesToFile(std::deque<Read>& singles, read_id_t& count) {
//    std::cout << "writesingles\n";
    std::string filename = PATH + "singles.fastq";
    std::string originals_file = PATH + "subreads.txt";
    std::ofstream outfile(filename, std::fstream::app);
    std::ofstream originals(originals_file, std::fstream::app);
    std::deque<Read>::iterator it;
//    for (it = single_SR_vec.begin(); it != single_SR_vec.end(); it++) {
    for (it = singles.begin(); it != singles.end(); it++) {
        read_id_t ID = count++;
        it->set_read_id(ID);
        std::string ori;
        int len1;
        outfile << "@" << ID << "\n";
        outfile << it->get_seq(0) << "\n";
        outfile << "+\n";
        outfile << it->get_phred(0) << "\n";

        std::unordered_map< read_id_t, OriginalIndex > subreads = it->get_original_reads();
        originals << ID;
        for (auto subread_it : subreads) {
            OriginalIndex original_index = subread_it.second;
            ori = original_index.forward ? "+" : "-";
            len1 = original_index.len1;
            originals << "\t" << subread_it.first << ":" << ori << ":" << original_index.index1;
            if (original_index.is_paired) {
                originals << "," << original_index.index2 << ":" << len1 << "," << original_index.len2;
            }
            else {
                originals << ":" << len1;
            }
        }
        originals << "\n";
    }
    outfile.close();
    originals.close();
}

void SRBuilder::writePairsToFile(std::deque<Read>& pairs, read_id_t& count) {
//    std::cout << "writepairs\n";
    std::string filename1 = PATH + "paired1.fastq";
    std::string filename2 = PATH + "paired2.fastq";
    std::string originals_file = PATH + "subreads.txt";
    std::ofstream outfile1(filename1, std::fstream::app);
    std::ofstream outfile2(filename2, std::fstream::app);
    std::ofstream originals(originals_file, std::fstream::app);
    std::deque<Read>::iterator it;
//    for (it = paired_SR_vec.begin(); it != paired_SR_vec.end(); it++) {
    for (it = pairs.begin(); it != pairs.end(); it++) {
        read_id_t ID = count++;
        it->set_read_id(ID);
        std::string ori;
        int len1;
        std::string seq1 = it->get_seq(1);
        std::string seq2 = it->get_seq(2);
        assert (seq1.size() > 0);
        assert (seq2.size() > 0);
        outfile1 << "@" << ID << "\n";
        outfile1 << it->get_seq(1) << "\n";
        outfile1 << "+\n";
        outfile1 << it->get_phred(1) << "\n";
        outfile2 << "@" << ID << "\n";
        outfile2 << it->get_seq(2) << "\n";
        outfile2 << "+\n";
        outfile2 << it->get_phred(2) << "\n";

        std::unordered_map< read_id_t, OriginalIndex > subreads = it->get_original_reads();
        originals << ID;
        for (auto subread_it : subreads) {
            OriginalIndex original_index = subread_it.second;
            ori = original_index.forward ? "+" : "-";
            len1 = original_index.len1;
            originals << "\t" << subread_it.first << ":" << ori << ":" << original_index.index1;
            if (original_index.is_paired) {
                originals << "," << original_index.index2 << ":" << len1 << "," << original_index.len2;
            }
            else {
                originals << ":" << len1;
            }
        }
        originals << "\n";
    }
    outfile1.close();
    outfile2.close();
    originals.close();
}
