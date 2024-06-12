//============================================================================
// Name        : EdgeCalculator.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Compute edges from overlaps file by computing overlap scores
//============================================================================

#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <algorithm> // std::swap

#include "EdgeCalculator.h"

// computes score for 2 specific nucleotides and their qualities
double EdgeCalculator::score(char nt1, char nt2, double p1, double p2, int & mismatch_count)
{
//    std::cout << "In score..\n";
    assert (nt1=='A' || nt1=='T' || nt1=='C' || nt1=='G' || nt1=='N');
    assert (nt2=='A' || nt2=='T' || nt2=='C' || nt2=='G' || nt2=='N');
//    if (!program_settings.error_correction && nt1 != nt2) {
//        return 1;
//    }
	double p;
	if (nt1=='N' || nt2=='N') {
//	    p = 0.25;
//        mismatch_count++;
        return 1;
	}
    else if (nt1==nt2) {
        p = (1-p1)*(1-p2) + (p1*p2)/3.0;
    }
    else {
        p = p1*(1-p2)/3.0 + p2*(1-p1)/3.0 + (2/9.0)*p1*p2;
        mismatch_count++;
    }
    assert (p > 0 && p <= 1);
    // only accept alignment if the probability that the reads come from the same strain is sufficiently high
    if (p < program_settings.mismatch) {
        return 2;
    }
	double lp = log(p);
	assert (lp == lp); // checks that lp != NaN
	assert (lp <= 0);
	return lp;
}


double EdgeCalculator::phred_to_prob(const int phred) {
    double P = pow(10, -phred/10.0);
    assert (P >= 0 && P <= 1);
    return P;
}


// Overlap score computation for 2 given sequences and phred scores, together with the overlap start position
double EdgeCalculator::overlap_score(std::string seq1, std::string seq2, std::string score1, std::string score2, const unsigned int pos, double & mismatch_rate)
{
//    std::cout << "In overlap_score..\n";
    assert (seq1.length() > 0);
    assert (seq2.length() > 0);
    assert (score1.length() > 0);
    assert (score2.length() > 0);
    mismatch_rate = 1.0; // will be updated after overlap computation
    assert (pos >= 0);
    if (pos >= seq1.length()) {
        std::cout << "pos: " << pos << ", seq length: " << seq1.length() << std::endl;
        return 0;
    }
    assert (pos < seq1.length());

    if (seq1.length() < program_settings.min_read_len || seq2.length() < program_settings.min_read_len) { // too short; TODO: make parameter setting
        return 0;
    }

    unsigned int L1 = seq1.length();
    unsigned int L2 = seq2.length();
    int L = std::min(L1-pos, L2); // overlap length

    std::vector<double> probs1;
    std::vector<double> probs2;
    for(int i=0; i<L; i++) {
        int Q1 = static_cast<int>(score1.at(i+pos)) - 33;
        int Q2 = static_cast<int>(score2.at(i)) - 33;
        double P1 = phred_to_prob(Q1);
        double P2 = phred_to_prob(Q2);
        assert (P1 == P1 && P1 >= 0 && P1 <= 1);
        assert (P2 == P2 && P2 >= 0 && P2 <= 1);
        probs1.push_back(P1);
        probs2.push_back(P2);
    }

    double total_score = 0.0;
    double total_len = 0.0; // count the number of overlapping A/C/T/G nucleotides
    int mismatch_count = 0;
    for (int i = 0; i < L; i++) {
//        if (seq1.at(i+pos) != 'N' && seq2.at(i) != 'N') {
//            double s = score(seq1.at(i+pos), seq2.at(i), probs1.at(i), probs2.at(i), mismatch_count);
//            if (s <= 0) {
//                total_score += s;
//                total_len += 1;
//            }
//            else { // there is an unacceptable mismatch in the overlap
//                return 0;
//            }
//        }
        double s = score(seq1.at(i+pos), seq2.at(i), probs1.at(i), probs2.at(i), mismatch_count);
        if (s <= 0) {
            total_score += s;
            total_len += 1;
        }
        else if (s == 1) { // N-base in overlap
            continue;
        }
        else { // there is an unacceptable mismatch in the overlap
            return 0;
        }
    }
    if (total_len == 0) {
        return 0;
    }
    mismatch_rate = float(mismatch_count)/total_len;
    if (mismatch_rate != mismatch_rate) { // check for NaN
        std::cerr << "mismatch rate NaN" << std::endl;
        exit(1);
    }
    total_score = (1.0/total_len)*total_score;
    return exp(total_score);
}


// Finds the overlap case for score computation of input overlap: collects the correct sequences, then calls overlap_score.
Edge EdgeCalculator::compute_overlap(const Overlap &overlap)
{
//    std::cout << "In compute_overlap\n";
	read_id_t id1 = overlap.get_id(1);
	read_id_t id2 = overlap.get_id(2);
	unsigned int pos1 = overlap.get_pos(1);
	unsigned int pos2 = overlap.get_pos(2);
	std::string ord = overlap.get_ord();
	bool ori1 = (overlap.get_ori(1) == "+") ? true : false;
	bool ori2 = (overlap.get_ori(2) == "+") ? true : false;
	if (!(program_settings.add_duplicates || program_settings.resolve_orientations)) { assert (ori1 && ori2); }
	int overlap_perc = overlap.get_perc();
	int overlap_len1 = overlap.get_len(1);
	int overlap_len2 = overlap.get_len(2);
    unsigned int index1;
    unsigned int index2;
    Read* read_1;
    Read* read_2;
    unsigned int single_count;
    single_count = fastq_storage->m_readcount_single;
    assert ((fastq_storage->m_ID_to_index).size() > 0);
    if ((fastq_storage->m_ID_to_index).find(id1) == (fastq_storage->m_ID_to_index).end()) {
        std::cout << id1 << "\n";
    }
    if ((fastq_storage->m_ID_to_index).find(id2) == (fastq_storage->m_ID_to_index).end()) {
        std::cout << id2 << "\n";
    }
    index1 = (fastq_storage->m_ID_to_index).at(id1);
    index2 = (fastq_storage->m_ID_to_index).at(id2);
//    std::cout << "check\n";
    read_1 = (fastq_storage->m_read_vec).at(index1);
    read_2 = (fastq_storage->m_read_vec).at(index2);
    node_id_t node1, node2;
    if (program_settings.add_duplicates) {
	    node1 = read_1->get_vertex_id(ori1);
	    node2 = read_2->get_vertex_id(ori2);
	}
	else {
	    node1 = read_1->get_vertex_id(true);
	    node2 = read_2->get_vertex_id(true);
	}
	assert (read_1 != read_2);
	assert (node1 != node2);
	std::string type1 = read_1->is_paired() ? "p" : "s";
	std::string type2 = read_2->is_paired() ? "p" : "s";
/*
    Edge nonedge(-1, 0, 0, 0, 0, "", read_1, read_2);
    if (id1 == id2) {
        self_overlap_count++;
        return Edge(0, 0, 0, 0, 0, "", read_1, read_2);
    }
*/
    double ov1 = 0;
    double ov2 = 0;
    if (single_count > 0) {
//        std::cout << "Considering overlaps with single-end reads\n";
        if (type1 == "s" && type2 == "s") {
//            std::cout << "S-S overlap\n";
            std::string seq1, seq2;
            std::string phred1, phred2;
            if (ori1) {
                seq1 = read_1->get_seq(0);
                phred1 = read_1->get_phred(0);
            }
            else {
                seq1 = read_1->get_rev_comp(0);
                phred1 = read_1->get_rev_phred(0);
            }
            if (ori2) {
                seq2 = read_2->get_seq(0);
                phred2 = read_2->get_phred(0);
            }
            else {
                seq2 = read_2->get_rev_comp(0);
                phred2 = read_2->get_rev_phred(0);
            }
            double score;
            double mismatch_rate;
            score = overlap_score(seq1, seq2, phred1, phred2, pos1, mismatch_rate);
            int pos3 = seq1.size() - pos1 - seq2.size();
            Edge edge(score, pos1, pos2, ori1, ori2, ord, read_1, read_2);
            edge.set_vertices(node1, node2);
            edge.set_extra_pos(pos3);
            edge.set_perc(overlap_perc);
            edge.set_len(overlap_len1, 0);
            if (!(mismatch_rate >= 0)) {
                std::cout << "mismatch rate: " << mismatch_rate << std::endl;
            }
            edge.set_mismatch(mismatch_rate);
            return edge;
        }
        else if (type1 == "s" && type2 == "p") {
//            std::cout << "S-P overlap\n";
//            double ov1, ov2;
            double mismatch1, mismatch2;
            if (ori1 && ori2) {
                ov1 = overlap_score(read_1->get_seq(0), read_2->get_seq(1), read_1->get_phred(0), read_2->get_phred(1), pos1, mismatch1);
                ov2 = overlap_score(read_1->get_seq(0), read_2->get_seq(2), read_1->get_phred(0), read_2->get_phred(2), pos2, mismatch2);
            }
            else if (!ori1 && ori2) {
                ov1 = overlap_score(read_1->get_rev_comp(0), read_2->get_seq(1), read_1->get_rev_phred(0), read_2->get_phred(1), pos1, mismatch1);
                ov2 = overlap_score(read_1->get_rev_comp(0), read_2->get_seq(2), read_1->get_rev_phred(0), read_2->get_phred(2), pos2, mismatch2);
            }
            else if (ori1 && !ori2) {
                ov1 = overlap_score(read_1->get_seq(0), read_2->get_rev_comp(2), read_1->get_phred(0), read_2->get_rev_phred(2), pos1, mismatch1);
                ov2 = overlap_score(read_1->get_seq(0), read_2->get_rev_comp(1), read_1->get_phred(0), read_2->get_rev_phred(1), pos2, mismatch2);
            }
            else {
                ov1 = overlap_score(read_1->get_rev_comp(0), read_2->get_rev_comp(2), read_1->get_rev_phred(0), read_2->get_rev_phred(2), pos1, mismatch1);
                ov2 = overlap_score(read_1->get_rev_comp(0), read_2->get_rev_comp(1), read_1->get_rev_phred(0), read_2->get_rev_phred(1), pos2, mismatch2);
            }
            double mismatch_rate = std::max(mismatch1, mismatch2);
            double score;
            if (ov1 > program_settings.edge_threshold && ov2 > program_settings.edge_threshold) {
                score = 0.5 * (ov1 + ov2);
            }
            else {
                score = std::min(ov1, ov2);
            }
            int pos3 = (read_1->get_seq(0)).size() - pos2 - (read_2->get_seq(2)).size();
            int pos4 = (read_1->get_seq(0)).size() - pos1 - (read_2->get_seq(1)).size();
            Edge edge(score, pos1, pos2, ori1, ori2, ord, read_1, read_2);
            edge.set_vertices(node1, node2);
            edge.set_extra_pos(pos3, pos4);
            edge.set_perc(overlap_perc);
            edge.set_len(overlap_len1, overlap_len2);
            edge.set_mismatch(mismatch_rate);
            return edge;
        }
        else if (type1 == "p" && type2 == "s") {
//            std::cout << "P-S overlap\n";
//            double ov1, ov2;
            double mismatch1, mismatch2;
            if (ori1 && ori2) {
                ov1 = overlap_score(read_1->get_seq(1), read_2->get_seq(0), read_1->get_phred(1), read_2->get_phred(0), pos1, mismatch1);
                ov2 = overlap_score(read_2->get_seq(0), read_1->get_seq(2), read_2->get_phred(0), read_1->get_phred(2), pos2, mismatch2);
            }
            else if (!ori1 && ori2) {
                ov1 = overlap_score(read_1->get_rev_comp(2), read_2->get_seq(0), read_1->get_rev_phred(2), read_2->get_phred(0), pos1, mismatch1);
                ov2 = overlap_score(read_2->get_seq(0), read_1->get_rev_comp(1), read_2->get_phred(0), read_1->get_rev_phred(1), pos2, mismatch2);
            }
            else if (ori1 && !ori2) {
                ov1 = overlap_score(read_1->get_seq(1), read_2->get_rev_comp(0), read_1->get_phred(1), read_2->get_rev_phred(0), pos1, mismatch1);
                ov2 = overlap_score(read_2->get_rev_comp(0), read_1->get_seq(2), read_2->get_rev_phred(0), read_1->get_phred(2), pos2, mismatch2);
            }
            else {
                ov1 = overlap_score(read_1->get_rev_comp(2), read_2->get_rev_comp(0), read_1->get_rev_phred(2), read_2->get_rev_phred(0), pos1, mismatch1);
                ov2 = overlap_score(read_2->get_rev_comp(0), read_1->get_rev_comp(1), read_2->get_rev_phred(0), read_1->get_rev_phred(1), pos2, mismatch2);
            }
            double mismatch_rate = std::max(mismatch1, mismatch2);
            double score;
            if (ov1 > program_settings.edge_threshold && ov2 > program_settings.edge_threshold) {
                score = 0.5 * (ov1 + ov2);
            }
            else {
                score = std::min(ov1, ov2);
            }
            int pos3 = (read_1->get_seq(2)).size() + pos2 - (read_2->get_seq(0)).size();
            int pos4 = (read_2->get_seq(0)).size() + pos1 - (read_1->get_seq(1)).size();
            Edge edge(score, pos1, pos2, ori1, ori2, ord, read_1, read_2);
            edge.set_vertices(node1, node2);
            edge.set_extra_pos(pos3, pos4);
            edge.set_perc(overlap_perc);
            edge.set_len(overlap_len1, overlap_len2);
            edge.set_mismatch(mismatch_rate);
            return edge;
        }
    }

    if (type1 == "p" && type2 == "p") {
//        std::cout << "P-P overlap\n";
//        double ov1, ov2;
        double mismatch1, mismatch2;
        if (!ori1 && ori2) {
            ov1 = overlap_score(read_1->get_rev_comp(2), read_2->get_seq(1), read_1->get_rev_phred(2), read_2->get_phred(1), pos1, mismatch1);
            if (ord == "1") {
                ov2 = overlap_score(read_1->get_rev_comp(1), read_2->get_seq(2), read_1->get_rev_phred(1), read_2->get_phred(2), pos2, mismatch2);
            }
            else if (ord == "2") {
                ov2 = overlap_score(read_2->get_seq(2), read_1->get_rev_comp(1), read_2->get_phred(2), read_1->get_rev_phred(1), pos2, mismatch2);
            }
        }
        else if (ori1 && !ori2) {
            ov1 = overlap_score(read_1->get_seq(1), read_2->get_rev_comp(2), read_1->get_phred(1), read_2->get_rev_phred(2), pos1, mismatch1);
            if (ord == "1") {
                ov2 = overlap_score(read_1->get_seq(2), read_2->get_rev_comp(1), read_1->get_phred(2), read_2->get_rev_phred(1), pos2, mismatch2);
            }
            else if (ord == "2") {
                ov2 = overlap_score(read_2->get_rev_comp(1), read_1->get_seq(2), read_2->get_rev_phred(1), read_1->get_phred(2), pos2, mismatch2);
            }
        }
        else if (ori1 && ori2) {
            ov1 = overlap_score(read_1->get_seq(1), read_2->get_seq(1), read_1->get_phred(1), read_2->get_phred(1), pos1, mismatch1);
            if (ord == "1") {
                ov2 = overlap_score(read_1->get_seq(2), read_2->get_seq(2), read_1->get_phred(2), read_2->get_phred(2), pos2, mismatch2);
            }
            else if (ord == "2") {
                ov2 = overlap_score(read_2->get_seq(2), read_1->get_seq(2), read_2->get_phred(2), read_1->get_phred(2), pos2, mismatch2);
            }
        }
        else {
            ov1 = overlap_score(read_1->get_rev_comp(2), read_2->get_rev_comp(2), read_1->get_rev_phred(2), read_2->get_rev_phred(2), pos1, mismatch1);
            if (ord == "1") {
                ov2 = overlap_score(read_1->get_rev_comp(1), read_2->get_rev_comp(1), read_1->get_rev_phred(1), read_2->get_rev_phred(1), pos2, mismatch2);
            }
            else if (ord == "2") {
                ov2 = overlap_score(read_2->get_rev_comp(1), read_1->get_rev_comp(1), read_2->get_rev_phred(1), read_1->get_rev_phred(1), pos2, mismatch2);
            }
        }
//        std::cout << "Scores: " << ov1 << " " << ov2 << std::endl;
        double mismatch_rate = std::max(mismatch1, mismatch2);
        double score;
        if (ov1 > program_settings.edge_threshold && ov2 > program_settings.edge_threshold) {
            score = 0.5 * (ov1 + ov2);
        }
        else {
            score = std::min(ov1, ov2);
        }
        int pos3;
        if (ord == "1") {
            pos3 = (read_1->get_seq(2)).size() - pos2 - (read_2->get_seq(2)).size();
        }
        else {
            if (ord != "2") {
                std::cout << overlap.get_overlap_line() << std::endl;
            }
            assert (ord == "2");
            pos3 = (read_1->get_seq(2)).size() + pos2 - (read_2->get_seq(2)).size();
        }
        int pos4 = (read_1->get_seq(1)).size() - pos1 - (read_2->get_seq(1)).size();
        Edge edge(score, pos1, pos2, ori1, ori2, ord, read_1, read_2);
        edge.set_vertices(node1, node2);
        edge.set_extra_pos(pos3, pos4);
        edge.set_perc(overlap_perc);
        edge.set_len(overlap_len1, overlap_len2);
        edge.set_mismatch(mismatch_rate);
        return edge;
    }
    else {
        std::cout << "Read types not recognized. Returning nonedge.\n";
        return Edge(0, 0, 0, 0, 0, "", read_1, read_2);
    }
}


// Initiates overlap score computations (in parallel) for overlaps from input vector
void EdgeCalculator::process_overlaps(std::vector<Overlap> overlaps_vec)
{
//    std::cout << "process_overlaps\n";
    unsigned int size = overlaps_vec.size();
    std::vector<Edge> edges;
    std::vector<Overlap> nonedge_overlaps;
	#pragma omp parallel num_threads(N_THREADS) shared(overlaps_vec, edges, nonedge_overlaps)
	{
	    std::vector<Edge> edges_this_thread;
        std::vector<Overlap> overlaps_this_thread;
	    #pragma omp for
        for (unsigned int i = 0; i < size; i++)
        {
            Overlap overlap = overlaps_vec.at(i);
            Edge edge = compute_overlap(overlap);
            if (edge.get_score() > program_settings.edge_threshold) {
                edges_this_thread.push_back(edge); // add edge to vector
            }
            else if (edge.get_mismatch_rate() != -1 && edge.get_mismatch_rate() <= program_settings.merge_contigs) {
                edges_this_thread.push_back(edge);
            }
            else if (edge.get_score() > program_settings.ov_threshold && edge.get_mismatch_rate() != -1) {
                assert(overlap.get_id(1) != overlap.get_id(2));
                overlaps_this_thread.push_back(overlap); // add overlap to nonedge vector
            }
        }
        #pragma omp critical(write_edges)
        {
            edges.insert(edges.end(), edges_this_thread.begin(), edges_this_thread.end());
        }
        #pragma omp critical(write_overlaps)
        {
            nonedge_overlaps.insert(nonedge_overlaps.end(), overlaps_this_thread.begin(), overlaps_this_thread.end());
        }
    }
    if (program_settings.verbose) {
        std::cout << "build edges / write overlaps to file\n";
    }
    #pragma omp parallel sections num_threads(N_THREADS)
    {
        #pragma omp section
        {
            unsigned int count = 0;
            unsigned int doubles = 0;
            for (std::vector<Edge>::iterator it1 = edges.begin(); it1 != edges.end(); it1++) {
//                // check for self-overlap:
//                if (it1->get_read(1) == it1->get_read(2)) {
//                    self_overlap_count++;
//                    continue;
//                }
                // check if edge already exists:
                double score;
                node_id_t v1 = it1->get_vertex(1);
                node_id_t v2 = it1->get_vertex(2);
                if (it1->get_pos(1) == 0) { // edge can be added in either direction; always direct them from small to large ID
                    if (v1 > v2) {
                        std::swap(v1, v2);
                        it1->swap_reads();
                    }
                }
                if (it1->get_perc() == 100) {
                    inclusion_count++;
                }
                // score = overlap_graph->checkEdge(v1, v2, /*reverse_allowed*/ true);
                bool opposite_orientations = (it1->get_ori(1) == it1->get_ori(2));
                score = overlap_graph->checkEdgeWithOri(v1, v2, opposite_orientations);
                if (score < 0) {
                    // edge does not yet exist, so add it now
                    overlap_graph->addEdge(*it1);
                    count++;
                    if (program_settings.ignore_inclusions && it1->get_perc() == 100 && it1->get_mismatch_rate() < 0.000001 && it1->get_mismatch_rate() >= 0) {
                        if (it1->get_extra_pos(1) < 0) {
                            if (it1->get_pos(1) == 0) { // otherwise not a true inclusion but simply the effect of rounding the overlap percentage
                                overlap_graph->inclusions[v1] = 1;
                            }
                        }
                        else {
                            overlap_graph->inclusions[v2] = 1;
                        }
                    }
                }
                else if (it1->get_score() >= score) {
                    // new edge scores better, so replace current edge in the graph
                    doubles++;
                    Edge* existing_edge = overlap_graph->getEdgeInfoWithOri(v1, v2, opposite_orientations, true);
                    if (score == it1->get_score()) {
                        // decide which edge is 'greater' to ensure deterministic behaviour
                        if (existing_edge->get_len(0) != it1->get_len(0)) {
                            if (existing_edge->get_len(0) > it1->get_len(0)) {
                                // existing edge has longer overlap
                                continue;
                            }
                        }
                        else if (existing_edge->get_mismatch_rate() != it1->get_mismatch_rate()) {
                            if (existing_edge->get_mismatch_rate() < it1->get_mismatch_rate()) {
                                // existing edge has lower mismatch rate
                                continue;
                            }
                        }
                        else if (existing_edge->get_vertex(1) != it1->get_vertex(1)) {
                            if (existing_edge->get_vertex(1) < it1->get_vertex(1)) {
                                // existing edge has different direction
                                continue;
                            }
                        }
                        else if (existing_edge->get_ori(1) != it1->get_ori(1)) {
                            if (existing_edge->get_ori(1)) {
                                // same direction but different orientation of v1
                                continue;
                            }
                        }
                        else if (existing_edge->get_ori(2) != it1->get_ori(2)) {
                            if (existing_edge->get_ori(2)) {
                                // same direction but different orientation of v2
                                continue;
                            }
                        }
                        else if (existing_edge->get_pos(1) != it1->get_pos(1)) {
                            if (existing_edge->get_pos(1) < it1->get_pos(1)) {
                                // existing edge has smaller overlap positions
                                continue;
                            }
                        }
                        else if (existing_edge->get_pos(2) != it1->get_pos(2)) {
                            if (existing_edge->get_pos(2) < it1->get_pos(2)) {
                                // existing edge has smaller overlap positions
                                continue;
                            }
                        }
                        else {
//                            std::cout << "Completely equivalent edge candidates: which to choose?" << std::endl;
                        }
                    }
                    int edgecount1 = overlap_graph->getEdgeCount();
                    if (existing_edge->get_vertex(1) == v1) {
                        overlap_graph->removeEdgeWithOri(v1, v2, opposite_orientations);
                    }
                    else {
                        overlap_graph->removeEdgeWithOri(v2, v1, opposite_orientations);
                    }
                    int edgecount2 = overlap_graph->getEdgeCount();
                    overlap_graph->addEdge(*it1);
                    int edgecount3 = overlap_graph->getEdgeCount();
                    assert (edgecount1 == edgecount3);
                    assert (edgecount1 == edgecount2 + 1);
                }
                else {
                    // new edge does not improve current edge, so do nothing
                    doubles++;
                }
            }
            if (program_settings.verbose) {
                std::cout << "Number of edges found: " << count << std::endl;
                std::cout << "Number of duplicates: " << doubles << std::endl;
            }
            dup_count += doubles;
        }
        #pragma omp section
        {
            std::ofstream new_overlapsfile;
            new_overlapsfile.open(program_settings.output_dir + "nonedge_overlaps.txt", std::fstream::out | std::fstream::app);
            for (std::vector<Overlap>::const_iterator it2 = nonedge_overlaps.begin(); it2 != nonedge_overlaps.end(); it2++) {
                // write nonedge to new overlapsfile
                new_overlapsfile << it2->get_overlap_line();
            }
            new_overlapsfile.close();
        }
    }
}


// Builds edges from overlaps file and add to the overlap graph
void EdgeCalculator::construct_edges()
{
//    std::cout << "In construct_edges... ";
//    const char* nonedge_file = (program_settings.output_dir + "nonedge_overlaps.txt").c_str();
//    std::cout << "Nonedge file: " << nonedge_file << std::endl;
    std::remove("nonedge_overlaps.txt");
    std::vector<Overlap> nonedge_overlaps;

    std::ifstream overlapsfile (program_settings.overlaps_file.c_str());
    unsigned int i = 0;
    const unsigned int overlaps_per_vec = 1000000;
    std::vector<Overlap> overlaps_vec;
    overlaps_vec.reserve(overlaps_per_vec);
    if (overlapsfile.is_open())     {
        if (program_settings.verbose) {
            std::cout << "reading overlaps file... \n";
        }
        std::stringstream ss;
        std::string tupleline;
        // read at most overlaps_per_vec overlaps into a vector and process this chunk; continue until all overlaps have been considered
        while (getline(overlapsfile, tupleline) && i < program_settings.max_overlaps)
        {
            i++;
            boost::trim_if(tupleline, boost::is_any_of("\t ")); // trim any outer spacing
            std::vector< std::string > tmp_vec;
            if (program_settings.allow_spaces) {
            	boost::algorithm::split(tmp_vec, tupleline, boost::is_any_of("\t "), boost::token_compress_on); // split by tabs OR spaces
            }
            else {
                ss << tupleline;
                std::string tmp;
                while (getline(ss, tmp, '\t')) { // split by tabs
                    tmp_vec.push_back(tmp);
                }
    			ss << "";
    			ss.clear();
            }
            unsigned int TUP_LEN = 13; // number of fields in overlap file according to format definition
//            assert (tmp_vec.size() == TUP_LEN);
            if (tmp_vec.size() != TUP_LEN) {
                std::cout << "incorrect overlap; skipping" << std::endl;
                continue;
            }
            Overlap overlap(tmp_vec);
            if (overlap.get_id(1) == overlap.get_id(2)) {
                continue;
            }
            if (program_settings.cliques && !program_settings.error_correction
                && overlap.get_perc() == 100) {
//                continue;
            }
            if (overlap.get_len(1) >= program_settings.min_overlap_len
                && overlap.get_type(1) == "s" && overlap.get_type(2) == "s") {
                if (overlap.get_perc() >= program_settings.min_overlap_perc) {
                    overlaps_vec.push_back(overlap);
                }
            }
            else if (overlap.get_len(1) >= 0.5*program_settings.min_overlap_len
                && overlap.get_len(2) >= 0.5*program_settings.min_overlap_len
                && (overlap.get_type(1) == "p" || overlap.get_type(2) == "p")) {
                if (overlap.get_perc() >= program_settings.min_overlap_perc) {
                    overlaps_vec.push_back(overlap);
                }
            }
            // relaxed edge condition: allow len1 + len2 >= threshold
            else if (program_settings.relax_PE_edges
                && overlap.get_len(1) + overlap.get_len(2) >= program_settings.min_overlap_len
                && (overlap.get_type(1) == "p" || overlap.get_type(2) == "p")) {
                if (overlap.get_perc() >= program_settings.min_overlap_perc) {
                    overlaps_vec.push_back(overlap);
                }
            }
            else { // store overlap, later on write it back to file (nonedge_overlaps)
                nonedge_overlaps.push_back(overlap);
            }
            if (overlaps_vec.size() == overlaps_per_vec) {
                process_overlaps(overlaps_vec); // process the currently collected overlaps
                overlaps_vec.clear(); // empty the vector
            }
        }
        if (overlaps_vec.size() > 0) { // process the remaining overlaps
            process_overlaps(overlaps_vec); // process the currently collected overlaps
            overlaps_vec.clear();
        }
        overlapsfile.close();
        if (program_settings.verbose) {
            std::cout << "Number of self-overlapping reads: " << self_overlap_count << "\n";
            std::cout << "Number of inclusion edges: " << inclusion_count << "\n";
        }
        if (program_settings.add_duplicates) {
            overlap_graph->addEquivalentEdges(); // for every ++edge also add the equivalent --edge, and similarly for +-/-+ edges
        }

        std::ofstream new_overlapsfile;
        new_overlapsfile.open(program_settings.output_dir + "nonedge_overlaps.txt", std::fstream::out | std::fstream::app);
        for (auto overlap_it : nonedge_overlaps) {
            // write discarded overlaps to new overlapsfile
            new_overlapsfile << overlap_it.get_overlap_line();
        }
        new_overlapsfile.close();
    }
    else {
        std::cerr << "Unable to open overlaps file";
        exit(1);
    }
}
