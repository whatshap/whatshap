//============================================================================
// Name        : OverlapGraph.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Construct an overlap graph for viral quasispecies assembly
//============================================================================

#include <fstream>
#include <assert.h>
#include <cstdlib>
#include <algorithm>
#include <stdio.h>
#include <boost/timer.hpp>
#include <map>

#include "OverlapGraph.h"


void OverlapGraph::getGraphStats() {
    int pp_count[5] = {0};
    int ps_count[5] = {0};
    int sp_count[5] = {0};
    int ss_count[5] = {0};

    std::map<std::string, int> ori_to_index;
    ori_to_index["++"] = 0;
    ori_to_index["+-"] = 1;
    ori_to_index["-+"] = 2;
    ori_to_index["--"] = 3;

    std::vector< std::list<Edge> >::const_iterator it;
    for (it = adj_out.begin(); it != adj_out.end(); it++) {
        std::list<Edge>::const_iterator it2;
        for (it2 = it->begin(); it2 != it->end(); it2++) {
            bool type1 = (it2->get_read(1))->is_paired();
            bool type2 = (it2->get_read(2))->is_paired();
            std::string ori;
            if (it2->get_ori(1) && it2->get_ori(2))
                ori = "++";
            else if (!it2->get_ori(1) && it2->get_ori(2))
                ori = "-+";
            else if (it2->get_ori(1) && !it2->get_ori(2))
                ori = "+-";
            else
                ori = "--";
            int index = ori_to_index[ori];

            if (type1 && type2) {
                pp_count[index]++;
                pp_count[4]++; // keep total pp-count
            }
            else if (type1 && !type2) {
                ps_count[index]++;
                ps_count[4]++; // keep total ps-count
            }
            else if (!type1 && type2) {
                sp_count[index]++;
                sp_count[4]++; // keep total sp-count
            }
            else if (!type1 && !type2) {
                ss_count[index]++;
                ss_count[4]++; // keep total ss-count
            }
            else {
                std::cout << "Overlap type not recognized...\n";
                return;
            }
        }
    }
    if (program_settings.verbose) {
        std::cout << "\nEdge statistics:\n";
        std::cout << "\t++\t+-\t-+\t--\tTotal\n";
        std::cout << "P-P\t" << pp_count[0] << "\t" << pp_count[1] << "\t" << pp_count[2] << "\t" << pp_count[3] << "\t" << pp_count[4] << "\n";
        std::cout << "P-S\t" << ps_count[0] << "\t" << ps_count[1] << "\t" << ps_count[2] << "\t" << ps_count[3] << "\t" << ps_count[4] << "\n";
        std::cout << "S-P\t" << sp_count[0] << "\t" << sp_count[1] << "\t" << sp_count[2] << "\t" << sp_count[3] << "\t" << sp_count[4] << "\n";
        std::cout << "S-S\t" << ss_count[0] << "\t" << ss_count[1] << "\t" << ss_count[2] << "\t" << ss_count[3] << "\t" << ss_count[4] << "\n\n";
    }
}


bool OverlapGraph::getOrientation(node_id_t v) {
    assert (v < vertex_orientations.size());
    return vertex_orientations[v]; // 1 if forward, 0 if reverse
}

node_id_t OverlapGraph::addVertex(read_id_t read_ID) {
    vertex_to_read.push_back(read_ID);
    vertex_count++;
    return vertex_to_read.size()-1;
}

void OverlapGraph::addEdge(Edge edge) {
//    std::cout << "In OverlapGraph::addEdge\n";
    node_id_t v = edge.get_vertex(1);
    node_id_t w = edge.get_vertex(2);
    adj_out[v].push_back(edge);
    adj_in[w].push_back(v);
    edge_count++;
}


Edge OverlapGraph::removeEdge(node_id_t v, node_id_t w) {
//    std::cout << "Removing edge " << v << " to " << w << "... ";
    bool found = false;
    std::list< Edge > adj_list_v = adj_out.at(v);
	std::list< Edge >::iterator it = adj_out.at(v).begin();
	assert (adj_list_v.size() > 0);
    Edge edge;
	for (node_id_t i = 0; i < adj_list_v.size(); i++) {
//	    std::cout << it->get_vertex(2) << " ";
		if (it->get_vertex(2) == w) {
            edge = *it;
		    int s1 = adj_out.at(v).size();
			it = adj_out.at(v).erase(it);
			int s2 = adj_out.at(v).size();
			assert (s1 == s2 + 1);
            edge_count--;
            found = true;
			break;
	    }
	    else {
	        it++;
	    }
    }
    if (!found) {
        std::cerr << "Edge to be removed not found...\n";
        std::cerr << v << " " << w << "\n";
        exit(1);
    }
    std::list< node_id_t > adj_list_w = adj_in.at(w);
    std::list< node_id_t >::iterator it2 = adj_in.at(w).begin();
    for (node_id_t i = 0; i < adj_list_w.size(); i++) {
		if (*it2 == v) {
		    int s1 = adj_in[w].size();
			it2 = adj_in[w].erase(it2);
			int s2 = adj_in[w].size();
			assert (s1 == s2 + 1);
			break;
	    }
	    else {
	        it2++;
	    }
    }
    return edge;
}

// Special edge removal, necessary when orientations have not yet been resolved
Edge OverlapGraph::removeEdgeWithOri(node_id_t v, node_id_t w, bool opposite_orientations) {
//    std::cout << "Removing edge " << v << " to " << w << "... ";
    bool found = false;
    std::list< Edge > adj_list_v = adj_out.at(v);
	std::list< Edge >::iterator it = adj_out.at(v).begin();
	assert (adj_list_v.size() > 0);
    Edge edge;
	for (node_id_t i = 0; i < adj_list_v.size(); i++) {
//	    std::cout << it->get_vertex(2) << " ";
        bool current_ori = (it->get_ori(1) == it->get_ori(2));
		if (it->get_vertex(2) == w && current_ori == opposite_orientations) {
            edge = *it;
		    int s1 = adj_out.at(v).size();
			it = adj_out.at(v).erase(it);
			int s2 = adj_out.at(v).size();
			assert (s1 == s2 + 1);
            edge_count--;
            found = true;
			break;
	    }
	    else {
	        it++;
	    }
    }
    if (!found) {
        std::cerr << "Edge to be removed not found...\n";
        std::cerr << v << " " << w << "\n";
        exit(1);
    }
    std::list< node_id_t > adj_list_w = adj_in.at(w);
    std::list< node_id_t >::iterator it2 = adj_in.at(w).begin();
    for (node_id_t i = 0; i < adj_list_w.size(); i++) {
		if (*it2 == v) {
		    int s1 = adj_in[w].size();
			it2 = adj_in[w].erase(it2);
			int s2 = adj_in[w].size();
			assert (s1 == s2 + 1);
			break;
	    }
	    else {
	        it2++;
	    }
    }
    return edge;
}


// check if edge exists; if yes, return overlap score, if no, return -1
double OverlapGraph::checkEdgeWithOri(node_id_t v, node_id_t w, bool opposite_orientations) {
    std::list< Edge >::const_iterator it;
    if (!adj_out.at(v).empty()) {
        for (it = adj_out[v].begin(); it != adj_out[v].end(); it++) {
            if (it->get_vertex(2) == w) {
                bool current_opp_ori = (it->get_ori(1) == it->get_ori(2));
                if (current_opp_ori == opposite_orientations) {
                    if (!(it->get_score() >= program_settings.edge_threshold || it->get_mismatch_rate() <= program_settings.merge_contigs)) {
                        std::cout << "score " << it->get_score() << std::endl;
                    }
                    assert (it->get_score() >= program_settings.edge_threshold || it->get_mismatch_rate() <= program_settings.merge_contigs);
                    // std::cout << "Edge found" << std::endl;
                    return it->get_score();
                }
	        }
	    }
    }
    if (!adj_out.at(w).empty()) {
        for (it = adj_out[w].begin(); it != adj_out[w].end(); it++) {
            if (it->get_vertex(2) == v) {
                bool current_opp_ori = (it->get_ori(1) == it->get_ori(2));
                if (current_opp_ori == opposite_orientations) {
                    assert (it->get_score() >= program_settings.edge_threshold || it->get_mismatch_rate() <= program_settings.merge_contigs);
                    // std::cout << "Edge found" << std::endl;
                    return it->get_score();
                }
            }
        }
    }
//    std::cout << "Edge not found" << std::endl;
    return -1;
}


// check if edge exists; if yes, return overlap score, if no, return -1
double OverlapGraph::checkEdge(node_id_t v, node_id_t w, bool reverse_allowed) {
//    std::cout << "In checkEdge()" << std::endl;
	std::list< Edge >::const_iterator it;
	if (!adj_out.at(v).empty()) {
	    for (it = adj_out[v].begin(); it != adj_out[v].end(); it++) {
		    if (it->get_vertex(2) == w) {
                if (!(it->get_score() >= program_settings.edge_threshold || it->get_mismatch_rate() <= program_settings.merge_contigs)) {
                    std::cout << "score " << it->get_score() << std::endl;
                }
		        assert (it->get_score() >= program_settings.edge_threshold || it->get_mismatch_rate() <= program_settings.merge_contigs);
//  		    std::cout << "Edge found" << std::endl;
			    return it->get_score();
	        }
	    }
    }
    if (reverse_allowed && !adj_out.at(w).empty()) {
	    for (it = adj_out[w].begin(); it != adj_out[w].end(); it++) {
		    if (it->get_vertex(2) == v) {
		        assert (it->get_score() >= program_settings.edge_threshold || it->get_mismatch_rate() <= program_settings.merge_contigs);
//  		    std::cout << "Edge found" << std::endl;
			    return it->get_score();
	        }
	    }
    }
//    std::cout << "Edge not found" << std::endl;
    return -1;
}


// return pointer to edge v->w if it exists, otherwise w->v
Edge* OverlapGraph::getEdgeInfo(node_id_t v, node_id_t w, bool reverse_allowed) {
    //std::cout << "getEdgeInfo " << v << ", " << w << std::endl;
	std::list< Edge >::iterator it;
	if (!adj_out.at(v).empty()) {
	    for (it = adj_out.at(v).begin(); it != adj_out.at(v).end(); it++) {
		    if (it->get_vertex(2) == w) {
			    return &(*it);
	        }
        }
    }
    if (reverse_allowed && !adj_out.at(w).empty()) {
        for (it = adj_out.at(w).begin(); it != adj_out.at(w).end(); it++) {
		    if (it->get_vertex(2) == v) {
			    return &(*it);
	        }
        }
    }
    std::cerr << v << " " << w << " Edge not found. Exiting.\n";
    exit(1);
}

// return pointer to edge v->w if it exists, otherwise w->v
Edge* OverlapGraph::getEdgeInfoWithOri(node_id_t v, node_id_t w,
        bool opposite_orientations, bool reverse_allowed) {
    //std::cout << "getEdgeInfo " << v << ", " << w << std::endl;
	std::list< Edge >::iterator it;
	if (!adj_out.at(v).empty()) {
	    for (it = adj_out.at(v).begin(); it != adj_out.at(v).end(); it++) {
            bool current_opp_ori = (it->get_ori(1) == it->get_ori(2));
		    if (it->get_vertex(2) == w && current_opp_ori == opposite_orientations) {
			    return &(*it);
	        }
        }
    }
    if (reverse_allowed && !adj_out.at(w).empty()) {
        for (it = adj_out.at(w).begin(); it != adj_out.at(w).end(); it++) {
            bool current_opp_ori = (it->get_ori(1) == it->get_ori(2));
		    if (it->get_vertex(2) == v && current_opp_ori == opposite_orientations) {
			    return &(*it);
	        }
        }
    }
    std::cerr << v << " " << w << " Edge not found. Exiting.\n";
    exit(1);
}

unsigned int OverlapGraph::getEdgeCount() {
    return edge_count;
}

unsigned int OverlapGraph::getBackEdgeCount() {
    return backedge_count;
}

unsigned int OverlapGraph::getVertexCount() {
    return vertex_count;
}


void OverlapGraph::writeGraphToFile() {
    std::string filename = PATH + "graph.txt";
    std::string tmpfile1 = PATH + "tmp_head.txt";
    std::string tmpfile2 = PATH + "tmp_graph.txt";
    if (program_settings.verbose) {
        std::cout << "WriteGraphToFile " << filename << "\n";
    }
    remove(filename.c_str());
    std::ofstream graph_file(tmpfile2);
//    graph_file << std::to_string(vertex_count) + "\n";
//    graph_file << std::to_string(2*(edge_count-twocycle_count)) + "\n";
	std::list< Edge >::const_iterator it1;
	std::vector< std::list< Edge > >::const_iterator it2;
	node_id_t i = 0;
	node_id_t j;
	unsigned int count = 0;
    if (program_settings.verbose) {
        std::cout << "count inclusion vertices: " << inclusions.count() << "\n";
    }
	for (it2 = adj_out.begin(); it2 != adj_out.end(); it2++) {
	    if (inclusions[i] == 1) {
            if (! it2->empty()) {
                std::cout << i << " " << (it2->front()).get_vertex(1) << " " << (it2->front()).get_vertex(2) << std::endl;
            }
            assert( it2->empty() );
            i++;
	        continue;
	    }
	    for (it1 = it2->begin(); it1 != it2->end(); it1++) {
	        j = it1->get_vertex(2);
	        if (inclusions[j] == 1) {
	            continue;
	        }
	        assert (j != i);
	        // check if reverse edge has been added
	        if (j < i && checkEdge(j, i, false) > 0) {
	            continue;
	        }
	        // write both forward and reverse edge to the graph file, since we
	        // consider the graph as undirected when enumerating maximal cliques.
		    std::string line1 = std::to_string(i) + "," + std::to_string(j) + "\n"; // forward edge
		    std::string line2 = std::to_string(j) + "," + std::to_string(i) + "\n"; // reverse edge
		    graph_file << line1;
		    graph_file << line2;
		    count++;
	    }
	    i++;
    }
    graph_file.close();
    // write head to final graph file
    std::ofstream head_file(tmpfile1);
    head_file << std::to_string(vertex_count) + "\n"; // number of vertices
    head_file << std::to_string(2*count) + "\n"; // number of edge lines
    head_file.close();
    std::string command = "cat " + tmpfile1 + " " + tmpfile2 + " > " + filename;
    int system_ret = system(command.c_str());
    if(system_ret != 0){
        // The system method failed
        std::cerr << "Adding head to file failed. Exiting..." << std::endl;
        exit(1);
    }
    remove(tmpfile1.c_str());
    remove(tmpfile2.c_str());
}


void OverlapGraph::writeDiGraphToFile() { // store directed edges
    std::string filename = PATH + "digraph.txt";
    if (program_settings.verbose) {
        std::cout << "WriteDiGraphToFile " << filename << "\n";
    }
    remove(filename.c_str());
    std::ofstream graph_file(filename);
	std::list< Edge >::const_iterator it1;
	std::vector< std::list< Edge > >::const_iterator it2;
	node_id_t i = 0;
	node_id_t j;
	for (it2 = adj_out.begin(); it2 != adj_out.end(); it2++) {
	    for (it1 = it2->begin(); it1 != it2->end(); it1++) {
	        j = it1->get_vertex(2);
	        assert (j != i);
		    std::string line1 = std::to_string(i) + "\t" + std::to_string(j) + "\n";
		    graph_file << line1;
	    }
	    i++;
    }
    graph_file.close();
}


void OverlapGraph::write2FASTG() { // store directed edges
    std::string filename = PATH + "graph.fastg";
    if (program_settings.verbose) {
        std::cout << "Write2FASTG " << filename << "\n";
    }
    remove(filename.c_str());
    std::ofstream graph_file(filename);
	std::list< Edge >::const_iterator it1;
	std::vector< std::list< Edge > >::const_iterator it2;
	graph_file << "#FASTG:begin;\n"; // write FASTG header
	graph_file << "#FASTG:version=1.0:assembly_name=\"test\";\n";
	node_id_t i = 0;
	node_id_t j;
	for (it2 = adj_out.begin(); it2 != adj_out.end(); it2++) {
	    unsigned int readcount = (fastq_storage->m_read_vec).size();
	    unsigned int singlecount = fastq_storage->m_readcount_single;
	    unsigned int read_pos;
	    std::string seq;
	    if (i < singlecount) {
	        read_pos = i;
	        Read* read = (fastq_storage->m_read_vec).at(read_pos);
	        assert (read->is_paired() == false);
	        seq = read->get_seq(0);
	    }
	    else if (i >= readcount && i < (readcount + singlecount)) {
	        read_pos = i - readcount;
	        Read* read = (fastq_storage->m_read_vec).at(read_pos);
	        assert (read->is_paired() == false);
	        seq = read->get_rev_comp(0);
	    }
	    else {
	        i++;
	        continue;
	    }
	    // write all adjacencies
		std::string line = ">" + std::to_string(i) + ":";
	    for (it1 = it2->begin(); it1 != it2->end(); it1++) {
	        j = it1->get_vertex(2);
	        assert (j != i);
	        if (j < singlecount || (j >= readcount && j < (readcount + singlecount))) {
	            // only add edges between single-end reads
		        line.append(std::to_string(j) + ",");
		    }
	    }
	    line.pop_back(); // remove ',' or ':'
	    line.append(";\n");
	    graph_file << line;
	    assert (seq.size() > 0);
	    graph_file << seq << "\n";
	    i++;
    }
    graph_file << "#FASTG:end;";
    graph_file.close();
}


void OverlapGraph::write2GFA(std::string filename) { // store directed edges
//    std::string filename = PATH + "graph.gfa";
    if (program_settings.verbose) {
        std::cout << "Write2GFA " << filename << "\n";
        std::cout << "Note: currently only S-S edges are written to GFA." << "\n";
    }
    remove(filename.c_str());
    std::ofstream graph_file(filename);
	std::list< Edge >::const_iterator it1;
	std::vector< std::list< Edge > >::const_iterator it2;
	graph_file << "H\tVN:Z:1.0\n"; // write GFA header
	// write a segment line for every node
	// write a link line for every edge without containment
	// write a containment line for every containment edge
	node_id_t i = 0;
	node_id_t j;
    unsigned int L_count = 0;
    unsigned int C_count = 0;
	for (it2 = adj_out.begin(); it2 != adj_out.end(); it2++) {
	    unsigned int readcount = (fastq_storage->m_read_vec).size();
	    unsigned int singlecount = fastq_storage->m_readcount_single;
	    unsigned int read_pos;
	    std::string seq;
	    if (i < singlecount) {
	        read_pos = i;
	        Read* read = (fastq_storage->m_read_vec).at(read_pos);
	        assert (read->is_paired() == false);
	        seq = read->get_seq(0);
	    }
	    else if (i >= readcount && i < (readcount + singlecount)) {
	        read_pos = i - readcount;
	        Read* read = (fastq_storage->m_read_vec).at(read_pos);
	        assert (read->is_paired() == false);
	        seq = read->get_rev_comp(0);
	    }
	    else {
	        i++;
	        continue;
	    }
	    // write all adjacencies
	    assert (seq.size() > 0);
		std::string line = "S\t" + std::to_string(i) + "\t" + seq + "\n";
	    graph_file << line;
	    for (it1 = it2->begin(); it1 != it2->end(); it1++) {
	        j = it1->get_vertex(2);
	        assert (j != i);
	        // only add edges between single-end reads:
	        if (j < singlecount || (j >= readcount && j < (readcount + singlecount))) {
	            int overlap_len = it1->get_len(0);
	            // check if the overlap is fully contained
	            std::string edge_line;
	            if (it1->get_perc() < 100) { // not fully contained
	                edge_line = "L\t" + std::to_string(i) + "\t+\t" + std::to_string(j) + "\t+\t";
	                edge_line.append(std::to_string(overlap_len) + "M\n");
	                L_count++;
	            }
	            else { // fully contained
                	edge_line = "L\t" + std::to_string(i) + "\t+\t" + std::to_string(j) + "\t+\t";
	                edge_line.append(std::to_string(overlap_len) + "M\n");
//	                edge_line = "C\t" + std::to_string(i) + "\t+\t" + std::to_string(j) + "\t+\t";
//                    int pos = it1->get_pos(1);
//	                edge_line.append(std::to_string(pos) + "\t" + std::to_string(overlap_len) + "M\n");
	                C_count++;
	            }
	            graph_file << edge_line;
		    }
	    }
	    i++;
    }
    graph_file.close();
    if (program_settings.verbose) {
        std::cout << "L_count = " << L_count << "\n";
        std::cout << "C_count = " << C_count << "\n";
        std::cout << "Total count = " << L_count + C_count << "\n";
    }
}


/* Write cycle edge found to file cycles.txt
 */
void OverlapGraph::reportCycle(node_id_t u, node_id_t v, bool remove) {
//    std::cout << "reporting cycle...\n" << u << " to " << v << "\n";
    std::ofstream cycle_file;
    cycle_file.open(PATH + "cycles.txt", std::fstream::out | std::fstream::app);
    if (!cycle_file) {
        std::cerr << "cycles.txt could not be opened for writing. Exiting reportCycle.\n";
        exit(1);
    }
    if (remove) {
        Edge edge = removeEdge(u, v);
        branching_edges.push_back(edge);
    }
    cycle_file << u << "\t" << v << std::endl;
    cycle_file.close();
}


void OverlapGraph::printAdjacencyLists() {
    std::cout << "Adjacency lists of outgoing edges: \n";
    std::list< Edge >::const_iterator it;
    for (node_id_t i=0; i<vertex_count; i++)
    {
        std::cout << i << " : ";
        for (it=adj_out[i].begin(); it!=adj_out[i].end(); it++)
            std::cout << it->get_vertex(2) << " ";
        std::cout << std::endl;
    }
}


void OverlapGraph::checkDuplicateEdges() {
    std::list< Edge >::const_iterator it;
    std::list< Edge >::const_iterator prev_it;
    node_id_t j;
    node_id_t j_prev;
    for (node_id_t i=0; i<vertex_count; i++)
    {
        for (it=adj_out[i].begin(); it!=adj_out[i].end(); it++) {
            j = it->get_vertex(2);
            if (it != adj_out[i].begin()) {
                if (j == j_prev) {
                    std::cout << i << " to " << j << " duplicate...\n";
                    std::cout << it->get_pos(1) << " " << it->get_pos(2) << "\n";
                    std::cout << it->get_ori(1) << " " << it->get_ori(2) << "\n";
                    std::cout << it->get_score() << "\n";
                    std::cout << prev_it->get_score() << "\n";
                    std::cout << (it->get_score() > prev_it->get_score()) << "\n";
                }
                assert (j != j_prev);
            }
            j_prev = j;
            prev_it = it;
        }
    }
    if (program_settings.verbose) {
        std::cout << "Duplicate edge check done.\n";
    }
}


void OverlapGraph::addEquivalentEdges() {
    /* for every edge in G, add the corresponding equivalent edge:
            ++ original -> -- equivalent
            +- original -> -+ equivalent
            -+ original -> +- equivalent
       where +- denotes that node1 is in normal orientation and node2 in reverse orientation
    */
    std::list< Edge > empty_list = {};
    std::vector< std::list< Edge > > extra_edges;
    extra_edges = std::vector< std::list< Edge > > (vertex_count, empty_list);
    for (node_id_t i=0; i < vertex_count; i++) {
        for (auto edge_it : adj_out.at(i)) {
            // compute equivalent edge
            double score = edge_it.get_score();
            bool ori1;
            bool ori2;
            int pos1 = edge_it.get_extra_pos(1);
            int pos2 = edge_it.get_extra_pos(2);
            std::string ord;
            Read* read_1;
            Read* read_2;
            if (pos1 < 0) {
                read_1 = edge_it.get_read(2);
                read_2 = edge_it.get_read(1);
                ori1 = !(edge_it.get_ori(2));
                ori2 = !(edge_it.get_ori(1));
                pos1 = -pos1;
                if (pos2 < 0) {
                    ord = "1";
                    pos2 = -pos2;
                }
                else {
                    if (edge_it.get_ord() == '-' || edge_it.get_ord() == '0') {
                        ord = "-";
                    }
                    else {
                        ord = "2";
                    }
                }
            }
            else {
                read_1 = edge_it.get_read(1);
                read_2 = edge_it.get_read(2);
                ori1 = !(edge_it.get_ori(1));
                ori2 = !(edge_it.get_ori(2));
                if (pos2 < 0) {
                    pos2 = -pos2;
                    ord = "2";
                }
                else {
                    if (edge_it.get_ord() == '-' || edge_it.get_ord() == '0') {
                        ord = "-";
                    }
                    else {
                        ord = "1";
                    }
                }
            }
            Edge edge(score, pos1, pos2, ori1, ori2, ord, read_1, read_2);
            // find the corresponding nodes
            node_id_t node1 = read_1->get_vertex_id(ori1);
            node_id_t node2 = read_2->get_vertex_id(ori2);
            edge.set_vertices(node1, node2);
            edge.set_len(edge_it.get_len(1), edge_it.get_len(2));
            edge.set_perc(edge_it.get_perc());
            // store in extra adjacency list
            extra_edges.at(node1).push_back(edge);
        }
    }
    // add all extra edges to G
    unsigned int count = 0;
    unsigned int doubles = 0;
    for (node_id_t i=0; i < vertex_count; i++) {
        for (auto it : extra_edges.at(i)) {
            // check if edge already exists:
            double score;
            node_id_t v1 = it.get_vertex(1);
            node_id_t v2 = it.get_vertex(2);
            if (it.get_pos(1) == 0 && v1 > v2) { // edge can be added in either direction; always direct them from small to large ID
                if (v1 > v2) {
                    std::swap(v1, v2);
                    it.swap_reads();
                }
            }
            score = checkEdge(v1, v2, /*reverse_allowed*/ false);
            if (score < 0) {
                // edge does not yet exist, so add it now
                addEdge(it);
                count++;
            }
            else if (it.get_score() > score) {
                // new edge scores better, so replace current edge in the graph
                int edgecount1 = getEdgeCount();
                removeEdge(v1, v2);
                int edgecount2 = getEdgeCount();
                addEdge(it);
                int edgecount3 = getEdgeCount();
                assert (edgecount1 == edgecount3);
                assert (edgecount1 == edgecount2 + 1);
                doubles++;
            }
            else {
                // new edge does not improve current edge, so do nothing
                doubles++;
            }
        }
    }
    if (program_settings.verbose) {
        std::cout << "Number of equivalent edges built: " << count << std::endl;
        std::cout << "Number of duplicates: " << doubles << std::endl;
    }
}


void OverlapGraph::sortEdges() {
    if (program_settings.verbose) {
        std::cout << "sortEdges: sort adjacency lists by non-overlap length" << std::endl;
    }
    std::vector< std::list< Edge > > new_adj_out;
    for (auto adj_list : adj_out) {
        // sort adjacency list by non-overlap length
        std::vector< std::pair<Edge, unsigned int> > pairs;
        for (auto edge_it : adj_list) {
            Edge edge = edge_it;
            pairs.push_back(std::make_pair(edge, edge_it.get_nonoverlap_len()));
        }
        std::sort(pairs.begin(), pairs.end(), [=](const std::pair<Edge, unsigned int>& a, const std::pair<Edge, unsigned int>& b)
        {
            if (a.second == b.second) {
                return a.first.get_vertex(2) < b.first.get_vertex(2);
            }
            else {
                return a.second < b.second;
            }
        }
        );
        assert (pairs.size() == adj_list.size());
        std::list< Edge > new_adj_list;
        for (auto pair_it : pairs) {
            Edge edge = pair_it.first;
            new_adj_list.push_back(edge);
        }
        new_adj_out.push_back(new_adj_list);
    }
    adj_out = new_adj_out;
    // rebuild adj_in accordingly
    std::list< node_id_t > empty_list2 = {};
    std::vector< std::list< node_id_t > > new_adj_in = std::vector< std::list< node_id_t > > (vertex_count, empty_list2);
    for (auto adj_list : new_adj_out) {
        for (auto edge_it : adj_list) {
            node_id_t node1 = edge_it.get_vertex(1);
            node_id_t node2 = edge_it.get_vertex(2);
            new_adj_in.at(node2).push_back(node1);
        }
    }
    adj_in = new_adj_in;
}


// Construct a dictionary that stores subread IDs (original fastq)
void OverlapGraph::buildOriginalsDict() {
    if (program_settings.verbose) {
        std::cout << "buildOriginalsDict... ";
    }
    if (program_settings.first_it) { // trivial originals dict
        read_id_t ID;
        std::pair< read_id_t, OriginalIndex > ID_idx_pair;
        OriginalIndex original_index;
        original_index.index1 = 0;
        original_index.forward = true;
        for (auto read_it : fastq_storage->m_singles_vec) {
            original_index.is_paired = 0;
            original_index.len1 = read_it.get_seq(0).size();
            ID = read_it.get_read_id();
            std::unordered_map< read_id_t, OriginalIndex > originals_map;
            ID_idx_pair = std::make_pair(ID, original_index);
            originals_map.insert(ID_idx_pair);
            original_ID_dict.insert(std::make_pair(ID, originals_map));
        }
        for (auto read : fastq_storage->m_paired_vec) {
            original_index.index2 = 0;
            original_index.is_paired = 1;
            original_index.len1 = read.get_seq(1).size();
            original_index.len2 = read.get_seq(2).size();
            ID = read.get_read_id();
            std::unordered_map< read_id_t, OriginalIndex > originals_map;
            ID_idx_pair = std::make_pair(ID, original_index);
            originals_map.insert(ID_idx_pair);
            original_ID_dict.insert(std::make_pair(ID, originals_map));
        }
    }
    else { // create originals dict from subreads file
        std::string line;
        std::stringstream ss;
        std::string filename = PATH + "subreads.txt";
        std::ifstream originals (filename.c_str());
        if (originals.is_open()) {
            std::string readID;
            std::string info;
            while (getline(originals, line)) {
                ss << line;
                getline(ss, readID, '\t');
                read_id_t ID = str_to_read_id(readID);
                std::unordered_map< read_id_t, OriginalIndex > originals_map;
                while (getline(ss, info, '\t')) {
                    std::pair< read_id_t, OriginalIndex > ID_idx_pair;
                    std::vector< std::string > tmp_vec;
                    boost::algorithm::split(tmp_vec, info, boost::is_any_of(":,"), boost::token_compress_on);
                    assert (tmp_vec.size() == 4 || tmp_vec.size() == 6);
                    read_id_t original_ID = str_to_read_id(tmp_vec[0]);
                    OriginalIndex original_index;
                    original_index.forward = (tmp_vec[1] == "+");
                    original_index.index1 = std::stol(tmp_vec[2]);
                    if (tmp_vec.size() == 6) {
                        original_index.index2 = std::stol(tmp_vec[3]);
                        original_index.len1 = std::stoi(tmp_vec[4]);
                        original_index.len2 = std::stoi(tmp_vec[5]);
                        original_index.is_paired = 1;
                    }
                    else {
                        original_index.len1 = std::stoi(tmp_vec[3]);
                        original_index.is_paired = 0;
                    }
                    ID_idx_pair = std::make_pair(original_ID, original_index);
                    originals_map.insert(ID_idx_pair);
                    info.clear();
                }
                original_ID_dict.insert(std::make_pair(ID, originals_map));
                ss << "";
                ss.clear();
            }
            originals.close();
        }
        else {
            std::cerr << "Unable to open subreads file";
            exit(1);
        }
    }
}
