//============================================================================
// Name        : SRBuilder.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Cluster reads in the overlap graph and construct super-reads
//============================================================================

#ifndef SRBUILDER_H_
#define SRBUILDER_H_

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include "Read.h"
#include "FastqStorage.h"
#include "OverlapGraph.h"
#include "EdgeCalculator.h"
#include "Types.h"
#include "Overlap.h"


struct OverlapCandidate {
    Read* SR1;
    Read* SR2;
    node_id_t common_node;
    read_id_t original_id;
};

// A class to build superreads from cliques and store them in a vector.
class SRBuilder
{
private:
    std::shared_ptr<FastqStorage> fastq_storage;
    std::shared_ptr<OverlapGraph> overlap_graph;
    unsigned int minCliqueSize;
    double minQual;
    int minOverlap;
    unsigned int N_THREADS;
    std::deque<Read> single_SR_vec;
    std::deque<Read> paired_SR_vec;
    std::deque<Read> trivial_SR_vec;
    std::deque<Read> tips_vec;
    std::map<node_id_t, read_id_t> nodes_to_new_IDs; // dictionary to get from nodes to new IDs
    std::string PATH;
    ProgramSettings program_settings;
    boost::dynamic_bitset<> visited;
    std::deque< std::vector< Read* > > nodes_to_SR;
    node_id_t new_read_count; // including trivial superreads

    double phred_to_prob(const int phred);
    int sort_vertices(std::vector< node_id_t > vertices, char type, node_id_t base_node, std::list<int> &pos_list, std::list<std::string> &seq_list, std::list<std::string> &qual_list, std::list<node_id_t> &sorted_vertices, int thread_id);
    int consensus(int total_len, std::list<int> &pos_list, std::list<std::string> &seq_list, std::list<std::string> &qual_list, std::string &cons_seq, std::string &cons_qual, bool subreads_needed, bool error_correction);
    bool consensus_pos(std::string nucleotides, std::string qualities, std::string &cons_seq, std::string& cons_qual);
    Read constructSuperread(std::vector<node_id_t> clique, read_id_t id, int thread_id);
    node_id_t process_cliques(const std::vector< std::vector<node_id_t> >& clique_vec, read_id_t& count);
    std::unordered_map< node_id_t, SubreadInfo > calcSubreadInfo(int trim_pos1, int trim_pos2, std::list<int> pos_list1, std::list<int> pos_list2, std::list<node_id_t> sorted_vertices1, std::list<node_id_t> sorted_vertices2);
    void filter_subreads(int num, node_id_t base_node, std::list< node_id_t > & sorted_vertices, std::list<int> & pos_list, std::list< std::string > & seq_list, std::list< std::string > & qual_list, std::list<int> & new_pos_list, std::list< std::string> & new_seq_list, std::list< std::string > & new_qual_list);
    std::vector< node_id_t > sortVerticesByEndpos(std::vector< std::pair<node_id_t, int> > pairs);
    Read merge_self_overlap(Read superread, EdgeCalculator & edge_calculator);

    // FindNextOverlaps.cpp: induce overlaps from current edges
    std::vector< std::set< read_id_t > > overlaps_found;
    void updateOverlap(Edge edge_info, edge_count_t& copied_count, edge_count_t& u2SR_count, edge_count_t& v2SR_count, edge_count_t& SR2SR_count, std::set< std::string >& overlap_set);
    int findCliqueIndex(node_id_t node, Read* superread, bool leftside, bool second_occ);
    bool computeOverlapData(Read* superread1, Read* superread2, int idx1l, int idx1r, int idx2l, int idx2r, Edge edge, int &new_pos1, int &new_pos2, char &ord1, char &ord2, std::string &type1, std::string &type2, int& overlap_perc, int& overlap_len1, int& overlap_len2);
    void processOverlaps(const std::vector<Edge>& edge_vec, edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set);
    void reconsiderEdgeOverlaps(edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set);
    void reconsiderNonedgeOverlaps(edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set);
    void findInclusionOverlaps(edge_count_t& total_copied_count, edge_count_t& total_u2SR_count, edge_count_t& total_v2SR_count, edge_count_t& total_SR2SR_count, std::set< std::string >& final_overlap_set);
    // FindNextOverlaps2.cpp: induce overlaps from common sub-reads
    Overlap deduceOverlap(OverlapCandidate candidate);
    boost::dynamic_bitset<> buildBitvec(Read* superread);
    void nodeDictApproach();
    void bitvecApproach();
    // FindNextOverlaps3.cpp: induce overlaps from common ORIGINAL reads
    Overlap getInducedOverlap(Read* SR1, Read* SR2);
    void nodeDictApproach(std::unordered_map< read_id_t, node_id_t > original_to_index);
    Overlap deduceOverlap(OverlapCandidate candidate, read_id_t original_id);

public:
	SRBuilder(std::shared_ptr<FastqStorage> fastq, std::shared_ptr<OverlapGraph> graph, ProgramSettings ps) {
        fastq_storage = fastq;
        overlap_graph = graph;
		minCliqueSize = ps.min_clique_size;
		minQual = ps.min_qual;
//		N_THREADS = ps.n_threads;
        N_THREADS = 1;
		PATH = ps.output_dir;
		program_settings = ps;
//		std::cout << "SRBuilder is being created.\n";
        node_id_t V = overlap_graph->getVertexCount();
        if (program_settings.fno == 3) {
            node_id_t s = program_settings.original_readcount;
            nodes_to_SR = std::deque< std::vector< Read* >> (s, std::vector< Read*>());
        }
        else {
            nodes_to_SR = std::deque< std::vector< Read* >> (V, std::vector< Read*>());
        }
        new_read_count = 0;
        clique_count = 0;
        next_overlaps_count = 0;
	}

	~SRBuilder(void) {
//        std::cout << "SRBuilder is being deleted" << std::endl;
    }

    void writeTrivialsToFile();
    void writeTipsToFile();
    void writeSinglesToFile(std::deque<Read>& singles, read_id_t& count);
    void writePairsToFile(std::deque<Read>& pairs, read_id_t& count);

    void cliquesToSuperreads();
    void mergeAlongEdges();

    // FindNextOverlaps.cpp:
    unsigned long findNextOverlaps();
    // FindNextOverlaps2.cpp:
    void findNextOverlaps2();
    // FindNextOverlaps3.cpp:
    void findNextOverlaps3();

    // statistics for log file
    node_id_t clique_count;
    node_id_t SR_singles_count;
    node_id_t SR_paired_count;
    node_id_t SR_trivials_count;
    edge_count_t next_overlaps_count;
};



#endif /* SRBUILDER_H_ */
