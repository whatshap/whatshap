//============================================================================
// Name        : Types.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Types used for ViralQuasispecies program
//============================================================================

#ifndef TYPES_H_
#define TYPES_H_

#include <string>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include <cstdlib>

struct ProgramSettings {
    std::string fastq_file;
    std::string singles_file;
    std::string paired1_file;
    std::string paired2_file;
    std::string overlaps_file;
    std::string output_dir;
    std::string id_correspondence;
    unsigned long max_overlaps;
    unsigned int n_threads;
    unsigned long max_reads;
    unsigned int min_clique_size;
    double min_qual;
    unsigned int min_overlap_perc;
    unsigned int min_overlap_len;
    double edge_threshold;
    double ov_threshold;
    bool allow_spaces;
    bool first_it;
    bool add_duplicates;
    bool resolve_orientations;
    unsigned int keep_singletons;
    bool error_correction;
    bool cliques;
    bool graph_only;
    int fno;
    unsigned long original_readcount;
    bool ignore_inclusions;
    double mismatch;
    bool optimize;
    bool no_inclusions;
    double merge_contigs;
    bool remove_multi_occ;
    unsigned int remove_trans;
    bool remove_branches;
    bool remove_tips;
    unsigned int min_read_len;
    std::string base_path;
    bool verbose;
    bool diploid;
    unsigned int max_tip_len;
    bool store_tips_separately;
    bool relax_PE_edges;
    std::string original_fastq;
    bool branch_reduction;
    unsigned int branch_SE_c;
    unsigned int branch_PE_c;
    bool careful;
};

struct IterationStats {
    unsigned long vertex_count;
    unsigned long edge_count;
    int max_read_len;
    int min_read_len;
    int average_read_len;
};

struct SubreadInfo {
    int index1;
    int index2;
    int startpos1;
    int startpos2;
};

struct OriginalIndex {
    long index1; // index of original read (/1) within contig (can be negative due to error correction)
    long index2; // index of original read (/2) within contig
    bool is_paired;
    bool forward; // orientation of original read within contig
    int len1; // length of original sequence (/1)
    int len2; // length of original sequence (/2)
};

typedef unsigned long read_id_t;
typedef unsigned long node_id_t;
typedef unsigned long edge_count_t;
typedef unsigned long long safe_edge_count_t;
typedef std::pair< node_id_t, node_id_t > node_pair_t;

inline read_id_t str_to_read_id(std::string input_id)
{
    return strtoul(input_id.c_str(), NULL, 0);
}

inline std::string read_id_to_str(read_id_t id)
{
    return std::to_string(id);
}

inline std::string build_rev_comp(std::string seq) {
    std::string::reverse_iterator it;
    std::string rev_comp = "";
    for (it = seq.rbegin(); it != seq.rend(); it++) {
        if (*it == 'A')
            rev_comp.append("T");
        else if (*it == 'T')
            rev_comp.append("A");
        else if (*it == 'C')
            rev_comp.append("G");
        else if (*it == 'G')
            rev_comp.append("C");
        else if (*it == 'N')
            rev_comp.append("N");
        else {
            std::cerr << "Invalid sequence character. Aborting.\n";
            exit(1);
        }
    }
    return rev_comp;
}

#endif /* TYPES_H_ */
