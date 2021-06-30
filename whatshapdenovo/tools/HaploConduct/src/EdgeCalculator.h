//============================================================================
// Name        : EdgeCalculator.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Compute edges from overlaps file by computing overlap scores
//============================================================================

#ifndef EDGECALCULATOR_H_
#define EDGECALCULATOR_H_


#include <cmath>
#include <iostream>

#include "FastqStorage.h"
#include "OverlapGraph.h"
#include "Edge.h"
#include "Overlap.h"
#include "Types.h"


// A class to read the overlaps from file and compute the corresponding edges;
// these edges are added to the overlap graph associated to the edge calculator.
class EdgeCalculator
{
private:
    int N_THREADS; // number of threads used for overlap computations
    ProgramSettings program_settings;
    std::shared_ptr<FastqStorage> fastq_storage;
    std::shared_ptr<OverlapGraph> overlap_graph;
    std::map< read_id_t, unsigned int > &ID_dict = fastq_storage->m_ID_to_index;

    double score(char nt1, char nt2, double p1, double p2, int & mismatch_count);
    Edge compute_overlap(const Overlap &overlap);
    void process_overlaps(std::vector<Overlap> overlaps_vec);

public:
    unsigned int self_overlap_count;
    unsigned int inclusion_count;
    unsigned int dup_count;

    EdgeCalculator(std::shared_ptr<FastqStorage> fastq, std::shared_ptr<OverlapGraph> graph, const ProgramSettings ps)
	{
//        std::cout << "EdgeCalculator is being created.\n";
        N_THREADS = ps.n_threads;
        program_settings = ps;
        fastq_storage = fastq;
        overlap_graph = graph;
        self_overlap_count = 0;
        inclusion_count = 0;
        dup_count = 0;
	}

	~EdgeCalculator(void)
    {
//        std::cout << "EdgeCalculator is being deleted" << std::endl;
    }

    double overlap_score(std::string seq1, std::string seq2, std::string score1, std::string score2, const unsigned int pos, double & mismatch_rate);
    void construct_edges();
    double phred_to_prob(const int phred);
};


#endif /* EDGECALCULATOR_H_ */
