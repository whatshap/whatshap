//============================================================================
// Name        : ViralQuasispecies.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : De novo viral quasispecies assembly using overlap graphs
//============================================================================

#include <iostream>
#include <string>
#include <vector>
#include <sys/time.h>
#include <cmath>
#include <map>
#include <omp.h>
#include <cstdlib>
#include <fstream>
#include <boost/program_options.hpp>
#include <unistd.h>
#include <limits.h>

#include "Types.h"
#include "EdgeCalculator.h"
#include "FastqStorage.h"
#include "OverlapGraph.h"
#include "SRBuilder.h"
#include "Edge.h"
#include "BranchReduction.h"


// define timestamp function to measure runtimes
typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp() {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t) now.tv_sec * 1000000;
}

int main(int argc, char *argv[])
{
    namespace po = boost::program_options;

    // Declare the supported options in a struct (see Types.h)
    ProgramSettings program_settings;

    // Gives the options of the input file.
    po::options_description desc("Program options");
    desc.add_options()
        ("help", "produce help message")
        ("fastq", po::value< std::string > (&program_settings.fastq_file), "path to fastq files: paired_1.fastq, paired_2.fastq and single.fastq")
        ("singles,s", po::value< std::string > (&program_settings.singles_file)->default_value(""), "path to single-end read fastq file")
        ("paired1", po::value< std::string > (&program_settings.paired1_file)->default_value(""), "path to paired-end read /1 fastq file")
        ("paired2", po::value< std::string > (&program_settings.paired2_file)->default_value(""), "path to paired-end read /2 fastq file")
        ("overlaps", po::value< std::string > (&program_settings.overlaps_file), "path to overlap file")
        ("output,O", po::value< std::string > (&program_settings.output_dir)->default_value(""), "path to output files")
        ("IDs", po::value< std::string > (&program_settings.id_correspondence)->default_value(""), "path to ID correspondence file")
        ("max_ov,MO", po::value< unsigned long > (&program_settings.max_overlaps)->default_value(100000000), "set the maximum number of overlaps considered")
        ("max_reads,MR", po::value< unsigned long > (&program_settings.max_reads)->default_value(100000000), "set the maximum number of reads used")
        ("threads,t", po::value< unsigned int > (&program_settings.n_threads)->default_value(1), "set the number of threads used")
        ("min_clique_size", po::value< unsigned int > (&program_settings.min_clique_size)->default_value(4), "set the minimum clique size for a superread")
        ("min_qual", po::value< double > (&program_settings.min_qual)->default_value(0.9), "set the minimum base quality (probability that it is correct) for a superread; otherwise an 'N' is inserted.")
        ("min_overlap_perc", po::value< unsigned int > (&program_settings.min_overlap_perc)->default_value(0), "set the minimum overlap percentage")
        ("min_overlap_len", po::value< unsigned int > (&program_settings.min_overlap_len)->default_value(150), "set the minimum overlap length (bp)")
        ("edge_threshold", po::value< double > (&program_settings.edge_threshold)->default_value(0.99), "set the minimal overlap score for creating an edge")
        ("ov_threshold", po::value< double > (&program_settings.ov_threshold)->default_value(0.9), "set the minimal overlap score for keeping non-edge overlap")
        ("allow_spaced_overlaps", po::value< bool > (&program_settings.allow_spaces)->default_value(false), "allow space-delimited overlaps instead of tabs")
        ("first_it", po::value< bool > (&program_settings.first_it)->default_value(true), "set to true when there is no subreads file")
        ("add_duplicates", po::value< bool > (&program_settings.add_duplicates)->default_value(false), "set to true when you prefer to deal with reverse complements by adding duplicate vertices")
        ("resolve_orientations", po::value< bool > (&program_settings.resolve_orientations)->default_value(true), "set to true when you prefer to deal with reverse complements by labelling vertices")
        ("keep_singletons", po::value< unsigned int > (&program_settings.keep_singletons)->default_value(0), "minimal read length for singletons not to be removed, should be 0 at pre-iterations")
        ("error_correction", po::value< bool > (&program_settings.error_correction)->default_value(false), "set to true when you only want to do error correction")
        ("cliques", po::value< bool > (&program_settings.cliques)->default_value(false), "set to true for clique-merging instead of edge-merging")
        ("ignore_inclusions", po::value< bool > (&program_settings.ignore_inclusions)->default_value(false), "set to true in order to ignore full inclusion edges in overlap graph")
        ("graph_only", po::value< bool > (&program_settings.graph_only)->default_value(false), "set to true when you only want to do the graph construction")
        ("FNO", po::value< int > (&program_settings.fno)->default_value(2), "set the FindNextOverlaps function desired")
        ("original_readcount", po::value< node_id_t > (&program_settings.original_readcount), "the number of original reads")
        ("mismatch", po::value< double > (&program_settings.mismatch)->default_value(0), "minimal score per position in overlap")
        ("optimize", po::value< bool > (&program_settings.optimize)->default_value(true), "optimize FNO by not reconsidering non-edge overlaps")
        ("no_inclusion_overlaps", po::value< bool > (&program_settings.no_inclusions)->default_value(false), "do not add full inclusion overlaps")
        ("merge_contigs", po::value< double > (&program_settings.merge_contigs)->default_value(0), "allow edge construction based on <merge_contigs> mismatch rate instead if overlap score is insufficient")
        ("remove_multi_occ", po::value< bool > (&program_settings.remove_multi_occ)->default_value(false), "remove clique nodes when used before, so use nodes at most once; to be used at merging iterations using cliques")
        ("remove_trans", po::value< unsigned int > (&program_settings.remove_trans)->default_value(0), "choose to (0) keep all edges, (1) remove transitive edges, (2) to remove double transitive edges, or (3) to remove triple transitive edges.")
        ("remove_branches", po::value< bool > (&program_settings.remove_branches)->default_value(false), "remove branches from overlap graph")
        ("remove_tips", po::value< bool > (&program_settings.remove_tips)->default_value(true), "remove tips from overlap graph to reduce branching")
        ("min_read_len", po::value< unsigned int > (&program_settings.min_read_len)->default_value(0), "set the minimum read length (bp) for allowing edges")
        ("max_tip_len", po::value< unsigned int > (&program_settings.max_tip_len)->default_value(150), "set the maximum extension length for a node to be considered a tip")
        ("separate_tips", po::value< bool > (&program_settings.store_tips_separately)->default_value(true), "store tip-sequences in a separate file, away from the contigs")
        ("base_path", po::value< std::string > (&program_settings.base_path)->default_value("."), "set path to SAVAGE directory containing quick-cliques-1.0")
        ("diploid", po::value< bool > (&program_settings.diploid)->default_value(false), "apply edge filtering for diploid genomes")
        ("relax_PE_edges", po::value< bool > (&program_settings.relax_PE_edges)->default_value(false), "relax edge restrictions for paired-end overlaps")
        ("original_fastq", po::value< std::string > (&program_settings.original_fastq)->default_value(""), "original reads for applying read-based branch reduction")
        // ("branch_min_ev", po::value< int > (&program_settings.branch_min_ev)->default_value(0), "minimum evidence required when applying read-based branch reduction")
        ("branch_reduction", po::value< bool > (&program_settings.branch_reduction)->default_value(false), "average coverage (sequencing depth) per haplotype, necessary information for read-based branch reduction (skipped otherwise)")
        ("branch_SE_c", po::value< unsigned int > (&program_settings.branch_SE_c)->default_value(0), "number of single-end input reads in original fastq")
        ("branch_PE_c", po::value< unsigned int > (&program_settings.branch_PE_c)->default_value(0), "number of paired-end input reads in original fastq")
        ("careful_diploid", po::value< bool > (&program_settings.careful)->default_value(true), "more careful merging by avoiding neighboring components")
        ("verbose,v", po::value< bool > (&program_settings.verbose)->default_value(false), "output additional information during assembly")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    po::notify(vm);

    if (!(vm.count("fastq") || vm.count("singles") || vm.count("paired1") || vm.count("paired2"))) {
        std::cerr << "No fastq file(s) provided.\n\n";
        std::cout << desc << "\n";
        return 1;
    }
    else if ((vm.count("fastq") && program_settings.singles_file.size() > 0) || (vm.count("fastq") && program_settings.paired1_file.size() > 0) || (vm.count("fastq") && program_settings.paired2_file.size() > 0)) {
        std::cerr << "Cannot combine --fastq option with --singles, --paired1 or --paired2. \n\n";
        std::cout << desc << "\n";
        return 1;
    }

    if ((vm.count("paired1") && !vm.count("paired2")) || (!vm.count("paired1") && vm.count("paired2"))) {
        std::cerr << "Only one paired-end fastq file provided.\n\n";
        std::cout << desc << "\n";
        return 1;
    }

    if (!vm.count("overlaps")) {
        std::cerr << "No overlaps file provided.\n\n";
        std::cout << desc << "\n";
        return 1;
    }

    if (!vm.count("original_readcount")) {
        std::cerr << "No original readcount provided.\n\n";
        std::cout << desc << "\n";
        return 1;
    }

    // if (program_settings.id_correspondence == "" && program_settings.first_it) {
    //     std::cout << "\nWARNING: no ID-correspondence provided at first iteration. Use the option --IDs <id_correspondence>.\n";
    // }

    if (program_settings.add_duplicates && program_settings.resolve_orientations) {
        std::cerr << "Add duplicates and resolve orientations are exclusive options, use at most 1.\n\n";
        std::cout << desc << "\n";
        return 1;
    }

    if (program_settings.error_correction && !program_settings.cliques) {
        std::cerr << "Error correction requires clique enumeration. Set --cliques=true.\n";
        std::cout << desc << "\n";
        return 1;
    }

    if (program_settings.remove_trans == 1 && program_settings.min_clique_size > 2 && program_settings.cliques) {
        std::cout << "WARNING: Removing transitive edges while minimum clique size is larger than 2; there cannot be any such cliques in the resulting overlap graph. Reduce --min_clique_size.\n";
    }
    else if (program_settings.remove_trans == 2 && program_settings.min_clique_size > 4 && program_settings.cliques) {
        std::cout << "WARNING: Removing double transitive edges while minimum clique size is larger than 4; there cannot be any such cliques in the resulting overlap graph. Reduce --min_clique_size.\n";
    }
    else if (program_settings.remove_trans == 3 && program_settings.min_clique_size > 8 && program_settings.cliques) {
        std::cout << "WARNING: Removing triple transitive edges while minimum clique size is larger than 8; there cannot be any such cliques in the resulting overlap graph. Reduce --min_clique_size.\n";
    }

    // create log file
    std::ofstream logfile;
    logfile.open(program_settings.output_dir + "viralquasispecies.log", std::fstream::out);
    time_t rawtime;
    time (&rawtime);
    logfile << ctime (&rawtime) << std::endl;

    // write parameter values to logfile
    logfile << "\nInput:\n";
    logfile << program_settings.singles_file << "\n";
    logfile << program_settings.paired1_file << "\n";
    logfile << program_settings.paired2_file << "\n";
    logfile << program_settings.overlaps_file << "\n\n";
    logfile << "Output directory: " << program_settings.output_dir << "\n";
    logfile << "Base path: " <<  program_settings.base_path << "\n";
    logfile << "Maximum number of overlaps: " <<  program_settings.max_overlaps << "\n";
    logfile << "Threads: " <<  program_settings.n_threads << "\n";
    logfile << "Minimum clique size: " <<  program_settings.min_clique_size << "\n";
    logfile << "Minimum base quality (superreads): " <<  program_settings.min_qual << "\n";
    logfile << "Minimal overlap percentage: " <<  program_settings.min_overlap_perc << "\n";
    logfile << "Minimal overlap length: " <<  program_settings.min_overlap_len << "\n";
    logfile << "Edge threshold: " <<  program_settings.edge_threshold << "\n";
    logfile << "Overlap threshold: " <<  program_settings.ov_threshold << "\n";
    logfile << "First it: " <<  program_settings.first_it << "\n";
    logfile << "Add duplicates: " <<  program_settings.add_duplicates << "\n";
    logfile << "Resolve read orientations: " <<  program_settings.resolve_orientations << "\n";
    logfile << "Keep singletons of length at least: " <<  program_settings.keep_singletons << "\n";
    logfile << "Enumerate cliques: " <<  program_settings.cliques << "\n";
    logfile << "Error correction: " <<  program_settings.error_correction << "\n";
    logfile << "Ignore inclusions: " <<  program_settings.ignore_inclusions << "\n";
    logfile << "Graph only: " << program_settings.graph_only << "\n";
    logfile << "FNO: " <<  program_settings.fno << "\n";
    logfile << "Mismatch prob: " <<  program_settings.mismatch << "\n";
    logfile << "Original readcount: " <<  program_settings.original_readcount << "\n";
    logfile << "Optimize: " << program_settings.optimize << "\n";
    logfile << "No inclusion overlaps: " << program_settings.no_inclusions << "\n";
    logfile << "Merge contigs: " << program_settings.merge_contigs << "\n";
    logfile << "Remove_multi_occ: " << program_settings.remove_multi_occ << "\n";
    logfile << "Remove transitive edges: " << program_settings.remove_trans << "\n";
    logfile << "Minimal read length: " <<  program_settings.min_read_len << "\n";
    logfile << "Remove branches: " <<  program_settings.remove_branches << "\n";
    logfile << "Remove tips: " <<  program_settings.remove_tips << "\n";
    logfile << "Maximal tip length: " <<  program_settings.max_tip_len << "\n";
    logfile << "Store tips separately: " <<  program_settings.store_tips_separately << "\n";
    logfile << "Relax PE edges: " <<  program_settings.relax_PE_edges << "\n";
    logfile << "Original fastq: " <<  program_settings.original_fastq << "\n";
    logfile << "Branch reduction: " <<  program_settings.branch_reduction << "\n";
    logfile << "Branch SE count: " <<  program_settings.branch_SE_c << "\n";
    logfile << "Branch PE count: " <<  program_settings.branch_PE_c << "\n";
    logfile << "Diploid: " <<  program_settings.diploid << "\n";
    logfile << "Careful diploid mode: " <<  program_settings.careful << "\n";
    logfile << "Verbose: " <<  program_settings.verbose << "\n";
    logfile.close();

    // define timestamps
    timestamp_t t0, t1;
    double time_s;

    // Read fastq file to vectors in fastq storage.
    t0 = get_timestamp();

    if (vm.count("fastq")) {
        program_settings.singles_file = program_settings.fastq_file + "/singles.fastq";
        program_settings.paired1_file = program_settings.fastq_file + "/paired1.fastq";
        program_settings.paired2_file = program_settings.fastq_file + "/paired2.fastq";
    }

    FastqStorage* input_fastq = new FastqStorage(program_settings);
    std::shared_ptr<FastqStorage> fastq_storage(input_fastq);
    t1 = get_timestamp();
    time_s = (t1 - t0) / 1000000.0L;
    if (program_settings.verbose) {
        std::cout << "FastqStorage ready! Construction took " << time_s << " seconds.\n";
    }
//	fastq_storage->writeIDsToFile("readIDs.txt"); // write IDs to file such that overlaps can be filtered
//    return 0;

    // Construct overlap graph and add a vertex for each read in fastq storage.
    t0 = get_timestamp();
    node_id_t graph_size;
    if (program_settings.add_duplicates) {
        graph_size = 2*(fastq_storage->get_readcount());
    }
    else {
        graph_size = fastq_storage->get_readcount();
    }
    OverlapGraph* graph = new OverlapGraph(graph_size, fastq_storage, program_settings);
    std::shared_ptr<OverlapGraph> overlap_graph(graph);
    if (program_settings.verbose) {
        std::cout << "Adding vertices...\n";
    }

    std::vector<Read *> read_vector = fastq_storage->m_read_vec; // contains all reads: first single-end, then paired-end
    for (auto read_iterator : read_vector){
        read_id_t id = read_iterator->get_read_id();
        node_id_t vertex_id = overlap_graph->addVertex(id);
        read_iterator->set_vertex_id(/*normal=*/ true, vertex_id);
    }

    if (program_settings.add_duplicates) { // also add a vertex for each reverse complementary read
        for (auto read_iterator : read_vector) {
            read_id_t id = read_iterator->get_read_id();
            node_id_t vertex_id = overlap_graph->addVertex(id);
            read_iterator->set_vertex_id(/*normal=*/ false, vertex_id);
        }
    }
    t1 = get_timestamp();
    time_s = (t1 - t0) / 1000000.0L;
    if (program_settings.verbose) {
        std::cout << "Overlap graph ready! Construction took " << time_s << " seconds.\n";
        std::cout << "Number of vertices: " << overlap_graph->getVertexCount() << std::endl;
    }
    // Construct edges by computing overlap scores for the overlaps given.
    EdgeCalculator edge_calculator(fastq_storage, overlap_graph, program_settings);
    t0 = get_timestamp();
    edge_calculator.construct_edges();
    t1 = get_timestamp();
    time_s = (t1 - t0) / 1000000.0L;
    if (overlap_graph->getEdgeCount() == 0) {
        if (program_settings.verbose) {
            std::cout << "There were no edges constructed, so there is nothing to be done." << std::endl;
        }
        std::string graphfile = program_settings.output_dir + "graph.txt";
        remove(graphfile.c_str()); // remove graph file to ensure pipeline termination
        return 0;
    }
    else if (program_settings.verbose) {
        std::cout << overlap_graph->getEdgeCount() << " edges have been constructed in " << time_s << " seconds.\n";
    }

    // Add vertex labels indicating read orientations
    overlap_graph->sortEdges(); // sort edges
    unsigned int conflict_count;
    overlap_graph->vertexLabellingHeuristic(conflict_count);
    // overlap_graph->printAdjacencyLists();
    overlap_graph->checkDuplicateEdges();

    // overlap_graph->getGraphStats();
    if (program_settings.ignore_inclusions) {
        overlap_graph->removeInclusions();
    }

    // Remove transitive edges as specified by program settings, if any
    overlap_graph->removeTransitiveEdges();

    // build a GFA file for analyzing the graph with Bandage
    overlap_graph->write2GFA(program_settings.output_dir + "graph.gfa");

    overlap_graph->buildOriginalsDict();

    // // Remove edges to obtain a diploid assembly
    // if (program_settings.diploid) {
    //     overlap_graph->reduceDiploidBranching();
    // }

    // Remove tips
    if (program_settings.remove_tips) {
        overlap_graph->removeTips();
    }
    // Reduce branches in the graph by evaluating read evidence
    if (program_settings.branch_reduction > 0) {
        // Adjust program settings to read from original input files
        ProgramSettings original_input = program_settings;
        original_input.singles_file = program_settings.original_fastq;
        original_input.paired1_file = "";
        original_input.paired2_file = "";
        // Now read and store original fastq file(s)
        FastqStorage* input_fastq = new FastqStorage(original_input);
        std::shared_ptr<FastqStorage> original_fastq(input_fastq);
        if (program_settings.branch_SE_c + 2*program_settings.branch_PE_c
            != original_fastq->get_readcount()) {
                std::cout << "branch_SE_c: " << program_settings.branch_SE_c << std::endl;
                std::cout << "branch_PE_c: " << program_settings.branch_PE_c << std::endl;
                std::cout << "original fastq: " << program_settings.original_fastq << std::endl;
            }
        assert (program_settings.branch_SE_c + 2*program_settings.branch_PE_c
            == original_fastq->get_readcount());
        // Use original fastq to resolve branches using read evidence
        BranchReduction branch_reduction(
            fastq_storage, original_fastq, program_settings, overlap_graph);
        branch_reduction.readBasedBranchReduction();
    }
    else if (program_settings.remove_branches) {
        // Find and remove all branches
        overlap_graph->removeBranches();
    }
    // else {
    //     overlap_graph->removeTips();
    //     overlap_graph->findBranches();
    // }

    // Find cycles in the overlap graph
    t0 = get_timestamp();
    overlap_graph->sortEdges(); // sort edges
    bool remove_backedges;
    if (program_settings.error_correction) {
        remove_backedges = false;
    }
    else {
        remove_backedges = true;
    }
    overlap_graph->cycleRemovalHeuristic(remove_backedges);
    t1 = get_timestamp();
    time_s = (t1 - t0) / 1000000.0L;
    if (program_settings.verbose) {
        std::cout << "Cycle detection (DFS) took " << time_s << " seconds.\n";
    }
    // write to log file
    logfile.open(program_settings.output_dir + "viralquasispecies.log", std::fstream::out | std::fstream::app);
    logfile << "Remove backedges: " << remove_backedges << std::endl;
    logfile << std::endl;
    logfile << "\nOutput:" << std::endl;
    logfile << "Vertex count: " << overlap_graph->getVertexCount() << std::endl;
    logfile << "Edge count: " << overlap_graph->getEdgeCount() << std::endl;
    logfile << "Duplicate overlaps: " << edge_calculator.dup_count << std::endl;
    logfile << "Self-overlaps: " << edge_calculator.self_overlap_count << std::endl;
    logfile << "Inclusions: " << edge_calculator.inclusion_count << std::endl;
    logfile << "Conflicting edges removed: " << conflict_count << std::endl;
    logfile << "Number of backedges (DFS): " << overlap_graph->getBackEdgeCount() << "\n" << std::endl;
    logfile.close();

//    if (program_settings.cliques) {
    overlap_graph->writeGraphToFile(); // write current graph to file for quick-cliques
    overlap_graph->write2GFA(program_settings.output_dir + "graph_trimmed.gfa"); // build a GFA file for analyzing the graph with Bandage
    if (program_settings.graph_only) {
        overlap_graph->writeDiGraphToFile(); // write digraph to file for comparison between runs
    }
    if (program_settings.graph_only) {
        return 0;
    }

    if (program_settings.cliques) {
        // Process the resulting overlap graph: run quick-cliques degeneracy for enumerating all maximal cliques
//        std::string command = program_settings.base_path + "/quick-cliques-1.0/bin/degeneracy <" + program_settings.output_dir + "graph.txt > " + program_settings.output_dir + "cliques.txt";
        std::string command = program_settings.base_path + "/quick-cliques/bin/qc --algorithm=degeneracy --input-file=" + program_settings.output_dir + "graph.txt > " + program_settings.output_dir + "cliques.txt";
        if (!program_settings.verbose) {
            command.append(" 2> /dev/null");
        }
        int system_ret = system(command.c_str());
        if(system_ret != 0){
            // The system method failed
            std::cerr << "ERROR: Quick cliques failed. Exiting..." << std::endl;
            exit(1);
        }
    }
//    return 0;

    // Compute consensus sequences from cliques.
    SRBuilder* SR_input = new SRBuilder(fastq_storage, overlap_graph, program_settings);
    std::shared_ptr<SRBuilder> superread_builder(SR_input);
//    superread_builder->buildOriginalsDict();
    if (program_settings.cliques) {
        if (program_settings.verbose) {
            std::cout << "Building superreads from clique file...\n";
        }
        t0 = get_timestamp();
        superread_builder->cliquesToSuperreads();
        t1 = get_timestamp();
        time_s = (t1 - t0) / 1000000.0L;
        if (program_settings.verbose) {
            std::cout << "Superread builder took " << time_s << " seconds.\n";
        }
    }
    else {
        if (program_settings.verbose) {
            std::cout << "Merging reads along edges...\n";
        }
        t0 = get_timestamp();
        overlap_graph->sortEdges(); // sort edges
        t1 = get_timestamp();
        time_s = (t1 - t0) / 1000000.0L;
        if (program_settings.verbose) {
            std::cout << "Sorting edges took " << time_s << " seconds.\n";
        }
        t0 = get_timestamp();
        superread_builder->mergeAlongEdges();
        t1 = get_timestamp();
        time_s = (t1 - t0) / 1000000.0L;
        if (program_settings.verbose) {
            std::cout << "Superread builder took " << time_s << " seconds.\n";
        }
    }

    // Find new overlaps for next iteration
    t0 = get_timestamp();
    std::string overlapcount;
    if (program_settings.fno == 1) {
        unsigned long count = superread_builder->findNextOverlaps();
        overlapcount = std::to_string(count);
    }
    else if (program_settings.fno == 2) {
        std::cout << "FNO2 no longer supported, use FNO1 or FNO3" << std::endl;
//        superread_builder->findNextOverlaps2();
        overlapcount = ".";
    }
    else {
        superread_builder->findNextOverlaps3();
        overlapcount = ".";
    }

    t1 = get_timestamp();
    time_s = (t1 - t0) / 1000000.0L;
    if (program_settings.verbose) {
        std::cout << "FindNextOverlaps took " << time_s << " seconds.\n";
    }
    // Write statistics to file
    std::string stats_filename = program_settings.output_dir + "stats.txt";
    std::ofstream stats_file(stats_filename, std::fstream::app);
    if (stats_file.is_open()) {
        stats_file << overlap_graph->getVertexCount()
            << "\t" << overlap_graph->getEdgeCount()
            << "\t" << overlapcount
            << "\n";
    }

    // write to log file
    logfile.open(program_settings.output_dir + "viralquasispecies.log", std::fstream::out | std::fstream::app);
    logfile << "Total clique count: " << superread_builder->clique_count << std::endl;
    logfile << "Superread (singles) count: " << superread_builder->SR_singles_count << std::endl;
    logfile << "Superread (paired) count: " << superread_builder->SR_paired_count << std::endl;
    logfile << "Superread (trivials) count: " << superread_builder->SR_trivials_count << std::endl;
    logfile << "Next overlaps count: " << superread_builder->next_overlaps_count << std::endl;
    logfile << "\n*****************************************\n" << std::endl;
    logfile.close();

    return 0;
}
