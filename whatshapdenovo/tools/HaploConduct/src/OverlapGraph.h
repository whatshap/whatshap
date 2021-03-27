//============================================================================
// Name        : OverlapGraph.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Construct an overlap graph for viral quasispecies assembly
//============================================================================

#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_

#include <list>
#include <stack>
#include <vector>
#include <set>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <sstream>

#include "Overlap.h"
#include "Types.h"
#include "Edge.h"
#include "FastqStorage.h"
//#include "SRBuilder.h"

class SRBuilder; // forward definition to enable friend class
class BranchReduction;

// Class to represent an overlap graph and the algorithms necessary to make it cycle-free.
// The maximal number of vertices must be specified at initialization to reserve the memory required.
class OverlapGraph
{
friend class SRBuilder;
friend class BranchReduction;

private:
    unsigned int vertex_count;
    unsigned int edge_count;
    unsigned int backedge_count;
    unsigned int twocycle_count;
    unsigned int threecycle_count;
    unsigned int graph_depth;
    std::vector< std::list< Edge > > adj_out;
    std::vector< std::list< node_id_t > > adj_in;
    std::vector< Edge > branching_edges;
    std::string PATH;
    std::shared_ptr<FastqStorage> fastq_storage;
    ProgramSettings program_settings;
    std::map<read_id_t, std::unordered_map< read_id_t, OriginalIndex > > original_ID_dict; // dict from current ID to original subread IDs and indexes

public:
    OverlapGraph(unsigned int V, std::shared_ptr<FastqStorage> fastq, ProgramSettings ps) {
        vertex_count = 0;
        edge_count = 0;
        backedge_count = 0;
        twocycle_count = 0;
        threecycle_count = 0;
        graph_depth = 0;
        PATH = ps.output_dir;
        fastq_storage = fastq;
        program_settings = ps;

//        std::cout << "OverlapGraph of size " << V << " is being created.\n";
        std::list< Edge > empty_list = {};
        std::list< node_id_t > empty_list2 = {};
        adj_out = std::vector< std::list< Edge > > (V, empty_list);
        adj_in = std::vector< std::list< node_id_t > > (V, empty_list2);
        inclusions = boost::dynamic_bitset<>(V);
    }

    ~OverlapGraph() {
//        std::cout << "OverlapGraph is being deleted.\n";
    }

    std::vector<read_id_t> vertex_to_read; // get read id from vertex id
    boost::dynamic_bitset<> vertex_orientations; // 1 if forward, 0 if reverse
    boost::dynamic_bitset<> inclusions; // 1 if involved in inclusion
    std::vector< std::vector< Edge > > inclusion_edges;

    // OverlapGraph.cpp: graph functions
	node_id_t addVertex(read_id_t read_ID);
    void addEdge(Edge edge);
    Edge removeEdge(node_id_t v, node_id_t w);
    Edge removeEdgeWithOri(node_id_t v, node_id_t w, bool opposite_orientations);
    double checkEdge(node_id_t v, node_id_t w, bool reverse_allowed=true);
    double checkEdgeWithOri(node_id_t v, node_id_t w, bool opposite_orientations);
    Edge* getEdgeInfo(node_id_t v, node_id_t w, bool reverse_allowed=true);
    Edge* getEdgeInfoWithOri(node_id_t v, node_id_t w,
            bool opposite_orientations, bool reverse_allowed=true);
    unsigned int getEdgeCount();
    unsigned int getBackEdgeCount();
    unsigned int getVertexCount();
    void writeGraphToFile();
    void writeDiGraphToFile();
    void write2FASTG();
    void write2GFA(std::string filename);
    void writeGraphToFile(std::vector< std::list< node_id_t >> &tmp_adj_out);
    void reportCycle(node_id_t u, node_id_t v, bool remove);
    void printAdjacencyLists();
    void getGraphStats();
    void checkDuplicateEdges();
    void addEquivalentEdges();
    bool getOrientation(node_id_t v);
    void sortEdges();
    void buildOriginalsDict();

    // GraphAlgos.cpp: algorithms
    std::vector< node_id_t > sortVerticesByIndegree();
    void labelVertices(std::list< Edge > & edges_to_be_moved, std::list< Edge > & edges_to_be_deleted, boost::dynamic_bitset<> & orientations, int rand_seed);
    void vertexLabellingHeuristic(unsigned int & conflict_count);
    void postprocessEdges();
//    void removeCycles();
    void cycleRemovalHeuristic(bool remove_edges);
    std::set< std::pair< node_id_t, node_id_t > > findCycles(int randomize);
    void dfs_helper(node_id_t parent, node_id_t node, boost::dynamic_bitset<> &marked, boost::dynamic_bitset<> &visited, std::vector< node_id_t >& path, std::set< std::pair< node_id_t, node_id_t > >& backedge_vec, int randomize);
    void findInduced3Cycles();
    std::vector< std::vector< node_id_t > > getEdgesForMerging();
    void findBranches();
    void findBranchfreeGraph(std::vector< std::list< node_id_t > > & cur_adj_in, std::vector< std::list< node_id_t > > & cur_adj_out, std::set< node_id_t > & remove_in, std::set< node_id_t > & remove_out);
    unsigned int findTransEdges(std::vector< std::list< node_id_t > > & cur_adj_in, std::vector< std::list< node_id_t > > & cur_adj_out, std::vector< std::list< node_id_t > > & new_adj_in, std::vector< std::list< node_id_t > > & new_adj_out, bool removeTrans);
    bool nonemptyIntersect(std::list< node_id_t > & list1, std::list< node_id_t > & list2);
    std::vector< std::list< node_id_t > > sortAdjLists(std::vector< std::list< node_id_t > > & input_lists);
    std::vector< std::list< node_id_t > > sortAdjOut(std::vector< std::list< Edge > > & input_lists);
    void removeBranches();
    void removeTransitiveEdges();
    void removeTips();
    void reduceDiploidBranching();
    void removeInclusions();
};

#endif /* OVERLAPGRAPH_H_ */
