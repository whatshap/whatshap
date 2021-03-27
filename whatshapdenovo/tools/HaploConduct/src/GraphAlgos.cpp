//============================================================================
// Name        : ViralQuasispecies.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Additional graph algorithms to extend OverlapGraph.cpp
//============================================================================

#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/timer.hpp>
#include <set>
#include <algorithm> // std::random_shuffle, std::count
#include <iterator> // std::next, std::back_inserter
#include <cstdlib> // std::rand, std::srand

#include "OverlapGraph.h"

void OverlapGraph::removeInclusions() {
    if (program_settings.verbose) {
        std::cout << "removeInclusions...\n";
    }
    // remove all in- and outgoing edges from nodes marked as inclusions
    std::set< std::pair< node_id_t, node_id_t > > edges_to_remove;
    for (node_id_t v=0; v < vertex_count; v++) {
        std::vector< Edge > edge_vec;
        if (inclusions[v] == 0) {
            continue;
        }
        for (auto edge : adj_out.at(v)) {
            node_id_t outneighbor = edge.get_vertex(2);
            edges_to_remove.insert( std::make_pair(v, outneighbor) );
            edge_vec.push_back(edge);
        }
        for (auto inneighbor : adj_in.at(v)) {
            edges_to_remove.insert( std::make_pair(inneighbor, v) );
            Edge* edge = getEdgeInfo(inneighbor, v, false);
            edge_vec.push_back(*edge);
        }
        inclusion_edges.push_back(edge_vec); // keep edges in separate vector for FNO1
    }
    for (auto node_pair : edges_to_remove) {
//        std::cout << "edge " << node_pair.first << " " << node_pair.second << std::endl;
        removeEdge(node_pair.first, node_pair.second);
    }
//    std::cout << "Inclusions removed = " << edges_to_remove.size() << std::endl;
}

void OverlapGraph::reduceDiploidBranching() {
    /* Reduce the number of branches in the overlap graph using the fact that
       the assembly should be diploid. */
    std::vector< std::list< node_id_t > > unique_out_extensions;
    int min_diploid_overlap = 30;
    // find all uniquely out-extending edges
    for (auto adj_list : adj_out) {
        // if out-degree = 1 -> keep edge
        if (adj_list.size() == 1
                && adj_list.front().get_len(0) >= min_diploid_overlap
                && adj_list.front().get_mismatch_rate() < 0.000001) {
            Edge edge = adj_list.front();
            node_id_t outneighbor = edge.get_vertex(2);
            unique_out_extensions.push_back(std::list< node_id_t >{outneighbor});
        }
        else {
            unique_out_extensions.push_back(std::list< node_id_t >{});
        }
    }
    std::vector< std::list< node_id_t > > unique_in_extensions;
    // find all uniquely in-extending edges
    for (auto adj_list : adj_in) {
        // if in-degree = 1 -> keep edge
        if (adj_list.size() == 1) {
            unique_in_extensions.push_back(adj_list);
        }
        else {
            unique_in_extensions.push_back(std::list< node_id_t >{});
        }
    }
    std::set< std::pair< node_id_t, node_id_t > > edges_to_be_deleted;
    // find all non-unique edges to delete
    for (auto u_out : unique_out_extensions) {
        if (!u_out.empty()) {
            assert (u_out.size() == 1);
            node_id_t outneighbor = u_out.front();
            for (auto inneighbor : adj_in.at(outneighbor)) {
                if (unique_out_extensions.at(inneighbor).empty() || unique_out_extensions.at(inneighbor).front() != outneighbor) {
                    edges_to_be_deleted.insert(std::make_pair(inneighbor, outneighbor));
                }
            }
        }
    }
    for (auto u_in : unique_in_extensions) {
        if (!u_in.empty()) {
            assert (u_in.size() == 1);
            node_id_t inneighbor = u_in.front();
            for (auto edge : adj_out.at(inneighbor)) {
                node_id_t outneighbor = edge.get_vertex(2);
                if (unique_in_extensions.at(outneighbor).empty() || unique_in_extensions.at(outneighbor).front() != inneighbor) {
                    edges_to_be_deleted.insert(std::make_pair(inneighbor, outneighbor));
                }
            }
        }
    }
    // remove edges
    for (auto node_pair : edges_to_be_deleted) {
        removeEdge(node_pair.first, node_pair.second);
    }
    std::cout << "Number of edges removed because of diploid argument: " << edges_to_be_deleted.size() << std::endl;
}

std::vector< std::vector< node_id_t > > OverlapGraph::getEdgesForMerging() {
    /* builds a list of node pairs corresponding to edges that are selected
       for merging (i.e. super-read construction) */
    boost::dynamic_bitset<> bitvec(vertex_count); // keeps track of vertices that have been processed
    std::vector< std::vector< node_id_t > > node_vec; // stores node pairs to be output
    node_id_t node = 0;
    for (auto adj_list : adj_out) {
        if (!bitvec[node] && adj_list.size() > 0) {
            std::vector< std::pair<node_id_t, int> > pairs;
            node_id_t outneighbor;
            for (auto edge_it : adj_list) {
                outneighbor = edge_it.get_vertex(2);
                assert (node == edge_it.get_vertex(1));
                int perc = edge_it.get_perc();
                pairs.push_back(std::make_pair(outneighbor, perc));
            }
//            // sort outgoing edges by decreasing overlap percentages
//            std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, int>& a, const std::pair<node_id_t, int>& b)
//            {
//                return a.second > b.second;
//            }
//            );
            for (auto pair_it : pairs) {
                outneighbor = pair_it.first;
                if (!bitvec[outneighbor]) {
                    std::vector< node_id_t > edge_nodes = {node, outneighbor};
                    node_vec.push_back(edge_nodes);
                    bitvec[node] = 1;
                    bitvec[outneighbor] = 1;
                    break;
                }
            }
        }
        node++;
    }
    return node_vec;
}

std::vector< node_id_t > OverlapGraph::sortVerticesByIndegree() {
    std::vector< std::pair<node_id_t, unsigned int> > pairs;
    node_id_t v = 0;
    for (auto it : adj_in) {
        pairs.push_back(std::make_pair(v, it.size()));
        v++;
    }
    std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, unsigned int>& a, const std::pair<node_id_t, unsigned int>& b)
    {
        if (a.second == b.second) {
            return a.first < b.first;
        }
        else {
            return a.second < b.second;
        }
    }
    );
    assert (pairs.size() == vertex_count);
    std::vector< node_id_t > sorted_vertices;
    for (auto pair_it : pairs) {
        node_id_t node = pair_it.first;
        assert (node >= 0 && node < vertex_count);
        sorted_vertices.push_back(node);
    }
    assert (sorted_vertices.size() == vertex_count);
    return sorted_vertices;
}

void OverlapGraph::vertexLabellingHeuristic(unsigned int & conflict_count) {
    if (program_settings.verbose) {
        std::cout << "Applying vertex labelling heuristic...\n";
    }
    std::list< Edge > min_edges_to_be_moved;
    std::list< Edge > min_edges_to_be_deleted;
    boost::dynamic_bitset<> opt_orientations;
    opt_orientations.resize(vertex_count, true); // all-ones vector: initially all labels forward
    if (program_settings.add_duplicates) { // second half of vertices are duplicates, hence reversed
        assert (!program_settings.resolve_orientations);
        unsigned int readcount = fastq_storage->get_readcount();
        for (node_id_t i=readcount; i<vertex_count; i++) { // change labels for duplicates
            opt_orientations[i] = 0;
        }
    }
    else if (program_settings.resolve_orientations) {
        if (program_settings.verbose) {
            std::cout << "resolving vertex orientations by BFS\n";
        }
        labelVertices(min_edges_to_be_moved, min_edges_to_be_deleted, opt_orientations, 1);
        // try k possible labellings and choose the best one
        std::list< Edge > edges_to_be_moved;
        std::list< Edge > edges_to_be_deleted;
        unsigned int delete_count = min_edges_to_be_deleted.size();
        boost::dynamic_bitset<> orientations;
        int count = 1;
        while (count < 100 && delete_count > 0) { // k = 100 tries
            count++;
            labelVertices(edges_to_be_moved, edges_to_be_deleted, orientations, count);
            if (edges_to_be_deleted.size() < delete_count) {
                min_edges_to_be_deleted = edges_to_be_deleted;
                min_edges_to_be_moved = edges_to_be_moved;
                delete_count = min_edges_to_be_deleted.size();
                opt_orientations = orientations;
            }
            edges_to_be_moved.clear();
            edges_to_be_deleted.clear();
        }
        // move edges
        if (program_settings.verbose) {
            std::cout << "moving edges where necessary..\n";
        }
        for (auto edge_it : min_edges_to_be_moved) {
            node_id_t u = edge_it.get_vertex(1);
            node_id_t v = edge_it.get_vertex(2);
            bool opposite_ori = (edge_it.get_ori(1) == edge_it.get_ori(2));
            removeEdgeWithOri(v, u, opposite_ori);
            addEdge(edge_it);
        }
        if (delete_count > 0) {
            // remove edges that were conflicting
            if (program_settings.verbose) {
                std::cout << "deleting " << delete_count << " conflicting edges..\n";
            }
            conflict_count = delete_count;
            for (auto edge_it : min_edges_to_be_deleted) {
                node_id_t u = edge_it.get_vertex(1);
                node_id_t v = edge_it.get_vertex(2);
                bool opposite_ori = (edge_it.get_ori(1) == edge_it.get_ori(2));
                removeEdgeWithOri(u, v, opposite_ori);
            }
        }
        else {
            if (program_settings.verbose) {
                std::cout << "A perfect labelling was found after " << count << " tries, no conflicting edges." << std::endl;
            }
            conflict_count = 0;
        }
    }
    vertex_orientations = opt_orientations;
}

void OverlapGraph::labelVertices(std::list< Edge > & edges_to_be_moved, std::list< Edge > & edges_to_be_deleted, boost::dynamic_bitset<> & orientations, int rand_seed) {
    // label vertices by BFS -> O(V+E)
    orientations.resize(vertex_count, true); // all-ones vector: initially all labels forward
    boost::dynamic_bitset<> visited(vertex_count);
    std::list<node_id_t> bfs;
    // sort vertices by increasing indegree
    std::vector< node_id_t > sorted_vertices = sortVerticesByIndegree();
    // process vertices in this order
    for (node_id_t i = 0; i < vertex_count; i++) {
        node_id_t start_node = sorted_vertices.at(i);
        if (!visited[start_node]) {
            bfs.push_back(start_node);
            visited[start_node] = 1; // orientation = 1 (i.e. forward) by default
        }
        while (!bfs.empty()) {
            node_id_t node = bfs.front();
            bfs.pop_front();
            // copy all neighbors (both in and out) to vector and order randomly
            std::vector< node_id_t > adj_vec(adj_in.at(node).begin(), adj_in.at(node).end());
            for (auto edge_it : adj_out.at(node)) {
                node_id_t neighbor = edge_it.get_vertex(2);
                adj_vec.push_back(neighbor);
            }
            std::srand( unsigned( rand_seed ) ); // set the 'random' seed (fixed per iteration)
            std::random_shuffle( adj_vec.begin(), adj_vec.end() );
            // recursively check all neighbors
            for (auto neighbor_it : adj_vec) {
                if (!visited[neighbor_it]) {
                    bfs.push_back(neighbor_it);
                    visited[neighbor_it] = 1;
                    Edge* edge = getEdgeInfo(node, neighbor_it);
                    if (edge->get_ori(1) == edge->get_ori(2)) {
                        orientations[neighbor_it] = orientations[node];
                    }
                    else {
                        orientations[neighbor_it] = !orientations[node];
                    }
                }
            }
        }
    }
    for (node_id_t i=0; i < vertex_count; i++) {
        assert (visited[i]);
    }
    // check all edges -> O(E)
    std::vector< std::list<Edge> >::iterator it1;
    std::list<Edge>::iterator it2;
    int switch_count = 0;
    int overall_count = 0;
    int move_count = 0;
    node_id_t i = 0;
    assert (edges_to_be_deleted.empty());
    assert (edges_to_be_moved.empty());
    for (it1 = adj_out.begin(); it1 != adj_out.end(); it1++) {
        std::list< Edge >::iterator it2 = it1->begin();
        unsigned int size = it1->size();
        for (unsigned int j=0; j < size; j++) {
            overall_count++;
            node_id_t u = it2->get_vertex(1);
            assert (u == i);
            node_id_t v = it2->get_vertex(2);
            bool ori1 = it2->get_ori(1);
            bool ori2 = it2->get_ori(2);
            bool type_v1 = orientations[u];
            bool type_v2 = orientations[v];
            if (ori1 == type_v1 && ori2 == type_v2) { // edge ok
                it2++;
            }
            else if ((ori1 == ori2 && type_v1 != type_v2) || (ori1 != ori2 && type_v1 == type_v2)) { // contradiction
//                    std::cout << "unsolvable vertex orientation due to conflicting edges... removing.\n";
                edges_to_be_deleted.push_back(*it2);
                it2++;
            }
            else { // edge agrees with labelling but the read orientations need to be flipped
                Edge edge = *it2;
                bool move = edge.switch_edge_orientation();
                assert (it2->get_vertex(1) == u);
                switch_count++;
                if (move) { // edge orientation was changed, so move edge to other adjacency list
                    move_count++;
                    edges_to_be_moved.push_back(edge);
//                        edges_to_be_inserted.push_back(*it2); // store edge that has to be moved
                    assert (edge.get_ori(1) == type_v2 && edge.get_ori(2) == type_v1);
                    assert (edge.get_vertex(1) == v);
                    it2++;
//                        it2 = adj_out.at(u).erase(it2);
//                        // also update adj_in because this is important for cycle check
//                        adj_in.at(v).remove(u);
//                        adj_in.at(u).push_back(v);
                }
                else {
                    it2->switch_edge_orientation();
                    assert (it2->get_ori(1) == type_v1 && it2->get_ori(2) == type_v2);
                    it2++;
                }
            }
        }
        i++;
    }
}


void OverlapGraph::dfs_helper(node_id_t parent, node_id_t node, boost::dynamic_bitset<> &marked, boost::dynamic_bitset<> &visited, std::vector<node_id_t>& path, std::set< std::pair< node_id_t, node_id_t > >& backedge_vec, int randomize) {
    if (program_settings.verbose) {
//        std::cout << "in dfs_helper" << std::endl;
    }
    if (marked[node]) {
        // node was seen before in this dfs-path, so there must be a cycle in G
        backedge_vec.insert(std::make_pair(parent, node));
        unsigned int len_path = 0;
        bool cycle_identified = false;
        for (auto it = path.rbegin(); it != path.rend(); it++) {
            len_path++;
            if (*it == node) {
                cycle_identified = true;
                break;
            }
        }
        assert (len_path > 0 && len_path <= path.size());
        assert (cycle_identified == true);
    }
    else if (!visited[node]) {
        marked[node] = 1;
        path.push_back(node);
        if (path.size() > graph_depth) {
            graph_depth = path.size(); // NOTE: this only the length of the longest DFS path that was traversed, the graph depth could be larger
        }
        // sort neighboring nodes
        std::vector< node_id_t > sorted_neighbors;
        if (randomize == 1) {
            // sort neighbors by increasing value of pos1
            std::vector< std::pair<node_id_t, int> > pairs;
            for (auto it : adj_out.at(node)) {
                pairs.push_back(std::make_pair(it.get_vertex(2), it.get_pos(1)));
            }
            std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, int>& a, const std::pair<node_id_t, int>& b)
            {
                if (a.second == b.second) {
                    return a.first < b.first;
                }
                else {
                    return a.second < b.second;
                }
            }
            );
            for (auto pair_it : pairs) {
                node_id_t v = pair_it.first;
                assert (v >= 0 && v < vertex_count);
                sorted_neighbors.push_back(v);
            }
        }
        else if (randomize == 2) {
            // sort neighbors by decreasing overlap score
            std::vector< std::pair<node_id_t, double> > pairs;
            for (auto it : adj_out.at(node)) {
                pairs.push_back(std::make_pair(it.get_vertex(2), it.get_score()));
            }
            std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, double>& a, const std::pair<node_id_t, double>& b)
            {
                if (a.second == b.second) {
                    return a.first < b.first;
                }
                else {
                    return a.second > b.second;
                }
            }
            );
            for (auto pair_it : pairs) {
                node_id_t v = pair_it.first;
                assert (v >= 0 && v < vertex_count);
                sorted_neighbors.push_back(v);
            }
        }
        else if (randomize == 3) {
            // sort neighbors by decreasing overlap length
            std::vector< std::pair<node_id_t, int> > pairs;
            for (auto it : adj_out.at(node)) {
                pairs.push_back(std::make_pair(it.get_vertex(2), it.get_len(0)));
            }
            std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, int>& a, const std::pair<node_id_t, int>& b)
            {
                if (a.second == b.second) {
                    return a.first < b.first;
                }
                else {
                    return a.second > b.second;
                }
            }
            );
            for (auto pair_it : pairs) {
                node_id_t v = pair_it.first;
                assert (v >= 0 && v < vertex_count);
                sorted_neighbors.push_back(v);
            }
        }
        else if (randomize == 4) {
            // sort neighbors by increasing mismatch rate
            std::vector< std::pair<node_id_t, double> > pairs;
            for (auto it : adj_out.at(node)) {
                pairs.push_back(std::make_pair(it.get_vertex(2), it.get_mismatch_rate()));
            }
            std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, double>& a, const std::pair<node_id_t, double>& b)
            {
                if (a.second == b.second) {
                    return a.first < b.first;
                }
                else {
                    return a.second < b.second;
                }
            }
            );
            for (auto pair_it : pairs) {
                node_id_t v = pair_it.first;
                assert (v >= 0 && v < vertex_count);
                sorted_neighbors.push_back(v);
            }
        }
        else {
            // sort neighbors randomly
            for (auto it : adj_out.at(node)) {
                sorted_neighbors.push_back(it.get_vertex(2));
            }
            std::srand( unsigned( randomize ) ); // set the random seed
            std::random_shuffle( sorted_neighbors.begin(), sorted_neighbors.end() );
        }
        // now iterate: process neighbors according to sorted list
        for (auto it_v : sorted_neighbors) {
            dfs_helper(node, it_v, marked, visited, path, backedge_vec, randomize);
        }
        // move up one node to find the next dfs path
        marked[node] = 0;
        visited[node] = 1;
        assert (path.back() == node);
        path.pop_back();
    }
}

std::set< std::pair< node_id_t, node_id_t > > OverlapGraph::findCycles(int randomize) {
    if (program_settings.verbose) {
        std::cout << "findCycles.." << std::endl;
    }
    std::string filename = PATH + "cycles.txt";
    remove(filename.c_str());
    // sort vertices by increasing order of indegree and process nodes in this order
    std::vector< node_id_t > sorted_vertices = sortVerticesByIndegree();
    // DFS to find cycles -> O(V+E)
    boost::dynamic_bitset<> visited(vertex_count);
    boost::dynamic_bitset<> marked(vertex_count);
    std::set< std::pair< node_id_t, node_id_t > > backedge_vec;
    std::vector<node_id_t> path;
    for (auto i : sorted_vertices) {
        if (!visited[i]) {
            dfs_helper(vertex_count, i, marked, visited, path, backedge_vec, randomize);
        }
    }
    return backedge_vec;
}

void OverlapGraph::cycleRemovalHeuristic(bool remove_edges) {
    if (program_settings.verbose) {
        std::cout << "in cycleRemovalHeuristic..." << std::endl;
    }
    std::set< std::pair< node_id_t, node_id_t > > cur_backedge_vec;
    std::set< std::pair< node_id_t, node_id_t > > opt_backedge_vec;
    opt_backedge_vec = findCycles(/*randomize*/ 1);
//    std::cout << "try 1: " << opt_backedge_vec.size() << std::endl;
    int count = 1;
    while (count < 20 && opt_backedge_vec.size() > 0) {
        // try k-1 (k=20) other dfs trees randomly to see if we can improve
        count++;
        cur_backedge_vec = findCycles(/*randomize*/ count);
//        std::cout << "try " << count << ": " << cur_backedge_vec.size() << std::endl;
        if (cur_backedge_vec.size() < opt_backedge_vec.size()) {
            opt_backedge_vec = cur_backedge_vec;
        }
    }
    backedge_count = opt_backedge_vec.size();
    if (program_settings.verbose) {
        if (opt_backedge_vec.empty()) {
            std::cout << "Overlap graph is cycle-free :)\n";
        }
        else {
            std::cout << "\nTHE GRAPH CONTAINS CYCLES!!! Number of back-edges: " << backedge_count << "\n";
        }
    }
    for (auto node_pair : opt_backedge_vec) {
        reportCycle(node_pair.first, node_pair.second, remove_edges);
    }
    if (remove_edges && program_settings.verbose) {
        std::cout << "New edge count " << edge_count << std::endl;
    }
}

void OverlapGraph::removeTips() {
    if (program_settings.verbose) {
        std::cout << "removeTips..." << std::endl;
    }
    unsigned int tip_count = 0;
    unsigned int max_tip_len = program_settings.max_tip_len;
    std::set< std::pair<node_id_t, node_id_t> > edges_to_remove;
    // find all outgoing tips
    for (node_id_t i = 0; i < vertex_count; i++) {
        bool cont = true;
        std::list< Edge > adj_list = adj_out.at(i);
        if (adj_list.size() <= 1) { // definitely no tips
            continue;
        }
        bool alltips = true;
        std::vector< std::pair<node_id_t, node_id_t> > local_tips;
        std::vector< Read* > local_tip_reads;
        for (auto edge1 = adj_list.begin(); edge1 != adj_list.end() && cont; edge1++) {
            node_id_t v1 = edge1->get_vertex(2);
            // check if i->v1 is a dead end ('tip')
            if (adj_out.at(v1).empty()) {
                if (edge1->ext_len(1) == 0) { // inclusion edge -> always a tip
                    tip_count += 1;
                    edges_to_remove.insert( std::make_pair(i, v1)); // make sure read is removed
                    edge1->get_read(2)->set_tip(); // mark read as tip
                }
                else if (edge1->ext_len(1) < max_tip_len) {
                    tip_count += 1;
                    local_tips.push_back( std::make_pair(i, v1) );
                    local_tip_reads.push_back( edge1->get_read(2) );
                }
            }
            else {
                alltips = false;
            }
        }
        if (!alltips) {
            edges_to_remove.insert(local_tips.begin(), local_tips.end());
            // mark sequences as tips
            for (auto read_ptr : local_tip_reads) {
                read_ptr->set_tip();
            }
        }
    }
    if (program_settings.verbose) {
        std::cout << "Number of out-tip edges: " << tip_count << std::endl;
    }
    // find all incoming tips
    for (node_id_t i = 0; i < vertex_count; i++) {
        bool cont = true;
        std::list< node_id_t > adj_list = adj_in.at(i);
        if (adj_list.size() <= 1) { // definitely no tips
            continue;
        }
        bool alltips = true;
        std::vector< std::pair<node_id_t, node_id_t> > local_tips;
        std::vector< Read* > local_tip_reads;
        for (auto v1 = adj_list.begin(); v1 != adj_list.end() && cont; v1++) {
            // check if i->v1 is a dead end ('tip')
            if (adj_in.at(*v1).empty()) {
                Edge* edge = getEdgeInfo(*v1, i, false);
                if (edge->ext_len(0) == 0) { // inclusion edge -> always a tip
                    tip_count += 1;
                    edges_to_remove.insert( std::make_pair(*v1, i)); // make sure read is removed
                    edge->get_read(1)->set_tip(); // mark read as tip
                }
                else if (edge->ext_len(0) < max_tip_len) {
                    tip_count += 1;
                    local_tips.push_back( std::make_pair(*v1, i) );
                    local_tip_reads.push_back( edge->get_read(1) );
                }
            }
            else {
                alltips = false;
            }
        }
        if (!alltips) {
            edges_to_remove.insert(local_tips.begin(), local_tips.end());
            // mark sequences as tips
            for (auto read_ptr : local_tip_reads) {
                read_ptr->set_tip();
            }
        }
    }
    if (program_settings.verbose) {
        std::cout << "Final number of tip edges: " << tip_count << std::endl;
    }
    for (auto node_pair : edges_to_remove) {
        // std::cout << node_pair.first << " " << node_pair.second << std::endl;
        // std::cout << "size adj_out: " << adj_out.at(node_pair.first).size() << " " << adj_out.at(node_pair.second).size() << std::endl;
        // std::cout << "size adj_in: " << adj_in.at(node_pair.first).size() << " " << adj_in.at(node_pair.second).size() << std::endl;
        Edge edge = removeEdge(node_pair.first, node_pair.second);
        branching_edges.push_back(edge);
    }
}

void OverlapGraph::findBranches() {
    if (program_settings.verbose) {
        std::cout << "findBranches..." << std::endl;
    }
    if (program_settings.min_overlap_perc > 0) {
        std::cout << "NOTE: removing branches while min_overlap_perc > 0" << std::endl;
//        return;
    }
    unsigned int removal_count_in = 0;
    unsigned int removal_count_out = 0;
    unsigned int tip_count = 0;
    unsigned int node_count_in = 0;
    unsigned int node_count_out = 0;
    // find all outgoing branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        bool cont = true;
        std::list< Edge > adj_list = adj_out.at(i);
        for (auto edge1 = adj_list.begin(); edge1 != adj_list.end() && cont; ++edge1) {
            node_id_t v1 = edge1->get_vertex(2);
            // check if i->v1 is a dead end ('tip')
            if (adj_out.at(v1).empty()) {
                tip_count += 1;
                continue;
            }
            for (auto edge2 = edge1; ++edge2 != adj_list.end() && cont; /**/) {
                node_id_t v2 = edge2->get_vertex(2);
                // check if i->v2 is a dead end ('tip')
                if (adj_out.at(v2).empty()) {
                    tip_count += 1;
                }
                else if (checkEdge(v1, v2) == -1) {
                    removal_count_out += adj_list.size();
                    node_count_out++;
                    cont = false;
                }
            }
        }
    }
    if (program_settings.verbose) {
        std::cout << "Number of out-tip edges: " << tip_count << std::endl;
    }
    // find all incoming branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        bool cont = true;
        std::list< node_id_t > adj_list = adj_in.at(i);
        for (auto v1 = adj_list.begin(); v1 != adj_list.end() && cont; ++v1) {
            // check if i->v1 is a dead end ('tip')
            if (adj_in.at(*v1).empty()) {
                tip_count += 1;
                continue;
            }
            for (auto v2 = v1; ++v2 != adj_list.end() && cont; /**/) {
                // check if i->v2 is a dead end ('tip')
                if (adj_in.at(*v2).empty()) {
                    tip_count += 1;
                }
                else if (checkEdge(*v1, *v2) == -1) {
                    removal_count_in += adj_list.size();
                    node_count_in++;
                    cont = false;
                }
            }
        }
    }
    if (program_settings.verbose) {
        std::cout << "Final number of tip edges: " << tip_count << std::endl;
        std::cout << "Number of edges to be removed in order to eliminate all out-branches: " << removal_count_out << std::endl;
        std::cout << "Number of edges to be removed in order to eliminate all in-branches: " << removal_count_in << std::endl;
        std::cout << "Number of nodes out-disconnected: " << node_count_out << std::endl;
        std::cout << "Number of nodes in-disconnected: " << node_count_in << std::endl;
        // unsigned int remaining_out = edge_count - removal_count;
        // std::cout << "Number of edges remaining: " << remaining << std::endl;
    }
}

void OverlapGraph::findBranchfreeGraph(std::vector< std::list< node_id_t > > & cur_adj_in, std::vector< std::list< node_id_t > > & cur_adj_out, std::set< node_id_t > & remove_in, std::set< node_id_t > & remove_out) {
    // NOTE: this function assumes that transitive edges have been removed first!!
    assert (program_settings.remove_trans == 1);
    if (program_settings.verbose) {
        std::cout << "findBranchfreeGraph..." << std::endl;
    }
    if (program_settings.min_overlap_perc > 0) {
        std::cout << "NOTE: removing branches while min_overlap_perc > 0" << std::endl;
    }
    // find all outgoing branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        std::list< node_id_t > adj_list = cur_adj_out.at(i);
        if (adj_list.size() > 1) {
            remove_out.insert(i);
//            remove_in.insert(adj_list.begin(), adj_list.end());
        }
    }
    // find all incoming branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        std::list< node_id_t > adj_list = cur_adj_in.at(i);
        if (adj_list.size() > 1) {
            remove_in.insert(i);
//            remove_out.insert(adj_list.begin(), adj_list.end());
        }
    }
    if (program_settings.verbose) {
        std::cout << "Number of nodes out-disconnected: " << remove_out.size() << std::endl;
        std::cout << "Number of nodes in-disconnected: " << remove_in.size() << std::endl;
    }
}


unsigned int OverlapGraph::findTransEdges(std::vector< std::list< node_id_t > > & cur_adj_in, std::vector< std::list< node_id_t > > & cur_adj_out, std::vector< std::list< node_id_t > > & new_adj_in, std::vector< std::list< node_id_t > > & new_adj_out, bool removeTrans) {
//    std::cout << "Find transitive edges..." << std::endl;
    // assumes adjacency lists are sorted
    unsigned int total_edges_kept = 0;
    unsigned int total_edges_removed = 0;
    node_id_t node1 = 0;
    node_id_t node2;
    std::list< node_id_t > list1;
    std::list< node_id_t > list2;
    for (auto adj_per_node : cur_adj_out) {
        list1 = adj_per_node;
        for (auto node_it : adj_per_node) {
            node2 = node_it;
            list2 = cur_adj_in.at(node2);
            bool transitive = nonemptyIntersect(list1, list2);
            if (transitive != removeTrans) {
                new_adj_out.at(node1).push_back(node2);
                new_adj_in.at(node2).push_back(node1);
                total_edges_kept++;
            }
            else {
                total_edges_removed++;
            }
        }
        node1++;
    }
    if (program_settings.verbose) {
        std::cout << total_edges_kept << " edges kept" << std::endl;
    }
    return total_edges_kept;
}


bool OverlapGraph::nonemptyIntersect(std::list< node_id_t > & list1, std::list< node_id_t > & list2) {
    // assumes list1 and list2 are sorted
    auto it1 = list1.begin();
    auto it2 = list2.begin();
    while (it1 != list1.end() && it2 != list2.end()) {
        if (*it1 == *it2) { // nonempty intersection
            return true;
        }
        else if (*it1 < *it2) {
            it1++;
        }
        else {
            it2++;
        }
    }
    return false;
}

std::vector< std::list< node_id_t > > OverlapGraph::sortAdjLists(std::vector< std::list< node_id_t > > & input_lists) {
    std::vector< std::list< node_id_t > > output_lists;
    for (auto list_it : input_lists) {
        list_it.sort();
        output_lists.push_back(list_it);
    }
    return output_lists;
}

std::vector< std::list< node_id_t > > OverlapGraph::sortAdjOut(std::vector< std::list< Edge > > & input_lists) {
    // sort adj_out by increasing outneighbor ID
    std::vector< std::list< node_id_t > > sorted_nodes_out;
    std::vector< std::list< Edge > > sorted_edges_out;
    for (auto adj_list : adj_out) {
        std::list< node_id_t > neighbors;
        std::list< Edge > edges;
        std::vector< std::pair< node_id_t, Edge > > pairs;
        node_id_t outneighbor;
        for (auto edge_it : adj_list) {
            outneighbor = edge_it.get_vertex(2);
            pairs.push_back(std::make_pair(outneighbor, edge_it));
        }
        std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, Edge>& a, const std::pair<node_id_t, Edge>& b)
        {
            return a.first < b.first;
        }
        );
        for (auto pair_it : pairs) {
            neighbors.push_back(pair_it.first);
            edges.push_back(pair_it.second);
        }
        sorted_nodes_out.push_back(neighbors);
        sorted_edges_out.push_back(edges);
    }
    adj_out = sorted_edges_out;
    return sorted_nodes_out;
}

void OverlapGraph::removeBranches() {
//    std::cout << "removeBranches..." << std::endl;
    // create sorted adjacency lists
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    sorted_adj_in = sortAdjLists(adj_in);
    sorted_adj_out = sortAdjOut(adj_out);
    // obtain graph without transitive edges
    std::vector< std::list< node_id_t > > new_adj_in;
    std::vector< std::list< node_id_t > > new_adj_out;
    new_adj_in = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    new_adj_out = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    findTransEdges(sorted_adj_in, sorted_adj_out, new_adj_in, new_adj_out, true);
    // DO WE NEED THIS???

    // find all branches
    std::set< node_id_t > remove_in;
    std::set< node_id_t > remove_out;
    findBranchfreeGraph(new_adj_in, new_adj_out, remove_in, remove_out);
    // remove branches from new graph
    for (auto node2 : remove_in) {
        // note: new_adj_out is not updated accordingly!!
        new_adj_in.at(node2).clear();
    }
    for (auto node1 : remove_out) {
        // note: new_adj_in is not updated accordingly!!
        new_adj_out.at(node1).clear();
    }
//    std::cout << "branches removed" << std::endl;
    // assign vertices to connected components in branch-free graph
    std::map< node_id_t, unsigned int > component_map;
    boost::dynamic_bitset<> visited(vertex_count);
    unsigned int current_component = 0;
    for (node_id_t i = 0; i < vertex_count; i++) {
        if (!visited[i]) {
            std::list< node_id_t > stack;
            stack.push_back(i);
            visited[i] = 1;
            while (!stack.empty()) {
                node_id_t node = stack.front();
                stack.pop_front();
                component_map.insert(std::make_pair(node, current_component));
                for (auto out_nb : new_adj_out.at(node)) {
                    int count_this_edge = std::count(
                        new_adj_in.at(out_nb).begin(),
                        new_adj_in.at(out_nb).end(),
                        node
                    );
                    if (count_this_edge == 0) {
                        // this edge was previously removed from adj_in but not
                        // yet from adj_out so we need to ignore it
                        continue;
                    }
                    if (!visited[out_nb]) {
                        stack.push_back(out_nb);
                        visited[out_nb] = 1;
                    }
                }
                for (auto in_nb : new_adj_in.at(node)) {
                    int count_this_edge = std::count(
                        new_adj_out.at(in_nb).begin(),
                        new_adj_out.at(in_nb).end(),
                        node
                    );
                    if (count_this_edge == 0) {
                        // this edge was previously removed from adj_out but not
                        // yet from adj_in so we need to ignore it
                        continue;
                    }
                    if (!visited[in_nb]) {
                        stack.push_back(in_nb);
                        visited[in_nb] = 1;
                    }
                }
            }
            current_component++;
        }
    }
    if (program_settings.verbose) {
        std::cout << "Total number of components " << current_component << std::endl;
    }
    // 2. remove all edges of overlap graph between different components of branch-free graph
    std::list< std::pair< node_id_t, node_id_t > > edges_to_remove;
    for (node_id_t i = 0; i < vertex_count; i++) {
        for (auto edge : adj_out.at(i)) {
            assert (edge.get_vertex(1) == i);
            node_id_t j = edge.get_vertex(2);
            if (component_map.at(i) != component_map.at(j)) {
                edges_to_remove.push_back(std::make_pair(i, j));
            }
        }
    }
    for (auto node_pair : edges_to_remove) {
        Edge edge = removeEdge(node_pair.first, node_pair.second);
//        std::cout << edge.get_vertex(1) << " " << edge.get_vertex(2) << std::endl;
        branching_edges.push_back(edge);
    }
    if (program_settings.verbose) {
        std::cout << "Number of edges removed: " << edges_to_remove.size() << std::endl;
        std::cout << "Number of edges remaining: " << edge_count << std::endl;
    }
}

void OverlapGraph::removeTransitiveEdges() {
    if (program_settings.remove_trans == 0) {
        return;
    }
    if (program_settings.verbose) {
        std::cout << "removeTransitiveEdges..." << std::endl;
    }
    // create sorted adjacency lists
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    sorted_adj_in = sortAdjLists(adj_in);
    sorted_adj_out = sortAdjOut(adj_out);
    // obtain graph of all transitive edges
    std::vector< std::list< node_id_t > > new_adj_in;
    std::vector< std::list< node_id_t > > new_adj_out;
    new_adj_in = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    new_adj_out = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    // new_adj_in, new_adj_out will contain all edges that are found to be transitive
    unsigned int transitive_count = findTransEdges(sorted_adj_in, sorted_adj_out, new_adj_in, new_adj_out, false);
    // now iterate on the graph of transitive edges we just found
    for (unsigned int it=1; it<program_settings.remove_trans; it++) {
        sorted_adj_in = sortAdjLists(new_adj_in);
        sorted_adj_out = sortAdjLists(new_adj_out);
        for (unsigned int j = 0; j < vertex_count; j++) {
            new_adj_out.at(j).clear();
            new_adj_in.at(j).clear();
        }
        transitive_count = findTransEdges(sorted_adj_in, sorted_adj_out, new_adj_in, new_adj_out, false);
    }
    // if doing branch reduction, remove branches lacking evidence based on 3-cliques
    std::set< node_pair_t > edges_to_be_deleted;
    if (program_settings.remove_trans == 1
                && program_settings.branch_reduction > 0) {
        // for every transitive edge, check branches
        node_id_t node1 = 0;
        for (auto neighbors : new_adj_out) {
            for (auto node2 : neighbors) {
                // node1->node2 is transitive
                int ovlen = getEdgeInfo(node1, node2, false)->get_len(0);
                // check out-branches
                for (auto edge_out : adj_out.at(node1)) {
                    if (edge_out.get_len(0) <= ovlen) {
                        node_id_t v = edge_out.get_vertex(2);
                        edges_to_be_deleted.insert( std::make_pair(node1, v) );
                    }
                }
                // check in-branches
                for (auto v_in : adj_in.at(node2)) {
                    if (getEdgeInfo(v_in, node2, false)->get_len(0) <= ovlen) {
                        edges_to_be_deleted.insert( std::make_pair(v_in, node2) );
                    }
                }
            }
            node1++;
        }
    }
    // finally we need to remove all single/double/triple transitive edges that we found
    if (1.0*transitive_count > 0.5*edge_count) {
        // better create new adjacency lists than remove edges from current lists
        std::vector< std::list< Edge > > final_adj_out;
        unsigned int remaining_edge_count = 0;
        unsigned int current_edge_count = 0;
        unsigned int del_count = 0;
        for (auto adj_list : adj_out) {
            std::list< Edge > final_adj_list;
            if (adj_list.size() > 0) {
                node_id_t node1 = adj_list.begin()->get_vertex(1);
                std::list< node_id_t > trans_adj_list = new_adj_out.at(node1);
                std::list< Edge >::const_iterator edge_it;
                std::list< node_id_t >::const_iterator trans_it = trans_adj_list.begin();
                for (edge_it = adj_list.begin(); edge_it != adj_list.end(); ) {
                    node_id_t node2 = edge_it->get_vertex(2);
                    Edge edge = *edge_it;
                    if (trans_it == trans_adj_list.end()) {
                        edge_it++;
                    }
                    else if (node2 < *trans_it) {
                        edge_it++;
                    }
                    else if (node2 == *trans_it) { // transitive edge
                        trans_it++;
                        edge_it++;
                        continue; // don't add edge
                    }
                    else {
                        trans_it++;
                    }
                    // check if edge is scheduled for deletion
                    if (edges_to_be_deleted.find(std::make_pair(node1, node2))
                            != edges_to_be_deleted.end()) {
                        // ignore edge
                        del_count++;
                    }
                    else {
                        // keep edge
                        final_adj_list.push_back(edge);
                        remaining_edge_count++;
                        current_edge_count++;
                    }
                }
            }
            final_adj_out.push_back(final_adj_list);
        }
        adj_out = final_adj_out;
        // also build new adj_in
        std::vector< std::list< node_id_t > > final_adj_in (vertex_count, std::list< node_id_t >());
        node_id_t node1 = 0;
        for (auto adj_list : adj_out) {
            for (auto edge_it : adj_list) {
                node_id_t node2 = edge_it.get_vertex(2);
                final_adj_in.at(node2).push_back(node1);
            }
            node1++;
        }
        adj_in = final_adj_in;
        // update edge count
        assert (transitive_count <= edge_count);
        edge_count = edge_count - transitive_count - del_count;
        if (program_settings.verbose) {
            std::cout << "transitive edge count: " << transitive_count << std::endl;
            std::cout << "deleted branching edges: " << del_count << std::endl;
            std::cout << "remaining edge count: " << remaining_edge_count << std::endl;
        }
        assert (edge_count == remaining_edge_count);
    }
    else {
        // remove transitive edges from current overlap graph
        node_id_t node1 = 0;
        for (auto neighbors : new_adj_out) {
            for (auto node2 : neighbors) {
                removeEdge(node1, node2);
            }
            node1++;
        }
        for (auto edge : edges_to_be_deleted) {
            if (checkEdge(edge.first, edge.second, false) >= 0) {
                removeEdge(edge.first, edge.second);
            }
        }
    }

    // findTransEdges(sorted_adj_in, sorted_adj_out, new_adj_in, new_adj_out, false);
    // if (program_settings.remove_trans == 1) { // remove transitive edges
    //     node_id_t node1 = 0;
    //     for (auto neighbors : new_adj_out) {
    //         for (auto node2 : neighbors) {
    //             removeEdge(node1, node2);
    //         }
    //         node1++;
    //     }
    // }
    // else {
    //     // find all double transitive edges
    //     std::vector< std::list< node_id_t > > new_adj_in2;
    //     std::vector< std::list< node_id_t > > new_adj_out2;
    //     new_adj_in2 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    //     new_adj_out2 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    //     findTransEdges(new_adj_in, new_adj_out, new_adj_in2, new_adj_out2, false);
    //     if (program_settings.remove_trans == 2) { // remove double transitive edges
    //         node_id_t node1 = 0;
    //         for (auto neighbors : new_adj_out2) {
    //             for (auto node2 : neighbors) {
    //                 removeEdge(node1, node2);
    //             }
    //             node1++;
    //         }
    //     }
    //     else {
    //         // find all triple transitive edges
    //         std::vector< std::list< node_id_t > > new_adj_in3;
    //         std::vector< std::list< node_id_t > > new_adj_out3;
    //         new_adj_in3 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    //         new_adj_out3 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    //         findTransEdges(new_adj_in2, new_adj_out2, new_adj_in3, new_adj_out3, false);
    //         // remove triple transitive edges
    //         node_id_t node1 = 0;
    //         for (auto neighbors : new_adj_out3) {
    //             for (auto node2 : neighbors) {
    //                 removeEdge(node1, node2);
    //             }
    //             node1++;
    //         }
    //     }
    // }
    // std::cout << "Number of edges remaining: " << edge_count << std::endl;
}


///**************************************************
//                 Test functions
//**************************************************/

//void OverlapGraph::testCycles() {
//    std::vector< std::list< unsigned int> > test_adj_list;
//    for (unsigned int i = 0; i < 5; i++) {
//        std::list<unsigned int> list;
//        list.push_back((i + 1) % 5);
//        test_adj_list.push_back(list);
//    }
//    removeCycles(test_adj_list, 5);
//}

//void OverlapGraph::removeCycles(std::vector< std::list< node_id_t >> &tmp_adj_out, unsigned int V) {
//    std::cout << "removeCycles.. testing\n";
//    // DFS to find cycles -> O(V+E)
//    boost::dynamic_bitset<> visited(V);
//    boost::dynamic_bitset<> marked(V);
//    for (unsigned int i = 0; i < V; i++) {
//        if (!visited[i]) {
//            dfs_helper(V, i, marked, visited, tmp_adj_out);
//        }
//    }
//    std::cout << "Number of back-edges: " << backedge_count << "\n";
//}


//void OverlapGraph::dfs_helper(node_id_t parent, node_id_t node, boost::dynamic_bitset<> &marked, boost::dynamic_bitset<> &visited, std::vector< std::list< node_id_t >> &tmp_adj_out) {
////    std::cout << "dfs_helper..\n";
//    if (marked[node]) {
//        reportCycle(parent, node, true);
//        backedge_count++;
//    }
//    else if (!visited[node]) {
//        marked[node] = 1;
//        std::list< node_id_t >::const_iterator it;
//        for (it = tmp_adj_out[node].begin(); it != tmp_adj_out[node].end(); it++) {
//            node_id_t neighbor = *it;
//            dfs_helper(node, neighbor, marked, visited, tmp_adj_out);
//        }
//        visited[node] = 1;
//        marked[node] = 0;
//    }
//}
