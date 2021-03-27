//============================================================================
// Name        : BranchReduction.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.1
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Use read evidence to reduce branches in the overlap graph
//============================================================================

#include <math.h> /* floor */
#include <fstream>

#include "BranchReduction.h"

/*
IMPORTANT:
When using initial paired-end input reads as single-end sequences, the input
must be ordered singles-paired1-paired2 with corresponding integer read IDs from
0 to #SE+2*#PE. The number of SE and PE input reads must be specified when
invoking readBasedBranchReduction, as well as the minimum evidence required for
keeping a branching edge, and the initial read file(s).
Otherwise, make sure PE_count = 0.

APPROACH:
Identify all branches in the graph and for every branch u->(v_0,...,v_k):
1. Compute a list of all FIRST k difference positions between any pair of
    branch sequences --> diff_list
2. For all i, find all common subreads between u and v_i, where PE read IDs from
    the same fragment are considered identical (i.e. modulo #PE). For all common
    subreads, check if the corresponding sequence is identical to the contig
    sequence at all positions of diff_list. If so, add to the "evidence set" of
    edge i.
3. Find all branching components by following in-branches connected to
    out-branches and vice versa.
4. For every component, take the intersection of evidence for edges which appear
    twice (once an in-branch, once an out-branch). Then compute the unique
    evidence set per edge by comparing evidence from all edges in the component.
5. Remove all branching edges that do not have sufficient unique evidence.
*/

void BranchReduction::readBasedBranchReduction() {
    if (program_settings.verbose) {
        std::cout << "readBasedBranchReduction" << std::endl;
    }
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    sorted_adj_in = overlap_graph->sortAdjLists(overlap_graph->adj_in);
    sorted_adj_out = overlap_graph->sortAdjOut(overlap_graph->adj_out);
    std::set< node_id_t > branch_in;
    std::set< node_id_t > branch_out;
    // find all branches in the graph and process one by one
    overlap_graph->findBranchfreeGraph(
        sorted_adj_in, sorted_adj_out, branch_in, branch_out
    );
    std::list< Edge > missing_edges;
    std::vector< std::pair< std::list< node_id_t >, int > > final_branch_in (
        overlap_graph->getVertexCount(), std::make_pair(std::list< node_id_t >(), 0)
    );
    std::pair< std::list< node_id_t >, int > branch;
    for (auto node : branch_in) {
        bool outbranch = false;
        branch = findBranchingEvidence(
            node, sorted_adj_in.at(node), missing_edges, outbranch
        );
        if (!branch.first.empty()) {
            final_branch_in.at(node) = branch;
        }
    }
    std::vector< std::pair< std::list< node_id_t >, int > > final_branch_out (
        overlap_graph->getVertexCount(), std::make_pair(std::list< node_id_t >(), 0)
    );
    for (auto node : branch_out) {
        bool outbranch = true;
        branch = findBranchingEvidence(
            node, sorted_adj_out.at(node), missing_edges, outbranch
        );
        if (!branch.first.empty()) {
            final_branch_out.at(node) = branch;
        }
    }
    // add missing edges to vector of removed branching edges, to ensure that
    // they will be found in the next iteration
    for (auto edge : missing_edges) {
        overlap_graph->branching_edges.push_back(edge);
    }
    // find branching components; components with false branches (i.e. missing edges)
    // are skipped automatically and all edges removed
    std::list< node_pair_t > edges_to_remove;
    findBranchingComponents(final_branch_in, final_branch_out, edges_to_remove);

    // keep track of neighboring components
    std::vector< std::set< unsigned int > > neighboring_components;
    if (program_settings.careful) {
        // build a map of nodes to components
        std::map< node_id_t, std::set< unsigned int > > nodes_to_components;
        for (node_id_t node=0; node < overlap_graph->getVertexCount(); node++) {
            nodes_to_components.insert(
                std::make_pair(node, std::set< unsigned int >{})
            );
        }
        unsigned int idx = 0;
        for (auto component_info : branching_components) {
            for (auto node_pair : component_info.first) {
                nodes_to_components.at(node_pair.first).emplace(idx);
                nodes_to_components.at(node_pair.second).emplace(idx);
            }
            idx++;
        }
        // now lookup nodes and keep track of neighboring components
        for (auto component_info : branching_components) {
            std::set< unsigned int > neighbors;
            for (auto node_pair : component_info.first) {
                std::set< unsigned int > s1, s2;
                s1 = nodes_to_components.at(node_pair.first);
                neighbors.insert(s1.begin(), s1.end());
                s2 = nodes_to_components.at(node_pair.second);
                neighbors.insert(s2.begin(), s2.end());
            }
            neighboring_components.push_back(neighbors);
        }
    }
    else {
        for (auto component : branching_components) {
            std::set< unsigned int > neighbors;
            neighboring_components.push_back(neighbors);
        }
    }

    // keep track of components for which edges were kept
    std::set< unsigned int > components_kept;

    // read evidence thresholds from file
    // std::cout << "Reading evidence threshold table" << std::endl;
    std::map< int, int > evidence_threshold_map;
    std::string filename = "evidence_threshold_table.tsv";
    std::ifstream thresholds_file(filename.c_str());
    if (thresholds_file.is_open()) {
        std::string line;
        std::stringstream ss;
        while (getline(thresholds_file, line)) {
            if (line.empty() || line.front() == '#') {
                continue;
            }
            ss << line;
            std::string tmp;
            getline(ss, tmp, '\t');
            int dist = std::stoi(tmp);
            getline(ss, tmp, '\t');
            getline(ss, tmp, '\t');
            int min_ev = std::stoi(tmp);
            evidence_threshold_map[dist] = min_ev;
            ss.clear();
        }
        thresholds_file.close();
    }
    else {
        std::cerr << "Unable to open 'evidence_threshold_table.tsv'.\n";
        exit(1);
    }

    // select UNIQUE evidence and add edges with insufficient evidence to edges_to_remove
    unsigned int idx = 0;
    for (auto component_info : branching_components) {
        std::vector< node_pair_t > component = component_info.first;
        int dist = component_info.second;
        std::set< unsigned int > neighbors = neighboring_components.at(idx);
        bool skip = false;
        for (auto comp_idx : neighbors) {
            if (comp_idx == idx) {
                continue;
            }
            else if (components_kept.find(comp_idx) != components_kept.end()) {
                edges_to_remove.insert(
                    edges_to_remove.end(), component.begin(), component.end()
                );
                skip = true;
                continue;
            }
        }
        if (!skip) {
            int min_evidence;
            auto min_ev_ptr = evidence_threshold_map.find(dist);
            if (min_ev_ptr != evidence_threshold_map.end()) {
                min_evidence = min_ev_ptr->second;
                if (program_settings.verbose) {
                    std::cout << "min_evidence = " << min_evidence << std::endl;
                }
                bool keep = countUniqueEvidence(component, min_evidence, edges_to_remove);
                if (keep) {
                    components_kept.insert(idx);
                }
            }
            else {
                // distance too large, remove this component
                if (program_settings.verbose) {
                    std::cout << "distance too large..." << dist << std::endl;
                }
                edges_to_remove.insert(
                    edges_to_remove.end(), component.begin(), component.end()
                );
            }
        }
        idx++;
    }
    // remove duplicates from edges_to_remove
    edges_to_remove.sort();
    edges_to_remove.unique();

//     if (false) {
// //    if (program_settings.diploid && edges_to_remove.size() == overlap_graph->getEdgeCount()) {
//         // no edges remaining, try with diploid branch reduction instead
//         overlap_graph->removeTips();
//         overlap_graph->reduceDiploidBranching();
//         overlap_graph->removeBranches();
//     }
    // now remove all selected edges from overlap graph
    for (auto node_pair : edges_to_remove) {
//        std::cout << "removing " << node_pair.first << " " << node_pair.second << std::endl;
        if (overlap_graph->checkEdge(node_pair.first, node_pair.second, false) < 0) {
            std::cout << "edge not found, reverse: " << overlap_graph->checkEdge(node_pair.second, node_pair.first, false) << std::endl;
            exit(1);
        }
        Edge edge = overlap_graph->removeEdge(node_pair.first, node_pair.second);
        overlap_graph->branching_edges.push_back(edge);
    }
//        std::cout << edges_to_remove.size() << " edges removed" << std::endl;
}

std::pair< std::list< node_id_t >, int > BranchReduction::findBranchingEvidence(
    node_id_t node1, std::list< node_id_t > neighbors,
    std::list< Edge > & missing_edges, bool outbranch
) {
    // std::cout << "findBranchingEvidence" << std::endl;
    assert (neighbors.size() > 1);
    std::list< node_id_t > final_branch (neighbors.begin(), neighbors.end());
    final_branch.push_front(node1);
    // build list of difference positions (first difference for every pair)
    std::vector< std::string > sequence_vec;
    std::vector< int > startpos_vec;
    std::pair< std::list< int >, int> diff_info;
    std::vector< node_pair_t > missing_inclusion_edges;
    std::vector< node_id_t > neighbors_vec(neighbors.begin(), neighbors.end());
    if (outbranch) {
        diff_info = buildDiffListOut(node1, neighbors_vec, sequence_vec, startpos_vec, missing_inclusion_edges, missing_edges);
    }
    else {
        diff_info = buildDiffListIn(node1, neighbors_vec, sequence_vec, startpos_vec, missing_edges);
    }
    std::list< int > diff_list = diff_info.first;
    int distance = diff_info.second;
    // build evidence list of subread IDs per neighbor
    std::unordered_map< read_id_t, OriginalIndex > subreads1 = overlap_graph->original_ID_dict.at(node1);
    std::unordered_map< node_id_t, std::list< read_id_t > > evidence_per_neighbor;
    std::vector< std::string >::const_iterator seq_it = sequence_vec.begin();
    std::vector< int >::const_iterator startpos_it = startpos_vec.begin();
    for (auto node2 : neighbors) {
        // find common subreads from originals dict
        std::list< read_id_t > evidence_list;
        std::string contig = *seq_it;
        int startpos = *startpos_it;
        seq_it++;
        startpos_it++;
        std::unordered_map< read_id_t, OriginalIndex > subreads2 = overlap_graph->original_ID_dict.at(node2);
        for (auto subread : subreads2) {
            read_id_t subread_id = subread.first;
            std::unordered_map< read_id_t, OriginalIndex >::const_iterator common_subread;
            common_subread = subreads1.find(subread_id);
            std::unordered_map< read_id_t, OriginalIndex >::const_iterator common_PE;
            read_id_t subread_id_PE;
            if (subread_id >= SE_count + PE_count) {
                // try its /1 mate
                common_PE = subreads1.find(subread_id - PE_count);
                subread_id_PE = subread_id - PE_count;
            }
            else if (subread_id >= SE_count) {
                // try its /2 mate
                common_PE = subreads1.find(subread_id + PE_count);
                subread_id_PE = subread_id + PE_count;
            }
            else {
                // single-end subread so no mate to try
                common_PE = subreads1.end();
                subread_id_PE = 0; // dummy
            }
            // first check the subread itself
            bool result1 = false;
            if (common_subread != subreads1.end()) {
                Read* original_read = original_fastq->get_read(subread_id);
                int index = subread.second.index1;
                std::string sequence;
                if (subread.second.forward) {
                    sequence = original_read->get_seq(0);
                }
                else {
                    sequence = build_rev_comp(original_read->get_seq(0));
                }
                result1 = checkReadEvidence(contig, startpos, sequence, index, diff_list);
            }
            if (result1) { // single-end evidence found
                assert (subread_id < 2*program_settings.original_readcount);
                evidence_list.push_back(subread_id);
            }
            // now also try its mate
            bool result2 = false;
            if (common_PE != subreads1.end()) {
                Read* original_read = original_fastq->get_read(subread_id);
                std::string sequence;
                int index = subread.second.index1;
                if (subread.second.forward) {
                    sequence = original_read->get_seq(0);
                }
                else {
                    sequence = build_rev_comp(original_read->get_seq(0));
                }
                result2 = checkReadEvidence(contig, startpos, sequence, index, diff_list);
            }
            if (result2) { // paired-end evidence found
                read_id_t joint_id = program_settings.original_readcount + std::min(subread_id, subread_id_PE);
                assert (joint_id < 2*program_settings.original_readcount);
                evidence_list.push_back(joint_id);
            }
        }
        evidence_list.sort();
        evidence_list.unique();
        evidence_per_neighbor.insert(std::make_pair(node2, evidence_list));
    }
    // take care of 'missing edges' due to inclusions
    for (auto node_pair : missing_inclusion_edges) {
        evidence_per_neighbor.at(node_pair.first).clear();
        if (neighbors.size() == 2) { // no actual branch
            final_branch.clear();
        }
        else { // real branch, only remove inclusion node
            final_branch.remove(node_pair.first);
        }
        // // identify newly created branches and remove corresponding edges;
        // // these branches will be resolved at a next iteration of ViralQuasispecies
        // node_id_t outnode = node_pair.first;
        // node_id_t innode = node_pair.second;
        // for (auto out_edge : overlap_graph->adj_out.at(outnode)) {
        //     edges_to_remove.push_back(std::make_pair(outnode, out_edge.get_vertex(2)));
        // }
        // for (auto in_edge : overlap_graph->adj_in.at(innode)) {
        //     edges_to_remove.push_back(std::make_pair(in_edge, innode));
        // }
    }
    // store evidence in evidence_per_edge vector (private class member)
    std::list< node_id_t >::const_iterator branch_it = final_branch.begin();
    branch_it++; // skip first entry, this contains the branching node
    for (auto neighbor : neighbors) {
        if (branch_it != final_branch.end() && neighbor == *branch_it) {
            // neighbor also in final_branch so we store its evidence
            std::unordered_map< safe_edge_count_t, std::list< read_id_t > >::iterator ev_it;
            std::list< read_id_t > current_evidence = evidence_per_neighbor.at(neighbor);
            // std::cout << "current evidence for " << node1 << " " << neighbor << std::endl;
            // for (auto id : current_evidence) {
            //     std::cout << id << " ";
            // }
            // std::cout << std::endl;
            safe_edge_count_t index;
            if (outbranch) {
                index = edgeToEvidenceIndex(node1, neighbor);
            }
            else {
                index = edgeToEvidenceIndex(neighbor, node1);
            }
            ev_it = evidence_per_edge.find(index);
            if (ev_it != evidence_per_edge.end()) {
                // intersect evidence sets
                std::list< read_id_t > existing_evidence = ev_it->second;
                std::list< read_id_t >::iterator it2 = existing_evidence.begin();
                std::list< read_id_t > intersection;
                while (it2 != existing_evidence.end()) {
                    auto find_ev = std::find(current_evidence.begin(), current_evidence.end(), *it2);
                    if (find_ev != current_evidence.end()) {
                        intersection.push_back(*it2);
                    }
                    // else {
                    //     it2 = ev_it->second.erase(it2);
                    // }
                    it2++;
                }
                evidence_per_edge.at(index) = intersection;
            }
            else {
                // insert new evidence set
                evidence_per_edge.insert(std::make_pair(index, current_evidence));
            }
            branch_it++;
        }
    }
    assert (branch_it == final_branch.end());
    return std::make_pair(final_branch, distance);
}

std::pair< std::list< int >, int > BranchReduction::buildDiffListOut(
    node_id_t node1,
    std::vector< node_id_t > neighbors,
    std::vector< std::string > & sequence_vec,
    std::vector< int > & startpos_vec,
    std::vector< node_pair_t > & missing_inclusion_edges,
    std::list< Edge > & missing_edges
) {
    // build a list of all FIRST difference positions between any pair of branch sequences
    // std::cout << "buildDiffListOut" << std::endl;
    std::list< int > diff_list;
    std::vector< Edge* > edge_vec;
    for (auto node : neighbors) {
        Edge* edge = overlap_graph->getEdgeInfo(node1, node, /*reverse_allowed*/ false);
        int pos = edge->get_pos(1);
        Read* read = edge->get_read(2);
        assert (!read->is_paired());
        std::string sequence;
        if (overlap_graph->getOrientation(node1)) {
            sequence = read->get_seq(0);
        }
        else {
            sequence = read->get_rev_comp(0);
        }
        sequence_vec.push_back(sequence);
        startpos_vec.push_back(pos);
        edge_vec.push_back(edge);
    }
    // do all pairwise sequence comparisons until a difference is found;
    /* note: could be done more efficiently by evaluating all sequences in one
        for loop, but this would take extensive bookkeeping. Since the number
        of neighbors at the same branch is not expected to be big, all pairwise
        comparisons should be fine computationally. */
    node_id_t node_i, node_j;
    std::string seq_i, seq_j;
    std::string subseq_i, subseq_j;
    int pos_i, pos_j;
    int relative_pos;
    std::vector< int > diff_pos;
    int len;
    int startpos;
    std::vector< int > distance_vec;
    for (unsigned int i=0; i < neighbors.size(); i++) {
        node_i = neighbors.at(i);
        seq_i = sequence_vec.at(i);
        pos_i = startpos_vec.at(i);
        for (unsigned int j=i+1; j < neighbors.size(); j++) {
            node_j = neighbors.at(j);
            seq_j = sequence_vec.at(j);
            pos_j = startpos_vec.at(j);
            assert (node_i != node_j);
            if (pos_i < pos_j) {
                relative_pos = pos_j - pos_i;
                if (relative_pos > int(seq_i.size() - program_settings.min_overlap_len)) {
                    missing_inclusion_edges.push_back(std::make_pair(node_i, node_j));
                    continue;
                }
                len = std::min(seq_i.size()-relative_pos, seq_j.size());
                subseq_i = seq_i.substr(relative_pos, len);
                subseq_j = seq_j.substr(0, len);
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_j;
            }
            else {
                relative_pos = pos_i - pos_j;
                if (relative_pos > int(seq_j.size() - program_settings.min_overlap_len)) {
                    missing_inclusion_edges.push_back(std::make_pair(node_j, node_i));
                    continue;
                }
                len = std::min(seq_j.size()-relative_pos, seq_i.size());
                subseq_i = seq_i.substr(0, len);
                subseq_j = seq_j.substr(relative_pos, len);
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_i;
            }
            assert (len > 0);
            for (auto pos : diff_pos) {
                diff_list.push_back(pos + startpos);
            }
            if (diff_pos.empty()) {
                // identical overlap -> add corresponding edge; NOTE: this should hardly occur!
                double score = program_settings.edge_threshold; // overlap score
                int pos1 = relative_pos;
                int pos2 = 0; // no PE-overlaps allowed
                bool ori1; // 1 if NORMAL, 0 if REVERSE (vertex 1)
                bool ori2; // similar for vertex 2
                Read* read1; // pointer to read corresponding to vertex1
                Read* read2; // pointer to read corresponding to vertex2
                std::string ord = "-";
                node_id_t vertex1; // out-vertex
                node_id_t vertex2; // in-vertex
                int overlap_perc = (int)floor(100 * len / std::min(seq_i.size(), seq_j.size())) ;

                if (pos_i < pos_j || (pos_i == pos_j && node_i < node_j)) {
                    // node_i first
                    ori1 = edge_vec.at(i)->get_ori(2);
                    ori2 = edge_vec.at(j)->get_ori(2);
                    read1 = edge_vec.at(i)->get_read(2);
                    read2 = edge_vec.at(j)->get_read(2);
                    vertex1 = node_i;
                    vertex2 = node_j;
                }
                else {
                    // node_j first
                    ori1 = edge_vec.at(j)->get_ori(2);
                    ori2 = edge_vec.at(i)->get_ori(2);
                    read1 = edge_vec.at(j)->get_read(2);
                    read2 = edge_vec.at(i)->get_read(2);
                    vertex1 = node_j;
                    vertex2 = node_i;
                }
                Edge new_edge(score, pos1, pos2, ori1, ori2, ord, read1, read2);
                new_edge.set_vertices(vertex1, vertex2);
                new_edge.set_perc(overlap_perc);
                new_edge.set_len(len, 0);
                missing_edges.push_back(new_edge);
                // mark node as false out branch
                false_out_branches.insert(node1);
            }
            else if (i==0) {
                distance_vec.push_back(diff_pos.at(0)+startpos);
            }
        }
    }
    // calculate variation distance
    int dist;
    if (!distance_vec.empty()) {
        dist = 0.5 * (
            *std::min_element(distance_vec.begin(), distance_vec.end())
            + *std::max_element(distance_vec.begin(), distance_vec.end())
        );
    }
    else {
        dist = 0;
    }
    // remove duplicate entries
    diff_list.sort();
    diff_list.unique();
    return std::make_pair(diff_list, dist);
}

std::pair< std::list< int >, int > BranchReduction::buildDiffListIn(
    node_id_t node1,
    std::vector< node_id_t > neighbors,
    std::vector< std::string > & sequence_vec,
    std::vector< int > & startpos_vec,
    std::list< Edge > & missing_edges
) {
    // build a list of all FIRST difference positions between any pair of branch sequences
//    std::cout << "buildDiffListIn" << std::endl;
    std::list< int > diff_list;
    std::vector< int > pos_vec;
    std::vector< Edge* > edge_vec;
    int node1_len = 0;
    for (auto node : neighbors) {
        Edge* edge = overlap_graph->getEdgeInfo(node, node1, /*reverse_allowed*/ false);
        int pos = edge->get_pos(1);
        Read* read = edge->get_read(1);
        assert (!read->is_paired());
        std::string sequence;
        if (overlap_graph->getOrientation(node1)) {
            sequence = read->get_seq(0);
        }
        else {
            sequence = read->get_rev_comp(0);
        }
        sequence_vec.push_back(sequence);
        pos_vec.push_back(pos);
        edge_vec.push_back(edge);
        if (node1_len == 0) {
            node1_len = edge->get_read(2)->get_len();
        }
    }
    /* since we are dealing wint an in-branch, we need to infer the startpos_vec entries
        once we have those positions, we can proceed similar to the out-branch case,
        except for the diff_pos computation where we need to reverse the sequences
        first */
    int max_pos = *std::max_element(pos_vec.begin(), pos_vec.end());
    for (auto pos : pos_vec) {
        startpos_vec.push_back(max_pos - pos);
    }
    // do all pairwise sequence comparisons until a difference is found;
    /* note: could be done more efficiently by evaluating all sequences in one
        for loop, but this would take extensive bookkeeping. Since the number
        of neighbors at the same branch is not expected to be big, all pairwise
        comparisons should be fine computationally. */
    node_id_t node_i, node_j;
    std::string seq_i, seq_j;
    std::string subseq_i, subseq_j;
    int pos_i, pos_j;
    int relative_pos;
    std::vector< int > diff_pos;
    int len;
    int startpos;
    int overlap_len;
    std::vector< int > distance_vec;
    for (unsigned int i=0; i < neighbors.size(); i++) {
        for (unsigned int j=i+1; j < neighbors.size(); j++) {
            node_i = neighbors.at(i);
            node_j = neighbors.at(j);
            seq_i = sequence_vec.at(i);
            seq_j = sequence_vec.at(j);
            pos_i = startpos_vec.at(i);
            pos_j = startpos_vec.at(j);
            overlap_len = std::min(
                seq_i.size()-pos_vec.at(i), seq_j.size()-pos_vec.at(j)
            );
            if (pos_i < pos_j) {
                relative_pos = pos_j - pos_i;
                assert (relative_pos <= int(seq_i.size() - program_settings.min_overlap_len)); // in-branches can't be the result of inclusions
                len = std::min(seq_i.size()-relative_pos, seq_j.size());
                subseq_i = seq_i.substr(relative_pos, len);
                subseq_j = seq_j.substr(0, len);
                std::reverse( subseq_i.begin(), subseq_i.end() );
                std::reverse( subseq_j.begin(), subseq_j.end() );
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_j;
            }
            else {
                relative_pos = pos_i - pos_j;
                assert (relative_pos <= int(seq_j.size() - program_settings.min_overlap_len)); // in-branches can't be the result of inclusions
                len = std::min(seq_j.size()-relative_pos, seq_i.size());
                subseq_i = seq_i.substr(0, len);
                subseq_j = seq_j.substr(relative_pos, len);
                std::reverse( subseq_i.begin(), subseq_i.end() );
                std::reverse( subseq_j.begin(), subseq_j.end() );
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_i;
            }
            assert (len > 0);
            for (auto pos : diff_pos) {
                diff_list.push_back(len - pos + startpos);
            }
            if (diff_pos.empty()) {
                // identical overlap -> add corresponding edge; NOTE: this should hardly occur!
                double score = program_settings.edge_threshold; // overlap score
                int pos1 = relative_pos;
                int pos2 = 0; // no PE-overlaps allowed
                bool ori1; // 1 if NORMAL, 0 if REVERSE (vertex 1)
                bool ori2; // similar for vertex 2
                Read* read1; // pointer to read corresponding to vertex1
                Read* read2; // pointer to read corresponding to vertex2
                std::string ord = "-";
                node_id_t vertex1; // out-vertex
                node_id_t vertex2; // in-vertex
                int overlap_perc = (int)floor(100 * len / std::min(seq_i.size(), seq_j.size()));

                if (pos_i < pos_j || (pos_i == pos_j && node_i < node_j)) {
                    // node_i first
                    ori1 = edge_vec.at(i)->get_ori(1);
                    ori2 = edge_vec.at(j)->get_ori(1);
                    read1 = edge_vec.at(i)->get_read(1);
                    read2 = edge_vec.at(j)->get_read(1);
                    vertex1 = node_i;
                    vertex2 = node_j;
                }
                else {
                    // node_j first
                    ori1 = edge_vec.at(j)->get_ori(1);
                    ori2 = edge_vec.at(i)->get_ori(1);
                    read1 = edge_vec.at(j)->get_read(1);
                    read2 = edge_vec.at(i)->get_read(1);
                    vertex1 = node_j;
                    vertex2 = node_i;
                }
                Edge new_edge(score, pos1, pos2, ori1, ori2, ord, read1, read2);
                new_edge.set_vertices(vertex1, vertex2);
                new_edge.set_perc(overlap_perc);
                new_edge.set_len(len, 0);
                missing_edges.push_back(new_edge);
                // also mark branch node as false
                false_in_branches.insert(node1);
            }
            else if (i==0) {
                distance_vec.push_back(diff_pos.at(0) + node1_len - overlap_len);
            }
        }
    }
    // calculate variation distance
    int dist;
    if (!distance_vec.empty()) {
        dist = 0.5 * (
            *std::min_element(distance_vec.begin(), distance_vec.end())
            + *std::max_element(distance_vec.begin(), distance_vec.end())
        );
    }
    else {
        dist = 0;
    }
    // remove duplicate entries
    diff_list.sort();
    diff_list.unique();
    return std::make_pair(diff_list, dist);
}



std::vector< int > BranchReduction::findDiffPos(std::string seq1, std::string seq2) {
    // given two input strings, find the first position where they disagree
    assert (seq1.size() == seq2.size());
    std::string::const_iterator it1 = seq1.begin();
    std::string::const_iterator it2 = seq2.begin();
    std::vector< int > diff_pos;
    int pos = 0;
    while (it1 != seq1.end()) {
        if (*it1 != *it2) {
            diff_pos.push_back(pos);
            if (diff_pos.size() == 100) {
                // store at most 100 positions per sequence pair
                return diff_pos;
            }
        }
        it1++;
        it2++;
        pos++;
    }
    return diff_pos;
}


bool BranchReduction::checkReadEvidence(std::string contig, int startpos, std::string read, int index, std::list< int > diff_list) {
    // check if subread agrees with contig on diff_list positions
    bool true_evidence = false;
    int read_start = startpos + index;
    int read_end = read_start + read.size();
    int contig_start = startpos;
    int contig_end = startpos + contig.size();
    for (auto diff_pos : diff_list) {
        if (diff_pos < read_start || diff_pos >= read_end) {
            // read does not overlap this diff_pos
            continue;
        }
        else if (diff_pos < contig_start || diff_pos >= contig_end) {
            // contig does not overlap this diff_pos
            continue;
        }
        // check if read agrees with contig base
        if (read.at(diff_pos - read_start) != contig.at(diff_pos - contig_start)) {
            // disagreement --> false evidence
            true_evidence = false;
            break; // no need to continue checking
        }
        else {
            true_evidence = true; // at least one diff_pos covered
        }
    }
    return true_evidence;
}

void BranchReduction::findBranchingComponents(
    std::vector< std::pair< std::list< node_id_t >, int > > final_branch_in,
    std::vector< std::pair< std::list< node_id_t >, int > > final_branch_out,
    std::list< node_pair_t > & edges_to_remove
) {
    // find all branching components by following in-branches connected to out-branches
    // and vice versa
    // std::cout << "findBranchingComponents" << std::endl;
    std::unordered_map< node_id_t, bool > visited_in_branches;
    std::unordered_map< node_id_t, bool > visited_out_branches;
    std::unordered_map< node_id_t, std::list< node_id_t > > branch_in_map;
    std::unordered_map< node_id_t, std::list< node_id_t > > branch_out_map;
    std::unordered_map< node_id_t, int > branch_in_dist_map;
    std::unordered_map< node_id_t, int > branch_out_dist_map;
    // store all branches in a dict for easy access
    for (auto branch_info : final_branch_in) {
        std::list< node_id_t > branch = branch_info.first;
        if (branch.empty()) {
            continue;
        }
        node_id_t node = branch.front();
        branch.pop_front();
        visited_in_branches.insert(std::make_pair(node, false));
        branch_in_map.insert(std::make_pair(node, branch));
        branch_in_dist_map.insert(std::make_pair(node, branch_info.second));
    }
    for (auto branch_info : final_branch_out) {
        std::list< node_id_t > branch = branch_info.first;
        if (branch.empty()) {
            continue;
        }
        node_id_t node = branch.front();
        branch.pop_front();
        visited_out_branches.insert(std::make_pair(node, false));
        branch_out_map.insert(std::make_pair(node, branch));
        branch_out_dist_map.insert(std::make_pair(node, branch_info.second));
    }
    // now build components
    std::vector< std::vector< node_pair_t > > size2_components;
    for (auto branch : branch_in_map) {
        node_id_t node = branch.first;
        if (visited_in_branches.at(node)) {
            continue;
        }
        std::list< node_id_t > neighbors = branch.second;
        std::vector< node_pair_t > new_component;
        bool has_false_branch;
        if (false_in_branches.find(node) != false_in_branches.end()) {
            has_false_branch = true;
        }
        else {
            has_false_branch = false;
        }
        for (auto in_neighbor : neighbors) {
            assert (node != in_neighbor);
            node_pair_t node_pair = std::make_pair(in_neighbor, node);
            new_component.push_back(node_pair);
        }
        visited_in_branches.at(node) = true;

        // compute variation distance
        int dist1 = branch_in_dist_map.at(node);
        std::pair<int, node_id_t> dist_node_pair;
        dist_node_pair = extendComponentOut(new_component, neighbors,
            has_false_branch, visited_in_branches, visited_out_branches,
            branch_in_map, branch_out_map, branch_out_dist_map
        );
        int dist2 = dist_node_pair.first;
        node_id_t outnode = dist_node_pair.second;
        Edge* edge = overlap_graph->getEdgeInfo(outnode, node, /*reverse_allowed*/ false);
        int len1 = edge->get_read(1)->get_len();
        int len2 = edge->get_read(2)->get_len();
        int overlap_len = edge->get_len(0);
        if (overlap_len < 100) {
            if (dist1 < len2 - overlap_len + 100) {
                dist1 = len2 - overlap_len + 100;
            }
            if (dist2 < len1 - overlap_len + 100) {
                dist2 = len1 - overlap_len + 100;
            }
            // dist1 = std::max(dist1, len2 - overlap_len + 100);
            // dist2 = std::max(dist2, len1 - overlap_len + 100);
        }
        else {
            if (dist1 < len2) { dist1 = len2; }
            if (dist2 < len1) { dist2 = len1; }
        }
        int dist = dist1 + dist2 - len1 - len2 + edge->get_len(0);
        assert (dist >= 100);
        // std::cout << "distance = " << dist << std::endl;

        // filter duplicate edges from component
        std::sort(new_component.begin(), new_component.end(),
            [](node_pair_t pair1, node_pair_t pair2) {
                if (pair1.first != pair2.first) {
                    return pair1.first < pair2.first;
                }
                else {
                    return pair1.second < pair2.second;
                }
            });
        auto end_it = std::unique(new_component.begin(), new_component.end());
        new_component.resize(std::distance(new_component.begin(), end_it));

        if (has_false_branch) {
            // component contains a false branch due to a missing edge;
            // remove entire component and add missing edge at next iteration
            for (auto node_pair : new_component) {
                edges_to_remove.push_back(node_pair);
            }
        }
        // else if (program_settings.careful && new_component.size() == 2) {
        //     size2_components.push_back(new_component);
        //     visited_in_branches.at(node) = false;
        // }
        else {
            branching_components.push_back(std::make_pair(new_component, dist));
        }
    }
    // for (auto comp : size2_components) {
    //     node_id_t node = comp.at(0).second;
    //     if (visited_out_branches.find(node) != visited_out_branches.end() &&
    //         visited_out_branches.at(node)) {
    //         // don't allow merging a node on both sides in the same iteration
    //         // TODO: only remove out-branch if in-branch is actually kept
    //         for (auto node_pair : comp) {
    //             edges_to_remove.push_back(node_pair);
    //         }
    //     }
    //     else {
    //         branching_components.push_back(comp);
    //         visited_in_branches.at(node) = true;
    //     }
    // }

    // process remaining out-branches; since we already did all in-branches, the
    // remaining components are trivial
    for (auto branch : branch_out_map) {
        node_id_t node = branch.first;
        // read_id_t read_id = overlap_graph->vertex_to_read.at(node);
        if (visited_out_branches.at(node)) {
            continue;
        }
        // else if (program_settings.careful &&
        //     visited_in_branches.find(node) != visited_in_branches.end() &&
        //     visited_in_branches.at(node)) {
        //     // don't allow merging a node on both sides in the same iteration
        //     // TODO: only remove out-branch if in-branch is actually kept
        //     for (auto innode : branch.second) {
        //         edges_to_remove.push_back(std::make_pair(branch.first, innode));
        //     }
        //     continue;
        // }
        std::list< node_id_t > neighbors = branch.second;
        std::vector< node_pair_t > new_component;
        for (auto out_neighbor : neighbors) {
            assert (node != out_neighbor);
            node_pair_t node_pair = std::make_pair(node, out_neighbor);
            new_component.push_back(node_pair);
        }
        // compute variation distance
        int dist1 = branch_out_dist_map.at(node);
        int dist2;
        node_id_t innode = neighbors.front();
        Edge* edge = overlap_graph->getEdgeInfo(node, innode, /*reverse_allowed*/ false);
        int len1 = edge->get_read(1)->get_len();
        int len2 = edge->get_read(2)->get_len();
        int overlap_len = edge->get_len(0);
        if (overlap_len < 100) {
            if (dist1 < len1 - overlap_len + 100) {
                dist1 = len1 - overlap_len + 100;
            }
            dist2 = len2 - overlap_len + 100;
            // dist1 = std::max(dist1, len2 - overlap_len + 100);
            // dist2 = std::max(dist2, len1 - overlap_len + 100);
        }
        else {
            if (dist1 < len1) { dist1 = len1; }
            dist2 = len2;
        }
        int dist = dist1 + dist2 - len1 - len2 + edge->get_len(0);
        assert (dist >= 100);
        // std::cout << "distance = " << dist << std::endl;
        if (false_out_branches.find(node) != false_out_branches.end()) {
            for (auto node_pair : new_component) {
                edges_to_remove.push_back(node_pair);
            }
        }
        else {
            branching_components.push_back(std::make_pair(new_component, dist));
        }
        visited_out_branches.at(node) = true;
    }
    return;
}

std::pair<int, node_id_t> BranchReduction::extendComponentOut(
    std::vector< node_pair_t > & component,
    std::list< node_id_t > neighbors, bool & has_false_branch,
    std::unordered_map< node_id_t, bool > & visited_in_branches,
    std::unordered_map< node_id_t, bool > & visited_out_branches,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_in_map,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_out_map,
    std::unordered_map< node_id_t, int > & branch_out_dist_map) {
    // extend component iteratively
//    std::cout << "extendComponentOut" << std::endl;
    std::pair<int, node_id_t> dist_node_pair;
    bool extended = false;
    for (auto node : neighbors) {
        auto visited = visited_out_branches.find(node);
        if (visited == visited_out_branches.end() || visited->second == true) {
            continue;
        }
        if (false_out_branches.find(node) != false_out_branches.end()) {
            has_false_branch = true;
        }
        std::list< node_id_t > branch = branch_out_map.at(node);
        dist_node_pair = std::make_pair(branch_out_dist_map.at(node), node);
        extended = true;
        for (auto out_neighbor : branch) {
            assert (node != out_neighbor);
            node_pair_t node_pair = std::make_pair(node, out_neighbor);
            component.push_back(node_pair);
        }
        visited_out_branches.at(node) = true;
        extendComponentIn(component, branch, has_false_branch, visited_in_branches,
            visited_out_branches, branch_in_map, branch_out_map, branch_out_dist_map);
    }
    if (!extended) {
        dist_node_pair = std::make_pair(0, neighbors.front());
    }
    return dist_node_pair;
}

void BranchReduction::extendComponentIn(std::vector< node_pair_t > & component,
    std::list< node_id_t > neighbors, bool & has_false_branch,
    std::unordered_map< node_id_t, bool > & visited_in_branches,
    std::unordered_map< node_id_t, bool > & visited_out_branches,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_in_map,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_out_map,
    std::unordered_map< node_id_t, int > & branch_out_dist_map) {
    // extend component iteratively
//    std::cout << "extendComponentIn" << std::endl;
    for (auto node : neighbors) {
        auto visited = visited_in_branches.find(node);
        if (visited == visited_in_branches.end() || visited->second == true) {
            continue;
        }
        if (false_in_branches.find(node) != false_in_branches.end()) {
            has_false_branch = true;
        }
        std::list< node_id_t > branch = branch_in_map.at(node);
        for (auto in_neighbor : branch) {
            assert (node != in_neighbor);
            node_pair_t node_pair = std::make_pair(in_neighbor, node);
            component.push_back(node_pair);
        }
        visited_in_branches.at(node) = true;
        extendComponentOut(component, branch, has_false_branch, visited_in_branches,
            visited_out_branches, branch_in_map, branch_out_map, branch_out_dist_map);
    }
    return;
}

bool BranchReduction::countUniqueEvidence( std::vector< node_pair_t > component,
        int min_evidence, std::list< node_pair_t > & edges_to_remove) {
    // compute effective evidence per edge: remove evidence that is shared
    // between other edges and count evidence load
//    std::cout << "countUniqueEvidence" << std::endl;
    std::unordered_map< safe_edge_count_t, std::list< read_id_t > > unique_evidence_per_edge;
    std::unordered_map< safe_edge_count_t, unsigned int > index_to_mapID;
    std::vector< bool > evidence_status;
    bool typical_double_branch; // are we in a typical double branch situation (4 edges, 2 in-branches, 2 out-branches)
    std::set< node_id_t > in_nodes;
    std::set< node_id_t > out_nodes;
    unsigned int idx = 0;
    bool keep_component; // return value: true if any edges are kept, false otherwise
    for (auto node_pair : component) {
        assert (node_pair.first != node_pair.second);
        in_nodes.insert(node_pair.second);
        out_nodes.insert(node_pair.first);
        safe_edge_count_t mapID = edgeToEvidenceIndex(node_pair.first, node_pair.second);
        index_to_mapID.insert(std::make_pair(idx, mapID));
        auto evidence = evidence_per_edge.find(mapID);
        if (evidence == evidence_per_edge.end()) {
            std::cout << "mapID not found for edge " << node_pair.first << " " << node_pair.second << std::endl;
        }
        else {
            if (evidence_per_edge.at(mapID).empty()) {
                evidence_status.push_back(0);
            }
            else {
                evidence_status.push_back(1);
            }
        }
        unique_evidence_per_edge.insert(std::make_pair(mapID, std::list< read_id_t >()));
        idx++;
    }
    // check if we have a typical diploid branch: 4 nodes of which 2 in-nodes
    // and 2 out-nodes, connected by 3 or 4 edges
    if ((component.size() == 3 || component.size() == 4) && in_nodes.size() == 2 && out_nodes.size() == 2) {
        typical_double_branch = true;
    }
    else {
        typical_double_branch = false;
    }
    unsigned int component_edge_count = evidence_status.size();
    // filter evidence such that we only count UNIQUE read support
    while (*std::max_element(evidence_status.begin(), evidence_status.end()) == 1) {
        // get maximum value over all evidence iterators
        std::vector< read_id_t > current_evidence;
        for (unsigned int idx = 0; idx < component_edge_count; idx++) {
            if (evidence_status.at(idx) == 1) {
                safe_edge_count_t mapID = index_to_mapID.at(idx);
                assert (!evidence_per_edge.at(mapID).empty());
                current_evidence.push_back(evidence_per_edge.at(mapID).front());
            }
        }
        std::sort(current_evidence.begin(), current_evidence.end());
        assert (!current_evidence.empty());
        read_id_t current_min = current_evidence.front();
        read_id_t current_max = current_evidence.back();
//        std::cout << "current_max " << current_max << " current_min " << current_min << std::endl;
        assert (current_max < 2*program_settings.original_readcount);
        assert (current_min < 2*program_settings.original_readcount);
        bool unique_min;
        if (current_evidence.size() == 1) {
            unique_min = true; // only evidence remaining so definitely unique
        }
        else if (current_min < current_evidence.at(1)) {
            unique_min = true;
        }
        else {
            unique_min = false;
        }
        for (unsigned int idx = 0; idx < component_edge_count; idx++) {
            if (evidence_status.at(idx) == 1) {
                safe_edge_count_t mapID = index_to_mapID.at(idx);
                std::list< read_id_t > evidence = evidence_per_edge.at(mapID);
                if (evidence.front() == current_min) {
                    if (unique_min) { // keep evidence
                        unique_evidence_per_edge.at(mapID).push_back(current_min);
                    }
                    evidence_per_edge.at(mapID).pop_front();
                    if (evidence.size() == 1) {
                        // reached end of evidence list, mark node as finished
                        evidence_status.at(idx) = 0;
                    }
                }
            }
        }
    }
    // try to resolve typical branches
    if (program_settings.diploid && typical_double_branch) {
        //std::cout << "typical double branch" << std::endl;
        std::vector< std::pair<safe_edge_count_t, int> > pairs;
        for (auto ev : unique_evidence_per_edge) {
            pairs.push_back(std::make_pair(ev.first, ev.second.size()));
        }
        std::sort(pairs.begin(), pairs.end(), [=](const std::pair<safe_edge_count_t, int>& a, const std::pair<safe_edge_count_t, int>& b)
        {
            return a.second < b.second;
        }
        );
        std::vector< node_pair_t > supported_edges;
        std::vector< node_pair_t > unsupported_edges;
        int max_count = 0;
        node_pair_t max_edge;
        for (auto ev_idx : pairs) {
            node_pair_t node_pair = evidenceIndexToEdge(ev_idx.first);
            std::list< read_id_t > evidence = unique_evidence_per_edge.at(ev_idx.first);
            evidence.sort();
            evidence.unique();
            int count = evidence.size();
            if (count > max_count) {
                max_count = count;
                max_edge = node_pair;
            }
            if (count > 0) {
                supported_edges.push_back(node_pair); // supported edges are now ordered by evidence load
            }
            else {
                unsupported_edges.push_back(node_pair);
            }
            // std::cout << "evidence load for edge " << node_pair.first << "," << node_pair.second << " : ";
            // for (auto id : evidence) {
            //     std::cout << id << " ";
            // }
            // std::cout << std::endl;
        }
        keep_component = supported_edges.size() > 0; // return value indicating whether any edges were kept
        if (supported_edges.size() == 1) {
            // keep this edge and it's non-conflicting buddy
            // std::cout << "1 supported edge" << std::endl;
            for (auto remove_pair : unsupported_edges) {
                if (remove_pair.first == max_edge.first ||
                    remove_pair.second == max_edge.second)
                {
                    edges_to_remove.push_back(remove_pair);
                }
            }
            return keep_component;
        }
        else if (supported_edges.size() == 2 &&
                 supported_edges.at(0).first != supported_edges.at(1).first &&
                 supported_edges.at(0).second != supported_edges.at(1).second)
        {
            // std::cout << "2 supported edges" << std::endl;
            for (auto remove_pair : unsupported_edges) {
                edges_to_remove.push_back(remove_pair);
            }
            return keep_component;
        }
        else if (supported_edges.size() == 2) {
            std::cout << "conflict - 2 supported edges" << std::endl;
            // keep both if their evidence loads differ by at most 0.5*min_evidence
            // otherwise keep the max edge
            assert (pairs.size() >= 2);
            bool keep_complement = false;
            if (pairs.at(0).second - pairs.at(1).second > 0.5*min_evidence) {
                edges_to_remove.push_back(supported_edges.at(1));
                keep_complement = true;
            }
            for (auto remove_pair : unsupported_edges) {
                if (!keep_complement || remove_pair.first == max_edge.first ||
                    remove_pair.second == max_edge.second)
                {
                    edges_to_remove.push_back(remove_pair);
                }
            }
            return keep_component;
        }
        else if (supported_edges.size() > 2) {
            std::cout << supported_edges.size() << " supported edges" << std::endl;
            // keep the non-conflicting pair of edges of highest support
            // first try max_edge + complement
            int load1 = 0;
            int load2 = 0;
            int i = 0;
            for (auto remove_pair : supported_edges) {
                if (remove_pair != max_edge && (
                    remove_pair.first == max_edge.first ||
                    remove_pair.second == max_edge.second))
                {
                    load2 += pairs.at(i).second;
                }
                else {
                    load1 += pairs.at(i).second;
                }
                i++;
            }
            if (load1 >= load2) {
                // keep max_edge and complement
                for (auto remove_pair : unsupported_edges) {
                    if (remove_pair != max_edge && (
                        remove_pair.first == max_edge.first ||
                        remove_pair.second == max_edge.second))
                    {
                        edges_to_remove.push_back(remove_pair);
                    }
                }
                for (auto remove_pair : supported_edges) {
                    if (remove_pair != max_edge && (
                        remove_pair.first == max_edge.first ||
                        remove_pair.second == max_edge.second))
                    {
                        edges_to_remove.push_back(remove_pair);
                    }
                }
            }
            else {
                // keep alternative edge pair
                for (auto remove_pair : unsupported_edges) {
                    if (remove_pair == max_edge || (
                        remove_pair.first != max_edge.first &&
                        remove_pair.second != max_edge.second))
                    {
                        edges_to_remove.push_back(remove_pair);
                    }
                }
                for (auto remove_pair : supported_edges) {
                    if (remove_pair == max_edge || (
                        remove_pair.first != max_edge.first &&
                        remove_pair.second != max_edge.second))
                    {
                        edges_to_remove.push_back(remove_pair);
                    }
                }
            }
            return keep_component;
        }
    }
    // finally analyze evidence and reduce branches by removing edges with
    // insufficient evidence
    //std::cout << "effective_evidence_counts" << std::endl;
    node_id_t node1, node2;
    node_pair_t node_pair;
    keep_component = false; // gets updated below, indicating whether any edges are kept
    for (auto ev : unique_evidence_per_edge) {
        node_pair_t node_pair = evidenceIndexToEdge(ev.first);
        node1 = node_pair.first;
        node2 = node_pair.second;
        std::list< read_id_t > evidence = ev.second;
        evidence.sort();
        evidence.unique();
        int count = evidence.size();
        if (count < min_evidence) {
            if (overlap_graph->checkEdge(node1, node2, false) < 0) {
                std::cout << "edge not found, reverse: "
                    << overlap_graph->checkEdge(node2, node1, false)
                    << std::endl;
                exit(1);
            }
            edges_to_remove.push_back(node_pair);
        }
        else {
            keep_component = true;
        }
        if (program_settings.verbose) {
            std::cout << "evidence load for edge " << node1 << "," << node2 << " : ";
            for (auto id : evidence) {
                std::cout << id << " ";
            }
            std::cout << std::endl;
        }
    }
    return keep_component;
}

safe_edge_count_t BranchReduction::edgeToEvidenceIndex(node_id_t node, node_id_t neighbor) {
    // map edge (node pair) to a unique index for storage of evidence in evidence_per_edge
    assert (node < overlap_graph->getVertexCount());
    assert (neighbor < overlap_graph->getVertexCount());
    safe_edge_count_t index = node * overlap_graph->getVertexCount() + neighbor;
    assert (index < (overlap_graph->getVertexCount())*(overlap_graph->getVertexCount()));
    return index;
}

node_pair_t BranchReduction::evidenceIndexToEdge(safe_edge_count_t index) {
    // map storage index back to original node pair
    assert (index < (overlap_graph->getVertexCount())*(overlap_graph->getVertexCount()));
    node_id_t node = floor(index / overlap_graph->getVertexCount());
    node_id_t neighbor = index - node * overlap_graph->getVertexCount();
    assert (node < overlap_graph->getVertexCount());
    assert (neighbor < overlap_graph->getVertexCount());
    return std::make_pair(node, neighbor);
}

bool BranchReduction::compareNodepairs (node_pair_t pair1, node_pair_t pair2) {
    if (pair1.first != pair2.first) {
        return pair1.first < pair2.first;
    }
    else {
        return pair1.second < pair2.second;
    }
}
