#ifndef CALLER_H
#define CALLER_H
#include<iostream>
#include <deque>
#include <string>
#include <vector>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <thread>
#include <utility>
class Caller{
public:
	Caller(std::string &reference, int k, int window);
	void add_read(int bam_alignment_pos,std::vector<std::vector<int>> &bam_alignment_cigartuples, std::string &bam_alignment_query, const std::string outfile);
	void all_variants(std::deque<std::pair<int,int>> &variant_list);
	void finish();
	std::pair<int, int> get_column(int pos);
	void pop_column();
	void process_complete_columns(int target_pos, const std::string outfile);
	void advance_to(int target_pos);
 	void enumerate_reference_kmers(std::string &reference, int k);
	void enumerate_kmers(int pos, std::string &query_string, int k, std::vector<std::vector<int>> &cigartuples);
	void final_pop(const std::string outfile);
	int kmer;
	int pos;
	int k;
	int window;
	int ref_pos;
	int target_pos_pre;
	int ref_kmer;
	std::deque<std::deque<std::pair<int,int>>::iterator> iterators;
	std::deque<std::deque<std::pair<int,int>>> kmer_generators;
	std::deque<int> kmer_generators_finished;
	std::deque<std::pair<int,int>> current_kmers;
	std::deque<std::pair<int,int>> ref_kmer_generator;
	std::deque<std::unordered_map<int,int>> pileup_columns;
	std::deque<int> ref_kmers;
	std::unordered_map<int,int> pileup_column;
	std::deque<std::pair<int,int>>::iterator i1;
	std::pair<int,int> temp_pair1;
	std::pair<int,int> temp_pair2;
	std::pair<int,int> temp_pair3;


};
#endif
