#include <fstream>
#include <cassert>
#include<cmath>
#include <unordered_map>
#include "caller.h"

using namespace std;
std::unordered_map<int,int> empty_dict;
std::unordered_map<int,int> none_dict= {{-1, -1}};
std::deque<std::pair<int,int>> variantslist;
std::deque<std::pair<int,int>> enum_refkmers;
std::deque<std::pair<int,int>> enum_kmers;
std::unordered_map<int,int> seen_comb;
std::unordered_map<int,int> comb_count;
ofstream writer;
Caller::Caller(std::string &reference, int k, int window){
	this->k=k;
	this->enumerate_reference_kmers(reference,k);
	this->ref_kmer_generator = enum_refkmers;
	this->i1= this->ref_kmer_generator.begin();
	kmer= (*this->i1).first;//first reference kmer that any read mapped to
	pos=(*this->i1).second;//position of the first reference kmer the read mapped to
	if (this->i1 != this->ref_kmer_generator.end()){
		this->i1++;
	}
	this->pileup_columns.push_back(empty_dict);
	this->ref_kmers.push_back(kmer);
	this->ref_pos = pos;
	this->window=window;
}
void Caller::all_variants(std::deque<std::pair<int,int>> &variant_list){
	variantslist =variant_list;
}
void Caller::final_pop(const std::string outfile){
	int final_ref_pos= (enum_refkmers[enum_refkmers.size()-1]).second;//final position in the reference
	this->process_complete_columns(final_ref_pos, outfile);
}
void Caller::add_read(int bam_alignment_pos, std::vector<std::vector<int>> &bam_alignment_cigartuples, std::string &bam_alignment_query, const std::string outfile){
	this->target_pos_pre= bam_alignment_pos;//initially the target position is where the read maps to reference
	this->enumerate_kmers(bam_alignment_pos,bam_alignment_query, this->k, bam_alignment_cigartuples);
	this->kmer_generators.push_back(enum_kmers);//all the read kmers enumerated and pushed as a deque to a deque of all the previously generated read kmers
	this->kmer_generators_finished.push_back(false);//not yet finished with kmer generators
	int iterator_index = this->kmer_generators.size()-1;//index of the last deque (just added) of enum read kmers
	this->iterators.push_back(this->kmer_generators[iterator_index].begin());//set an iterator to the beginning of the last added enum kmers deque in kmer_generators
	temp_pair2= *this->iterators[this->iterators.size()-1];//get the read kmer and its position corresponding to the iterator pushed in the previous step
	kmer= temp_pair2.first;
	pos= temp_pair2.second;
	this->current_kmers.push_back(std::pair<int,int> (kmer,pos));//add to the list of kmers currently being dealt with
	if(this->iterators[this->iterators.size()-1]!= this->kmer_generators[iterator_index].end()-1){
			this->iterators[this->iterators.size()-1]++;
	}//if the respesctive iterator is not already at the end move it one step forward
	std::pair<int, int> pair_getcolumn = this->get_column(pos);//this would give the reference kmer that this read kmer just read correspond to
	//and also its index in the dictionary where all the columns are being stored. SO, at this index we have a record of what read kmers have
	//been mapped to this reference kmer till now and their count
	ref_kmer= pair_getcolumn.first;
	signed index_getcolumn= pair_getcolumn.second;
	if (index_getcolumn>=0){
		this->pileup_columns[index_getcolumn][temp_pair2.first]+=1;//increment the count for this readkmer mapping to the ref kmer we got above
	}
	int target_pos= this->target_pos_pre + this->k - 1;
	this->process_complete_columns(target_pos, outfile);
}

void Caller::finish(){
	;
}
std::pair<int, int> Caller::get_column(int pos){
	//
	int index = pos - this->ref_pos;//how far the read kmer is away from the reference start, it tells how much of the reference we have already covered
	if (index >=0){//if we are at or to the right of the reference start
		while (int (this->pileup_columns.size()) <= index){
			kmer= (*this->i1).first;
			pos= (*this->i1).second;
			if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}
			this->ref_kmers.push_back(kmer);
			this->pileup_columns.push_back(empty_dict);
		}//for all positions traversed to the right starting from the ref pos, get the ref kmers, push then to ref_kmers deque
		//and push an empty dictionary for each of them to the pileup columns deque
		return std::pair<int,int> (this->ref_kmers[index], index);//the ref_kmer at the position we are right now, and the actual position
	} else {//if we are somewhere on the left of ref pos(which I don't really know how would happen)
		kmer= (*this->i1).first;
		pos= (*this->i1).second;
		if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}

		this->ref_kmers.push_back(kmer);

		this->pileup_columns.push_back(empty_dict);

		return std::pair<int,int> (kmer,-1);
	}

}
void Caller::pop_column(){
	//finally output the ref kmer , read kmer, counts for each position that is complete
	int result_ref_pos;
	int result_ref_kmer;
	std::unordered_map<int,int> result_pileup_kmers;
	int result_kmer;
	int result_count;
	if (int (this->pileup_columns.size()) > 0){
		result_ref_pos= this->ref_pos;
		result_ref_kmer= this->ref_kmers.front();
		result_pileup_kmers= this->pileup_columns.front();
		this->ref_kmers.pop_front();
		this->pileup_columns.pop_front();
	}
	else {//if there is no ref kmer dictionary in pile up columns,
		// add one for the first non read ref kmer from the ref_kmer_generator and move iterator to next position
		kmer= (*this->i1).first;
		pos=(*this->i1).second;
		if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}
		assert (pos == this->ref_pos);
		result_ref_pos= this->ref_pos;
		result_ref_kmer=kmer;
		result_pileup_kmers= none_dict;
	}
	this->ref_pos += 1;
	std::pair<int,int> var=variantslist.front();
	int variantposition=var.first;
	int var_length= var.second-1;
	int varstart= variantposition-window;
	int varend= variantposition+var_length+window+this->k-1;
	//as kmer end positions are being recorded, we skip first k-1 kmers after the window as they would be starting inside
	std::pair<int,int> next_var= variantslist[1];
	int next_variantposition= next_var.first;
	int next_var_length= next_var.second-1;
	if (result_ref_pos>=varstart and result_ref_pos<=varend){
		//if position inside variant window don't output anything
		;
	}
	else if (int(variantslist.size())>0 and result_ref_pos>= (next_variantposition-window) and result_ref_pos<= (next_variantposition+next_var_length+window)){
		//if reference position such that it is not inside the current variant window but is inside the next
		 //it is safe to remove the first variant as we know all readkmers would be further to the right
			variantslist.pop_front();
	}
	else{
	if (result_pileup_kmers!=none_dict and int (result_pileup_kmers.size()) > 0){
		for(std::unordered_map<int,int> ::iterator it = result_pileup_kmers.begin(); it != result_pileup_kmers.end(); ++it){
			result_kmer=it->first;
			result_count=it->second;
		  writer<<result_ref_pos<<"\t"<<result_ref_kmer<<"\t"<<result_kmer<<"\t"<<result_count<<endl;

	}}}}
 void Caller::process_complete_columns(int target_pos, const std::string outfile){
	/*
	Perform calling of columns that are complete, i.e. they cannot receive
	more reads because subsequent reads are further to the right'''
	# compute largest position for which k-mer pileup can safely be generated*/

	std::vector<std::pair<std::pair<int,int>,std::unordered_map<int,int>>> complete_columns;
	this->advance_to(target_pos);//get all readkmers from all reads read so far upto target pos and add their counts
	writer.open(outfile, ios::app);
	while (this->ref_pos <= target_pos){
		this->pop_column();
		}
	writer.close();
	}
void Caller::advance_to(int target_pos){
	/*
	Add all k-mer from all reads up to target_pos to pileup_columns.
	*/
	int i;
	for (i=0; i< int (this->kmer_generators.size()) ;i++){
		std::deque<std::pair<int,int>> kmer_generator= this->kmer_generators[i];//getting all the enumerated kmer deques
		try{
			temp_pair3 = this->current_kmers[i];//for the enum kmer deque we got above what was the kmer, pos we have read last
			kmer = temp_pair3.first;
			pos= temp_pair3.second;
			while (pos <= target_pos){//this would keep on popping read kmers until they are towards the left of our target position
				if(this->iterators[i]!= this->kmer_generators[i].end()){
					temp_pair3= *this->iterators[i];//the read kmer and its position based on iterator index in each deque of enum kmers
					kmer=temp_pair3.first;
					pos=temp_pair3.second;
					std::pair<int, int> pair_getcolumn = this->get_column(pos);//check to which ref kmer this read kmer maps to
					ref_kmer= pair_getcolumn.first;
					signed index_getcolumn= pair_getcolumn.second;
					if (index_getcolumn>=0){
						this->pileup_columns[index_getcolumn][temp_pair3.first]+=1;//increase the count for this read kmer
					}
					this->iterators[i]++;//move the iterator in this part to one position forward
				}
				else{

					throw(1);
				}
			}
			this->current_kmers[i] = std::pair<int,int> (kmer,pos);//update the read kmer just read
		}
		catch(int e){
			this->kmer_generators_finished[i] = true;//the read kmers have all been dealt with
		}
	}
	while ((this->kmer_generators.size() > 0) && this->kmer_generators_finished[0]){
		//removing everything for the reads that have been dealt completely
		this->current_kmers.pop_front();
		this->kmer_generators.pop_front();
		this->iterators.pop_front();
		this->kmer_generators_finished.pop_front();
	}
}
void Caller::enumerate_reference_kmers(std::string &reference, int k){
	int A = 'A';
  int C = 'C';
  int G = 'G';
  int T = 'T';
  int c = 0;
  int h = 0;
  int mask = (1 << (2*k)) - 1;
  int i=0;
	enum_refkmers.clear();
	int length_reference = reference.length();

  for (i=0; i<length_reference; i++){
		c = reference[i];
		if (c == A){
			h = ((h << 2) | 0) & mask;
				}
		else if (c == C){
			h = ((h << 2) | 1) & mask;
				}
		else if (c == G){
			h = ((h << 2) | 2) & mask;
				}
		else if (c == T){
			h = ((h << 2) | 3) & mask;
				}
		//else{
		//	h = ((h << 2) | 0) & mask;
		//		}
		if (i >= k-1){
			enum_refkmers.push_back(std::pair<int,int> (h,i+1));
				}
		}
}

void Caller::enumerate_kmers(int pos, std::string &query_string, int k, std::vector<std::vector<int>> &cigartuples){
	//Generates all kmers in the read and yields pairs (kmer_hash, position),
	//where kmer_hash is a binary representation of the kmer and postion is
	//the reference position this kmer has been aligned to.
        int h = 0;
        //int n_kmer=0;
        int mask = (1 << (2*k)) - 1;
        int cigar_index = 0;
        int cigar_op = cigartuples[cigar_index][0];
        int cigar_length = cigartuples[cigar_index][1];
        int BAM_CMATCH = 0;     // M
        int BAM_CINS = 1 ;      // I
        int BAM_CDEL = 2 ;      // D
        int BAM_CREF_SKIP = 3;  // N
        int BAM_CSOFT_CLIP = 4; // S
        int BAM_CHARD_CLIP = 5; // H
        int BAM_CPAD = 6 ;     // P
        int BAM_CEQUAL = 7 ;   //=
        int BAM_CDIFF = 8  ;    // X
        int BAM_CBACK = 9 ;     // B
        int A = 'A';
        int C = 'C';
        int G = 'G';
        int T = 'T';
        int i =0;
        int consecutive = 0;
				enum_kmers.clear();
				int length_query_string= query_string.length();
        while (i < length_query_string){
					//process cigar entries that don't consume a character from the read
					while (true){
						if ((cigar_op == BAM_CDEL) or (cigar_op == BAM_CREF_SKIP)){
							pos += cigar_length;
						}
						else if (cigar_op == BAM_CSOFT_CLIP){
							consecutive = 0;
						}
						else if ((cigar_length == 0) or (cigar_op == BAM_CHARD_CLIP)){

						}
						else{
							break;
						}
            cigar_index += 1;
						cigar_op = cigartuples[cigar_index][0];
						cigar_length = cigartuples[cigar_index][1];
				}
				if (i >= length_query_string){
					break;
				}
				// update hash
				if (query_string[i] == A){
					h = ((h << 2) | 0) & mask;
					consecutive += 1;
				}
				else if (query_string[i] == C){
					h = ((h << 2) | 1) & mask;
					consecutive += 1;
				}
				else if (query_string[i] == G){
					h = ((h << 2) | 2) & mask;
					consecutive += 1;
				}
				else if (query_string[i] == T){
					h = ((h << 2) | 3) & mask;
					consecutive += 1;
				}
				else {
				//	n_kmer=1;
					consecutive += 1;
				}
				if (consecutive >= k){
					//if (n_kmer ==0){
						enum_kmers.push_back(std::pair<int,int> (h,pos+1));
						//}
					//else if (n_kmer==1){
					//	n_kmer=0;
					//	}
				}
				//consume one character of read
				assert (cigar_length > 0);
				if ((cigar_op == BAM_CMATCH) or (cigar_op == BAM_CEQUAL) or (cigar_op == BAM_CDIFF)){
					cigar_length -= 1;
					pos += 1;
				}
				else if (cigar_op == BAM_CINS){
					cigar_length -= 1;
				}
				else {
					assert (false);
				}
				i += 1;
				}
		}
