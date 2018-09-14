/*

  Copyright (C) 2017-2018 Stefano Beretta, Marco Dell'Acqua, Murray
  Patterson, Yuri Pirola, Simone Zaccaria

  Distributed under the MIT license.

  You should have received a copy of the MIT license along with this
  program.

*/

#include <stdlib.h>
#include <string.h>
#include <limits>
#include <math.h>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <stdexcept>

#include "basictypes.h"
#include "binomialcoefficient.h"
#include "combinations.h"
#include "balancedcombinations.h"
#include "hapchatcolumniterator.cpp"

// Log messages with DEBUG priority and higher
#define LOG_MSG
#define LOG_THRESHOLD LOG_LEVEL_FATAL
// Include log facilities. It should the last include!!
#include "log.h"

class HapChatCore{

private:

  ReadSet* readset_;
  int sample_id_;
  vector<pair<Read*,Read*>> superreads_;
  unsigned int optimal_;

  bool unweighted_;

  double errorrate_;
  double alpha_;

  bool balancing_;
  unsigned int balancecov_;
  double balanceratio_;

public:	

  HapChatCore(ReadSet* read_set)
    : readset_(read_set),
      sample_id_(0),
      optimal_(0u),
      unweighted_(false),
      errorrate_(0.05),
      alpha_(0.01),
      balancing_(false),
      balancecov_(20),
      balanceratio_(0.5) {

    DEBUG(read_set->toString());

    readset_->reassignReadIds();
    if(readset_->size() > 0) {
      // since hapchat should only receive the reads corresponding to
      // a single individual, we can suppose that the sample_id is
      // uniform, so we take said from the first read
      sample_id_ = readset_->get(0)->getSampleID();
    }

    superreads_.emplace_back(
			     new Read("superread_0_0",-1,-1,sample_id_),
			     new Read("superread_1_0",-1,-1,sample_id_)
			     );

    run_core();
  }


void get_super_reads(vector<ReadSet*>* output_read_set) {

  for(unsigned int k=0;k<superreads_.size();k++) {

    output_read_set->at(k)->add(superreads_[k].first);
    output_read_set->at(k)->add(superreads_[k].second);
  }
}


unsigned int get_optimal_cost() {

  return optimal_;
}


vector<bool> * get_optimal_partitioning() {
  vector<bool> * partitioning = new vector<bool>(readset_->size(),false);
  // TODO: implement this
  return partitioning;
}


int get_length() {

  return superreads_.size();
}


void run_core()
{

  INFO("HapCHAT Starting...");

  //Counter threshold_coverage = 30;

  const constants_t constants;
  
  INFO("Parameters:");
  INFO("Discard weights? " << (unweighted_?"True":"False"));
  INFO("Error rate: " << errorrate_);
  INFO("Alpha: " << alpha_);
  INFO("Balancing? " << (balancing_?"True":"False"));
  INFO("Balancing activates after coverage: " << balancecov_);
  INFO("Balance ratio: " << balanceratio_);

  //Readset section
  HapChatColumnIterator hap=HapChatColumnIterator(readset_);
  if(!hap.has_next_block()) {
    return;
  }

  //Initializing the starting parameters: no competitive section

  //Pre-compute binomial values
  BinomialCoefficient::initialize_binomial_coefficients(MAX_COVERAGE, MAX_COVERAGE);
  computeK(MAX_COVERAGE, alpha_, errorrate_);

  Counter step = 0;
  Cost OPT = 0;
  Counter counter_columns;

  Counter MAX_COV = 0;
  Counter MAX_L = 0;
  Counter MAX_K = 0;
  Counter MAX_GAPS = 0;

  vector<vector<char> > haplotype_blocks1;
  vector<vector<char> > haplotype_blocks2;

  hap.reset();

  vector<bool> haplotype1(hap.column_count());
  vector<bool> haplotype2(hap.column_count());

  if(hap.column_count() > 0)
    dp(constants, haplotype1, haplotype2, step, OPT, MAX_COV, MAX_L, MAX_K, MAX_GAPS, hap);

  counter_columns = hap.column_count();

  vector<char> output_block1(hap.column_count());
  vector<char> output_block2(hap.column_count());

  fill_haplotypes(haplotype1, haplotype2, output_block1, output_block2, hap);
  DEBUG("Filled haplotypes");

  haplotype_blocks1.push_back(output_block1);
  haplotype_blocks2.push_back(output_block2);
	
  INFO("");
  INFO("OPTIMUM:  " << OPT);
  INFO("");

  INFO("MAX_COV:  " << MAX_COV);
  INFO("MAX_L:  " << MAX_L);
  INFO("MAX_K:  " << MAX_K);
  INFO("MAX_GAPS:  " << MAX_GAPS);
  INFO("# of columns:  "  << counter_columns);
  INFO("");

  DEBUG("<<>> Writing haplotypes...");
  //write_haplotypes(haplotype_blocks1, haplotype_blocks2, std::cout);

  make_super_reads(hap.get_positions(), haplotype_blocks1[0], haplotype_blocks2[0]);
  optimal_ = OPT.get_cost();

}


Pointer next(const Pointer &indexer_pointer, const int &total_size, const int &shift)
{
  return (indexer_pointer + shift) % total_size;
}


Pointer prev(const Pointer &indexer_pointer, const int &total_size, const int &shift)
{
  return (indexer_pointer + total_size - shift) % total_size;
}

 
bool check_end(HapChatColumnIterator hap, const vector<Column> &input, const Pointer &pointer, const bool &re_run)
{
  return (!re_run && (hap.is_ended() && (input[pointer][0].get_read_id() == -1)));
}

 
void complement_mask(BitColumn &mask, const Counter &length, const constants_t &constants)
{
  mask ^= (constants.ones<<length).flip();
}

 
string column_to_string(const BitColumn &mask, const unsigned int &len) {
  string str = mask.to_string();
  reverse(str.begin(), str.end());
  return str.substr(0, len);
}


template <typename T>
 
void replace_if_less(T& a, const T& b) {
  if(a > b) { a = b;}
}


void fill_haplotypes(const vector<bool> &haplotype1, const vector<bool> &haplotype2, vector<bool> &complete_haplo1, vector<bool> &complete_haplo2, HapChatColumnIterator hap)
{

  hap.reset();
  vector<bool>::const_iterator ihap1 = haplotype1.begin();
  vector<bool>::const_iterator ihap2 = haplotype2.begin();

  vector<bool>::iterator iout1 = complete_haplo1.begin();
  vector<bool>::iterator iout2 = complete_haplo2.begin();

  while(hap.has_next()) { 

    hap.get_column();
    *iout1 = *ihap1;
    *iout2 = *ihap2;

    ++ihap1;
    ++ihap2;
    ++iout1;
    ++iout2;
  }
}


void fill_haplotypes(const vector<bool> &haplotype1, const vector<bool> &haplotype2, vector<char> &output_block1, vector<char> &output_block2, HapChatColumnIterator hap)
{

  hap.reset();
  vector<bool>::const_iterator ihap1 = haplotype1.begin();
  vector<bool>::const_iterator ihap2 = haplotype2.begin();

  vector<char>::iterator iout1 = output_block1.begin();
  vector<char>::iterator iout2 = output_block2.begin();

  while(hap.has_next()) {

    hap.get_column();
    *iout1 = (*ihap1)? '1' : '0';
    *iout2 = (*ihap2)? '1' : '0';

    ++ihap1;
    ++ihap2;
    ++iout1;
    ++iout2;
  }
}


void write_haplotypes(const vector<vector<char> > &haplotype_blocks1, const vector<vector<char> > &haplotype_blocks2, ostream &ofs)
{
  vector<vector<char> >::const_iterator ivv = haplotype_blocks1.begin();
  ofs << *ivv;
  for(ivv = haplotype_blocks1.begin() + 1;
      ivv != haplotype_blocks1.end();
      ++ivv) {
    ofs << '|' << *ivv;
  }
  ofs << endl;

  ivv = haplotype_blocks2.begin();
  ofs << *ivv;
  for(ivv = haplotype_blocks2.begin() + 1;
      ivv != haplotype_blocks2.end();
      ++ivv) {
    ofs << '|' << *ivv;
  }
  ofs << endl;
}


void dp(const constants_t &constants,
        vector<bool> &haplotype1, vector<bool> &haplotype2, Counter &step_global, Cost &OPT_global,
        Counter &MAX_COV_global, Counter &MAX_L_global, Counter &MAX_K_global,
        Counter &MAX_GAPS_global, HapChatColumnIterator hap)
{
	
  Counter MAX_COV = 0;
  Counter MAX_K = 0;
  Counter MAX_L = 0;
  Counter MAX_GAPS = 0;
  Counter num_col = 0;
  vector<Counter> sum_successive_L;
  vector<vector<Counter> > scheme_backtrace;
	hap.reset();
  computeInputParams(num_col, MAX_COV, MAX_L, MAX_K, MAX_GAPS, sum_successive_L,
                     scheme_backtrace, hap);

  MAX_COV_global = max(MAX_COV_global, MAX_COV);
  MAX_K_global = max(MAX_K_global, MAX_K);
  MAX_L_global = max(MAX_L_global, MAX_L);
  MAX_GAPS_global = max(MAX_GAPS_global, MAX_GAPS);


  DEBUG(">> Initialized starting parameters");
  INFO("::== Starting parameters:  MAX_COV = " << MAX_COV << " // MAX_L = " << MAX_L << " // MAX_K = " << MAX_K << " // MAX_GAPS = " << MAX_GAPS);
  DEBUG("::== no of columns:     " << num_col);
  //DEBUG("-->> sum_successive_L:  " << sum_successive_L);

  //.:: ALLOCATION MEMORY

  DEBUG(">> Starting allocation of memory");
  //Allocation of memory for input window
  vector<Column> input(2 * (MAX_L - 1) + 1,
                       Column(MAX_COV,
                              Entry(-1, Entry::BLANK, 0)));
  Pointer input_pointer = 0;
  TRACE("-->> input allocated");


  //Allocation of memory for backward indexer
  //The index in j of the shared elements between p and j
  vector<vector<Pointer> > back_indexer(2 * (MAX_L - 1) + 1,
                                      vector<Pointer>(MAX_COV,
                                                      -1));
  //Equal to indexer_pointer
  TRACE("-->> back indexer allocated");

  //Allocation of memory for forward indexer
  //The index in p of the shared elements between p and j
  vector<vector<Pointer> > forw_indexer(2 * (MAX_L - 1) + 1,
                                      vector<Pointer>(MAX_COV,
                                                      -1));
  const Pointer indexer_pointer = MAX_L - 1;
  TRACE("-->> forw indexer allocated");

  //Allocation of memory for pos_gaps
  //The considered gaps are the ones in the column with the lower index
  vector<BitColumn> pos_gaps(2 * (MAX_L - 1) + 1);
  //Equal to indexer_pointer
  TRACE("-->> pos_gaps allocated");

  //Allocation of memory for num_pos_gaps
  //The considered gaps are the ones in the column with the lower index
  vector<Counter> num_pos_gaps(2 * (MAX_L - 1) + 1);
  //Equal to indexer_pointer
  TRACE("-->> num_pos_gaps allocated");


  //Allocation of memory for vector of k_j
  vector<Counter> k_j(2 * (MAX_L - 1) + 1,
                      MAX_K);
  //its pointer is equal to input_pointer
  TRACE("-->> k_j allocated");

  //Allocation of memory for homozigous costs
  vector<Counter> homo_cost(2 * (MAX_L - 1) + 1,
                            MAX_COUNTER);
  //its pointer is equal to input_pointer
  TRACE("-->> homo_cost allocated");

  //Allocation of memory for homozigous weights
  vector<Cost> homo_weight(2 * (MAX_L - 1) + 1,
                           Cost::INFTY);
  //its pointer is equal to input_pointer
  TRACE("-->> homo_weight allocated");

  //Allocation of memory for prevision matrix
  //[Destinatary of prevision][Who make the prevision][Indexof(mask of who makes prevision on common fragments)]

  vector<vector<vector<Cost> > > prevision(MAX_L,
                                           vector<vector<Cost> > (MAX_L,
                                                                  vector<Cost>(0)));

  for(unsigned int j = 0; j < MAX_L; j++) {
    for(unsigned int q = 0; q < MAX_L; q++) {
      prevision[j][q].resize(sum_successive_L[q],
                             Cost::INFTY);
    }
  }
  Pointer prevision_pointer = 0;
  TRACE("-->> prevision allocated");

  //Allocation of memory for OPT vector
  vector<Cost> OPT(MAX_L + 1,
                   Cost::INFTY);        //+ 1 since I need OPT[j - L]
  Pointer OPT_pointer = 0;
  TRACE("-->> OPT allocated");

  //Allocation of memory for backtrace column
  vector<vector<vector<Backtrace1> > > backtrace_table1(num_col);
  vector<vector<vector<bool> > > backtrace_table2_haplotypes(num_col);
  vector<vector<vector<bool> > > backtrace_table2_new_block(num_col);

  for(unsigned int j = 0; j < num_col; j++) {
      backtrace_table1[j].resize(scheme_backtrace[j].size());
      DEBUG("j: " << j << " scheme[j] size: " << scheme_backtrace[j].size());
      backtrace_table2_haplotypes[j].resize(scheme_backtrace[j].size());
      backtrace_table2_new_block[j].resize(scheme_backtrace[j].size());
      for(unsigned int q = 0; q < backtrace_table1[j].size(); q++) {
          backtrace_table1[j][q].resize(scheme_backtrace[j][q]);
          DEBUG("j: " << j << " q: " << q << " scheme[j] size: " << scheme_backtrace[j][q]);
          backtrace_table2_haplotypes[j][q].resize(scheme_backtrace[j][q]);
          backtrace_table2_new_block[j][q].resize(scheme_backtrace[j][q]);
      }
  }
  TRACE("-->> Backtrace table allocated");


  vector<bool> is_homozygous(num_col);
  vector<bool> homo_haplotypes(num_col);
  vector<Backtrace1> best_heterozygous1(num_col);
  vector<bool> best_heterozygous2_haplotypes(num_col);
  vector<bool> best_heterozygous2_new_block(num_col);


  DEBUG(">> Completed allocation of memory");

  //INITIALIZATION

  Combinations generator;
  BalancedCombinations balanced_generator;
  hap.reset();
  Counter step = 0;
  Column column;

  //Place the first and the next L columns in the positions of input data structure

  input_pointer = 0;

  Counter l = 0;
  bool flag = true;
	
  do {
    
    Pointer new_l_pointer = next(input_pointer, input.size(), l);
		
    if(l == 0) {
      column.clear();
    } else {
			flag = hap.has_next();
      column=hap.get_column();
      
    }
    insert_col_and_update(input, k_j, homo_cost, homo_weight, new_l_pointer, column, homo_haplotypes, step + l);

    l++;
    //The short circuit && is fundamental to avoid an unexpected has_next that read a column
  } while(l < MAX_L && flag);

  DEBUG(">> Initialization completed");


  //  .::: BASE CASE :::.
  DEBUG(".:: Basic Step: " << step);

  BitColumn colj;
  BitColumn gaps_mask;
  BitColumn proj;
  BitColumn mask;
  BitColumn comb_no_gaps;
  BitColumn comb_gaps;
  Cost current_cost(Cost::INFTY);
  Cost current_best(Cost::INFTY);
  Counter cov_j(0);
  Counter num_gaps(0);
  bool feasibility;
  bool has_successive;
  bool solution_existence(true);
  Counter temp_jump(-1);
  Counter temp_index(0);
  bool temp_haplotypes(false);
  bool temp_new_block(false);

  BitColumn cut_mask;
  BitColumn mask_qj;

  //Base case for OPT
  OPT[OPT_pointer] = 0;

  //k_j, homo_weight and homo_cost for base case
  k_j[input_pointer] = 0;
  homo_weight[input_pointer] = 0;
  homo_cost[input_pointer] = 0;

  current_cost = 0;

  //Make a prevision for all the successive column
  has_successive = true;
  Counter p = 1;
  //XXX: Is it necessary?
  do {
    //All the condition that have to be satisfied for the next column
    const Pointer new_p_pointer = next(input_pointer, input.size(), p - 1);
    feasibility = (p - 1 == 0) || (homo_cost[new_p_pointer] <= k_j[new_p_pointer]);

    if (p >= MAX_L || forw_indexer[indexer_pointer + p][0] == -1 || !feasibility) {
      has_successive = false;
    } else {
      //The number of elements shared between p and j
      Pointer new_prevision_pointer = next(prevision_pointer, prevision.size(), p);

      prevision[new_prevision_pointer][p].resize(1);
      prevision[new_prevision_pointer][p][0] = current_cost;

      p++;
    }
  } while (has_successive);

  // INC-K
  double k_j_inc = k_j[input_pointer];
  bool re_run_k = false;

  DEBUG("-->> Basic case completed  -- current_cost: " << current_cost);

  //DP
  
  //For all the columns

  while(!check_end(hap, input, next(input_pointer, input.size(), 1), re_run_k))// INC-K && solution_existence)
    {
      current_best = Cost::INFTY;
      solution_existence = false;
      temp_jump = -1;
      temp_index = 0;
      temp_haplotypes = false;
      temp_new_block = false;
      step++;
      step_global++;
      DEBUG("STARTING STEP:  " << step);

      // >>>>>>>>>>>>>>>>>>>>>> UPDATE DATA STRUCTURE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //.:: Update input
			
      //.:: Update common pointer
      input_pointer = next(input_pointer, input.size(), 1);

      Pointer new_input_pointer = next(input_pointer, input.size(), MAX_L - 1);

      // INC-K 
      if(!re_run_k) {
          //.:: Read Column
          column = hap.get_column();
          
          
          insert_col_and_update(input, k_j, homo_cost, homo_weight, new_input_pointer,
                                column, homo_haplotypes, step + (MAX_L - 1));
          k_j_inc = k_j[input_pointer];
      }

      //.:: Update indexers


      //For all the q successive columns
      for(unsigned int q = 1; q < MAX_L; q++)
        {
          intersect(input[next(input_pointer, input.size(), q)],
                    input[input_pointer],
                    indexer_pointer + q,
                    forw_indexer, back_indexer,
                    pos_gaps, num_pos_gaps, false);

          //If the just considered column q did not have any common element, we can fill as
          //  empty all the remaining columns and terminate.
          if(forw_indexer[indexer_pointer + q][0] == -1)
            {
              for(unsigned int p = q + 1; p < MAX_L; p++)
                {
                  forw_indexer[indexer_pointer + p][0] = -1;
                  back_indexer[indexer_pointer + p][0] = -1;
                  //pos_gaps[indexer_pointer + p].reset();
                }
              q = MAX_L;
            }
        }

      //For all the previous q columns
      for(unsigned int q = 1; q < MAX_L; q++)
        {
          intersect(input[prev(input_pointer, input.size(), q)],
                    input[input_pointer],
                    indexer_pointer - q,
                    forw_indexer, back_indexer,
                    pos_gaps, num_pos_gaps, true);

          //If the just considered column q did not have any common element, we can fill as
          //  empty all the remaining columns and terminate.
          if(forw_indexer[indexer_pointer -  q][0] == -1)
            {
              for(unsigned int p = q + 1; p < MAX_L; p++)
                {
                  forw_indexer[indexer_pointer - p][0] = -1;
                  back_indexer[indexer_pointer - p][0] = -1;
                  //pos_gaps[indexer_pointer - p].reset();
                }
              q = MAX_L;
            }
        }


      //.:: Update prevision
      //XXX: Ricontrollare
      prevision_pointer = next(prevision_pointer, prevision.size(), 1);
      Pointer new_prevision_pointer = next(prevision_pointer, prevision.size(), MAX_L - 1);
      Pointer last_input_pointer = next(input_pointer, input.size(), MAX_L - 1);

      //Notice: prevision[new_prevision_pointer].size() == MAX_L
      for(unsigned int i = 1; i < prevision[new_prevision_pointer].size(); i++)
        {
          Pointer prec_input_pointer = prev(last_input_pointer, input.size(), i);
          unsigned int common_gaps = 0;
          Counter active_common = compute_active_common(input[prec_input_pointer], input[last_input_pointer], common_gaps);

          unsigned int combinations = 0;

          combinations = BinomialCoefficient::cumulative_binomial_coefficient(active_common - common_gaps, k_j[prec_input_pointer]) << common_gaps;

          // if(active_common != common_gaps) {
          //   combinations = binom_coeff::cumulative_binomial_coefficient(active_common - common_gaps, k_j[prec_input_pointer]) << common_gaps;
          // } else {
          //   //combinations = binom_coeff::cumulative_binomial_coefficient(common_gaps, common_gaps);
          //   combinations = 1 << common_gaps;
          // }

            // INC-K
            if(re_run_k && prevision[new_prevision_pointer][i].size() < combinations) {
                DEBUG("Resize new_prevision_pointer[" << new_prevision_pointer << "][" << i << "] to " << combinations);
                prevision[new_prevision_pointer][i].resize(combinations);
            }
	    if(re_run_k && i < backtrace_table1[step].size() && backtrace_table1[step][i].size() < combinations) {
                DEBUG("Resize backtrace tables[" << step << "][" << i << "] to " << combinations);
                backtrace_table1[step][i].resize(combinations);
                backtrace_table2_haplotypes[step][i].resize(combinations);
                backtrace_table2_new_block[step][i].resize(combinations);
            }

          fill(prevision[new_prevision_pointer][i].begin(),
               (prevision[new_prevision_pointer][i].begin() + combinations),
               Cost::INFTY);
        }


      //.:: Update OPT

      OPT_pointer = next(OPT_pointer, OPT.size(), 1);
      OPT[OPT_pointer] = Cost::INFTY;


      DEBUG(">> Update data structure completed");

      //>>>>>>>>>>>>>>>>>>>>> END UPDATE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //>>>>>>>>>>>>>>>>>>>>> ITERATIVE STEP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //Binary repesentation of the column and compute of coverage
      represent_column(input[input_pointer], colj, cov_j, gaps_mask, num_gaps);
      project(proj, colj, gaps_mask, cov_j);

      DEBUG("...| Column: " <<  column_to_string(colj, cov_j) << " -- current coverage: " << cov_j << " and current k: " << k_j[input_pointer]);
      DEBUG("...| #of gaps: " << num_gaps << "  and their positions: " << column_to_string(gaps_mask, cov_j));
      DEBUG("...| Column with gaps removed: " << column_to_string(proj, cov_j - num_gaps));

      //Initializing OPT[j] = infinite
      //XXX: Is it redundant??
      //OPT[OPT_pointer] = Cost::INFTY;

      //We have already computed k, homo_weights and homo_cost for current column

      //First option for the value of OPT[j]
      if(homo_cost[input_pointer] <= k_j[input_pointer])
        {
          //XXX: Can I remove this check?
          Cost temp = homo_weight[input_pointer] + OPT[prev(OPT_pointer, OPT.size(), 1)];
          if(temp < OPT[OPT_pointer]) {
            OPT[OPT_pointer] =  temp;
            solution_existence = true;
            is_homozygous[step] = true;
            DEBUG(".:: Column: " << step << " can be homozygous with a cost: " << OPT[OPT_pointer]);
          }
        }

      //Enumerate all the combinations

      bool loop_var;
      if(balancing_ and (cov_j > balancecov_)) {
        INFO("STEP " << step_global << " column has coverage : " << cov_j << " -- balancing with ratio : " << balanceratio_);
	balanced_generator.initialize(cov_j - num_gaps, k_j[input_pointer], proj, balanceratio_);
	loop_var = balanced_generator.has_next();
      }
      else {
	generator.initialize_cumulative(cov_j - num_gaps, k_j[input_pointer]);
	loop_var = generator.has_next();
      }

      while(loop_var)
        {
	  if(balancing_ and (cov_j > balancecov_)) {
	    balanced_generator.next();
	    balanced_generator.get_combination(comb_no_gaps);
	  }
	  else {
	    generator.next();
	    generator.get_combination(comb_no_gaps);
	  }

          TRACE("Combination of not gaps: " << column_to_string(comb_no_gaps, cov_j - num_gaps));

          Counter comb_gaps_int = 0;
          do {
            BitColumn comb_gaps(comb_gaps_int);

            TRACE("Combination of gaps: " << column_to_string(comb_gaps, num_gaps));

            make_mask(mask, gaps_mask, cov_j, comb_gaps, comb_no_gaps);

            TRACE("|--------");
            TRACE("|== Mask: " << column_to_string(mask, cov_j));

            //Initialize D[j, C'j] to infinite
            current_cost = Cost::INFTY;

            //Compute C'j
            //corrected_colj = colj ^ mask;
            TRACE("-->> corrected column: " << column_to_string(colj ^ mask, cov_j));

            //The column cannot be transformed into an homozygous column
            //if(corrected_colj.any() && (corrected_colj.count() != cov_j) )
            //{
            //Compute the weight of the mask
            Cost weight_mask = 0;

            if (unweighted_) {
              weight_mask = Cost((Cost::cost_t)mask.count());
            } else {
              compute_weight_mask(mask, input[input_pointer], weight_mask);
            }

            //Compute current_cost that corresponds to D[j, Bj]

            Counter q = 1;
            Pointer new_homo_pointer = prev(input_pointer, input.size(), q - 1);
            bool has_previous = true;
            Cost cumulative_homo = 0;

            do {
              //All the condition that have to be satisfied for the next column
              feasibility = (q - 1 == 0) || (homo_cost[new_homo_pointer] <= k_j[new_homo_pointer]);

              if (q >= MAX_L || forw_indexer[indexer_pointer - q][0] == -1 || !feasibility) {
                has_previous = false;
              } else {
                Counter active_qj(0);
                Cost temp(0);

                Pointer new_q_pointer = prev(input_pointer, input.size(), q);

                //First Mask
                extract_common_mask(input[new_q_pointer], indexer_pointer - q, input[input_pointer],
                                    mask, back_indexer, forw_indexer, mask_qj, active_qj);

                const Counter ungaps_q_corrected = (mask_qj&((pos_gaps[indexer_pointer - q]|constants.zeroes).flip())).count();

                if(ungaps_q_corrected <= k_j[new_q_pointer])
                  {
                    //Counter index = generator.cumulative_indexof(mask_qj, active_qj);
                    //Leave .count()
                    Counter index = compute_index_of(mask_qj, active_qj, num_pos_gaps[indexer_pointer - q],
                                                     pos_gaps[indexer_pointer - q]);
                    temp = prevision[prevision_pointer][q][index] + weight_mask + cumulative_homo;
                    if(temp < current_cost) {
                      current_cost = temp;
                      solution_existence = true;

                      temp_jump = q;
                      temp_index = index;
                      temp_haplotypes = backtrace_table2_haplotypes[step - q][q][index];
                      temp_new_block = false;
                    }
                    TRACE("-->> Temporary current cost: " << current_cost);
                    TRACE("---->> the previous equal heterozigous is " << (step - q)
                          << "  -- its mask: " << column_to_string(mask_qj, active_qj));
                  }

                //Complement
                complement_mask(mask_qj, active_qj, constants);

                if((active_qj - num_pos_gaps[indexer_pointer - q] - ungaps_q_corrected) <= k_j[new_q_pointer])
                  {
                    //Counter index = generator.cumulative_indexof(mask_qj, active_qj);
                    Counter index = compute_index_of(mask_qj, active_qj, num_pos_gaps[indexer_pointer - q],
                                                     pos_gaps[indexer_pointer - q]);
                    temp = prevision[prevision_pointer][q][index] + weight_mask + cumulative_homo;
                    if(temp < current_cost) {
                      current_cost = temp;
                      solution_existence = true;

                      temp_jump = q;
                      temp_index = index;
                      temp_haplotypes = !backtrace_table2_haplotypes[step - q][q][index];
                      temp_new_block = false;
                    }
                    TRACE("-->> Temporary current cost: " << current_cost);
                    TRACE("---->> the previous equal heterozigous is " << (step - q)
                          << "  -- its mask: " << column_to_string(mask_qj, active_qj));
                  }

                q++;

                new_homo_pointer = prev(input_pointer, input.size(), q - 1);

                cumulative_homo += homo_weight[new_homo_pointer];
                TRACE("----> Cumulative homo: " << cumulative_homo << "  with q:  " << q);
              }
            } while(has_previous);

            TRACE("-->> Best current cost (D[j, C'j]): "<< current_cost);

            //Third case of the recursion for D[j, C'j]
            //XXX: Check carefully!
            if(q <= MAX_L && feasibility) {
              Cost temp = OPT[prev(OPT_pointer, OPT.size(), q)] + weight_mask + cumulative_homo;
              if(temp < current_cost) {
                current_cost = temp;
                solution_existence = true;

                temp_jump = q;
                temp_index = 0;
                temp_haplotypes = false;
                temp_new_block = true;
                TRACE("<<>> Third case of recursion - First heterozigous of new block");
                TRACE("..OPT[previous] = " << OPT[prev(OPT_pointer, OPT.size(), q)] << " - weight:  " << weight_mask << " - cumulative_homo: " << cumulative_homo);
                //TRACE(".:: Column: " << step << " can be heterozigous with a cost: " << current_cost);
                //TRACE("====> Best correction:  " << column_to_string(mask, cov_j));
              }
            }

            //Make a prevision for all the seccessive column
            has_successive = true;
            Counter p = 1;

            //All the condition that have to be satisfied for the next column

            do {
              Pointer new_homo_pointer = next(input_pointer, input.size(), p - 1);
              feasibility = (p - 1 == 0) || (homo_cost[new_homo_pointer] <= k_j[new_homo_pointer]);

              if (p >= MAX_L || forw_indexer[indexer_pointer + p][0] == -1 || !feasibility) {
                has_successive = false;
              } else {
                //The number of elements shared between p and j
                Counter active_pj = 0;

                cut(mask, cut_mask, back_indexer[indexer_pointer + p], active_pj);
                TRACE("-->> Successive column: " << (step + p)
                      << " -- Prevision cost: " << current_cost
                      << " -- Common elements:  " << active_pj << " -- Cut mask: " << cut_mask
                      << "---" << column_to_string(cut_mask, active_pj));

                //Counter index = generator.cumulative_indexof(cut_mask, active_pj);
                Counter index = compute_index_of(cut_mask, active_pj, num_pos_gaps[indexer_pointer + p],
                                                 pos_gaps[indexer_pointer + p]);
                Pointer new_prevision_pointer = next(prevision_pointer, prevision.size(), p);
                Cost& temp = prevision[new_prevision_pointer][p][index];
                if(current_cost < temp) {
                    temp = current_cost;
                    backtrace_table1[step][p][index].jump = temp_jump;
                    backtrace_table1[step][p][index].index = temp_index;
                    backtrace_table2_haplotypes[step][p][index] = temp_haplotypes;
                    backtrace_table2_new_block[step][p][index] = temp_new_block;
                }
                TRACE("USCITO:  ");
                p++;
              }
            } while(has_successive);


            if(current_cost < current_best) {
              current_best = current_cost;

              best_heterozygous1[step].jump = temp_jump;
              best_heterozygous1[step].index = temp_index;
              best_heterozygous2_haplotypes[step] = temp_haplotypes;
              best_heterozygous2_new_block[step] = temp_new_block;
            }


            //Update value of OPT for the current column
            //XXX: Try <= to advantage heterozigosity
            if(current_cost < OPT[OPT_pointer]) {
              OPT[OPT_pointer] = current_cost;
              is_homozygous[step] = false;
              DEBUG(".:: Column: " << step << " can be heterozigous with a cost: " << OPT[OPT_pointer]);
              DEBUG("====> Best correction:  " << column_to_string(mask, cov_j));
            }
            TRACE("-->> OPT: " << OPT[OPT_pointer]);
            //}
            ++comb_gaps_int;
          } while (comb_gaps_int < (unsigned int)(1 << num_gaps));

	  if(balancing_ and (cov_j > balancecov_)) {
	    loop_var = balanced_generator.has_next();
	  }
	  else {
	    loop_var = generator.has_next();
	  }
        }

      if (step_global % 500 == 0) {
        INFO(".:: Step: " << step_global << "  ==>  OPT: " << OPT[OPT_pointer] + OPT_global);
      } else {
        DEBUG(".:: Step: " << step_global << "  ==>  OPT: " << OPT[OPT_pointer]);
      }
      //End of DP cycle for all the columns

      // INC-K
      if(solution_existence) {
          re_run_k = false;
      } else {
	  int old_k = k_j[input_pointer];
	  if(k_j_inc <= 0) {
	      k_j_inc = 1;
	  } else {
	      k_j_inc = k_j_inc + log2(k_j_inc) + 1;
	  }
          k_j[input_pointer] = floor(k_j_inc);
          INFO("STEP " << step_global << " INCREMENT k from " << old_k << " to " << k_j[input_pointer]);
          input_pointer = prev(input_pointer, input.size(), 1);
          prevision_pointer = prev(prevision_pointer, prevision.size(), 1);
          OPT_pointer = prev(OPT_pointer, OPT.size(), 1);
          step--;
          step_global--;
          re_run_k = true;
      }
    }

	
  if(solution_existence) {
    DEBUG("*** SUCCESS ***");
    DEBUG("> Optimal block cost:  " << OPT[OPT_pointer]);
    OPT_global += OPT[OPT_pointer];
    DEBUG("===> Optimal global cost:  " << OPT_global);

    reconstruct_haplotypes(backtrace_table1, backtrace_table2_haplotypes, backtrace_table2_new_block,
                           is_homozygous, homo_haplotypes,
                           best_heterozygous1, best_heterozygous2_haplotypes, best_heterozygous2_new_block,
                           haplotype1, haplotype2);

  } else {
    INFO("*** NO SOLUTION ***");
    INFO("<<>> No feasible solution exists with these parameters -- alpha = " << alpha_ << " and error rate = " << errorrate_);
    INFO("<<>> The last not feasible column is:  " << step << "  with coverage = " << cov_j << " and k = " << k_j[input_pointer]);
    exit(EXIT_FAILURE);
  }
}


//Change such that we do not take all the input in memory!!
void computeInputParams(Counter &num_cols, Counter &MAX_COV, Counter &MAX_L,
                        Counter &MAX_K, Counter &MAX_GAPS, vector<Counter> &sum_successive_L,
                        vector<vector<Counter> > &scheme_backtrace,
                        HapChatColumnIterator hap)
{
 
	hap.reset();
  num_cols = hap.column_count() + 1; //We add a starting dummy empty column

  vector<Column> input(num_cols);
  vector<unsigned int> homo_cost(num_cols);
  vector<unsigned int> k_j(num_cols);
  vector<Column>::iterator input_iterator(input.begin());
  vector<Counter> rows(num_cols * MAX_COVERAGE, 0);
  Column prova;
  do {
      Counter count_major = 0;
      Counter count_minor = 0;
      Counter count_gaps = 0;
      Column &current_column = *input_iterator;
      Column read_column;

      //XXX: Use current_column = column_read.get_next()
      if(input_iterator == input.begin()) {
        read_column.clear();
      } else {
        read_column=hap.get_column();
      }
      current_column.resize(read_column.size(), Entry(-1, Entry::BLANK, 0));

      for(unsigned int i = 0; i < read_column.size(); ++i) {
        current_column[i].set_read_id(read_column[i].get_read_id());
        current_column[i].set_allele_type(read_column[i].get_allele_type());
        current_column[i].set_phred_score(read_column[i].get_phred_score());
        //current_column[i].set_gap(read_column[i].is_gap());
				
        if(!(current_column[i].get_allele_type()==Entry::BLANK)) {
          if(current_column[i].get_allele_type() == Entry::REF_ALLELE) {
            ++count_major;
          } else if (current_column[i].get_allele_type() == Entry::ALT_ALLELE) {
            ++count_minor;
          } else {
            cerr << "ERROR: read invalid entry of type: " << current_column[i].get_allele_type() << endl;
            exit(EXIT_FAILURE);
          }
        } else {
          ++count_gaps;
        }

        ++rows[current_column[i].get_read_id()];
      }

      //sufficient condition to check the feasibility for the homozygous transformation
      homo_cost[input_iterator - input.begin()] = std::min(count_major, count_minor);
      homo_cost[input_iterator - input.begin()] = MAX_COVERAGE + 1;

      k_j[input_iterator - input.begin()] = computeK(count_minor + count_major);

      MAX_COV = std::max(static_cast<Counter>((*input_iterator).size()), MAX_COV);
      MAX_K = std::max(k_j[input_iterator - input.begin()], MAX_K);
      MAX_GAPS = std::max(count_gaps, MAX_GAPS);

      ++input_iterator;
  } while((input_iterator != input.end()) & hap.has_next());
  MAX_L = *max_element(rows.begin(), rows.end());
  MAX_L = std::max(MAX_L, static_cast<Counter>(2));

  Counter MAX_CONS_HOMO = 0;    //The maximum number of consecutive homozigous columns

  sum_successive_L.resize(MAX_L, 0);
  scheme_backtrace.clear();
  scheme_backtrace.resize(num_cols);
  for(unsigned int i = 0; i < input.size(); i++)
    {
      unsigned int y = 1;
      unsigned int k_temp = k_j[i];
      unsigned int current_cons_homo = 0;   //The maximum number of consecutive homozigous columns assuming i the first
      bool flag = true;
      scheme_backtrace[i].push_back(0);

      while(y < MAX_L && (i + y) < input.size())
        {
          unsigned int common_gaps = 0;
          Counter active_common = compute_active_common(input[i], input[i + y], common_gaps);

          //XXX: Add MAX_COMB_K and MAX_COMB_GAPS??

          Counter result = 0;
          result = BinomialCoefficient::cumulative_binomial_coefficient(active_common - common_gaps, k_temp) << common_gaps;

          // if(active_common != common_gaps) {
          //   result = binom_coeff::cumulative_binomial_coefficient(active_common - common_gaps, k_temp) << common_gaps;
          // } else {
          //   result = binom_coeff::cumulative_binomial_coefficient(common_gaps, common_gaps);
          // }

          sum_successive_L[y] = max(sum_successive_L[y], result);

          if(flag) {
            //XXX: Can I add && active_common != 0?
            if( (homo_cost[i + y] <= k_j[i + y]) && active_common != 0) {
              ++current_cons_homo;
              scheme_backtrace[i].push_back(result);
            } else {
              flag = false;
              scheme_backtrace[i].push_back(result);
            }
          }

          y++;
        }

      MAX_CONS_HOMO = std::max(MAX_CONS_HOMO, current_cons_homo);
    }

  //+1 is necessary to count the first heterozygous column before the longest sequence of homozugouses
  //and another +1 to count the heterozygous column after that

  //min it is necessary since MAX_CONS_HOMO + 2 can be larger than the real MAX_L
  MAX_L = std::min(MAX_L, MAX_CONS_HOMO + 2);
}


//XXX: Can I leave parameter q and pass as parameter just one column of forw_indexer and back_indexer??????????
void intersect(const Column &colQ, const Column &colJ, const Pointer &q,
               vector<vector<Pointer> > &forw_indexer, vector<vector<Pointer> > &back_indexer,
               vector<BitColumn> &pos_gaps, vector<Counter> &num_pos_gaps, const bool &q_is_back)
{
  size_t i = 0;
  size_t j = 0;
  size_t count = 0;
  pos_gaps[q].reset();
  num_pos_gaps[q] = 0;

  while ((i < colQ.size()) &&
         (j < colJ.size()) &&
         (colJ[j].get_read_id() != -1) &&
         (colQ[i].get_read_id() != -1)) {
    if (colQ[i].get_read_id() == colJ[j].get_read_id()) {
      forw_indexer[q][count] = i;
      back_indexer[q][count] = j;

      // pos_gaps[q].set(count, (q_is_back && colQ[i].is_gap()) || (!q_is_back && colJ[j].is_gap()));
      // num_pos_gaps[q] += pos_gaps[q].test(count);

      if(q_is_back) {
        if(colQ[i].get_allele_type()==Entry::BLANK) {
          pos_gaps[q].set(count);
          ++num_pos_gaps[q];
        }
      } else {
        if(colJ[j].get_allele_type()==Entry::BLANK) {
          pos_gaps[q].set(count);
          ++num_pos_gaps[q];
        }
      }

      ++i;
      ++j;
      ++count;
    } else if (colQ[i].get_read_id() < colJ[j].get_read_id()) {
      ++i;
    } else {
      ++j;
    }
  }

  if(count < forw_indexer[q].size()) {
    forw_indexer[q][count] = -1;
    back_indexer[q][count] = -1;
  }
}


void represent_column(const Column &column, BitColumn &result, Counter &cov,
                      BitColumn &gaps_mask, Counter &num_gaps)
{
  result.reset();
  gaps_mask.reset();
  cov = 0;
  num_gaps = 0;
  while(cov < column.size() && column[cov].get_read_id() != -1) {
    if(column[cov].get_allele_type()==Entry::BLANK) {
      gaps_mask.set(cov);
      ++num_gaps;
    }
    result.set(cov, (column[cov].get_allele_type() == Entry::ALT_ALLELE));
    ++cov;
  }
}


void make_mask(BitColumn &mask, const BitColumn &mask_gaps, const unsigned int &cov,
               const BitColumn &comb_gaps, const BitColumn &comb_no_gaps)
{
  mask.reset();
  unsigned int i_no_gaps = 0;
  unsigned int i_gaps = 0;

  for(unsigned int i = 0; i < cov; ++i) {
    // mask.set(i, (mask_gaps.test(i) && comb_gaps.test(i_gaps)) || (!mask_gaps.test(i) && comb_no_gaps.test(i_no_gaps)) );
    // i_gaps += mask_gaps.test(i);
    // i_no_gaps += !mask_gaps.test(i);
    if(mask_gaps[i]) {
      mask.set(i, comb_gaps[i_gaps++]);
    } else {
      mask.set(i, comb_no_gaps[i_no_gaps++]);
    }
  }
}


/*
  given some column col, i.e., 1010, project out those positions in
  mask_out, i.e., 0110, resulting in projection 10 (note that coverage
  cov is 4 here)

  somewhat the inverse of 'make_mask' above
*/
void project(BitColumn & projection, const BitColumn & col, const BitColumn & mask_out, const unsigned int & cov) {

  projection.reset();
  unsigned int j = 0;

  for(unsigned int i = 0; i < cov; ++i) {

    if(!mask_out[i]) {
      projection.set(j, col[i]);
      ++j;
    }
  }
}


unsigned int compute_index_of(const BitColumn &mask, const unsigned int &cov, const unsigned int &num_gaps,
                              const BitColumn &pos_gaps)
{
  BitColumn comb_gaps;
  BitColumn comb_no_gaps;

  unsigned int i_gaps = 0;
  unsigned int i_no_gaps = 0;

  for(unsigned int i = 0; i < cov; ++i) {
    // comb_gaps.set(i_gaps, pos_gaps.test(i) && mask.test(i));
    // comb_no_gaps.set(i_no_gaps, !pos_gaps.test(i) && mask.test(i));
    // i_gaps += pos_gaps.test(i);
    // i_no_gaps += !pos_gaps.test(i);
    if(pos_gaps[i]) {
      comb_gaps.set(i_gaps++, mask[i]);
    } else {
      comb_no_gaps.set(i_no_gaps++, mask[i]);
    }
  }

  return (BinomialCoefficient::cumulative_indexof(comb_no_gaps, cov - num_gaps) << num_gaps) | ((unsigned int) comb_gaps.to_ulong());
  //  return generator.cumulative_indexof(comb_no_gaps, cov - num_gaps) +
  //  ((unsigned int) comb_gaps.to_ulong()) * binom_coeff::binomial_coefficient(cov - num_gaps, k);
}


void cut(const BitColumn &in_col, BitColumn &cut_mask, const vector<Pointer> &indexer, Counter &active_pj)
{
  active_pj = 0;
  cut_mask.reset();

  const size_t isize= indexer.size();
  const Pointer* iactive= &indexer[0];
  while(active_pj < isize && *iactive != -1) {
    cut_mask.set(active_pj, in_col[*iactive]);
    ++active_pj;
    ++iactive;
  }
}


void extract_common_mask(const Column &column_q, const Pointer &q_pointer,
                         const Column &column_j, const BitColumn &mask_colj,
                         const vector<vector<Pointer> > &back_indexer,
                         const vector<vector<Pointer> > &forw_indexer,
                         BitColumn &mask_qj, Counter &active_qj)
{
  mask_qj.reset();
  active_qj = 0;
  //ungaps_q_corrected = 0;

  const Pointer* forw_indexer_q= &forw_indexer[q_pointer][0];
  const Pointer* back_indexer_q= &back_indexer[q_pointer][0];
  const size_t back_indexer_size= back_indexer[q_pointer].size();
  while (active_qj < back_indexer_size && *back_indexer_q != -1) {
    //XXX: leave if
    if ((column_q[*forw_indexer_q].get_allele_type() !=
         column_j[*back_indexer_q].get_allele_type())
        != mask_colj[*back_indexer_q]) {
      mask_qj.set(active_qj, 1);
    }
    //ungaps_q_corrected += !column_q[*forw_indexer_q].is_gap() && mask_qj.test(active_qj);
    // mask_qj.set(active_qj, ((column_q[*forw_indexer_q].get_allele_type() !=
    //                          column_j[*back_indexer_q].get_allele_type())
    //                         != mask_colj.test(*back_indexer_q)));
    ++active_qj;
    ++forw_indexer_q;
    ++back_indexer_q;
  }
}


int compute_active_common(const Column &colJ, const Column &colQ, unsigned int &common_gaps)
{
  int active_common = 0;
  common_gaps = 0;
  unsigned int i = 0;
  unsigned int j = 0;

  while ((i < colQ.size()) &&
         (j < colJ.size()) &&
         (colJ[j].get_read_id() != -1) &&
         (colQ[i].get_read_id() != -1)) {
    if (colQ[i].get_read_id() == colJ[j].get_read_id()) {
      if(colJ[j].get_allele_type()==Entry::BLANK) {
        ++common_gaps;
      }
      ++i;
      ++j;
      ++active_common;
    } else if (colQ[i].get_read_id() < colJ[j].get_read_id()) {
      ++i;
    } else {
      ++j;
    }
  }
  return active_common;
}


void insert_col_and_update(vector<Column> &input, vector<Counter> &k_j, vector <Counter> &homo_cost,
                           vector<Cost> &homo_weight, const Pointer &pointer, const Column &column,
                           vector<bool> &kind_homozygous, const Counter &step)
{
  Counter count_major = 0;
  Cost weight_major = 0;
  Counter count_minor = 0;
  Cost weight_minor = 0;
  unsigned int i = 0;
  bool gap;
  for(i = 0; i < column.size(); i++)
    {
    	gap=false;
      const long int column_read_id = column[i].get_read_id();
      const Entry::allele_t column_allele_type = column[i].get_allele_type();
      const unsigned int column_phred_score = (unweighted_)? 1 : column[i].get_phred_score();
			
      input[pointer][i].set_read_id(column_read_id);
      input[pointer][i].set_allele_type(column_allele_type);
      input[pointer][i].set_phred_score(column_phred_score);
      if(column[i].get_allele_type()==Entry::BLANK){
			gap=true;
			input[pointer][i].set_allele_type(Entry::REF_ALLELE);
			}
      //input[pointer][i].set_gap(column[i].is_gap());
      if(!gap) {
        if(column_allele_type == Entry::ALT_ALLELE) {
          ++count_minor;
          weight_minor += column_phred_score;
          //XXX: togliere questo if
        } else if (column_allele_type == Entry::REF_ALLELE) {
          ++count_major;
          weight_major += column_phred_score;
        }
      }
    }
		//cout<<count_major<<"||"<<count_minor;
  if(i < input[pointer].size() && input[pointer][i].get_read_id() != -1)
    {
      input[pointer][i].set_read_id(-1);
      input[pointer][i].set_allele_type(Entry::BLANK);
      input[pointer][i].set_phred_score(0);
    }

  //.:: Update k_j for current column

  k_j[pointer] = computeK(count_minor + count_major);


  //.:: Update homozygous cost

  homo_cost[pointer] = MAX_COUNTER;
  homo_weight[pointer] = Cost::INFTY;

  if(count_minor <= k_j[pointer] && weight_minor < homo_weight[pointer]) {
    homo_cost[pointer] = count_minor;
    homo_weight[pointer] = weight_minor;
    if(step < kind_homozygous.size()) {
      kind_homozygous[step] = true;
    }
  }

  if(count_major <= k_j[pointer] && weight_major < homo_weight[pointer]) {
    homo_cost[pointer] = count_major;
    homo_weight[pointer] = weight_major;
    if(step < kind_homozygous.size()) {
      kind_homozygous[step] = false;
    }
  }

  homo_cost[pointer] = MAX_COVERAGE + 1;

  //Add for the all-heterozygous assumption
  //homo_weight[pointer] =  homo_weight[pointer] + 100000;column.size()==1 ? 0 : 1000000;
}


void compute_weight_mask(const BitColumn &mask, const Column &column, Cost &weight_mask) {
  BitColumn comb(mask);
  weight_mask = 0;
  int pos = 0;
  //Column::const_iterator i_column = column.begin();
  int temp = 0;

  while(comb.any())
    {
      temp = ffsl(comb.to_ulong());
      pos += temp;
      //i_column += temp - 1;
      weight_mask += column[pos - 1].get_phred_score();
      //weight_mask += *i_column.get_phred_score();
      comb>>=(temp);
      //++i_column;
    }
}


void reconstruct_haplotypes(const vector<vector<vector<Backtrace1> > > &backtrace_table1,
                            const vector<vector<vector<bool> > > &backtrace_table2_haplotypes,
                            const vector<vector<vector<bool> > > &backtrace_table2_new_block,
                            const vector<bool> &is_homozygous,
                            const vector<bool> &homo_haplotypes,
                            const vector<Backtrace1> &best_heterozygous1,
                            const vector<bool> &best_heterozygous2_haplotypes,
                            const vector<bool> &best_heterozygous2_new_block,
                            vector<bool> &haplotype1, vector<bool> &haplotype2) {

  Counter col = backtrace_table1.size() - 1;
  haplotype1.resize(col);
  haplotype2.resize(col);
	
  while(col > 0) {
    while(is_homozygous[col]) {
      if(homo_haplotypes[col]) {
        haplotype1[col - 1] = false;
        haplotype2[col - 1] = false;
      } else {
        haplotype1[col - 1] = true;
        haplotype2[col - 1] = true;
      }

      --col;
    }
	
    Backtrace1 back1 = best_heterozygous1[col];
    
    bool back2_haplotypes = best_heterozygous2_haplotypes[col];
    bool back2_new_block = best_heterozygous2_new_block[col];
    bool flag = col > 0;
    while (flag) {
      if(back2_haplotypes) {
        haplotype1[col - 1] = false;
        haplotype2[col - 1] = true;
      } else {
        haplotype1[col - 1] = true;
        haplotype2[col - 1] = false;
      }

      for(int i = 0; i < (back1.jump - 1); i++) {
        --col;
        if(homo_haplotypes[col]) {
          haplotype1[col - 1] = false;
          haplotype2[col - 1] = false;
        } else {
          haplotype1[col - 1] = true;
          haplotype2[col - 1] = true;
        }
      }

      --col;

      if(back2_new_block || col == 0) {
        flag = false;
      } else {
        flag = true;
        back2_haplotypes = backtrace_table2_haplotypes[col][back1.jump][back1.index];
        back2_new_block = backtrace_table2_new_block[col][back1.jump][back1.index];
        back1 = backtrace_table1[col][back1.jump][back1.index];
      }

    }
  }
}


Counter computeK(const Counter &cov, const double &alpha=0.0, const double &error_rate=0.0)
{
  static bool computed = false;
  static vector<Counter> ks(cov + 1, 0);

  if (!computed) {
    for(Counter i = 1; i < ks.size(); ++i) {
      Counter k = 0;

      double cumulative =  pow(1.0 - error_rate, i);

      while(!(1.0 - cumulative <= alpha) && (k < i)) {
        ++k;
        cumulative += (double)BinomialCoefficient::binomial_coefficient(i, k) * pow(error_rate, k) * pow(1.0 - error_rate, i - k);
      }

      ks[i] = k;
    }
    computed = true;
  }

  return ks[cov];
}


int map_fragment(const vector<char> &read, const vector<unsigned int> &weights, const unsigned int offset,
                 const vector<bool> &haplo1_in, const vector<bool> &haplo2_in, unsigned int &total_errors)
{
  unsigned int distance1 = 0;
  unsigned int distance2 = 0;

  for(unsigned int col = 0; col < read.size(); ++col) {
    if(haplo1_in[col + offset]) {
      if(read[col] != '1') {
        distance1 += weights[col];
      }
    } else {
      if(read[col] != '0') {
        distance1 += weights[col];
      }
    }

    if(haplo2_in[col + offset]) {
      if(read[col] != '1') {
        distance2 += weights[col];
      }
    } else {
      if(read[col] != '0') {
        distance2 += weights[col];
      }
    }
  }

  if(distance1 <= distance2) {
    total_errors += distance1;
    return 1;
  } else {
    total_errors += distance2;
    return 2;
  }
}


void make_super_reads(vector <unsigned int> positions, vector <char> haplo1, vector<char> haplo2){

  int allele;
  size_t k = superreads_.size() - 1; // TODO: make this a bit nicer
  for(unsigned int i=0;i<positions.size();i++){
    if(haplo1[i]=='1') allele=1;
    else allele=0;
    superreads_[k].first->addVariant(positions[i],allele,30);
    if(haplo2[i]=='1') allele=1;
    else allele=0;
    superreads_[k].second->addVariant(positions[i],allele,30);
  }
}

};
