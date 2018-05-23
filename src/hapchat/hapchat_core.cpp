/**
 *
 * HapCHAT: Adaptive haplotype assembly for efficiently leveraging
 * high coverage in long read HapCHAT
 *
 * Copyright (C) 2017  Stefano Beretta, Murray Patterson and Simone Zaccaria
 *
 * Distributed under the terms of the GNU General Public License (GPL)
 *
 * This file is part of HapCHAT.
 *
 * HapCHAT is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * HapCHAT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HapCHAT.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#include <stdlib.h>
#include <string.h>
#include <limits>
#include <math.h>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <stdexcept>


#include "basic_types.h"
#include "binomial.h"
#include "combinations.h"
#include "balanced_combinations.h"
#ifdef LOAD_REVISION
#include "revision.h"
#endif

// Log messages with DEBUG priority and higher
#define LOG_MSG
#define LOG_THRESHOLD LOG_LEVEL_INFO
// Include log facilities. It should the last include!!
#include "log.cpp"
#include "fondamental.cpp"


class hapchatcore{
private:
Fondamental fonda;
options_t options;
public:	
	hapchatcore(ReadSet* read_set){
	options= whatshap_options(read_set);
	runCore();
	};


int runCore()
{
#if defined(VCS_DATE) && defined(VCS_SHORT_HASH) && defined(VCS_WC_MODIFIED)
  INFO("HapCHAT (" VCS_BRANCH "@" VCS_SHORT_HASH << (VCS_WC_MODIFIED ? "-dirty" : "-clean") << ")");
#else
  INFO("HapCHAT");
#endif
  INFO("Starting...");

  //Counter threshold_coverage = 30;

  const constants_t constants;

  
  INFO("Arguments:");
  INFO("Initialized? " << (options.options_initialized?"True":"False"));
  INFO("Input filename: '" << options.input_filename << '\'');
  INFO("Haplotype filename: '" << options.haplotype_filename << '\'');
  INFO("Discard weights? " << (options.unweighted?"True":"False"));
  INFO("Do not add X's? " << (options.no_xs?"True":"False"));
  INFO("All-heterozygous assumption? " << (options.all_heterozygous?"True":"False"));
  INFO("Input as unique block? " << (options.unique?"True":"False"));
  INFO("Error rate: " << options.error_rate);
  INFO("Alpha: " << options.alpha);
  INFO("Balancing? " << (options.balancing?"True":"False"));
  INFO("Balancing activates after coverage: " << options.balance_cov);
  INFO("Balance ratio: " << options.balance_ratio);
	
  if (!options.options_initialized) {
    FATAL("Arguments not correctly initialized! Exiting..");
    exit(EXIT_FAILURE);
  }
  //Readset section
	HapCHATcore hap=HapCHATcore(options.readset,options.unique);
  //Initializing the starting parameters: no competitive section
	
  //Pre-compute binomial values
  binom_coeff::initialize_binomial_coefficients(MAX_COVERAGE, MAX_COVERAGE);
  fonda.computeK(MAX_COVERAGE, options.alpha, options.error_rate);

  Counter step = 0;
  Cost OPT = 0;
  Counter counter_block = 0;
  Counter counter_columns = 0;
  Counter counter_inhomo = 0;
  //BlockReader blockreader(options.input_filename, MAX_COVERAGE, options.unweighted, options.unique);

  Counter MAX_COV = 0;
  Counter MAX_L = 0;
  Counter MAX_K = 0;
  Counter MAX_GAPS = 0;
  Counter XS1 = 0;
  Counter XS2 = 0;
  Counter TOTAL_MISMATCHES = 0;

  vector<vector<char> > haplotype_blocks1;
  vector<vector<char> > haplotype_blocks2;
	hap.reset();
  while(hap.hasNextBlock()) {
  	
    //Block block = blockreader.get_block();
    DEBUG("BLOCK: "<< counter_block);

    //ColumnReader1 columnreader_jump(block, !options.all_heterozygous);

    vector<bool> haplotype1(hap.columnCount());
    vector<bool> haplotype2(hap.columnCount());

    if(hap.columnCount()> 0) {
      fonda.dp(constants, options, haplotype1, haplotype2, step, OPT,
         MAX_COV, MAX_L, MAX_K, MAX_GAPS, counter_block++,hap);
    } else {
      DEBUG("jumped");
      ++counter_block;
    }

    //ColumnReader1 columnreader_nojump(block, false);

    counter_columns += hap.columnCount();
   // counter_inhomo += (columnreader_nojump.num_cols() - columnreader_jump.num_cols());

    if(!options.no_xs) {
      vector<bool> filled_haplo1(hap.columnCount());
      vector<bool> filled_haplo2(hap.columnCount());

      DEBUG("Starting fill");

      fonda.fill_haplotypes( haplotype1, haplotype2, filled_haplo1, filled_haplo2, options,hap);

      DEBUG("Filled haplotypes");

      vector<char> xs_haplotype1(hap.columnCount());
      vector<char> xs_haplotype2(hap.columnCount());

      fonda.add_xs(filled_haplo1, filled_haplo2, xs_haplotype1, xs_haplotype2, options,
             XS1, XS2, TOTAL_MISMATCHES,hap);

      DEBUG("Added X's");

      haplotype_blocks1.push_back(xs_haplotype1);
      haplotype_blocks2.push_back(xs_haplotype2);

    } else {
      vector<char> output_block1(hap.columnCount());
      vector<char> output_block2(hap.columnCount());

      fonda.fill_haplotypes(haplotype1, haplotype2, output_block1, output_block2, options,hap);

      DEBUG("Filled haplotypes");

      haplotype_blocks1.push_back(output_block1);
      haplotype_blocks2.push_back(output_block2);
    }

  }

  INFO("");

  INFO("OPTIMUM:  " << OPT);

  INFO("");

  INFO("MAX_COV:  " << MAX_COV);
  INFO("MAX_L:  " << MAX_L);
  INFO("MAX_K:  " << MAX_K);
  INFO("MAX_GAPS:  " << MAX_GAPS);
  INFO("# of blocks:  " << counter_block);
  INFO("# of columns:  "  << counter_columns);
  INFO("# of homozygous in input:  " << counter_inhomo);

  INFO("");

  INFO("X's INSERTED IN THE FIRST HAPLOTYPE:  " << XS1);
  INFO("X's INSERTED IN THE SECOND HAPLOTYPE:  " << XS2);
  INFO("TOTAL MISMATCHES:  " << TOTAL_MISMATCHES);

  DEBUG("<<>> Writing haplotypes...");
  ofstream ofs;
  try {
    ofs.open(options.haplotype_filename.c_str(), ios::out);
    fonda.write_haplotypes(haplotype_blocks1, haplotype_blocks2, ofs);
  } catch(exception & e) {
    ERROR("::::::: Error writing haplotype to \"" << options.haplotype_filename << "\": " << e.what());
    //write_haplotypes(haplotype_blocks1, haplotype_blocks2, cout);
    return EXIT_FAILURE;
  }
}
};




