#include <cassert>
#include <limits>
#include "columncostcomputer.h"

using namespace std;

ColumnCostComputer::ColumnCostComputer(const std::vector<const Entry*>& column, const std::vector<unsigned int>& read_marks, unsigned int inheritance_val, std::vector<std::vector<int>>& triples)
: column(column), read_marks(read_marks), inheritance_val(inheritance_val), partitioning(0), num_of_triples(triples.size()), triples(triples) {
  this->roots = compute_roots(triples);
  this->num_of_roots=roots.size();
  this->cost_partitions(num_of_roots);
  
  for (auto& triple : triples){
    unsigned int triple_inheritance = inheritance_val & 3;
    // mother and father
    for(int i = 0; i < 2; i++) {
      int m = triple[i];
      if (haps.find(m) == haps.end()) {
	int pos = -1;
	for(unsigned j = 0; j < roots.size(); j++){
	  if(roots[j] == m) {
	    pos = j;
	    break;
	  }
	}
	if (pos >= 0) {
	  haps.insert(m, std::make_pair<unsigned>(2*pos, 2*pos + 1));
	}
      }
    }
    //child
    int c = triple[2];
    int haps1 = haps[triple[0]][triple_inheritance / 2];
    int haps2 = haps[triple[1]][triple_inheritance % 2];
    haps.insert(c, std::make_pair<unsigned>(haps1, haps2);

    inheritance_val = (inheritance_val >> 2);
  }
}

   std::vector<unsigned int> compute_roots(std::vector<this->triple>& triples){
   std::set<unsigned int> parents;
   std::set<unsigned int> children;
   std::set<unsigned int> intersect;
   std::vector<unsigned int> roots;
   
   for (auto& triple : triples){
      parents.insert(triple[0]);
      parents.insert(triple[1]);
      children.insert(triple[2]);
   }
      std::set_intersection(parents.begin(), parents.end(), children.begin(), children.end(), std::inserter(intersect, intersect.end()));
      std::set_difference(parents.begin(), parents.end(), intersect.begin(), intersect.end(), std::inserter(roots, roots.begin()));
      return roots;
 }

void ColumnCostComputer::set_partitioning(unsigned int partitioning) {
  // compute cost from scratch
  unsigned int border = 2*this->num_of_roots-1;
  for(unsigned int i = 0; i <= border; i++) {
    cost_partition[i]=std::make_pair(0,0);
  }

  this->partitioning = partitioning;
  for (vector<const Entry*>::const_iterator it = column.begin(); it != column.end(); ++it) {
    auto& entry = **it;
    bool entry_in_partition1 = (partitioning & ((unsigned int)1)) == 0;
    unsigned int ind_id= read_marks[entry.get_read_id()];
   switch (entry.get_allele_type()) {
      case Entry::MAJOR_ALLELE:
        (entry_in_partition1?cost_partition[this->haps[ind_id].first]:cost_partition[this->haps[ind_id].second])[1] += entry.get_phred_score();
        break;
      case Entry::MINOR_ALLELE:
        (entry_in_partition1?cost_partition[this->haps[ind_id].first]:cost_partition[this->haps[ind_id].second])[0] += entry.get_phred_score();
        break;
      case Entry::BLANK:
        break;
      default:
        assert(false);
    }

   
    partitioning = partitioning >> 1;
  }
}

void ColumnCostComputer::update_partitioning(int bit_to_flip) {
  // update cost based on the changed bit
  const Entry& entry = *column[bit_to_flip];
  partitioning = partitioning ^ (((unsigned int)1) << bit_to_flip);
  bool entry_in_partition1 = (partitioning & (((unsigned int)1) << bit_to_flip)) == 0;
  unsigned int ind_id= read_marks[entry.get_read_id()];
  switch (entry.get_allele_type()) {
    case Entry::MAJOR_ALLELE:
      (entry_in_partition1?cost_partition[this->haps[ind_id].second]:cost_partition[this->haps[ind_id].first])[1] -= entry.get_phred_score();
      (entry_in_partition1?cost_partition[this->haps[ind_id].first]:cost_partition[this->haps[ind_id].second])[1] += entry.get_phred_score();
      break;
    case Entry::MINOR_ALLELE:
      (entry_in_partition1?cost_partition[this->haps[ind_id].second]:cost_partition[this->haps[ind_id].first])[0] -= entry.get_phred_score();
      (entry_in_partition1?cost_partition[this->haps[ind_id].first]:cost_partition[this->haps[ind_id].second])[0] += entry.get_phred_score();
      break;
    case Entry::BLANK:
      break;
    default:
      assert(false);
  }
}

unsigned int ColumnCostComputer::get_cost( std::vector<std::vector<int>>& genotypes) {
  unsigned int best_cost = numeric_limits<unsigned int>::max();
  // Enumerate all possible assignments of alleles to haplotypes and 
  // compute costs for those which are compatible with genotypes.
  // TODO: This can be done more efficiently.
  unsigned int border = std::pow(4,this->num_of_roots);
  std::vector<unsigned int> duplicate_roots = this->roots;
  
  for (unsigned int i=0; i<border; ++i) {
    std::map<unsigned, std::pair<unsigned int>> enumerate_haps;
    unsigned int counter = 0;
    for (unsigned int j = 0; j <= 2; j++){
      unsigned int triple_inheritance = inheritance_val & 3;
      for (unsigned int k = 0; k <= 2; k++){
	if (std::find(duplicate_roots.begin(), duplicate_roots.end(), triples[j][k]) != duplicate_roots.end()){
	  duplicate_roots.erase(std::remove(duplicate_roots.begin(), duplicate_roots.end(), triples[j][k]), duplicate_roots.end());
	  unsigned int allele1= ((i >> counter) & 1);
	  unsigned int allele2 = ((i >> (counter+1)) & 1);
	  enumerate_haps[triples[j][k]] = std::make_pair<unsigned int,unsigned int>(allele1, allele2);
	  counter++;
	  if (allele1+ allele2 != genotypes[j][k]) continue;
	} else {
		if (enumerate_haps.find(triples[j][k]) == enumerate_haps.end()) {
	        unsigned int allele_c1 = ((triple_inheritance & 1) == 0)?enumerate_haps[triples[j][0]].first:enumerate_haps[triples[j][0]].second;
		unsigned int allele_c2 = (((triple_inheritance>>1) & 1) == 0)?enumerate_haps[triples[j][1]].first:enumerate_haps[triples[j][1]].second;
		 enumerate_haps[triples[j][k]] = std::make_pair<unsigned int,unsigned int>(allele_c1, allele_c2);
		if (allele_c1 + allele_c2 != genotypes[j][k]) continue;
		}
	}
      } inheritance_val = (inheritance_val >> 2);
    } 
    unsigned int cost;
    
    for (unsigned int p=0;p< this->roots.size(); p++)
    {
      unsigned int temp1= enumerate_haps[this->roots[p]].first;
      unsigned int temp2= enumerate_haps[this->roots[p]].second;
      cost+= cost_partition[2*p][temp1];
       cost+= cost_partition[2*p+1][temp2];
    }
    if (cost < best_cost) {
      best_cost = cost;
    }
  }
  
  if (best_cost == numeric_limits<unsigned int>::max()) {
    throw std::runtime_error("Error: Mendelian conflict");
  }
  
  return best_cost;
}

ColumnCostComputer::trio_alleles_t ColumnCostComputer::get_alleles(std::vector<std::vector<int>>& genotypes) {
  // TODO: avoid code duplication
  unsigned int best_cost = numeric_limits<unsigned int>::max();
  unsigned int second_best_cost = numeric_limits<unsigned int>::max();
  trio_alleles_t result(Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES);
  std::map<unsigned int, ind_alleles_t> population_alleles;
  for (unsigned int i=0; i<16; ++i) {
    unsigned int allele_m1 = i & 1;
    unsigned int allele_m2 = (i >> 1) & 1;
    if (allele_m1 + allele_m2 != genotypem) continue;
    
    unsigned int allele_f1 = (i >> 2) & 1;
    unsigned int allele_f2 = (i >> 3) & 1;
    if (allele_f1 + allele_f2 != genotypef) continue;
    
    unsigned int allele_c1 = ((inheritance_val & 1) == 0)?allele_m1:allele_m2;
    unsigned int allele_c2 = (((inheritance_val>>1) & 1) == 0)?allele_f1:allele_f2;
    if (allele_c1 + allele_c2 != genotypec) continue;
    
    unsigned int cost = cost_partition_m1[allele_m1] + cost_partition_m2[allele_m2] + cost_partition_f1[allele_f1] + cost_partition_f2[allele_f2]; 
    if (cost <= best_cost) {
      second_best_cost = best_cost;
      best_cost = cost;
      result = trio_alleles_t(
        (allele_m1==0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE,
        (allele_m2==0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE,
        (allele_f1==0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE,
        (allele_f2==0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE,
        (allele_c1==0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE,
        (allele_c2==0)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE
      );
    }
  }
  
  if (best_cost == numeric_limits<unsigned int>::max()) {
    throw std::runtime_error("Error: Mendelian conflict");
  }
  
  if (second_best_cost == best_cost) {
    // TODO: Having too equal scores does not necissarily imply that the case is
    //       undecidable for all three individuals. Be smarter here.
    return trio_alleles_t(Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES);
  }
  
  return result;
}

unsigned int ColumnCostComputer::get_weight(bool second_haplotype) {
  throw std::runtime_error("Not yet implemented for trios.");
}
