#include <cassert>
#include <limits>
#include "columncostcomputer.h"

using namespace std;

ColumnCostComputer::ColumnCostComputer(const std::vector<const Entry*>& column, const std::vector<unsigned int>& read_marks, unsigned int inheritance_val)
: column(column), read_marks(read_marks), inheritance_val(inheritance_val) {
  cost_partition_m1[0] = 0;
  cost_partition_m1[1] = 0;
  cost_partition_m2[0] = 0;
  cost_partition_m2[1] = 0;
  cost_partition_f1[0] = 0;
  cost_partition_f1[1] = 0;
  cost_partition_f2[0] = 0;
  cost_partition_f2[1] = 0;
  cost_partition1[0] = 0;
  cost_partition1[1] = 0;
  cost_partition2[0] = 0;
  cost_partition2[1] = 0;
  partitioning = 0;
}
//used in superreads
void ColumnCostComputer::set_partitioning_m(unsigned int partitioning) {
  // compute cost from scratch
 // cout<<"i am in set"<<endl;
  cost_partition1[0] = 0;
  cost_partition1[1] = 0;
  cost_partition2[0] = 0;
  cost_partition2[1] = 0;
  this->partitioning = partitioning;
  for (vector<const Entry*>::const_iterator it = column.begin(); it != column.end(); ++it) {
    auto& entry = **it;
    if(read_marks[entry.get_read_id()]==1){
      bool entry_in_partition1 = (partitioning & ((unsigned int)1)) == 0;
     // cout<< "partitioning"<< partitioning<< "binary"<< (*it)->get_allele_type() << (*it)->get_phred_score();
      switch ((*it)->get_allele_type()) {
        case Entry::MAJOR_ALLELE:
          (entry_in_partition1?cost_partition1:cost_partition2)[1] += (*it)->get_phred_score();
        // cout<<"major***\t"<<cost_partition1[0]<<"\t"<<cost_partition2[0]<<"\t"<<cost_partition1[1]<<"\t"<<cost_partition2[1]<<endl;
          break;
        case Entry::MINOR_ALLELE:
          (entry_in_partition1?cost_partition1:cost_partition2)[0] += (*it)->get_phred_score();
        // cout<<"minor****\t"<<cost_partition1[0]<<"\t"<<cost_partition2[0]<<"\t"<<cost_partition1[1]<<"\t"<<cost_partition2[1]<<endl;
          break;
        case Entry::BLANK:
          break;
        default:
          assert(false);
    }
    }
    partitioning = partitioning >> 1;
  }
}
//used in superreads
void ColumnCostComputer::set_partitioning_f(unsigned int partitioning) {
  // compute cost from scratch
  cost_partition1[0] = 0;
  cost_partition1[1] = 0;
  cost_partition2[0] = 0;
  cost_partition2[1] = 0;
  this->partitioning = partitioning;
  for (vector<const Entry*>::const_iterator it = column.begin(); it != column.end(); ++it) {
    auto& entry = **it;
    bool entry_in_partition1 = (partitioning & ((unsigned int)1)) == 0;
    if(read_marks[entry.get_read_id()]==2){    
     // cout<< "partitioning"<< partitioning<< "binary"<< (*it)->get_allele_type() << (*it)->get_phred_score();
      switch ((*it)->get_allele_type()) {
        case Entry::MAJOR_ALLELE:
          (entry_in_partition1?cost_partition1:cost_partition2)[1] += (*it)->get_phred_score();
        // cout<<"major***\t"<<cost_partition1[0]<<"\t"<<cost_partition2[0]<<"\t"<<cost_partition1[1]<<"\t"<<cost_partition2[1]<<endl;
          break;
        case Entry::MINOR_ALLELE:
          (entry_in_partition1?cost_partition1:cost_partition2)[0] += (*it)->get_phred_score();
        // cout<<"minor****\t"<<cost_partition1[0]<<"\t"<<cost_partition2[0]<<"\t"<<cost_partition1[1]<<"\t"<<cost_partition2[1]<<endl;
          break;
        case Entry::BLANK:
          break;
        default:
          assert(false);
    }
    }
    partitioning = partitioning >> 1;
  }
}
//used in superreads
void ColumnCostComputer::set_partitioning_c(unsigned int partitioning) {
  // compute cost from scratch
  cost_partition1[0] = 0;
  cost_partition1[1] = 0;
  cost_partition2[0] = 0;
  cost_partition2[1] = 0;
  this->partitioning = partitioning;
  for (vector<const Entry*>::const_iterator it = column.begin(); it != column.end(); ++it) {
    auto& entry = **it;
    bool entry_in_partition1 = (partitioning & ((unsigned int)1)) == 0;
    if(read_marks[entry.get_read_id()]==0){    
     // cout<< "partitioning"<< partitioning<< "binary"<< (*it)->get_allele_type() << (*it)->get_phred_score();
      switch ((*it)->get_allele_type()) {
        case Entry::MAJOR_ALLELE:
          (entry_in_partition1?cost_partition1:cost_partition2)[1] += (*it)->get_phred_score();
        // cout<<"major***\t"<<cost_partition1[0]<<"\t"<<cost_partition2[0]<<"\t"<<cost_partition1[1]<<"\t"<<cost_partition2[1]<<endl;
          break;
        case Entry::MINOR_ALLELE:
          (entry_in_partition1?cost_partition1:cost_partition2)[0] += (*it)->get_phred_score();
        // cout<<"minor****\t"<<cost_partition1[0]<<"\t"<<cost_partition2[0]<<"\t"<<cost_partition1[1]<<"\t"<<cost_partition2[1]<<endl;
          break;
        case Entry::BLANK:
          break;
        default:
          assert(false);
    }
    }
    partitioning = partitioning >> 1;
  }
}

void ColumnCostComputer::set_partitioning(unsigned int partitioning) {
  // compute cost from scratch
  cost_partition_m1[0] = 0;
  cost_partition_m1[1] = 0;
  cost_partition_m2[0] = 0;
  cost_partition_m2[1] = 0;
  cost_partition_f1[0] = 0;
  cost_partition_f1[1] = 0;
  cost_partition_f2[0] = 0;
  cost_partition_f2[1] = 0;
  this->partitioning = partitioning;
  for (vector<const Entry*>::const_iterator it = column.begin(); it != column.end(); ++it) {
    auto& entry = **it;
    bool entry_in_partition1 = (partitioning & ((unsigned int)1)) == 0;
    switch(read_marks[entry.get_read_id()]) {
      case 0u: // child
        switch (entry.get_allele_type()) {
          case Entry::MAJOR_ALLELE:
            switch (inheritance_val) {
              case 0:
                (entry_in_partition1?cost_partition_m1:cost_partition_f1)[1] += entry.get_phred_score();
                break;
              case 1:
                (entry_in_partition1?cost_partition_m2:cost_partition_f1)[1] += entry.get_phred_score();
                break;
              case 2:
                (entry_in_partition1?cost_partition_m1:cost_partition_f2)[1] += entry.get_phred_score();
                break;
              case 3:
                (entry_in_partition1?cost_partition_m2:cost_partition_f2)[1] += entry.get_phred_score();
                break;
              default:
                assert(false);
            }
            break;
          case Entry::MINOR_ALLELE:
            switch (inheritance_val) {
              case 0:
                (entry_in_partition1?cost_partition_m1:cost_partition_f1)[0] += entry.get_phred_score();
                break;
              case 1:
                (entry_in_partition1?cost_partition_m2:cost_partition_f1)[0] += entry.get_phred_score();
                break;
              case 2:
                (entry_in_partition1?cost_partition_m1:cost_partition_f2)[0] += entry.get_phred_score();
                break;
              case 3:
                (entry_in_partition1?cost_partition_m2:cost_partition_f2)[0] += entry.get_phred_score();
                break;
              default:
                assert(false);
            }
            break;
          case Entry::BLANK:
            break;
          default:
            assert(false);
        }        
        break;
      case 1u: // mother
        switch (entry.get_allele_type()) {
          case Entry::MAJOR_ALLELE:
            (entry_in_partition1?cost_partition_m1:cost_partition_m2)[1] += entry.get_phred_score();
            break;
          case Entry::MINOR_ALLELE:
            (entry_in_partition1?cost_partition_m1:cost_partition_m2)[0] += entry.get_phred_score();
            break;
          case Entry::BLANK:
            break;
          default:
            assert(false);
        }
        break;
      case 2u: // father
        switch (entry.get_allele_type()) {
          case Entry::MAJOR_ALLELE:
            (entry_in_partition1?cost_partition_f1:cost_partition_f2)[1] += entry.get_phred_score();
            break;
          case Entry::MINOR_ALLELE:
            (entry_in_partition1?cost_partition_f1:cost_partition_f2)[0] += entry.get_phred_score();
            break;
          case Entry::BLANK:
            break;
          default:
            assert(false);
        }
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
  switch(read_marks[entry.get_read_id()]) {
    case 0u: // child
      switch (entry.get_allele_type()) {
        case Entry::MAJOR_ALLELE:
          switch (inheritance_val) {
            case 0:
              (!entry_in_partition1?cost_partition_m1:cost_partition_f1)[1] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m1:cost_partition_f1)[1] += entry.get_phred_score();
              break;
            case 1:
              (!entry_in_partition1?cost_partition_m2:cost_partition_f1)[1] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m2:cost_partition_f1)[1] += entry.get_phred_score();
              break;
            case 2:
              (!entry_in_partition1?cost_partition_m1:cost_partition_f2)[1] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m1:cost_partition_f2)[1] += entry.get_phred_score();
              break;
            case 3:
              (!entry_in_partition1?cost_partition_m2:cost_partition_f2)[1] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m2:cost_partition_f2)[1] += entry.get_phred_score();
              break;
            default:
              assert(false);
          }
          break;
        case Entry::MINOR_ALLELE:
          switch (inheritance_val) {
            case 0:
              (!entry_in_partition1?cost_partition_m1:cost_partition_f1)[0] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m1:cost_partition_f1)[0] += entry.get_phred_score();
              break;
            case 1:
              (!entry_in_partition1?cost_partition_m2:cost_partition_f1)[0] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m2:cost_partition_f1)[0] += entry.get_phred_score();
              break;
            case 2:
              (!entry_in_partition1?cost_partition_m1:cost_partition_f2)[0] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m1:cost_partition_f2)[0] += entry.get_phred_score();
              break;
            case 3:
              (!entry_in_partition1?cost_partition_m2:cost_partition_f2)[0] -= entry.get_phred_score();
              (entry_in_partition1?cost_partition_m2:cost_partition_f2)[0] += entry.get_phred_score();
              break;
            default:
              assert(false);
          }
          break;
        case Entry::BLANK:
          break;
        default:
          assert(false);
      }        
      break;
    case 1u: // mother
      switch (entry.get_allele_type()) {
        case Entry::MAJOR_ALLELE:
          (!entry_in_partition1?cost_partition_m1:cost_partition_m2)[1] -= entry.get_phred_score();
          (entry_in_partition1?cost_partition_m1:cost_partition_m2)[1] += entry.get_phred_score();
          break;
        case Entry::MINOR_ALLELE:
          (!entry_in_partition1?cost_partition_m1:cost_partition_m2)[0] -= entry.get_phred_score();
          (entry_in_partition1?cost_partition_m1:cost_partition_m2)[0] += entry.get_phred_score();
          break;
        case Entry::BLANK:
          break;
        default:
          assert(false);
      }
      break;
    case 2u: // father
      switch (entry.get_allele_type()) {
        case Entry::MAJOR_ALLELE:
          (!entry_in_partition1?cost_partition_f1:cost_partition_f2)[1] -= entry.get_phred_score();
          (entry_in_partition1?cost_partition_f1:cost_partition_f2)[1] += entry.get_phred_score();
          break;
        case Entry::MINOR_ALLELE:
          (!entry_in_partition1?cost_partition_f1:cost_partition_f2)[0] -= entry.get_phred_score();
          (entry_in_partition1?cost_partition_f1:cost_partition_f2)[0] += entry.get_phred_score();
          break;
        case Entry::BLANK:
          break;
        default:
          assert(false);
      }
      break;
    default:
      assert(false);
    }
}

unsigned int ColumnCostComputer::get_cost(unsigned int genotypem, unsigned int genotypef, unsigned int genotypec) {
  unsigned int best_cost = numeric_limits<unsigned int>::max();
  // Enumerate all possible assignments of alleles to haplotypes and 
  // compute costs for those which are compatible with genotypes.
  // TODO: This can be done more efficiently.
  for (unsigned int i=0; i<16; ++i) {
    unsigned int allele_m1 = i & 1;
    unsigned int allele_m2 = (i >> 1) & 1;
    if (allele_m1 + allele_m2 != genotypem) continue;
    
    unsigned int allele_f1 = (i >> 2) & 1;
    unsigned int allele_f2 = (i >> 3) & 1;
    if (allele_f1 + allele_f2 != genotypef) continue;
    
    unsigned int allele_c1 = (inheritance_val & 1 == 0)?allele_m1:allele_m2;
    unsigned int allele_c2 = ((inheritance_val>>1) & 1 == 0)?allele_f1:allele_f2;
    if (allele_c1 + allele_c2 != genotypec) continue;
    
    unsigned int cost = cost_partition_m1[allele_m1] + cost_partition_m2[allele_m2] + cost_partition_f1[allele_f1] + cost_partition_f2[allele_f2]; 
    if (cost < best_cost) {
      best_cost = cost;
    }
  }
  
  if (best_cost == numeric_limits<unsigned int>::max()) {
    throw std::runtime_error("Error: Mendelian conflict");
  }
  
  return best_cost;
}

Entry::allele_t ColumnCostComputer::get_allele(bool second_haplotype) {
//   if (all_heterozygous) {
//     unsigned int cost0 = cost_partition1[0] + cost_partition2[1];
//     unsigned int cost1 = cost_partition1[1] + cost_partition2[0];
//     if (cost0 == cost1) {
//       return Entry::EQUAL_SCORES;
//     } else {
//       return ((cost0 < cost1) ^ second_haplotype)?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE;
//     }
//   } else {
    const unsigned int* cost = second_haplotype?cost_partition2:cost_partition1;
    if(cost[0] == cost[1]) {	
      return Entry::EQUAL_SCORES;
    } else {
      return (cost[0]<cost[1]?Entry::MAJOR_ALLELE:Entry::MINOR_ALLELE);
    }
  //}
}

unsigned int ColumnCostComputer::get_weight(bool s) {

  if(s) {
    // if cost_partition2[0] is less than cost_partition2[1], for
    // example, i.e., the haplotype allele will be 0, then the price
    // of flipping this 0 to a 1 -- on the corresponding super-read --
    // should be cost_partition2[1] - cost_partition2[0], so this is
    // the phred score on this new super-read
    return (cost_partition2[0]<=cost_partition2[1]?cost_partition2[1]-cost_partition2[0]:cost_partition2[0]-cost_partition2[1]);
  }
  else {
    return (cost_partition1[0]<=cost_partition1[1]?cost_partition1[1]-cost_partition1[0]:cost_partition1[0]-cost_partition1[1]);
  }
}
