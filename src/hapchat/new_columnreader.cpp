#include "new_columnreader.h"

void ColumnReader1::read_block() 
{
  unsigned int count_minor = 0;
  unsigned int count_major = 0;
  num_col = 0;

  std::vector<bool>::iterator ik = kind.begin();
  std::vector<bool>::iterator ii = is_homozygous.begin();
  for(Block::const_iterator ib = block.begin();
      ib != block.end();
      ++ib) {
    
    ++num_col;
    
    count_minor = 0;
    count_major = 0;
    for(Column::const_iterator ic = (*ib).begin();
        ic != (*ib).end();
        ++ic) {
            
      if(!(*ic).is_gap() & ((*ic).get_phred_score() != 0)) {
        if ((*ic).get_allele_type() == Entry::MAJOR_ALLELE) {
          ++count_major;
        } else if ((*ic).get_allele_type() == Entry::MINOR_ALLELE) {
          ++count_minor;
        } else {
          std::cerr << "ERROR: allele not equal to 0 or not 1 is in the input!" << std::endl;
          exit(-1);
        }
      }
    }

    if(count_minor == 0 || count_major == 0) {
      *ii = true;
      *ik = (count_major != 0)? false : true; 
      if(jump_homozygous) {
        --num_col;
      }
    }

    ++ii;
    ++ik;
  }
  
  first = block.begin();


  if(jump_homozygous) {
    while ((first != block.end()) && is_homozygous[first - block.begin()]) {
      ++first;
    }
  }
}



bool ColumnReader1::has_next()
{
  if(!next) {
    return false;
  } else {
    if(started) {
      ++iblock;
      if(jump_homozygous) {
        while ((iblock != block.end()) && is_homozygous[iblock - block.begin()]) {
          ++iblock;
        }
      } 
    } else {
      started = true;
    }
  
    next = iblock != block.end();

    return next;
  }
}



Column ColumnReader1::get_next()
{
  if(next) {
    if((*iblock).size() == 0) {
      std::cout << "ERRORRRRRRRRRRRRRRR" << std::endl;
      exit(-1);
    }
    return *iblock;
  } else {
    return Column(0, Entry(-1, Entry::BLANK, 0));
  }
}
