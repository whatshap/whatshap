#ifndef COLUMN_READER1_H
#define COLUMN_READER1_H

#include "basic_types.h"
#include "entry.h"

class ColumnReader1 {

public:

  ColumnReader1(const Block &b, const bool &jump) 
    {
      block = b;
      jump_homozygous = jump;
      
      is_homozygous = std::vector<bool> (block.size(), false);
      kind = std::vector<bool> (block.size(), false);
      
      read_block();
      
      restart();
    }

  ~ColumnReader1() { }

  bool has_next();
  Column get_next();

  void restart() { 
    iblock = first;
    started = false;
    next = first != block.end();
  }

  Counter num_cols() { return num_col; }

  bool was_homozygous() {return is_homozygous[iblock - block.begin()]; }
  bool homozigosity() {return kind[iblock - block.begin()]; }
    

private:
  
  Block block;
  Counter num_col;
  bool jump_homozygous;

  Block::const_iterator iblock;
  Block::const_iterator first;
  bool started;
  bool next;

  std::vector<bool> is_homozygous;
  std::vector<bool> kind;

  void read_block();
};

#endif

