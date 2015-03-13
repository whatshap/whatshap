#include <cassert>
#include "columnindexingscheme.h"
#include "columnindexingiterator.h"
#include <math.h>

using namespace std;
static int flag=0;

ColumnIndexingIterator::ColumnIndexingIterator(const ColumnIndexingScheme* parent) {
  flag=1;
  assert(parent != 0);
  this->parent = parent;
  this->graycodes = new GrayCodes(parent->read_ids.size());
  this->index = -1;
  this->forward_projection = -1;
  this->forward_dual_projection = -1; // used for HALF_TABLE
}
ColumnIndexingIterator::ColumnIndexingIterator(const ColumnIndexingScheme* parent, int decimalindx) {
  flag=2;
  assert(parent != 0);
  this->parent = parent;
  //int decimalindx22=parent->decimalindx;
  //cout<<"iterator..."<<decimalindx<<"\n";
  this->whitecodes = new WhiteCodes(parent->read_ids.size(), decimalindx);
  this->index = -1;
  this->forward_projection = -1;
  this->forward_dual_projection = -1; // used for HALF_TABLE
}


ColumnIndexingIterator::~ColumnIndexingIterator() {
  if(flag==1)
      delete graycodes;
  if(flag==2)
    delete whitecodes;
}

bool ColumnIndexingIterator::has_next() {
  return graycodes->has_next();
}

bool ColumnIndexingIterator::has_next_2() {
  return whitecodes->has_next();
}

void ColumnIndexingIterator::advance(int* bit_changed) {
  assert(graycodes->has_next());

  int graycode_bit_changed = -1;
  index = graycodes->get_next(&graycode_bit_changed);
  // first iteration?
  if (graycode_bit_changed == -1) {
    assert(index == 0);
    if (parent->forward_projection_mask != 0) {
      forward_projection = 0;
      //forward_dual_projection = (((unsigned int)1) << parent->forward_projection_width)-1;
    }
  } else {
    if (parent->forward_projection_mask != 0) {
      // index of bit in the forward_projection
      int bit_index = parent->forward_projection_mask->at(graycode_bit_changed);
      if (bit_index >= 0) {
	forward_projection = forward_projection ^ (((unsigned int)1) << bit_index);
	//forward_dual_projection = forward_dual_projection ^ (((unsigned int)1) << bit_index);
      }
    }
  }
  if (bit_changed != 0) {
    *bit_changed = graycode_bit_changed;
  }
}

void ColumnIndexingIterator::advance_2(int** bit_changed, int decimalindx) {
  assert(whitecodes->has_next());
 int* whitecode_bit_changed =(int*)malloc(sizeof(int)*30);
  for(int i=0; i<30; i++)
    whitecode_bit_changed[i] = -1;   
   
  index = whitecodes->get_next(&whitecode_bit_changed);
  if(whitecodes->has_next()==0)
  {
    index = whitecodes->get_next(&whitecode_bit_changed);
  }
 
  int white_flag = 0;
  for(int i=0; i<30; i++)
    if(whitecode_bit_changed[i]!=-1)
      white_flag=1;
  
  
  if (white_flag == 0) {
    assert(index == decimalindx);
    if (parent->forward_projection_mask != 0) {
          int len= parent->forward_projection_width-1;
          int forwardprojarray[len];
          int length= parent->read_ids.size();
          int indxarr[length];
          for (int i=0; i < length; ++i)  // assuming a 32 bit int
             indxarr[i]=decimalindx & (1 << i) ? 1 : 0;
          
          for(int i=0;i<length;i++)
          {
          
           int maskindx=parent->forward_projection_mask->at(i);
         
           if(maskindx>-1)
           forwardprojarray[maskindx]=indxarr[i];
          }
        
          int decforward=0;
          for (int j=0;j<len;j++)  
           decforward= (pow(2,j)*(forwardprojarray[j])) + decforward;
            forward_projection=decforward;
     
        }
        } else {
            if (parent->forward_projection_mask != 0) {
      // index of bit in the forward_projection
            int bit_index;
            for(int i=0; i<30; i++){
             if(whitecode_bit_changed[i]!=-1){
            bit_index = parent->forward_projection_mask->at(whitecode_bit_changed[i]);
            if (bit_index >= 0) {
              forward_projection = forward_projection ^ (((unsigned int)1) << bit_index);
            //forward_dual_projection = forward_dual_projection ^ (((unsigned int)1) << bit_index);
         }
        }
      }
   }
  }
   
   for(int i=0; i<30; i++){
        (*bit_changed)[i] = whitecode_bit_changed[i];
    }
    
  
  std::free(whitecode_bit_changed);
 
}

unsigned int ColumnIndexingIterator::get_forward_projection() {
  assert(index >= 0);
  return forward_projection;
}

unsigned int ColumnIndexingIterator::get_forward_dual_projection() {
  assert(index >= 0);
  return forward_dual_projection;
}

unsigned int ColumnIndexingIterator::get_backward_projection() {
  assert(index >= 0);
  return index & ((((unsigned int)1)<<parent->backward_projection_width) - 1);
}

unsigned int ColumnIndexingIterator::get_index() {
  assert(index >= 0);
  return index;
}

unsigned int ColumnIndexingIterator::get_partition() {
  assert(index >= 0);
  return index;
}

// for backtracking

unsigned int ColumnIndexingIterator::index_backward_projection(unsigned int i) {
  assert(i >= 0); // assert the proper boundaries
  assert(i < (((unsigned int)1) << parent->read_ids.size()));

  return i & ((((unsigned int)1) << parent->backward_projection_width) -1);
}

// in the HALF_TABLE case

unsigned int ColumnIndexingIterator::dual_index_backward_projection(unsigned int i) {
  assert(i >= 0); // assert the proper boundaries
  assert(i < (((unsigned int)1) << parent->read_ids.size()));
  
  // return backward projection of the dual
  return dual_index(i) & ((((unsigned int)1) << parent->backward_projection_width) -1);
}

unsigned int ColumnIndexingIterator::dual_index(unsigned int i) {
  assert(i >= 0);
  assert(i < (((unsigned int)1) << parent->read_ids.size()));

  unsigned int dual_i = i;
  unsigned int s = 1;
  for(int j=0; j< parent->read_ids.size(); ++j) { // compute the dual of i
    dual_i = dual_i ^ s;
    s = s << 1;
  }

  return dual_i;
}

unsigned int ColumnIndexingIterator::index_forward_projection(unsigned int i) {
  assert(i >= 0);
  assert(i < (((unsigned int)1) << parent->read_ids.size()));

  unsigned int i_forward_projection = 0;
  unsigned int s = 1;
  for(int j=0; j< parent->read_ids.size(); ++j) {
    unsigned int m = parent->forward_projection_mask->at(j);
    if(m != -1) {
      unsigned int s = (((unsigned int)1) << m);
      i_forward_projection += (s&i);
    }
  }
  
  return i_forward_projection;
}
