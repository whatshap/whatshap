#include <bitset>
#include <limits>
#include <cassert>
#include <cstdio>

#include "graycodes.h"

using namespace std;

#if 0  // ~original

GrayCodes::GrayCodes(int length) {
  //cout << "length is : " << length << endl;
  //cout << "max is : " << numeric_limits<GrayCodes::int_t>::digits << endl;
  assert(length <= numeric_limits<GrayCodes::int_t>::digits);
  this->length = length;
  this->s = ~((int_t)0);
  this->c = 0;
  this->i = -1;
  this->changed_bit = -1;
}

// mt
int GrayCodes::get_length() const { return length; }

bool GrayCodes::has_next() {
  return i < length;
}


GrayCodes::int_t GrayCodes::get_next(int* changed_bit) {
  GrayCodes::int_t result = c;

  //printf("[%d  bit_changed=%d  s=%u    c=%d ]  ", result, *changed_bit, s, c ); 

  if (changed_bit != 0) {
    *changed_bit = this->changed_bit;
  }
  i = 0;
  while (i<length) {
    int_t mask = ((GrayCodes::int_t)1) << i;
    if (((c&mask) ^ (s&mask)) != 0) {
      c = c ^ mask;
      this->changed_bit = i;
      break;
    }
    s = s ^ mask; 
    i += 1;
  }
 
  //printf("%d  bit_changed=%d  s=%u    c=%d  \n", result, *changed_bit, s, c ); 
  
  return result;
}
#else  // mt


GrayCodes::GrayCodes(int length) {
  //cout << "length is : " << length << endl;
  //cout << "max is : " << numeric_limits<GrayCodes::int_t>::digits << endl;
  assert(length <= numeric_limits<GrayCodes::int_t>::digits);
  this->length = (1<<length);   
  this->s = ~((int_t)0);
  this->c = 0;
  this->i = -1;
  this->changed_bit = -1;
}

// mt
int GrayCodes::get_length() const { return length; }
bool GrayCodes::has_next() {  
    return i < length;
}
void GrayCodes::reset(int length) {
  this->length = (1<<length);   
  this->s = ~((int_t)0);
  this->c = 0;
  this->i = -1;
  this->changed_bit = -1;
}


static GrayCodes::int_t binaryToGray(GrayCodes::int_t x) { return (x>>1)^x;}

#if defined(HAVE_BUILTIN_CLZ)
static inline int log2i(int x) {
    if (x==0) return 0;
    return sizeof(int) * 8 - __builtin_clz(x) - 1;
}
#else
static inline int log2i(int x) {
    if (x==0) return 0;
    int c=0;
    while((x = (x >> 1))) ++c;
    return c;
}
#endif


GrayCodes::int_t GrayCodes::get_current() {
    return  binaryToGray(i-1);
}

GrayCodes::int_t GrayCodes::get_next(int* changed_bit) {
    if (i == -1) { 
        ++i;
        *changed_bit = -1;
        c= binaryToGray(i);
        //printf("%d  bit_changed=-1\n", c);
        ++i;
        return c;
    }

    GrayCodes::int_t result = binaryToGray(i);
    *changed_bit = log2i(result ^ c);
    c=result;
    ++i;
    //printf("%d  bit_changed=%d \n", result, *changed_bit); 
    return result;
}
GrayCodes::int_t GrayCodes::get_next(int* changed_bit, const int idx) {
    if (i == -1) {
        i=idx;
        *changed_bit = -1;
        c= binaryToGray(idx);
        //	printf("%d   bit_changed=-1\n", c);
        ++i;
        return c;
    }
    GrayCodes::int_t result = binaryToGray(idx);
    *changed_bit = log2i(result ^ c);
    c=result;
    i=idx+1;
    //printf("%d  bit_changed=%d \n", result, *changed_bit); 
    return result;
}
#endif  // mt
