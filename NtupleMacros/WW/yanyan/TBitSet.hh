#ifndef _TBITSET_HH_
#define _TBITSET_HH_

#include <vector>
#include <ostream>
#include <iostream>
#include <string.h>

#include "TObject.h"

class TBitSet;

class TBitSet_index {
public:
  TBitSet_index(TBitSet *array, unsigned int index);

  /* assignment */
  void operator=(bool src);

private:
  TBitSet*     m_BitArray;        /* array index applies to */
  unsigned int m_Index;           /* index of bit in array */
};

class TBitSet : public TObject {
  
public:

  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TBitSet();
  TBitSet(int numBits);
  TBitSet(const std::vector<unsigned char> vect, int numBits);
  TBitSet(int word, int numBits);

  ~TBitSet();

  void Dump      (std::ostream &outStream);  // LSB starting from left!!
  void DumpBinary(std::ostream &outStream);  // LSB starting from right!!

  const unsigned int Size() { return m_NumBits; };

  /* plain english accessors */
  bool IsTrue (int bit)  { return (*this)[bit]; }
  bool IsFalse(int bit)  { return !IsTrue(bit); }
  void SetTrue(int bit)  { SetBit(bit); } 
  void SetFalse(int bit) { ClearBit(bit); } 
  void SetAllBitsFalse() { ClearAll();  }
  bool AllBitsAreTrue();
  bool AllOtherBitsAreTrue(unsigned int bit, unsigned int bit2=0xffffffff, unsigned int bit3=0xffffffff);
  int  GetIntegerWord();
  int  GetNumBits()      { return m_NumBits; }


  /* set/clear functions */
  void SetAll  (void);
  void ClearAll(void);
  void SetBit  (unsigned int bit);
  void ClearBit(unsigned int bit);
  
  TBitSet_index operator()(unsigned int bit);
  
  /* boolean operator */
  bool operator[] (unsigned int bit);
  bool operator== (TBitSet other);
  bool operator!= (TBitSet other);
  bool operator<  (TBitSet other);
  bool operator<= (TBitSet other);
  bool operator>  (TBitSet other);
  bool operator>= (TBitSet other);
  
  /* bitwise operators */
  TBitSet operator& (TBitSet other);
  TBitSet operator^ (TBitSet other);
  TBitSet operator| (TBitSet other);
  TBitSet operator~ (void);
  
  /* increment/decrement */
  void operator++ (void);          /* prefix */
  void operator++ (int dummy);     /* postfix */
  void operator-- (void);          /* prefix */
  void operator-- (int dummy);     /* postfix */
  
  /* assignments */
  void operator=  (TBitSet src);
  void operator&= (TBitSet src);
  void operator^= (TBitSet src);
  void operator|= (TBitSet src);
  void Not(void);                 /* negate (~=) */
  
  void operator<<= (unsigned int shifts);
  void operator>>= (unsigned int shifts);
  
protected:
  unsigned int m_NumBits;                 /* number of bits in the array */
  std::vector<unsigned char> m_Array;     /* vector of characters */
  
private:
  
  ClassDef(TBitSet,1)
};


/* printing operators */
ostream& operator<<(ostream& os,TBitSet* p);
ostream& operator<<(ostream& os,TBitSet& p);


#endif

