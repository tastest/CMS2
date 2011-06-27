
#include "TBitSet.hh"

ClassImp(TBitSet)

using namespace std;

/***************************************************************************
*                                 MACROS
***************************************************************************/

/* make CHAR_BIT 8 if it's not defined in limits.h */
#ifndef CHAR_BIT
#warning CHAR_BIT not defined.  Assuming 8 bits.
#define CHAR_BIT 8
#endif

/* Bits in int's and long int's */
#define BITS_IN_INT 32

/* position of bit within character */
#define BIT_CHAR(bit)         ((bit) / CHAR_BIT)

/* array index for character containing bit */
#define BIT_IN_CHAR(bit)      (1 << (CHAR_BIT - 1 - ((bit)  % CHAR_BIT)))

/* number of characters required to contain number of bits */
#define BITS_TO_CHARS(bits)   ((((bits) - 1) / CHAR_BIT) + 1)

/* most significant bit in a character */
#define MS_BIT                (1 << (CHAR_BIT - 1))

//-----------------------------------------------------------------------------
// Constructors and Destructor
//-----------------------------------------------------------------------------
/***************************************************************************
*   Method     : TBitSet - default constructor
*   Description: This is the TBitSet constructor.  It reserves memory
*                for the vector storing the array.
*   Parameters : numBits - number of bits in the array
*   Effects    : Allocates vectory for array bits
*   Returned   : None
***************************************************************************/
TBitSet::TBitSet() {

  /* allocate space for bit array */
  m_Array.reserve(BITS_TO_CHARS(0));
  
  /* insert unsigned chars of 0 */
  for(int i = 0; i < BITS_TO_CHARS(0); ++i) {
    m_Array.push_back(0);
  }
  m_NumBits = 0;

}

/***************************************************************************
*   Method     : TBitSet - constructor
*   Description: This is the TBitSet constructor.  It reserves memory
*                for the vector storing the array.
*   Parameters : numBits - number of bits in the array
*   Effects    : Allocates vectory for array bits
*   Returned   : None
***************************************************************************/
TBitSet::TBitSet(int numBits) {

  /* allocate space for bit array */
  m_Array.reserve(BITS_TO_CHARS(numBits));
  
  /* insert unsigned chars of 0 */
  for(int i = 0; i < BITS_TO_CHARS(numBits); ++i) {
    m_Array.push_back(0);
  }
  m_NumBits = numBits;

}

/***************************************************************************
*   Method     : TBitSet - constructor
*   Description: This is the TBitSet constructor.  It copies the
*                for contents of a vector of unsigned char into the
*                bitarray.
*   Parameters : vect - vector to be copied
*                numBits - number of bits in the array
*   Effects    : Allocates vectory for array bits
*   Returned   : None
***************************************************************************/
TBitSet::TBitSet(const std::vector<unsigned char> vect, int numBits) {

  m_Array = vect;
  m_NumBits = numBits;

}

/***************************************************************************
*   Method     : TBitSet - constructor
*   Description: This is the TBitSet constructor.  It reserves memory
*                for the vector storing the array.
*   Parameters : numBits - number of bits in the array
*   Effects    : Allocates vectory for array bits
*   Returned   : None
***************************************************************************/
TBitSet::TBitSet(int word, int numBits) {

  /* allocate space for bit array */
  m_Array.reserve(BITS_TO_CHARS(numBits));
  
  /* insert unsigned chars of 0 */
  for(int i = 0; i < BITS_TO_CHARS(numBits); ++i) {
    m_Array.push_back(0);
  }
  m_NumBits = numBits;

  if (numBits > BITS_IN_INT) {
    printf("Warning in TBitSet Constructor: # of bits exceeds 32!\n");
  }

  for (int i = numBits-1; i >= 0; i--) {
    if   ((word & (0x1 << i)) != 0) SetBit  (i);
    else                            ClearBit(i);
  }
}

/***************************************************************************
*   Method     : ~TBitSet - destructor
*   Description: This is the TBitSet destructor.  At this point it's
*                just a place holder.
*   Parameters : None
*   Effects    : None
*   Returned   : None
***************************************************************************/
TBitSet::~TBitSet() {
//    /* I think std::vector destructor will automatically delete m_Array */
//  std::vector<unsigned char> dummy(0);
//  m_Array.swap(dummy); /* this deletes the internal memory in the array */
}

//-----------------------------------------------------------------------------
// Methods
//-----------------------------------------------------------------------------

/***************************************************************************
*   Method     : Dump
*   Description: This method dumps the conents of a bit array to stdout.
*                The format of the dump is a series of bytes represented in
*                hexadecimal. LSB starting from left!!
*   Parameters : outStream - stream to write to
*   Effects    : Array contents are dumped to stdout
*   Returned   : None
***************************************************************************/
void TBitSet::Dump(std::ostream &outStream) {

  outStream.width(2);
  outStream.fill('0');
  outStream << uppercase << hex << (int)(m_Array[0]);  /* first byte */
  
  for (unsigned int i = 1; i < m_Array.size(); i++) {

    /* remaining bytes with a leading space */
    outStream << " ";
    outStream.width(2);
    outStream.fill('0');
    outStream << (int)(m_Array[i]);
  }
  
  outStream << dec;
}

/***************************************************************************
*   Method     : DumpBinary
*   Description: This method dumps the conents of a bit array to stdout.
*                The format of the dump is a series of bytes represented in
*                binary. LSB starting from right!!
*   Parameters : outStream - stream to write to
*   Effects    : Array contents are dumped to stdout
*   Returned   : None
***************************************************************************/
void TBitSet::DumpBinary(std::ostream &outStream) {
  for (int i = m_NumBits-1; i >= 0; i--) {
    if ((*this)[i]) outStream << "1";
    else            outStream << "0";
  }
}
/***************************************************************************
*   Method     : AllBitsAreTrue
*   Description: Checks that all bits are true 
*   Parameters : None
*   Effects    : None
*   Returned   : True/False of all bits are true
***************************************************************************/
bool TBitSet::AllBitsAreTrue(void) {

  for (int i = m_NumBits-1; i >= 0; i--) {
    if (!((*this)[i])) return false;
  }
  return true;
}

/***************************************************************************
*   Method     : AllOtherBitsAreTrue
*   Description: Checks that all bits are true, ignoring supplied bit
*   Parameters : bit - to be ignored
*   Effects    : None
*   Returned   : True/False of all bits are true
***************************************************************************/
bool TBitSet::AllOtherBitsAreTrue(unsigned int bit, unsigned int bit2, unsigned int bit3) {

  for (int i = m_NumBits-1; i >= 0; i--) {
    if (i == bit) continue;
    if (i == bit2) continue; 
    if (i == bit3) continue; // if there are 4 billion bits the 0xffffffff defaule will screw up (oh well).
    if (!((*this)[i])) return false;
  }
  return true;
}

/***************************************************************************
*   Method     : GetIntegerWord
*   Description: Return an integer for bit set. True == 1, False == 0
*   Parameters : None
*   Effects    : None
*   Returned   : Integer
***************************************************************************/
int TBitSet::GetIntegerWord(void) {

  if (m_NumBits > BITS_IN_INT) {
    printf("Warning in TBitSet::GetIntegerWord: # of bits exceeds 32!\n");
  }

  int intWord = 0;
  for (int i = m_NumBits-1; i >= 0; i--) {
    if ((*this)[i]) intWord |= (0x1 << i);
  }
  return intWord;
}

/***************************************************************************
*   Method     : SetAll
*   Description: This method sets every bit in the bit array to 1.  This
*                method uses UCHAR_MAX to determine what it means to set
*                all bits in an unsigned char, so it is crucial that the
*                machine implementation of unsigned char utilizes all of
*                the bits in the memory allocated for an unsigned char.
*   Parameters : None
*   Effects    : Each of the bits used in the bit array are set to 1.
*                Unused (spare) bits are set to 0.
*   Returned   : None
***************************************************************************/
void TBitSet::SetAll(void) {

  int bits;
  unsigned char mask;
  
  /* set bits in all bytes to 1 */
  m_Array.assign(m_Array.size(), UCHAR_MAX);
  
  /* zero any spare bits so increment and decrement are consistent */
  bits = m_NumBits % CHAR_BIT;
  if (bits != 0) {
    mask = UCHAR_MAX << (CHAR_BIT - bits);
    m_Array[BIT_CHAR(m_NumBits - 1)] = mask;
  }
}

/***************************************************************************
*   Method     : ClearAll
*   Description: This method sets every bit in the bit array to 0.
*   Parameters : None
*   Effects    : Each of the bits in the bit array are set to 0.
*   Returned   : None
***************************************************************************/
void TBitSet::ClearAll(void) {

  /* set bits in all bytes to 0 */
  m_Array.assign(m_Array.size(), 0);
}

/***************************************************************************
*   Method     : SetBit
*   Description: This method sets a bit in the bit array to 1.
*   Parameters : bit - the number of the bit to set
*   Effects    : The specified bit will be set to 1.
*   Returned   : None
***************************************************************************/
void TBitSet::SetBit(unsigned int bit) {

  if (m_NumBits <= bit) {
    return;         /* bit out of range */
  }
  
  m_Array[BIT_CHAR(bit)] |= BIT_IN_CHAR(bit);
}

/***************************************************************************
*   Method     : ClearBit
*   Description: This method sets a bit in the bit array to 0.
*   Parameters : bit - the number of the bit to clear
*   Effects    : The specified bit will be set to 0.
*   Returned   : None
***************************************************************************/
void TBitSet::ClearBit(unsigned int bit) {
  
  unsigned char mask;

  if (m_NumBits <= bit) {
    return;         /* bit out of range */
  }

  /* create a mask to zero out desired bit */
  mask =  BIT_IN_CHAR(bit);
  mask = ~mask;

  m_Array[BIT_CHAR(bit)] &= mask;
}

/***************************************************************************
*   Method     : operator()
*   Description: Overload of the () operator.  This method approximates
*                array indices used for assignment.  It returns a
*                TBitSet_index which includes an = method used to
*                set bit values.
*   Parameters : bit - index of array bit
*   Effects    : None
*   Returned   : TBitSet_index (pointer to bit)
***************************************************************************/
TBitSet_index TBitSet::operator()(unsigned int bit) {

  TBitSet_index result(this, bit);
  
  return result;
}

/***************************************************************************
*   Method     : operator[]
*   Description: Overload of the [] operator.  This method returns the
*                value of a bit in the bit array.
*   Parameters : bit - index of array bit
*   Effects    : None
*   Returned   : The value of the specified bit.
***************************************************************************/
bool TBitSet::operator[](unsigned int bit) {
  
  return((m_Array[BIT_CHAR(bit)] & BIT_IN_CHAR(bit)) != 0);
}

/***************************************************************************
*   Method     : operator==
*   Description: overload of the == operator
*   Parameters : other - bit array to compare
*   Effects    : None
*   Returned   : True if this == other.  Otherwise false.
***************************************************************************/
bool TBitSet::operator==(TBitSet other) {

  if (m_NumBits != other.m_NumBits) {
    /* unequal sizes */
    return false;
  }
  
  return (this->m_Array == other.m_Array);
}

/***************************************************************************
*   Method     : operator!=
*   Description: overload of the != operator
*   Parameters : other - bit array to compare
*   Effects    : None
*   Returned   : True if this != other.  Otherwise false.
***************************************************************************/
bool TBitSet::operator!=(TBitSet other) {

  if (m_NumBits != other.m_NumBits) {
    /* unequal sizes */
    return true;
  }

  return (this->m_Array != other.m_Array);
}

/***************************************************************************
*   Method     : operator<
*   Description: overload of the < operator
*   Parameters : other - bit array to compare
*   Effects    : None
*   Returned   : True if this < other.  Otherwise false.
***************************************************************************/
bool TBitSet::operator<(TBitSet other) {
  
  if (m_NumBits != other.m_NumBits) {
    /* unequal sizes */
    return false;
  }

  return (this->m_Array < other.m_Array);
}

/***************************************************************************
*   Method     : operator<=
*   Description: overload of the <= operator
*   Parameters : other - bit array to compare
*   Effects    : None
*   Returned   : True if this <= other.  Otherwise false.
***************************************************************************/
bool TBitSet::operator<=(TBitSet other) {

  if (m_NumBits != other.m_NumBits) {
    /* unequal sizes */
    return false;
  }

  return (this->m_Array <= other.m_Array);
}

/***************************************************************************
*   Method     : operator>
*   Description: overload of the > operator
*   Parameters : other - bit array to compare
*   Effects    : None
*   Returned   : True if this > other.  Otherwise false.
***************************************************************************/
bool TBitSet::operator>(TBitSet other) {
  if (m_NumBits != other.m_NumBits) {
    /* unequal sizes */
    return false;
  }

  return (this->m_Array > other.m_Array);
}

/***************************************************************************
*   Method     : operator>=
*   Description: overload of the >= operator
*   Parameters : other - bit array to compare
*   Effects    : None
*   Returned   : True if this >= other.  Otherwise false.
***************************************************************************/
bool TBitSet::operator>=(TBitSet other) {
  
  if (m_NumBits != other.m_NumBits) {
    /* unequal sizes */
    return false;
  }

  return (this->m_Array >= other.m_Array);
}

/***************************************************************************
*   Method     : operator~
*   Description: overload of the ~ operator.  Negates all non-spare bits in
*                bit array
*   Parameters : None
*   Effects    : None
*   Returned   : value of this after bitwise not
***************************************************************************/
TBitSet TBitSet::operator~(void) {
  
  TBitSet result(this->m_Array, this->m_NumBits);
  result.Not();

  return result;
}

/***************************************************************************
*   Method     : operator&
*   Description: overload of the & operator.  Performs a bitwise and
*                between the source array and this bit array.
*   Parameters : other - bit array on righthand side of &
*   Effects    : None
*   Returned   : value of bitwise and of this and other.
***************************************************************************/
TBitSet TBitSet::operator&(TBitSet other) {
  
  TBitSet result(this->m_Array, this->m_NumBits);
  result &= other;
  
  return result;
}

/***************************************************************************
*   Method     : operator^
*   Description: overload of the ^ operator.  Performs a bitwise xor
*                between the source array and this bit array.
*   Parameters : other - bit array on righthand side of ^
*   Effects    : None
*   Returned   : value of bitwise xor of this and other.
***************************************************************************/
TBitSet TBitSet::operator^(TBitSet other) {

  TBitSet result(this->m_Array, this->m_NumBits);
  result ^= other;

  return result;
}

/***************************************************************************
*   Method     : operator|
*   Description: overload of the | operator.  Performs a bitwise or
*                between the source array and this bit array.
*   Parameters : other - bit array on righthand side of |
*   Effects    : None
*   Returned   : value of bitwise or of this and other.
***************************************************************************/
TBitSet TBitSet::operator|(TBitSet other) {

  TBitSet result(this->m_Array, this->m_NumBits);
  result |= other;
  
  return result;
}

/***************************************************************************
*   Method     : operator++ (prefix)
*   Description: overload of the ++ operator.  Increments the contents of
*                a bit array.  Overflows cause rollover.
*   Parameters : None
*   Effects    : Bit array contents are incremented
*   Returned   : None
***************************************************************************/
void TBitSet::operator++(void) {

  int i;
  unsigned char maxValue;     /* maximum value for current char */
  unsigned char one;          /* least significant bit in current char */

  if (m_Array.size() == 0) {
    return;         /* nothing to increment */
  }

  /* handle arrays that don't use every bit in the last character */
  i = (m_NumBits % CHAR_BIT);
  if (i != 0) {
    maxValue = UCHAR_MAX << (CHAR_BIT - i);
    one = 1 << (CHAR_BIT - i);
  }
  else {
    maxValue = UCHAR_MAX;
    one = 1;
  }

  for (i = BIT_CHAR(m_NumBits - 1); i >= 0; i--) {
    if (m_Array[i] != maxValue) {
      m_Array[i] = m_Array[i] + one;
      return;
    }
    else {
      /* need to carry to next byte */
      m_Array[i] = 0;
      
      /* remaining characters must use all bits */
      maxValue = UCHAR_MAX;
      one = 1;
    }
  }
}

/***************************************************************************
*   Method     : operator++ (postfix)
*   Description: overload of the ++ operator.  Increments the contents of
*                a bit array.  Overflows cause rollover.
*   Parameters : dumy - needed for postfix increment
*   Effects    : Bit array contents are incremented
*   Returned   : None
***************************************************************************/
void TBitSet::operator++(int dummy) {

  ++(*this);
}

/***************************************************************************
*   Method     : operator-- (prefix)
*   Description: overload of the -- operator.  Decrements the contents of
*                a bit array.  Underflows cause rollover.
*   Parameters : None
*   Effects    : Bit array contents are decremented
*   Returned   : None
***************************************************************************/
void TBitSet::operator--(void) {
  
  int i;
  unsigned char maxValue;     /* maximum value for current char */
  unsigned char one;          /* least significant bit in current char */
  
  if (m_Array.size() == 0) {
    return;         /* nothing to decrement */
  }

  /* handle arrays that don't use every bit in the last character */
  i = (m_NumBits % CHAR_BIT);
  if (i != 0) {
    maxValue = UCHAR_MAX << (CHAR_BIT - i);
    one = 1 << (CHAR_BIT - i);
  }
  else {
    maxValue = UCHAR_MAX;
    one = 1;
  }

  for (i = BIT_CHAR(m_NumBits - 1); i >= 0; i--) {
    if (m_Array[i] >= one) {
      m_Array[i] = m_Array[i] - one;
      return;
    }
    else {
      /* need to borrow from the next byte */
      m_Array[i] = maxValue;
      
      /* remaining characters must use all bits */
      maxValue = UCHAR_MAX;
      one = 1;
    }
  }
}

/***************************************************************************
*   Method     : operator-- (postfix)
*   Description: overload of the -- operator.  Decrements the contents of
*                a bit array.  Underflows cause rollover.
*   Parameters : dumy - needed for postfix decrement
*   Effects    : Bit array contents are decremented
*   Returned   : None
***************************************************************************/
void TBitSet::operator--(int dummy) {

  --(*this);
}

/***************************************************************************
*   Method     : operator=
*   Description: overload of the = operator.  Copies source contents into
*                this bit array.
*   Parameters : src - Source bit array
*   Effects    : Source bit array contents are copied into this array
*   Returned   : None
***************************************************************************/
void TBitSet::operator=(TBitSet src) {

  if (m_NumBits != src.m_NumBits) {

    /* don't do assignment with different array sizes */
    return;
  }

  if ((m_Array.size() == 0) || (src.m_Array.size() == 0)) {
    /* don't do assignment with unallocated array */
    return;
  }

  /* copy bits from source */
  this->m_Array = src.m_Array;
}

/***************************************************************************
*   Method     : operator&=
*   Description: overload of the &= operator.  Performs a bitwise and
*                between the source array and this bit array.  Contents of
*                this bit array will be the result.
*   Parameters : src - Source bit array
*   Effects    : Results of bitwise and are stored in this array
*   Returned   : None
***************************************************************************/
void TBitSet::operator&=(TBitSet src) {

  if (m_NumBits != src.m_NumBits) {
    /* don't do assignment with different array sizes */
    return;
  }

  /* AND array one unsigned char at a time */
  for(unsigned int i = 0; i < m_Array.size(); i++) {
    m_Array[i] = m_Array[i] & src.m_Array[i];
  }
}

/***************************************************************************
*   Method     : operator^=
*   Description: overload of the ^= operator.  Performs a bitwise not
*                between the source array and this bit array.  Contents of
*                this bit array will be the result.
*   Parameters : src - Source bit array
*   Effects    : Results of bitwise not are stored in this array
*   Returned   : None
***************************************************************************/
void TBitSet::operator^=(TBitSet src) {

  if (m_NumBits != src.m_NumBits) {
    /* don't do assignment with different array sizes */
    return;
  }

  /* XOR array one unsigned char at a time */
  for(unsigned int i = 0; i < m_Array.size(); i++) {
    m_Array[i] = m_Array[i] ^ src.m_Array[i];
  }
}

/***************************************************************************
*   Method     : operator|=
*   Description: overload of the |= operator.  Performs a bitwise or
*                between the source array and this bit array.  Contents of
*                this bit array will be the result.
*   Parameters : src - Source bit array
*   Effects    : Results of bitwise or are stored in this array
*   Returned   : None
***************************************************************************/
void TBitSet::operator|=(TBitSet src) {

  if (m_NumBits != src.m_NumBits) {
    /* don't do assignment with different array sizes */
    return;
  }

  /* OR array one unsigned char at a time */
  for(unsigned int i = 0; i < m_Array.size(); i++) {
    m_Array[i] = m_Array[i] | src.m_Array[i];
  }
}

/***************************************************************************
*   Method     : Not
*   Description: Negates all non-spare bits in bit array
*   Parameters : None
*   Effects    : Contents of bit array are negated.  Any spare bits are
*                left at 0.
*   Returned   : None
***************************************************************************/
void TBitSet::Not(void) {

  int bits;
  unsigned char mask;
  
  if (m_Array.size() == 0) {
    /* don't do not with unallocated array */
    return;
  }

  /* NOT array one unsigned char at a time */
  for(unsigned int i = 0; i < m_Array.size(); i++) {
    m_Array[i] = ~m_Array[i];
  }

  /* zero any spare bits so increment and decrement are consistent */
  bits = m_NumBits % CHAR_BIT;
  if (bits != 0) {
    mask = UCHAR_MAX << (CHAR_BIT - bits);
    m_Array[BIT_CHAR(m_NumBits - 1)] &= mask;
  }
}

/***************************************************************************
*   Method     : operator<<=
*   Description: overload of the <<= operator.  Performs a left shift on
*                this bit array.  Contents of this bit array will be the
*                result.
*   Parameters : shifts - number of bit positions to shift
*   Effects    : Results of the shifts are stored in this array
*   Returned   : None
***************************************************************************/
void TBitSet::operator<<=(unsigned int shifts) {

  unsigned int i;
  int chars = shifts / CHAR_BIT; /* number of whole byte shifts */
  
  shifts = shifts % CHAR_BIT;    /* number of bit shifts remaining */
  
  if (shifts >= m_NumBits) {
    /* all bits have been shifted off */
    this->ClearAll();
    return;
  }

  /* first handle big jumps of bytes */
  if (chars > 0) {
    for (i = 0; (i + chars) < m_Array.size(); i++) {
      m_Array[i] = m_Array[i + chars];
    }

    /* now zero out new bytes on the right */
    for (i = m_Array.size(); chars > 0; chars--) {
      m_Array[i - chars] = 0;
    }
  }

  /* now we have at most CHAR_BIT - 1 bit shifts across the whole array */
  for (i = 0; i < shifts; i++) {
    for (unsigned int j = 0; j < BIT_CHAR(m_NumBits - 1); j++) {
      m_Array[j] <<= 1;
      
      /* handle shifts across byte bounds */
      if (m_Array[j + 1] & MS_BIT) {
        m_Array[j] |= 0x01;
      }
    }
    
    m_Array[BIT_CHAR(m_NumBits - 1)] <<= 1;
  }
}

/***************************************************************************
*   Method     : operator>>=
*   Description: overload of the >>= operator.  Performs a right shift on
*                this bit array.  Contents of this bit array will be the
*                result.
*   Parameters : shifts - number of bit positions to shift
*   Effects    : Results of the shifts are stored in this array
*   Returned   : None
***************************************************************************/
void TBitSet::operator>>=(unsigned int shifts) {

  int i;
  char mask;
  int chars = shifts / CHAR_BIT;  /* number of whole byte shifts */
  
  shifts = shifts % CHAR_BIT;     /* number of bit shifts remaining */

  if (shifts >= m_NumBits) {
    /* all bits have been shifted off */
    this->ClearAll();
    return;
  }

  /* first handle big jumps of bytes */
  if (chars > 0) {
    for (i = BIT_CHAR(m_NumBits - 1); (i - chars) >= 0; i--) {
      m_Array[i] = m_Array[i - chars];
    }

    /* now zero out new bytes on the right */
    for (; chars > 0; chars--) {
      m_Array[chars - 1] = 0;
    }
  }

  /* now we have at most CHAR_BIT - 1 bit shifts across the whole array */
  for (i = 0; i < (int)shifts; i++) {
    for (unsigned int j = BIT_CHAR(m_NumBits - 1); j > 0; j--) {
      m_Array[j] >>= 1;
      
      /* handle shifts across byte bounds */
      if (m_Array[j - 1] & 0x01) {
        m_Array[j] |= MS_BIT;
      }
    }
    m_Array[0] >>= 1;
  }
  
  /***********************************************************************
   * zero any spare bits that are beyond the end of the bit array so
   * increment and decrement are consistent.
   ***********************************************************************/
  i = m_NumBits % CHAR_BIT;
  if (i != 0) {
    mask = UCHAR_MAX << (CHAR_BIT - i);
    m_Array[BIT_CHAR(m_NumBits - 1)] &= mask;
  }
}

/***************************************************************************
*   Method     : TBitSet_index - constructor
*   Description: This is the TBitSet_index constructor.  It stores a
*                pointer to the bit array and the bit index.
*   Parameters : array - pointer to bit array
*                index - index of bit in array
*   Effects    : Pointer to bit array and bit index are stored.
*   Returned   : None
***************************************************************************/
TBitSet_index::TBitSet_index(TBitSet *array, unsigned int index) {

  m_BitArray = array;
  m_Index = index;
}

/***************************************************************************
*   Method     : operator=
*   Description: overload of the = operator.  Sets the bit array bit to
*                the value of src.
*   Parameters : src - bit value
*   Effects    : Bit pointed to by this object is set to the value of
*                source.
*   Returned   : None
***************************************************************************/
void TBitSet_index::operator=(bool src) {

  if (m_BitArray == NULL) {
    return;     /* no array */
  }

  if (m_BitArray->Size() <= m_Index) {
    return;     /* index is out of bounds */
  }

  if (src) {
    m_BitArray->SetBit(m_Index);
  }
  else {
    m_BitArray->ClearBit(m_Index);
  }
}

/***************************************************************************
*   Method     : operator<<
*   Description: overload of the << operator for printing
*   Parameters : an ostream
*   Effects    : prints the contents to the ostream
*   Returned   : the original ostream
***************************************************************************/
ostream& operator<<(ostream& os,TBitSet* p)
{  
  p->DumpBinary(os);
  return os;
}

/***************************************************************************
*   Method     : operator<<
*   Description: overload of the << operator for printing
*   Parameters : an ostream
*   Effects    : prints the contents to the ostream
*   Returned   : the original ostream
***************************************************************************/
ostream& operator<<(ostream& os,TBitSet& p)
{  
  p.DumpBinary(os);
  return os;
}
