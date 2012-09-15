#ifndef SETUP_H
#define SETUP_H
#include <limits>
#include <cmath>
typedef unsigned int idx;	// INDEX
typedef double real;		// VALUES

inline real abs(const real& a)
{
  /*! 
    Oevrloads abs(real) to handle absolute values of floating point numbers without
    truncating to integer.
  */
  return fabs(a); 
} 

struct IndVal{
/*!
  The IndVal struct represents a non-zero position in a sparse matrix at column or row postion ind.
  The nonzero is on row/col IndVal::ind, and its value is IndVal::val
*/
    idx  ind;	// POSITION INDEX, CAN BE EITHER ROW OR COLUMN INDEX
    real val;	// VALUE
    IndVal(const idx ind, const real val):ind(ind), val(val){}
    IndVal():ind(0),val(0){}
};


static const idx MAX_INDEX= std::numeric_limits<idx>::max();

#endif

