#ifndef SETUP_H
#define SETUP_H
#include <limits>
typedef unsigned int idx;	// INDEX
typedef double real;		// VALUES

/*!
  The IndVal struct represents a non-zero value in a sparse matrix at column or row postion ind.
  The nonzero is on row/col IndVal::ind, and its value is IndVal::val
*/
struct IndVal{
    idx  ind;	// POSITION INDEX, CAN BE EITHER ROW OR COLUMN INDEX
    real val;	// VALUE
    IndVal(const idx ind, const real val):ind(ind), val(val){}
    IndVal():ind(0),val(0){}
};


static const idx MAX_INDEX= std::numeric_limits<idx>::max();

#endif

