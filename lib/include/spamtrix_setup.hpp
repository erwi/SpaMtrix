#ifndef SETUP_H
#define SETUP_H
#include <limits>
#include <cmath>
typedef unsigned int idx;	// INDEX
typedef double real;		// VALUES


//inline real abs(const real& a)
//{
//  /*! 
//    Oevrloads abs(real) to handle absolute values of floating point numbers without
//    truncating to integer.
//  */
//  return fabs(a); 
//} 

struct IndVal{
/*!
  The IndVal struct represents a non-zero position in a sparse matrix at column or row position ind.
  The nonzero is on row/col IndVal::ind, and its value is IndVal::val
*/
    idx  ind;	// POSITION INDEX, CAN BE EITHER ROW OR COLUMN INDEX
    real val;	// VALUE
    IndVal(const idx ind, const real val):ind(ind), val(val){}
    IndVal():ind(0),val(0){}
};


//inline bool compare_columns(const IndVal &iv1, const IndVal &iv2)
//{
//    /*!
//      Implements the < operation between a IndVal and a position index.
//      Returns iv.ind < ind
//      */
//    return (iv1.ind < iv2.ind); // IF FIRST INDEX IS SMALLER THAN SECOND
//}







static const idx MAX_INDEX= std::numeric_limits<idx>::max();

#endif

