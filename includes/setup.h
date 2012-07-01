#ifndef SETUP_H
#define SETUP_H
#include <limits>
typedef unsigned int idx;	// INDEX
typedef double real;		// VALUES

/*!
  The ColVal struct represents a non-zero value in a sparse matrix along column.
  The nonzero is on column ColVal::col, and its value is ColVal::val
*/
struct ColVal{
    idx col;	// COLUMN INDEX
    real val;	// VALUE
};


static const idx MAX_INDEX= std::numeric_limits<idx>::max();

#endif

