#ifndef FLEXIMATRIX_H
#define FLEXIMATRIX_H
#include <setup.h>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <stdio.h>
using std::vector;

class FlexiMatrix
{
/*!
  FlexiMatrix is a Flexible sparse matrix storage scheme consisting of rows of ColVal pairs
  stored using std::vectors of std::vectors of ColVals. The intended use is for Cholesky and
  LU decompositions.
  */



public:
    idx numRows, numCols;
    vector<vector<ColVal> > nonZeros;
    FlexiMatrix():numRows(0),numCols(0){}
    FlexiMatrix(const idx numRows, const idx numCols);


    // NOTE! THESE TWO HAVE NOT BEEN THOROUGHLY TESTED
    void addNonZero(const idx row, const idx col, const real val = 0.0 );
    void addNonZero(const idx row, const ColVal &cv);

    // GET VALUE CURERNTLY USES LINEAR SEARCH. CREATE A BINARY SEARCH VERSION!
    real getValue(const idx row, const idx col) const;
    void print() const; // PRINTS TO STDOUT
};











#endif // FLEXIMATRIX_H

