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
  FlexiMatrix is a Flexible sparse matrix storage scheme consisting of rows of IndVal pairs
  stored using std::vectors of std::vectors of IndVals.
  */



public:
    idx numDim1, numDim2;   // MATRIX DIMENSIONS.
    // IN CASE OF ROW COMPRESSED MATRIX numDim1 IS NUMBER OF ROWS, numDim2 IS NUMBER OF COLUMNS
    // IN CASE OF COLUMN COMPRESSED MATRIX numDim1 IS NUMBER OF COLUMNS, numDim2 IS NUMBER OF ROWS
    vector<vector<IndVal> > nonZeros;
    FlexiMatrix():numDim1(0),numDim2(0){}
    FlexiMatrix(const idx numDim1, const idx numDim2);


    // NOTE! THESE TWO HAVE NOT BEEN THOROUGHLY TESTED
    void addNonZero(const idx dim1, const idx dim2, const real val = 0.0 );
    void addNonZero(const idx dim1, const IndVal &iv);

    // GET VALUE CURERNTLY USES LINEAR SEARCH. CREATE A BINARY SEARCH VERSION!
    real getValue(const idx dim1, const idx dim2) const;
    void print() const; // PRINTS TO STDOUT
};











#endif // FLEXIMATRIX_H
