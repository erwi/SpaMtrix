#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <assert.h>

#include <spamtrix_densematrix.hpp>
DenseMatrix::DenseMatrix(const idx numRows, const idx numCols):
    numRows(numRows),
    numCols(numCols)
{

    values = std::vector<real>(numRows*numCols,0.0); //new real[numRows*numCols];

    if (values.size() !=  numRows*numCols )
    {
        std::cout << "error in " <<__func__ << " could not allocate dense matrix" << std::endl;
        exit(1);
    }

    //setAllValuesTo(0);

}

DenseMatrix::~DenseMatrix()
{
   // delete [] values;
}


void DenseMatrix::setAllValuesTo(const real v)
{
    for (idx c = 0 ; c < numCols ; ++c)
    {
        for (idx r = 0 ; r < numRows ; ++r)
        {
            idx ind = c*numRows + r;
            values[ind] = v;
        }
    }
}

real& DenseMatrix::operator()(const idx row, const idx col)
{
    /*!
      Row/column access operator
      */
#ifdef DEBUG
    assert(row<numRows);
    assert(col<numCols);
#endif

    const idx ind = col*numRows + row;
    return values[ind];
}

real DenseMatrix::operator ()(const idx row, const idx col) const
{
    /*!
      Row/column access operator
      */
#ifdef DEBUG
    assert(row<numRows);
    assert(col<numCols);
#endif

    const idx ind = col*numRows + row;
    return values[ind];
}


void DenseMatrix::print() const
{
    for (idx r = 0 ; r < numRows ; r++)
    {
        for (idx c = 0  ; c < numCols ; c++)
        {
            //std::cout << values[c*numRows + r] << " ";
            printf("%1.2f ", values[c*numRows + r] );
        }
        printf("\n"); fflush(stdout);
        //std::cout << std::endl;
    }
}
