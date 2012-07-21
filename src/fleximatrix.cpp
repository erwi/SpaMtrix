#include <fleximatrix.h>

bool compare_columns(const IndVal &iv1, const IndVal &iv2)
{
    /*!
      Implements the < operation between a IndVal and a position index.
      Returns iv.ind < ind
      */
    return (iv1.ind < iv2.ind); // IF FIRST INDEX IS SMALLER THAN SECOND

}

FlexiMatrix::FlexiMatrix(const idx numDim1, const idx numDim2):
    numDim1(numDim1),
    numDim2(numDim2),
    nonZeros(numDim1)
{

}

void FlexiMatrix::addNonZero(const idx dim1, const idx dim2, const real val)
{
/*!
    Adds storage for a non-zero at location dim1, dim1.
    The value will be initialised to val
*/
    // CHECK DIMENSION1 VECTOR SIZE
#ifdef DEBUG
    assert( dim1 < (idx) nonZeros.size() );
#endif
    this->addNonZero( dim1, IndVal(dim2,val) );

}

void FlexiMatrix::addNonZero(const idx dim1, const IndVal &iv)
{
/*!
  Adds IndVal pair iv to this matrix, if the same row/column is not already allocated
*/
    assert( dim1 < (idx) nonZeros.size() );

    // IF EMPTY ROW OR NEW VALUE PLACED AT END, JUST PUSH IT BACK AND EXIT
    if  ( (nonZeros[dim1].size() == 0 ) ||
          (nonZeros[dim1].back().ind < iv.ind) )
    {
        nonZeros[dim1].push_back(iv);
        return;
    }


    // FIND CORRECT POSITION BY SEARCHING FOR COLUMN POSITIONS.
    // ALL COLUMN VALUES MUST INCREASE, I.E MAINTAINING AN ASCENDIGLY
    // SORTED VECTOR OF NON-ZERO COLUMNS

    std::vector<IndVal>::iterator itr =
    std::upper_bound(nonZeros[dim1].begin(),
                     nonZeros[dim1].end(),
                     iv,
                     compare_columns
                     );

    // CHECK THAT AN ENTRY FOR THIS COLUMN DOES NOT ALREADY EXIST.
    // IF IT DOES, IT WOULD EXIST IMMEDIATELY BEFORE THE RETURNED ITERATOR
    idx i = itr - nonZeros[dim1].begin(); // ITERATOR IS AT i'th POSITION

    if ( ( nonZeros[dim1].begin() + (i-1) )->ind == iv.ind ) // CHECK ( i-1 )th POSITION
        return;

    nonZeros[dim1].insert(itr, iv); // INSERT TO LIST

}

real FlexiMatrix::getValue(const idx dim1, const idx dim2) const
{
/*!
  returns the value held at row, col. If id does not exist, return 0.0 by default
  */

    // LINEAR SEARCH.
    std::vector<IndVal>::const_iterator itr1 = nonZeros[dim1].begin();
    std::vector<IndVal>::const_iterator itr2 = nonZeros[dim1].end();
    for ( ; itr1 != itr2 ; itr1++)
    {
        if (itr1->ind == dim2)
            return itr1->val;
    }
    return 0.0;
}

void FlexiMatrix::print() const
{
/*!
  Prints the values held in the matrix to stdout
*/
    std::cout << "FlexiMatrix "<< numDim1 <<", "<<numDim2 << std::endl;
    for (idx r = 0 ; r < nonZeros.size() ; r++)
    {
        for (idx c = 0 ; c < numDim2 ; c++ )
        {
            real val = this->getValue(r, c);
            printf("%1.3f\t", val );
        }
        printf("\n");
    }
}
