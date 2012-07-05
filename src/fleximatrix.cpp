#include <fleximatrix.h>

bool compare_columns(const ColVal &cv1, const ColVal &cv2)
{
    /*!
      Implements the < operation between a ColVal and a column index.
      Returns cv.col < col
      */
    return (cv1.col < cv2.col);

}

FlexiMatrix::FlexiMatrix(const idx numRows, const idx numCols):
    numRows(numRows),
    numCols(numCols),
    nonZeros(numRows)
{

}

void FlexiMatrix::addNonZero(const idx row, const idx col, const real val)
{
/*!
    Adds storage for a non-zero at location row, column.
    The value will be initialised to val
*/
    // FIRST MAKE SURE ENOUGH ROWS EXIST
    assert( row < (idx) nonZeros.size() );
    this->addNonZero( row, ColVal(col,val) );

}

void FlexiMatrix::addNonZero(const idx row, const ColVal &cv)
{
/*!
  Adds ColVal pair cv to this matrix, if the same row/column is not already allocated
*/
    assert( row < (idx) nonZeros.size() );

    // IF EMPTY ROW OR NEW VALUE PLACED AT END, JUST PUSH IT BACK AND EXIT
    if  ( (nonZeros[row].size() == 0 ) ||
          (nonZeros[row].back().col < cv.col) )
    {
        nonZeros[row].push_back(cv);
        return;
    }


    // FIND CORRECT POSITION BY SEARCHING FOR COLUMN POSITIONS.
    // ALL COLUMN VALUES MUST INCREASE, I.E MAINTAINING AN ASCENDIGLY
    // SORTED VECTOR OF NON-ZERO COLUMNS

    std::vector<ColVal>::iterator itr =
    std::upper_bound(nonZeros[row].begin(),
                     nonZeros[row].end(),
                     cv,
                     compare_columns
                     );

    // CHECK THAT AN ENTRY FOR THIS COLUMN DOES NOT ALREADY EXIST.
    // IF IT DOES, IT WOULD EXIST IMMEDIATELY BEFORE THE RETURNED ITERATOR
    idx i = itr - nonZeros[row].begin(); // ITERATOR IS AT i'th POSITION

    if ( ( nonZeros[row].begin() + (i-1) )->col == cv.col ) // CHECK ( i-1 )th POSITION
        return;

    nonZeros[row].insert(itr, cv); // INSERT TO LIST

}

real FlexiMatrix::getValue(const idx row, const idx col) const
{
/*!
  returns the value held at row, col. If id does not exist, return 0.0 by default
  */

    // LINEAR SEARCH.
    std::vector<ColVal>::const_iterator itr1 = nonZeros[row].begin();
    std::vector<ColVal>::const_iterator itr2 = nonZeros[row].end();
    for ( ; itr1 != itr2 ; itr1++)
    {
        if (itr1->col == col)
            return itr1->val;
    }
    return 0.0;
}

void FlexiMatrix::print() const
{
/*!
  Prints the values held in the matrix to stdout
*/
    std::cout << "FlexiMatrix "<< numRows <<", "<<numCols << std::endl;
    for (idx r = 0 ; r < nonZeros.size() ; r++)
    {
        for (idx c = 0 ; c < numCols ; c++ )
        {
            real val = this->getValue(r, c);
            printf("%1.3f\t", val );
        }
        printf("\n");
    }
}
