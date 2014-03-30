#include <spamtrix_fleximatrix.hpp>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <stdio.h>

namespace SpaMtrix {
FlexiMatrix::FlexiMatrix(const idx numDim1, const idx numDim2):
    numDim1(numDim1),
    numDim2(numDim2),
    nonZeros(numDim1) { }

FlexiMatrix::FlexiMatrix(const IRCMatrix &A):
    numDim1(A.getNumRows()),
    numDim2(A.getNumCols()) {
    /*!
    Constructor that makes a copy of an IRCMatrix.
    If USES_OMP is defined, uses OpenMP to parallelise construction.*/
    nonZeros = std::vector<std::vector<IndVal> >(numDim1);
#ifdef USES_OMP
    #pragma omp parallel for
#endif
    for (idx r = 0; r < numDim1; r++) {
        idx rowStart = A.rows[r];
        idx rowEnd   = A.rows[r + 1];
        nonZeros[r] = std::vector<IndVal>(&A.cvPairs[rowStart] , &A.cvPairs[rowEnd]);
    }
}

FlexiMatrix::~FlexiMatrix() { }

idx FlexiMatrix::calcNumNonZeros() const {
    /*!
     * Returns number of non-zeros allocated. Does this by calculating it, so this gets
     * linearly slower with larger matrices
     */
    idx nnz(0);
    for (idx i = 0 ; i < nonZeros.size() ; ++i)
        nnz += nonZeros[i].size();
    return nnz;
}

void FlexiMatrix::addNonZero(const idx dim1, const idx dim2, const real val) {
    /*!
    Adds storage for a non-zero at location dim1, dim1.
    The value will be initialised to val
    */
    // CHECK DIMENSION1 VECTOR SIZE
#ifdef DEBUG
    if (dim1 >= (idx) nonZeros.size())
        std::cerr << "dim1 = " << dim1 << "nonZeros.size() = " << nonZeros.size() << std::endl;
    assert(dim1 < (idx) nonZeros.size());
#endif
    this->addNonZero(dim1, IndVal(dim2, val));
}

void FlexiMatrix::addNonZero(const idx dim1, const IndVal &iv) {
    /*!
    Adds IndVal pair iv to this matrix, if the same row/column is not already allocated
    */
    assert(dim1 < (idx) nonZeros.size());
    // IF EMPTY ROW OR NEW VALUE PLACED AT END, JUST APPEND TO ROW AND EXIT FUNCTION
    if ((nonZeros[dim1].size() == 0) ||
            (nonZeros[dim1].back().ind < iv.ind)) {
        nonZeros[dim1].push_back(iv);
        return;
    }
    // FIND CORRECT POSITION BY SEARCHING FOR COLUMN POSITIONS.
    // ALL COLUMN VALUES MUST INCREASE, I.E MAINTAINING AN ASCENDIGLY
    // SORTED VECTOR OF NON-ZERO COLUMNS
    // LOWER BOUND RETURNS ITERATOR TO FIRST ELEMENT THAT DOES NOT
    // COMPARE TO LESS THAN dim2
    std::vector<IndVal>::iterator itr =
        std::lower_bound(nonZeros[dim1].begin(),
                         nonZeros[dim1].end(),
                         iv,
    [](const IndVal & iv1, const IndVal & iv2) {
        return iv1.ind < iv2.ind;   // lamda column comparison
    }
                        );
    // IF ADDING DUPLICATE NONZER, ONLY WRITE VALUE TO EXISTING MEMORY
    if (itr->ind == iv.ind) {
        itr->val = iv.ind;
    }
    // OTEHRWISE ADD NEW VALUE
    else {
        nonZeros[dim1].insert(itr, iv);
    }
    // TODO - THIS SHOULD BE MAINTAINED INTERNALLY
    // MAINTAIN DIMENSION SIZES
    //this->numDim1 = numDim1 > nonZeros.size() ? numDim1 : nonZeros.size();
    //this->numDim2 = numDim2 > (iv.ind+1) ? numDim2 : (iv.ind+1);
}

real FlexiMatrix::getValue(const idx dim1, const idx dim2) const {
    /*!
    returns the value held at row, col. If id does not exist, return 0.0 by default
    */
    //if (nonZeros[dim1].empty() )
    IndVal temp(dim2, 0.0);
    // GET ITERATOR TO FIRST ELEMENT THAT DOES NOT COMPARE TO
    // LESS THAN dim2
    std::vector<IndVal>::const_iterator itr =
        std::lower_bound(nonZeros[dim1].begin() ,
                         nonZeros[dim1].end(),
                         temp,
    [](const IndVal & iv1, const IndVal & iv2) {
        return iv1.ind < iv2.ind;   // lamda column comparison
    }
                        );
    if (itr == nonZeros[dim1].end())
        return 0.0;
    else if (itr->ind == dim2)
        return itr->val;
    else
        return 0.0;
// LINEAR SEARCH.
    /*
        std::vector<IndVal>::const_iterator itr1 = nonZeros[dim1].begin();
        std::vector<IndVal>::const_iterator itr2 = nonZeros[dim1].end();
        for ( ; itr1 != itr2 ; itr1++)
        {
            if (itr1->ind == dim2)
                return itr1->val;
        }
        */
    return 0.0;
}

void FlexiMatrix::setValue(const idx dim1, const idx dim2, const real val) {
    /*!
      Sets value at (dim1,dim2) to value val. If a non-zero does not already
      exist at (dim1,dim2), a new non-zero is inserted
      */
#ifdef DEBUG
    assert(dim1 < numDim1);
    assert(dim2 < numDim2);
#endif
    real *nnzval;
    if (isNonZero(dim1, dim2, nnzval)) {    // IF NON-ZERO STORAGE ALREADY EXISTS, SET ITS VALUE
        *nnzval = val;
    } else                                  // OTHERWISE INSERT NEW NONZERO LOCATION
        addNonZero(dim1, dim2, val);
}


void FlexiMatrix::print() const {
    /*!
    Prints the values held in the matrix to stdout
    */
    std::cout << "FlexiMatrix " << numDim1 << ", " << numDim2 << std::endl;
    for (idx r = 0 ; r < nonZeros.size() ; r++) {
        for (idx c = 0 ; c < numDim2 ; c++) {
            real val = this->getValue(r, c);
            printf("%1.3f\t", val);
        }
        printf("\n");
    }
}

bool FlexiMatrix::isNonZero(const idx dim1, const idx dim2, real *&val) {
    /*!
    Returns true if memory is allocated for a non-zero value at index dim1,dim2.

    Input argument val is a pointer to location of addres of stored value if non-zero
    (or NULL pointer otherwise).
    */
    val = NULL; // ASSUME STORAGE IS NOT ALLOCATED AT dim1,dim2
    if (dim1 >= (idx) nonZeros.size())
        return false;
    // IF VECTOR HAS NOT INITIALISED
    if (nonZeros[dim1].empty())
        return false;
    // SET ITERATOR TO FIRST ELEMENT THAT DOES NOT COMPARE TO LESS
    // THAN dim2.
    IndVal iv(dim2, 0.0);
    std::vector<IndVal>::iterator itr =
        std::lower_bound(nonZeros[dim1].begin(),
                         nonZeros[dim1].end(),
                         iv,
    [](const IndVal & iv1, const IndVal & iv2) {
        return iv1.ind < iv2.ind;   // lamda column comparison
    }
                        );
// itr WILL POINT TO END IF ALL INDEXES ARE LESS THAN dim2
    if (itr == nonZeros[dim1].end())
        return false;
// itr NOW POINTS EITHER TO dim2 OR THE NEXT AFTER IT
    if ((*itr).ind == dim2) {
        val = &(itr->val); // SET POINTER TO VALUE AT INDEX dim1,dim2
        return true;
    } else
        return false; // NON-ZERO AT INDEX dim1,dim2 NOT FOUND
}
} // end namespace SpaMtrix
