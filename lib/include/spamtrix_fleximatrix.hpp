#ifndef FLEXIMATRIX_H
#define FLEXIMATRIX_H
#include <vector>
#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>

namespace SpaMtrix {
  class IRCMatrix;

  /*!
    FlexiMatrix is a Flexible sparse matrix storage scheme consisting of rows of IndVal pairs
    stored using std::vectors of std::vectors of IndVals.
  */
  class FlexiMatrix {
      size_t numCols_ = 0;
      std::vector<std::vector<IndVal> > nonZeros;
    public:

    FlexiMatrix(): numCols_(0) {}

    FlexiMatrix(const IRCMatrix &A);
    virtual ~FlexiMatrix();

    /** Add storage for non-zero at (row, col) */
    void addNonZero(size_t row, size_t col, real value = 0.0);

    [[nodiscard]] idx calcNumNonZeros() const;

    /** Return the value at (row,col) or 0 by default */
    [[nodiscard]] real getValue(size_t row, size_t col) const;
    [[nodiscard]] size_t getNumRows() const { return nonZeros.size(); };
    [[nodiscard]] size_t getNumCols() const { return numCols_; };

    [[nodiscard]] std::vector<IndVal>& row(size_t row);
    [[nodiscard]] const std::vector<IndVal>& row(size_t row) const;

    // SETS VALUE, IF IT DOES NOT EXIST, INSERTS IT
    void setValue(const idx dim1, const idx dim2, const real val);
    void print() const; // PRINTS TO STDOUT

    [[nodiscard]] bool isNonZero(const idx dim1, const idx dim2, real *&val);
  };
} // end namespace SpaMtrix










#endif // FLEXIMATRIX_H

