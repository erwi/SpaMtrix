#ifndef SETUP_H
#define SETUP_H

#include <cmath>
#include <limits>

/** @brief Integer type used for indices and dimensions. */
typedef unsigned int idx;

/** @brief Floating-point type used for matrix and vector values. */
typedef double real;

/**
 * @brief Sparse matrix entry storing an index-value pair.
 */
struct IndVal {
    /** @brief Column or row index. */
    idx ind;
    /** @brief Stored value. */
    real val;

    /**
     * @brief Construct an index-value pair.
     *
     * @param ind Index position.
     * @param val Stored value.
     */
    IndVal(const idx ind, const real val) : ind(ind), val(val) {}

    /** @brief Construct a zero-valued entry at index 0. */
    IndVal() : ind(0), val(0) {}
};

/** @brief Maximum representable index value. */
static const idx MAX_INDEX = std::numeric_limits<idx>::max();

#endif
