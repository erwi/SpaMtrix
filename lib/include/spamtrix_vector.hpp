#ifndef VECTOR_H
#define VECTOR_H


#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include <spamtrix_setup.hpp>

namespace SpaMtrix
{
/**
 * @brief Dynamically sized dense vector of scalar values.
 */
class Vector{
    std::vector<real> values;
public:
    /** @brief Construct an empty vector. */
    Vector(){}
    /**
     * @brief Construct a vector with the given length, initialized to zero.
     *
     * @param length Number of entries.
     */
    Vector(const idx length);
    /**
     * @brief Copy-construct a vector.
     *
     * @param v Source vector.
     */
    Vector(const Vector& v);
    /**
     * @brief Construct a vector from a raw array.
     *
     * @param val Source values.
     * @param length Number of entries to copy.
     */
    Vector(const real *val, const idx &length);
    /** @brief Destroy the vector. */
    virtual ~Vector();

    /**
     * @brief Access a mutable element by index.
     *
     * @param i Zero-based index.
     * @return Reference to the selected element.
     */
    real& operator[](const idx i){
#ifdef DEBUG
        assert( i < this->getLength() );
#endif
        return this->values[i];
    }

    /**
     * @brief Access a mutable element by index using function-call syntax.
     *
     * @param i Zero-based index.
     * @return Reference to the selected element.
     */
    real& operator()(const idx i){
        return (*this)[i];
    }
    /**
     * @brief Access a const element by index.
     *
     * @param i Zero-based index.
     * @return Const reference to the selected element.
     */
    const real& operator[](const idx i ) const{
#ifdef DEBUG
        assert( i < this->getLength() );
#endif
        return this->values[i];
    }

    /**
     * @brief Set every element to the same value.
     *
     * @param val Value to write into every entry.
     */
    void setAllValuesTo(const real val){
        for (idx i = 0 ; i < getLength() ; i++){
            values[i] = val;
        }
    }

    /** @brief Copy-assign from another vector. */
    Vector& operator=(const Vector& v);
    /** @brief Fill every element with a scalar value. */
    Vector& operator=(const real& a);
    /** @brief Subtract another vector in place. */
    Vector& operator-=(const Vector& v);
    /** @brief Add another vector in place. */
    Vector& operator+=(const Vector& v);
    /** @brief Add a scalar to every element in place. */
    Vector& operator+=(const real a);
    /** @brief Subtract a scalar from every element in place. */
    Vector& operator-=(const real a);
    /** @brief Scale the vector in place. */
    Vector& operator*=(const real a);
    /** @brief Return the vector difference. */
    const Vector operator-(const Vector& rhs) const;
    /** @brief Return the vector sum. */
    const Vector operator+(const Vector& rhs) const;
    /** @brief Return the scaled vector. */
    const Vector operator*(const real a) const;

    /**
     * @brief Resize the vector.
     *
     * @param length New length.
     */
    void resize(const idx length);

    /**
     * @brief Get the number of entries.
     *
     * @return Vector length.
     */
    idx getLength()const{return (idx) values.size();}

    /**
     * @brief Print the vector to stdout.
     *
     * @param name Optional label.
     */
    void print(const char* name = NULL)const;

    /**
     * @brief Compute the Euclidean norm.
     *
     * @return The 2-norm of the vector.
     */
    real getNorm() const;
    /** @brief Scale the vector to unit norm. */
    void normalise();

    //
    // FUNCTIONS REQUIRED BY PYTHON WRAPPERS
    //
    /** @brief Python wrapper length hook. */
    int __len__(){return (int) values.size();}
    /** @brief Python wrapper element access hook. */
    real __getitem__(int i){return values[i];}
    /** @brief Python wrapper element assignment hook. */
    void __setitem__(int i, real v){(*this)[i] = v;}
};


inline const Vector operator*(const real a, const Vector &v)
{
    return Vector(v)*=a;
}
} // end namespace SpaMtrix


#endif

