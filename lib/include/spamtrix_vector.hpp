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
class Vector{
    std::vector<real> values;
public:
    Vector(){} //:length(0),values(NULL){}
    Vector(const idx length);
    Vector(const Vector& v); // COPY CONSTRUCTOR
    Vector(const real *val, const idx &length);
    virtual ~Vector();

    real& operator[](const idx i){
#ifdef DEBUG
        assert( i < this->getLength() );
#endif
        return this->values[i];
    }

    real& operator()(const idx i){
        return (*this)[i];
    }
    const real& operator[](const idx i ) const{
#ifdef DEBUG
        assert( i < this->getLength() );
#endif
        return this->values[i];
    }

    void setAllValuesTo(const real val){
        for (idx i = 0 ; i < getLength() ; i++){
            values[i] = val;
        }
    }

    Vector& operator=(const Vector& v);
    Vector& operator=(const real& a);
    Vector& operator-=(const Vector& v);
    Vector& operator+=(const Vector& v);
    Vector& operator+=(const real a);
    Vector& operator-=(const real a);
    Vector& operator*=(const real a);
    const Vector operator-(const Vector& rhs) const;
    const Vector operator+(const Vector& rhs) const;
    const Vector operator*(const real a) const;

    void resize(const idx length);

    idx getLength()const{return (idx) values.size();}

    void print(const char* name = NULL)const;

    real getNorm() const;   // returns vector length
    void normalise();       // normalises to unit vector

    //
    // FUNCTIONS REQUIRED BY PYTHON WRAPPERS
    //
    int __len__(){return (int) values.size();}
    real __getitem__(int i){return values[i];}
    void __setitem__(int i, real v){(*this)[i] = v;}
};


inline const Vector operator*(const real a, const Vector &v)
{
    return Vector(v)*=a;
}
} // end namespace SpaMtrix


#endif

