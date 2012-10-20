#ifndef VECTOR_H
#define VECTOR_H

#include <setup.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
using std::cout;
using std::endl;
class Vector{

    //idx length;
    //real* values;
    std::vector<real> values;

public:

    Vector(){} //:length(0),values(NULL){}
    Vector(const idx length);
    Vector(const Vector& v); // COPY CONSTRUCTOR
    Vector(const real *val, const idx &length);
    ~Vector();

    real& operator[](const idx i)
    {
#ifdef DEBUG
        assert( i < this->getLength() );
#endif

        return this->values[i];
    }

    real& operator()(const idx i)
    {
        return (*this)[i];
    }


    const real& operator[](const idx i ) const
    {
#ifdef DEBUG
        assert( i < this->getLength() );
#endif
        return this->values[i];
    }

    void setAllValuesTo(const real val)
    {
        for (idx i = 0 ; i < getLength() ; i++)
            values[i] = val;
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

};


inline const Vector operator*(const real a, const Vector &v)
{
    return Vector(v)*=a;
}



#endif

