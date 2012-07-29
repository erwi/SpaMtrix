#ifndef VECTOR_H
#define VECTOR_H

#include <setup.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <stdlib.h>
using std::cout;
using std::endl;
class Vector{
  
  idx length;
  real* values;
public:
  
  Vector():length(0),values(NULL){};
  Vector(const idx length);
  Vector(const Vector& v); // COPY CONSTRUCTOR
  
  ~Vector();
   
  real& operator[](const idx i)
  {
    #ifdef DEBUG
    assert( i<this->length );
    assert( this->values );
    #endif
  
    return this->values[i];
  }
  const real& operator[](const idx i ) const
  {
    #ifdef DEBUG
    assert( i<this->length );
    assert( this->values );
    #endif
    return this->values[i]; 
  }
  
  void setAllValuesTo(const real val)
  {
      for (idx i = 0 ; i < length ; i++)
      values[i] = val;
  }

  
  
  Vector& operator=(const Vector& v);
  const Vector& operator-=(const Vector& v);
  const Vector& operator+=(const Vector& v);
  
  void resize(const idx length);
  
  idx getLength()const{return this->length;}

  void print(const char* name = NULL)const;

  
};


#endif

