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
   
  real& operator[](const idx i);
  const real& operator[](const idx i) const;
  
  Vector& operator=(const Vector& v);
  const Vector& operator-=(const Vector& v);
  const Vector& operator+=(const Vector& v);
  
  void resize(const idx length);
  void setAllValuesTo(const real val);
  idx getLength()const{return this->length;}
#ifdef DEBUG
  void print()const;
#endif
  
};


#endif

