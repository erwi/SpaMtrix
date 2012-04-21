#include <vector.h>

Vector::Vector(const idx length ):
length(length),
values(NULL)
{
    values = new real[length];
  
    if (!values)
    {
      cout << "error in "<< __func__ << "could not allocate " << length << "elements" << endl; 
    }
    setAllValuesTo(0.0);
}

Vector::~Vector()
{
  this->length = 0;
  delete [] values;
}

Vector::Vector(const Vector& v):
  length(v.length),
  values(NULL)
{
  values = new real[length];
  if (!values)
  {
    cout << "error in " <<__func__<< " could not allocate "<< length <<" elements "<< endl;
    exit(1);
  }
  *this = v;
    
}

Vector& Vector::operator=(const Vector& v)
{
  if ( &v == this)
    return *this;
  
  if (length != v.length )
  {
    length = v.length;
    delete [] values;
    values = NULL;
    values = new real[length];
  }
  
  for (idx i = 0 ; i < length ; ++i )
    this->values[i] = v.values[i];
  
  return *this;
}

void Vector::resize(const idx length)
{
/*!
 *RESIZES THIS VECTOR TO GIVEN length 
 */
  if (this->values)
  {
    delete [] values;
    values = NULL;
  }
  this->length = length;
  values = new real[length];
  
  if (!values)
  {
    cout << "error in "<< __func__ << " could not allocate for " << length << " elements "<< endl;
    exit(1);
  }
  this->setAllValuesTo(0.0);
}

void Vector::setAllValuesTo(const real val)
{
  for (idx i = 0 ; i < length ; i++)
    values[i] = val;
}


real& Vector::operator[](const idx i)
{
#ifdef DEBUG
  assert( i<this->length );
  assert( this->values );
#endif
  
  return this->values[i];
}

const real& Vector::operator[](const idx i ) const
{
 #ifdef DEBUG
  assert( i<this->length );
  assert( this->values );
#endif
  
  return this->values[i]; 
}

const Vector& Vector::operator+=(const Vector& v)
{
#ifdef DEBUG
  assert(this->length == v.length);
#endif

  for (idx i = 0 ; i < length ; ++i)
    this->values[i]+=v.values[i];
  
  return *this;
}

const Vector& Vector::operator-=(const Vector& v)
{
#ifdef DEBUG
  assert(this->length == v.length);
#endif
  for (idx i = 0 ; i < length ; ++i)
    this->values[i]-=v.values[i];
  
  return *this;
}


//============================================
// 	DEBUG FUNCTIONS
//============================================
#ifdef DEBUG
void Vector::print() const
{
  cout << "Vector length " << this->length << endl;
  for (idx i = 0 ; i < length ; i++ )
  {
    cout <<"Vector["<<i<<"] = " << values[i] << endl;
  }
}
#endif



