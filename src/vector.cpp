#include <vector.h>

Vector::Vector(const idx length )
{
    //values = new real[length];
    //values.reserve( length );
    values = std::vector<real>(length,0.0);
    if (values.size() != length )
    {
      cout << "error in "<< __func__ << "could not allocate " << length << "elements" << endl; 
    }
    //setAllValuesTo(0.0);
}

Vector::~Vector()
{

}

Vector::Vector(const Vector& v)
{
    values = v.values;
    if (values.size() != v.getLength() )
    {
        cout << "error in " <<__func__<< " could not allocate "<< v.getLength() <<" elements "<< endl;
        exit(1);
    }
}

Vector& Vector::operator=(const Vector& v)
{
    if ( &v == this)
        return *this;

    values = v.values;
    return *this;
}

Vector& Vector::operator=(const real& a)
{
    idx len = this->getLength();
    
#ifdef USES_OMP
#pragma omp parallel for
#endif
    for (idx i = 0 ; i < len ; ++i)
        values[i] = a;
    return *this;
}



void Vector::resize(const idx length)
{
/*!
 *RESIZES THIS VECTOR TO GIVEN length 
 */

    if (this->getLength() == length )
        return;
    else
        values = std::vector<real>(length,0.0);

}







Vector& Vector::operator+=(const Vector& v)
{
#ifdef DEBUG
    assert(this->getLength() == v.getLength());
#endif
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i)
        this->values[i]+=v.values[i];

    return *this;
}
Vector& Vector::operator-=(const Vector& v)
{
#ifdef DEBUG
    assert(this->getLength() == v.getLength() );
#endif
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i)
        this->values[i]-=v.values[i];

    return *this;
}

Vector& Vector::operator+=(const real a)
{
    idx len = this->getLength();
    for (idx i = 0; i < len ; ++i)
        values[i]+= a;
    return *this;
}

Vector& Vector::operator-=(const real a)
{
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i)
        values[i]-= a;
    return *this;
}
Vector& Vector::operator *=(const real a)
{
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i)
        values[i]*= a;
    return *this;
}


const Vector Vector::operator+(const Vector& rhs) const
{
    idx len = this->getLength();
#ifdef DEBUG
    assert(len == rhs.getLength() );
#endif
    return Vector(*this)+= rhs;
}
const Vector Vector::operator-(const Vector& rhs) const
{
    idx len = this->getLength();
#ifdef DEBUG
    assert(len == rhs.getLength() );
#endif
    return Vector(*this)-= rhs;
}
const Vector Vector::operator*(const real a) const
{
    return Vector(*this)*=a;
}

real Vector::getNorm() const
{
/*!
  returns scalar length of this vector.
  */

    real sum(0);
#ifdef USES_OMP
#pragma omp parallel for reduction(+:sum)
#endif
    for (idx i = 0 ; i < this->getLength() ; i++ )
        sum += values[i]*values[i];

    return sqrt(sum);


}

void Vector::normalise()
{
    /*!
      scales each component to make this a unit vector.
      If length is intially zero, returns without doing anything
      */

    real norm = this->getNorm();

    if (norm==0.0) // do nothing is zero-vector
        return;

    real k = 1.0 / norm;
#ifdef USES_OMP
    #pragma omp parallel for
#endif
    for (idx i = 0 ; i < this->getLength() ; i++)
        values[i]*= k;

}


//============================================
// 	DEBUG FUNCTIONS
//============================================

void Vector::print(const char* name) const
{
    cout << "Vector length " << this->getLength() << endl;
  for (idx i = 0 ; i < getLength() ; i++ )
  {
    if (name)
      cout << name;
    else
      cout << "Vector";
    cout <<"["<<i<<"] = " << values[i] << endl;
  }
}




