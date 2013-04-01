#include <spamtrix_vector.hpp>

namespace SpaMtrix
{

Vector::Vector(const idx length ){
    values = std::vector<real>(length,0.0);
    if (values.size() != length ){
      std::cout << "error in "<< __func__ << "could not allocate " << length << "elements" << std::endl; 
    }
}

Vector::~Vector(){
}

Vector::Vector(const Vector& v){
    values = v.values;
    if (values.size() != v.getLength() ){
        std::cout << "error in " <<__func__<< " could not allocate "<< v.getLength() <<" elements "<< std::endl;
        exit(1);
    }
}

Vector::Vector(const real *val, const idx &length){
  /*! Constructor that copies values from an existing array of reals*/
  values = std::vector<real>(val, val+length );
}


Vector& Vector::operator=(const Vector& v){
    if ( &v == this)
        return *this;
    values = v.values;
    return *this;
}

Vector& Vector::operator=(const real& a){
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i){
        values[i] = a;
    }
    return *this;
}

void Vector::resize(const idx length){
/*!
 *RESIZES THIS VECTOR TO GIVEN length 
 */
    if (this->getLength() == length ){
        return;
    }
    else{
        values = std::vector<real>(length,0.0);
    }

}

Vector& Vector::operator+=(const Vector& v){
#ifdef DEBUG
    assert(this->getLength() == v.getLength());
#endif
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i){
        this->values[i]+=v.values[i];
    }
    return *this;
}
Vector& Vector::operator-=(const Vector& v){
#ifdef DEBUG
    assert(this->getLength() == v.getLength() );
#endif
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i){
        this->values[i]-=v.values[i];
    }
    return *this;
}

Vector& Vector::operator+=(const real a){
    idx len = this->getLength();
    for (idx i = 0; i < len ; ++i){
        values[i]+= a;
    }
    return *this;
}

Vector& Vector::operator-=(const real a){
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i){
        values[i]-= a;
    }
    return *this;
}
Vector& Vector::operator *=(const real a){
    idx len = this->getLength();
    for (idx i = 0 ; i < len ; ++i){
        values[i]*= a;
    }
    return *this;
}


const Vector Vector::operator+(const Vector& rhs) const{
#ifdef DEBUG
    assert(this->getLength() == rhs.getLength() );
#endif
    return Vector(*this)+= rhs;
}
const Vector Vector::operator-(const Vector& rhs) const{
#ifdef DEBUG
    assert(this->getLength() == rhs.getLength() );
#endif
    return Vector(*this)-= rhs;
}
const Vector Vector::operator*(const real a) const{
    /*!Scales by scalar*/
    return Vector(*this)*=a;
}

real Vector::getNorm() const{
/*!
  returns scalar length of this vector.
  */
    real sum(0);
    for (idx i = 0 ; i < this->getLength() ; i++ ){
        sum += values[i]*values[i];
    }
    return sqrt(sum);
}

void Vector::normalise(){
    /*!
      scales each component to make this a unit vector.
      If length is intially zero, returns without doing anything
      */
    real norm = this->getNorm();
    if (norm==0.0){ // do nothing is zero-vector
        return;
    }
    real k = 1.0 / norm;
    for (idx i = 0 ; i < this->getLength() ; i++){
        values[i]*= k;
    }
}
//============================================
// 	DEBUG FUNCTIONS
//============================================
void Vector::print(const char* name) const{
/*!
 * Prints vector values to stdout.
*/
    std::cout << "Vector length " << this->getLength() << std::endl;
    for (idx i = 0 ; i < getLength() ; i++ ){
        if (name){
            std::cout << name;
        }
        else{
            std::cout << "Vector";
        }
        std::cout <<"["<<i<<"] = " << values[i] << std::endl;
    }
}
} // end namespace SpaMtrix
