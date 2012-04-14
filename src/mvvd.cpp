
/*! DOUBLE PRECISION VECTOR */

#include <iostream>                                 
#include <mvvd.h>

/*! Constructor, empty vector*/
Vector::Vector():
                    p_(0),
                    dim_(0) ,
                    ref_(0){}
/*! Constructor of vector of length n. Manages its own memory space.*/
Vector::Vector( int n) :
                    p_(new double[n]),
                    dim_(n),
                    ref_(0)
{
    /*! \param n = length of vector to be created.*/

    if (p_ == NULL)
    {
        std::cerr << "Error: NULL pointer in Vector(int) constructor " << "\n";
        std::cerr << "       Most likely out of memory... " << "\n";
        exit(1);
    }
}

 
/*! Constructor of vector of length n, where all components are set to value v*/
Vector::Vector( int n, const double& v) :
                    p_(new double[n]),
                    dim_(n), ref_(0)
{
    /*! \param n = length of vector to be created.*/
    /*! \param v = value to which all vector components are set to*/
    if (p_ == NULL)
    {
        std::cerr << "Error: NULL pointer in Vector(int) constructor " << "\n";
        std::cerr << "       Most likely out of memory... " << "\n";
        exit(1);
    }
    for (int i=0; i<n; i++)
        p_[i] = v;
}

// operators and member functions
/*! Makes this vector a pointer to some other vectors data. Added by Eero*/
void Vector::point_to(double* p , int len)
{
    /*! \param *p = pointer to external data*/
    /*! \param len = length of external data*/

    // CLEAR OLD DATA FIRST, IF OWNED BY THIS VECTOR
    if (( this->p_ ) &&     // IF POINTS TO DATA AND
        ( !this->ref_ ) )   // IF OWNS THE DATA
    {
        delete [] this->p_;
    }

    this->p_    = p;
    this->dim_  = len;
    this->ref_  = 1;
}

/*! Assignment operator, assigns value to all vector components*/
Vector& Vector::operator=(const double & m)
{
/*!
    \param m = value assigend to all vector components.
*/

    // TODO : SURELY MANUAL LOOP UNROLLING DOES NOTHING IN THIS DAY AND AGE...

    // unroll loops to depth of length 4

    int N = size();

    int Nminus4 = N-4;
    int i;

    for (i=0; i<Nminus4; )
    {
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
    }

    for (; i<N; p_[i++] = m);   // finish off last piece...


    return *this;
}

/*! Resizes vector */
Vector& Vector::newsize( int n)
{
    /*! \param n = new length of vector. */

    if (ref_ )                  // is this structure just a pointer?
    {
        {
            std::cerr << "MV_Vector::newsize can't operator on references.\n";
            exit(1);
        }
    }
    else
    if (dim_ != n )                     // only delete and new if
    {                                   // the size of memory is really
        if (p_) delete [] p_;           // changing, otherwise just
        p_ = new double[n];              // copy in place.
        if (p_ == NULL)
        {
            std::cerr << "Error : NULL pointer when resizing vector in " << __func__ <<std::endl;// "\n";
            exit(1);
        }
        dim_ = n;
    }


    return *this;
}

    

/*! Assigment operator */
Vector& Vector::operator=(const Vector & m)
{

    // CHECK FOR SELF-REFERENCE, ADDED BY EERO
    if (this == &m)
        return *this;

    int N = m.dim_;
    int i;

    if (ref_ )                  // is this structure just a pointer?
    {
        if (dim_ != m.dim_)     // check conformance,
        {
            std::cerr << "MV_VectorRef::operator=  non-conformant assignment.\n";
            exit(1);
        }

        // handle overlapping matrix references
        if ((m.p_ + m.dim_) >= p_)
        {
            // overlap case, copy backwards to avoid overwriting results
            for (i= N-1; i>=0; i--)
                p_[i] = m.p_[i];
        }
        else
        {
            for (i=0; i<N; i++)
                p_[i] = m.p_[i];
        }
                
    }
    else
    {
        newsize(N);

        // no need to test for overlap, since this region is new
        for (i =0; i< N; i++)       // careful not to use bcopy()
            p_[i] = m.p_[i];        // here, but double::operator= double.
    }
    return *this;   
}

/*! Copy-constructor*/
Vector::Vector(const Vector & m) :
                p_(new double[m.dim_]),
                dim_(m.dim_) , ref_(0)
{
    if (p_ == NULL)
    {
        std::cerr << "Error:  Null pointer in Vector(const MV_Vector&); " << "\n";
        exit(1);
    }

    int N = m.dim_;

    for (int i=0; i<N; i++)
        p_[i] = m.p_[i];
}

/*! Constructor of vector of length with values assigned from raw array.
Makes copy of raw data and manages its own memory*/
Vector::Vector(double* d,  int n) :
                p_(new double[n]),
                dim_(n) ,
                ref_(0)
{
    /*!
    \param d = pointer to raw data array with values
    \param n = length of raw data
*/

    if (p_ == NULL)
    {
        std::cerr << "Error: Null pointer in Vector(double*, int) " << "\n";
        exit(1);
    }
    for (int i=0; i<n; i++)
        p_[i] = d[i];

}



Vector::Vector(const double* d,  int n) : p_(new double[n]),
      dim_(n) , ref_(0)
{
    if (p_ == NULL)
    {
        std::cerr << "Error: Null pointer in Vector(double*, int) " << "\n";
        exit(1);
    }
    for (int i=0; i<n; i++)
        p_[i] = d[i];

}


//Vector Vector::operator()(void)
//{
//    return Vector(p_, dim_, ref_);
//}


//const Vector Vector::operator()(void) const
//{
//    return Vector(p_, dim_, MV_Vector_::ref);
//}


Vector Vector::operator()(const MV_VecIndex &I)
{
    // default parameters
    if (I.all())
        return Vector(p_, dim_, ref_); // CHECK THAT ref_ DOESNT DO NASTIES
    else
    {
    // check that index is not out of bounds
    //
        if ( I.end() >= dim_)
        {
            std::cerr << "MV_VecIndex: (" << I.start() << ":" << I.end() << 
                ") too big for matrix (0:" << dim_ - 1 << ") " << "\n";
            exit(1);
        }
        return Vector(p_+ I.start(), I.end() - I.start() + 1, ref_); // CHECK THAT ref_ DOESNT DO NASTIES
    }
}


const Vector Vector::operator()(const MV_VecIndex &I) const
{
    // check that index is not out of bounds
    //
    if (I.all())
        return Vector(p_, dim_, ref_); // CHECK THAT ref_ DOESNT DO NASTIES
    else
    {
      if ( I.end() >= dim_)
      {
        std::cerr << "MV_VecIndex: (" << I.start() << ":" << I.end() << 
                ") too big for matrix (0:" << dim_ - 1 << ") " << "\n";
        exit(1);
      }
      return Vector(p_+ I.start(), I.end() - I.start() + 1, ref_);  // CHECK THAT ref_ DOESNT DO NASTIES
    }
}


Vector::~Vector()
{
        if (p_ && !ref_ ) delete [] p_;
}


std::ostream&   operator<<(std::ostream& s, const Vector& V)
{
    int N = V.size();

    for (int i=0; i< N; i++)
        s << V(i) << "\n";
    
    return s;
}


