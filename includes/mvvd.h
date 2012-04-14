
#ifndef VECTOR_H
#define VECTOR_H



#include <stdlib.h>

// for formatted printing of matrices
//#include <sstream>


#ifdef MV_VECTOR_BOUNDS_CHECK
#   include <assert.h>
#endif

#include "mvvind.h" // INDEX VECTOR - THIS SEEMS SUPRFLUOUS REMOVE IF NOT NEEDED

// this is really used as a sort of global constant. The reason
// for creating its own type is that so it can be overloaded to perform
// a deep or shallow assignement.  (Any variable of type MV_Vector_::ref_type
// has only one possible value: one.)
//   It is included as a seperate file to avoid multiple definitions.

//#include "mvvrf.h"

class Vector
{                                                                      
    private:
           double *p_;
           int dim_;
           bool ref_;  // 0 or 1; does this own its own memory space?
    public:                                                            


        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    Vector();
    Vector( int);
    Vector( int, const double&);
                                                     
                                                    
    Vector(double*,  int);
    Vector(const double*,  int);
    Vector(const Vector &);
    // Eero's changes
    void point_to(double* p, int len);
    
    
    // reference of an exisiting data structure
    //
    // note that ref() is initalized with i rather than 1.
    // this is so compilers will not generate a warning that i was
    // not used in the construction.  (MV_Vector::ref_type is an enum that
    // can *only* have the value of 1.
    //
    
    Vector(double* d,  int N, bool ref) :
                            p_(d), 
                            dim_(N), 
                            ref_( ref ) {}

    Vector(const Vector &V, bool ref)   :
                            p_(V.p_), 
                            dim_(V.dim_), 
                            ref_(ref ) {}

    ~Vector();
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       

    double& operator()( int i)
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }
    const  double& operator()( int i) const 
    {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
    }

    double&      operator[]( int i)
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }
    const  double&       operator[]( int i) const 
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }


    Vector operator()(const MV_VecIndex &I) ;
    Vector operator()(void);
    const Vector operator()(void) const;
    const Vector operator()(const MV_VecIndex &I) const;

    inline  int             size() const { return dim_;}
    inline  int             dim() const { return dim_;}
    inline int                      ref() const { return  ref_;}
    inline int                      null() const {return dim_== 0;}
            //
            // Create a new *uninitalized* vector of size N
            Vector & newsize( int );
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    Vector & operator=(const Vector&);
    Vector & operator=(const double&);


    friend std::ostream& operator<<(std::ostream &s, const Vector &A);

};                                                                     



#endif  
