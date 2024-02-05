#ifndef LUINCPRECONDITIONER_H
#define LUINCPRECONDITIONER_H

#include <spamtrix_preconditioner.hpp>
#include <spamtrix_fleximatrix.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>

namespace SpaMtrix {
    class LUIncPreconditioner: public Preconditioner {
        FlexiMatrix M;  // stores both L and U. L is unit lower diagonal, so its
                    // diagonal values are not stored
        LUIncPreconditioner(){}

        /**
            Performs forward substitution on the lower diagonal part of the incomplete LU factorisation (L).<br>

            L is a lower unit diagonal matrix, stored in the lower part of matrix 'a'.
            The diagonal entries of L are not explicitly stored, since they are all ones.<br>
            <br>

            x[i] = { b[i] - sum( L[i,j]*x[j] ) } / L[i,j]; <br>
            simplifies to<br>
            x[i] = b[i] - sum( L[i,j]*x[j] ) <br>
        */
        void forwardSubstitution(Vector &x, const Vector &b) const;

        /**
            Performs back-substitution of upper diagonal part of the incomplete LU factorisation (U).
            U is the upper diagonal matrix, stored in the upper diagonal part of 'a'.
        */
        void backwardSubstitution(Vector &x, const Vector &b) const;
    public:
        /**
        * Creates an incomplere LU preconditioner with zero drop tolerance, LU(0).
        * The implementation is based on the IKJ version of the ILU factorisation
        * in Saad's book (ALGORITHM 10.3: General ILU Factorization, IKJ Version)
        */
        LUIncPreconditioner(const IRCMatrix &A);
        virtual ~LUIncPreconditioner();
        void print() const;
        void solveMxb(Vector &x, const Vector &b) const;
    };
}
#endif
