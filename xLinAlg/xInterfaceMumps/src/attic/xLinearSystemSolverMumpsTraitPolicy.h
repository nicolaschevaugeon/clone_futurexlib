/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XTRAITPOLICYLINEARSYSTEMSOLVERMUMPS_H
#define _XTRAITPOLICYLINEARSYSTEMSOLVERMUMPS_H


#include <cstring>
#include "xCSRVector.h"
#include "xDenseMatrix.h"
#include "xTraitsMatrix.h"


namespace lalg
{

// traits for xLinearSystemSolverMumps class
template < typename PATTERN, typename DEFINED >
class xTraitsLinearSystemSolverMumpsMatrixType;

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixUnSymPSym,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 0;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixUnSymPSym,xTraitMatrixDefinitePositive>
{
    public:
        static const int sym = 0;
};


template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixUnSym,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 0;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixUnSym,xTraitMatrixDefinitePositive >
{
    public:
        static const int sym = 0;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixDiagonal,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 2;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixDiagonal,xTraitMatrixDefinitePositive >
{
    public:
        static const int sym = 1;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixLowerSym,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 2;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixLowerSym,xTraitMatrixDefinitePositive >
{
    public:
        static const int sym = 1;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixUpperSym,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 2;
};

template < >
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixUpperSym,xTraitMatrixDefinitePositive >
{
    public:
        static const int sym = 1;
};

// Policies for xLinearSystemSolverMumps class
template < typename T, typename V >
class xPolicyLinearSystemSolverMumpsRHS;

template < >
class xPolicyLinearSystemSolverMumpsRHS < double,xCSRVector >
{
    public:
        static void * copyBInX(const xCSRVector & B, xCSRVector & X)
        {
            std::copy(B.begin(), B.end(), X.begin());
            return ( (void *) X.GetArray() );
        }
        static int nbColRhs(const xCSRVector & X)
        {
            return( 1 );
        }
        static int nbRowRhs(const xCSRVector & X)
        {
            return( X.size() );
        }
        static void * beginRhs(const xCSRVector & X)
        {
            return( (void *) X.GetArray() );
        }
};

template < >
class xPolicyLinearSystemSolverMumpsRHS < double,xDenseMatrix >
{
    public:
        static void * copyBInX(const xDenseMatrix & B, xDenseMatrix & X)
        {
            const int nbl = B.nline();
            const int nbc = B.ncolumn();
            const double *b = B.ReturnValPointer();
            double *x = X.ReturnValPointer();
            memcpy ( (void *) ( x ), (void *) ( b ), sizeof( double )*nbl*nbc );
            return ( (void *) x );
        }
        static int nbColRhs(const xDenseMatrix & X)
        {
            return( X.ncolumn());
        }
        static int nbRowRhs(const xDenseMatrix & X)
        {
            return( X.nline());
        }
        static void * beginRhs(const xDenseMatrix & X)
        {
            return( (void *) X.ReturnValPointer() );
        }
};

} // end namespace

#endif
