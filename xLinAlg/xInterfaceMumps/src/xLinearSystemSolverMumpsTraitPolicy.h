/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XTRAITPOLICYLINEARSYSTEMSOLVERMUMPS_H
#define _XTRAITPOLICYLINEARSYSTEMSOLVERMUMPS_H


#include <cstring>
#include <complex>
#include "xCSRVector.h"
#include "xDistVector.h"
#include "xDenseMatrix.h"
#include "xTraitsMatrix.h"


namespace xlinalg
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
class xTraitsLinearSystemSolverMumpsMatrixType < xTraitMatrixUnSymPSym,xTraitMatrixDefinitePositive >
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
class xPolicyLinearSystemSolverMumpsRHSOneType;

template < >
class xPolicyLinearSystemSolverMumpsRHSOneType < double,xCSRVector >
{
    public:
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

template < typename T >
class xPolicyLinearSystemSolverMumpsRHSOneType < T,std::vector<T> >
{
    public:
        static int nbColRhs(const std::vector<T> & X)
        {
            return( 1 );
        }
        static int nbRowRhs(const std::vector<T> & X)
        {
            return( X.size() );
        }
        static void * beginRhs(const std::vector<T> & X)
        {
            return( (void *) X.data() );
        }
};


template < typename T >
class xPolicyLinearSystemSolverMumpsRHSOneType < T,xDistVector<T> >
{
    public:
        static int nbColRhs(const xDistVector<T> & X)
        {
            return( 1 );
        }
        // no shur stuff for now, maybe never
        static int nbRowRhs(const xDistVector<T> & X)
        {
            throw -1;
            return 0;
        }
        static void * beginRhs(const xDistVector<T> & X)
        {
            throw -1;
            return nullptr;
        }
};



template < >
class xPolicyLinearSystemSolverMumpsRHSOneType < double,xDenseMatrix >
{
    public:
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

// -----

template < typename T, typename V, typename U >
class xPolicyLinearSystemSolverMumpsRHSTwoType;

template < >
class xPolicyLinearSystemSolverMumpsRHSTwoType < double,xCSRVector, xCSRVector >
{
    public:
        static void * copyBInX(const xCSRVector & B, xCSRVector & X, std::vector < double > & rhs_gather)
        {
            std::copy(B.begin(), B.end(), X.begin());
            return ( (void *) X.GetArray() );
        }
        static void  brodcast( xCSRVector & X, std::vector < double > & rhs_gather, int n, bool do_brd, MPI_Datatype data_type, MPI_Comm univ)
        {
            if (do_brd)
                MPI_Bcast( (void *) X.GetArray(), n, data_type, 0, univ);
        }
};
template < typename T >
class xPolicyLinearSystemSolverMumpsRHSTwoType < T,xDistVector<T>, std::vector < T > >
{
    public:
        static void * copyBInX(const xDistVector<T> & B, std::vector < T > & X, std::vector < T > & rhs_gather)
        {
            B.gather(X,0);
            return ( (void *) X.data() );
        }
        static void  brodcast( std::vector < T > & X, std::vector < T > & rhs_gather, int n, bool do_brd, MPI_Datatype data_type, MPI_Comm univ)
        {
            if (do_brd)
                MPI_Bcast( (void *) X.data(), n, data_type, 0, univ);
        }
};

template < typename T >
class xPolicyLinearSystemSolverMumpsRHSTwoType < T,xDistVector<T>, xDistVector<T> >
{
    public:
        static void * copyBInX(const xDistVector<T> & B, xDistVector<T> & X, std::vector < T > & rhs_gather)
        {
            B.gather(rhs_gather,0);
            return ( (void *) rhs_gather.data() );
        }
        static void  brodcast( xDistVector<T> & X, std::vector < T > & rhs_gather, int n, bool do_brd, MPI_Datatype data_type, MPI_Comm univ)
        {
            if (do_brd)
                X.scatter(rhs_gather,0);
            else
                throw -1; // As X is distributed if you do not brodcast you simply lose the result stored in rhs_gather !
        }
};




template < >
class xPolicyLinearSystemSolverMumpsRHSTwoType < double,xDenseMatrix, xDenseMatrix >
{
    public:
        static void * copyBInX(const xDenseMatrix & B, xDenseMatrix & X, std::vector < double > & rhs_gather)
        {
            const int nbl = B.nline();
            const int nbc = B.ncolumn();
            const double *b = B.ReturnValPointer();
            double *x = X.ReturnValPointer();
            memcpy ( (void *) ( x ), (void *) ( b ), sizeof( double )*nbl*nbc );
            return ( (void *) x );
        }
        static void  brodcast( xDenseMatrix & X, std::vector < double > & rhs_gather, int n, MPI_Datatype data_type, MPI_Comm univ)
        {
            MPI_Bcast( (void *) X.ReturnValPointer(), n*X.ncolumn(), data_type, 0, univ);
        }
};
} // end namespace

#endif
