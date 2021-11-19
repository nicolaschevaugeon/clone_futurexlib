/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _XTRAITPOLICYLINEARSYSTEMSOLVERSUPERLU_H
#define _XTRAITPOLICYLINEARSYSTEMSOLVERSUPERLU_H


#include "xLinearSystemSolverSuperLUDataType.h"
#include "xCSRVector.h"
#include "xDenseMatrix.h"
#include "xTraitsMatrix.h"

#ifndef  _SOLVERSUPERLU_HEADER_SET 
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=                                                                                                                                         ==
//=   Begin dangerous part (see main header for coment)                                                                                     ==
//=                                                                                                                                         ==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
#include  "xLinearSystemSolverSuperLUcplx.h"
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=                                                                                                                                         ==
//=   End dangerous part                                                                                                                    ==
//=                                                                                                                                         ==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
#endif

namespace xlinalg
{

// traits for xLinearSystemSolverSuperLU class
template < typename T >
class xLinearSystemSolverSuperLUSubType ;

template < >
class xLinearSystemSolverSuperLUSubType<float>
{
    public:
        typedef float sub_type;
};

template < >
class xLinearSystemSolverSuperLUSubType<double>
{
    public:
        typedef double sub_type;
};

#ifdef CSUPERLU
template < >
class xLinearSystemSolverSuperLUSubType<superLuCmplx::complex>
{
    public:
        typedef float sub_type;
};
#endif
#ifdef ZSUPERLU
template < >
class xLinearSystemSolverSuperLUSubType<doublecomplex>
{
    public:
        typedef double sub_type;
};
#endif


template < typename PATTERN, typename DEFINED >
class xTraitsLinearSystemSolverSuperLUMatrixType;

template < >
class xTraitsLinearSystemSolverSuperLUMatrixType < xTraitMatrixUnSymPSym,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 0;
};

template < >
class xTraitsLinearSystemSolverSuperLUMatrixType < xTraitMatrixUnSymPSym,xTraitMatrixDefinitePositive>
{
    public:
        static const int sym = 0;
};


template < >
class xTraitsLinearSystemSolverSuperLUMatrixType < xTraitMatrixUnSym,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 0;
};

template < >
class xTraitsLinearSystemSolverSuperLUMatrixType < xTraitMatrixUnSym,xTraitMatrixDefinitePositive >
{
    public:
        static const int sym = 0;
};

template < >
class xTraitsLinearSystemSolverSuperLUMatrixType < xTraitMatrixDiagonal,xTraitMatrixNonSingular >
{
    public:
        static const int sym = 0;
};

template < >
class xTraitsLinearSystemSolverSuperLUMatrixType < xTraitMatrixDiagonal,xTraitMatrixDefinitePositive >
{
    public:
        static const int sym = 0;
};

// Policies for xLinearSystemSolverSuperLU class
template < typename T, typename V >
class xPolicyLinearSystemSolverSuperLURHS;

template < >
class xPolicyLinearSystemSolverSuperLURHS < double,xCSRVector >
{
    public:
        static void * copyBInX(const xCSRVector & B, xCSRVector & X)
        {
            X=B;
            return ( (void *) X.GetArray() );
        }
        static double * getPointer( xCSRVector & X)
        {
            return ( X.GetArray() );
        }
        static int nbColRhs(const xCSRVector & X)
        {
            return( 1 );
        }
};

template < >
class xPolicyLinearSystemSolverSuperLURHS < double,xDenseMatrix >
{
    public:
        static void * copyBInX(const xDenseMatrix & B, xDenseMatrix & X)
        {
            X=B;
            double *x = X.ReturnValPointer();
            return ( (void *) x );
        }
        static double * getPointer( xDenseMatrix & X)
        {
            return (  X.ReturnValPointer() );
        }
        static int nbColRhs(const xDenseMatrix & X)
        {
            return( X.ncolumn());
        }
};

} // end namespace

#endif
