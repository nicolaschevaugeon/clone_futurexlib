/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _TRAITSMATRIX_H
#define _TRAITSMATRIX_H

// pattern type
struct xTraitMatrixLowerSym { };
struct xTraitMatrixUpperSym { };
struct xTraitMatrixUnSym { };
struct xTraitMatrixUnSymPSym { };
struct xTraitMatrixDiagonal { };
struct xTraitMatrixBande { };
// fixed or variable structure : zero saved or not
struct xTraitMatrixForcedAssemblyOnZero { };
struct xTraitMatrixNoAssemblyOnZero { };
// Storage type
struct xTraitMatrixSparceCSC { };
struct xTraitMatrixSparceDCSC { };
struct xTraitMatrixSparceCSR { };
struct xTraitMatrixSparceDCSR { };
struct xTraitMatrixSparceCOO { };
struct xTraitMatrixSparceDCOO { };
struct xTraitMatrixDense { };
// Defined type
struct xTraitMatrixDefinitePositive { };
struct xTraitMatrixNonSingular { };
// c/fortran indexation
struct xTraitMatrixCindex { };
struct xTraitMatrixFindex { };


#endif
