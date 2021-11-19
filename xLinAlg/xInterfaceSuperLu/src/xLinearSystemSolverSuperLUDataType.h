/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _SOLVERSUPERLUDATATYPE_H
#define _SOLVERSUPERLUDATATYPE_H

/// chose your type here
//! for now include of SuperLU pacakage are mutualy exclusive as library might be ! :-(
//! only one to be uncoment
//
//#define SSUPERLU 1
#define DSUPERLU 1
//#define CSUPERLU 1
//#define ZSUPERLU 1

/*
 * equivalence :
 * SSUPERLU => T=float
 * DSUPERLU => T=double
 * CSUPERLU => T=superLuCmplx::complex with following definition of complex : namespace superLuCmplx { typedef struct { float r, i; } complex; } 
 * ZSUPERLU => T=doublecomplex with following definition of doublecomplex : typedef struct { double r, i; } doublecomplex;
 *
 * those complex type are explicitely defined in xLinearSystemSolverSuperLU.h
 *
 * use of namespace superLuCmplx avoid ambiguous ‘complex’ reference with 'template<class _Tp> struct std::complex'
 *
 * should work, only tested with T=double for now
 * interface compile in all case but never try to make a executable
 *
 */

#endif
