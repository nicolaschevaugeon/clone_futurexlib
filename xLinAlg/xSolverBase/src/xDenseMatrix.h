/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _DENSE_MATRIX__H
#define _DENSE_MATRIX__H

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include "xCSRVector.h"
#include "xTraitsMatrix.h"

// used blas functions

namespace xlinalg
{
/// xDenseMatrix (Column Major Storage)
class xDenseMatrix
{
  public:
   /// Fixed traits for assemble/solver
   typedef xTraitMatrixUnSym matrix_pattern;
   typedef xTraitMatrixDense matrix_storage;
   typedef xTraitMatrixNonSingular matrix_defined;
   typedef xTraitMatrixNoAssemblyOnZero matrix_assembly_on_zero;
   typedef xTraitMatrixCindex matrix_indexing;

   /// Construct an empty  matrix (nl=nc =0)
   xDenseMatrix();
   /// Construct a dense nl*nc matrix
   xDenseMatrix(const int nl, const int nc);
   /// Construct a dense nl*nl matrix
   xDenseMatrix(const int nl);
   /// Copy Constructor
   xDenseMatrix(const xDenseMatrix &in);
   /// Constructor from stream
   xDenseMatrix(std::istream &, int options = 0);
   /// Destructor
   ~xDenseMatrix();
   /// resize the matrix (old data are lost, initialization with v_init defaulted to zero)
   void resize(const int nl, const int nc, double v_init = 0.);
   /// add the Matrix A to the current matrix.
   /*! throw an exception if current and A don't have the same dimension !*/
   void operator+=(const xDenseMatrix &A);
   /// substract the Matrix A to the current matrix.
   /*! throw an exception if current and A don't have the same dimension !*/
   void operator-=(const xDenseMatrix &A);
   xDenseMatrix &operator=(const xDenseMatrix &A);

   /// return a read only reference to element i, j of the matrix
   /*!
     (i, j indices start at 0 at go up to nl-1 et nc-1 )
     No check are done on the boundary !
   !*/
   const double &operator()(const int i, const int j) const;
   /// return a reference to element i, j of the matrix
   /*!
     (i, j indices start at 0 at go up to nl-1 et nc-1 )
     No check are done on the boundary !
   !*/
   double &operator()(const int i, const int j);
   /// return the Number of lines of the matrix
   int nline() const { return nl; };
   /// return the Number of column of the matrix
   int ncolumn() const { return nc; };
   /// add val to element i j of the matrix (follow fortran convention : first element is i =1, j=1)
   /*! This function is there for compatibility with xfem Assemblers. !*/
   void AddMatrix(int i, int j, double val);

   /// add val to element i,col_ass of the matrix
   /*! This function is there for compatibility with xfem Assemblers for rhs  and give the opportunity
    *  to use xDenseMatrix as a real container for muti-right hand side
    *  It follow fortran convention for row  : first element is i =1
    *  It follow C convention for colonne  : first element is j =0
    *  colone is drived by menber col_ass
    *  !*/
   void AddVal(const int &i, const double &val);

   /// set colone where xfem Assemblers will store values
   void setAssCol(const int col_num);

   /// return a copy of element i, j of the matrix. (follow fortran convention : first element is i =1, j=1)
   double GetMatrix(int i, int j) const;
   /// return a pointer to the array of element.
   double *ReturnValPointer() { return a; };
   /// const version of ReturnValPointer. User can't change the value in the returned array
   const double *ReturnValPointer() const { return a; };
   /// print to screen the matrix
   void print2screen() const;
   /// export the matrix to matlab
   void export2Matlab(std::string fileName, std::string variableName = "M") const;
   void export2Matlab(std::ostream &out, std::string variableName = "M") const;
   /// export the matrix to a text file
   /*! precision set the number of significant digit in the output.
     option =0 : just output the components of the matrix.
     option = 1: first output the number of line and columns, so that the file is
     readeable by the constructor that take an input strema
    */
   void export2txt(std::string fileName, unsigned precision = 5, int option = 0) const;

   void EndOfAssembly();
   //// return nb lines of the matrix
   int GetNbUnknown() const { return nl; };
   void operator*=(const double &alpha);
   /// return the maxnorm of the matrix : return the absolute value of the element having the max absolute value
   double maxnorm() const;
   friend xDenseMatrix operator*(const xDenseMatrix &A, const xDenseMatrix &B);
   friend xCSRVector operator*(const xDenseMatrix &A, const xCSRVector &X);
   friend std::ostream &operator<<(std::ostream &output, const xDenseMatrix &A);

  private:
   int nl, nc;
   double *a;
   int col_ass;
};

/// Compute y = alpha*A.x +beta*y
void gemv(bool TRANSA, double ALPHA, const xDenseMatrix &A, const xCSRVector &x, double BETA, xCSRVector &y);
/// Compute C= alpha*A.B + beta*C
void gemm(bool TRANSA, bool TRANSB, double ALPHA, const xDenseMatrix &A, const xDenseMatrix &B, double BETA, xDenseMatrix &C);

void syrk(const char &UPLO, const bool &TRANS, const double &ALPHA, const xDenseMatrix &A, const double &BETA, xDenseMatrix &C);

/// return C = A*B , A, B et C xDenseMatrix
/*! Here for convenience ... use gemm in priority : this version will add the cost of a temporary for C, unless optimized out by
 * the compiler !*/
xDenseMatrix operator*(const xDenseMatrix &A, const xDenseMatrix &B);
// compute B = alpha*A +B;
void axpy(const double &alpha, const xDenseMatrix &A, xDenseMatrix &B);
// compute B = transpose(A);
void transpose(const xDenseMatrix &A, xDenseMatrix &B);

/// return y = A*x A xDenseMatrix, x, y xCSRVector.
/*! Here for convenience ... use gemv in priority : this version will add the cost of a temporary for y, unless optimized out by
 * the compiler !*/
xCSRVector operator*(const xDenseMatrix &A, const xCSRVector &x);

std::ostream &operator<<(std::ostream &output, const xDenseMatrix &A);

}  // namespace xlinalg
#endif
