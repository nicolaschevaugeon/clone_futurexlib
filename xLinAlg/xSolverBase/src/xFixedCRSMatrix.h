/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xFixedCRSMatrix_
#define _xFixedCRSMatrix_
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <set>
#include <vector>

#include "xDiagonalMatrix.h"
#include "xTraitsMatrix.h"

#ifdef HAVE_MKL
#include "mkl_spblas.h"
#endif

namespace xlinalg
{
//! An implementation of the compressed row storage format.
template <class T>
class xFixedCRSMatrix
{
  public:
   //! fixed traits for assemble/solver
   // arbitrary chosen as a unsymetric container with non zero termes all presents (fixed graph)
   // naturaly as it's name says  it, the Storage type is compressed row storage format= Compressed Sparce Row
   typedef xTraitMatrixUnSym matrix_pattern;
   typedef xTraitMatrixSparceCSR matrix_storage;
   typedef xTraitMatrixNonSingular matrix_defined;
   typedef xTraitMatrixNoAssemblyOnZero matrix_assembly_on_zero;

   //! When calling this constructor, memory is alocated for the user.
   xFixedCRSMatrix(const int &_nrows, const int &_ncolumns, const int &_nnz)
       : nrows(_nrows),
         ncolumns(_ncolumns),
         nnz(_nnz),
         values(new T[nnz]),
         columns(new int[nnz]),
         rowindex(new int[nrows + 1]),
         owner(true)
   {
      rowindex[0] = 0;
      rowindex[nrows] = nnz;
   };

   //! For this constructor, the data are supposed to be allocated by the user, and directly given, to the constructor
   xFixedCRSMatrix(const int &_nrows, const int &_ncolumns, const int &_nnz, T *_values, int *_columns, int *_rowindex)
       : nrows(_nrows), ncolumns(_ncolumns), nnz(_nnz), values(_values), columns(_columns), rowindex(_rowindex), owner(false)
   {
      rowindex[0] = 0;
      rowindex[nrows] = nnz;
   };

   ~xFixedCRSMatrix()
   {
      if (owner)
      {
         delete[] values;
         delete[] columns;
         delete[] rowindex;
      }
   }

   //! indexing the matrix, const version
   const T &operator()(const int &i, const int &j) const { return (const_cast<xFixedCRSMatrix<T> *>(this))->operator()(i, j); }

   //! indexing of the matrix
   T &operator()(const int &i, const int &j)
   {
      const int istart = rowindex[i];
      const int iend = rowindex[i + 1];
      int *pos = std::lower_bound(&columns[istart], &columns[iend], j);
      assert(pos != &columns[iend] && (*pos == j));
      int k = std::distance(columns, pos);
      return values[k];
   }

   //! xfem style AddMatrix operation (i and j are given starting at one (fortran convention))
   void AddMatrix(int i, int j, T val) { this->operator()(i - 1, j - 1) += val; }

   //! output the non zero term of the matrix on screen. Row by row, along with the index of each term
   void print()
   {
      for (int i = 0; i < nrows; ++i)
      {
         std::cout << i << std::endl;
         int starti = rowindex[i];
         int endi = rowindex[i + 1];
         for (int k = starti; k < endi; ++k)
         {
            std::cout << "A(" << i << ", " << columns[k] << ")= " << values[k] << std::endl;
         }
      }
   };

  public:
   int nrows;
   int ncolumns;
   int nnz;
   T *values;
   int *columns;
   int *rowindex;
   bool owner;
};

//! implementation of a blas style Matrix vector multiplication for xFixedCRSMatrix, simplifyed interface
template <class T>
void gemv(const T &alpha, const xFixedCRSMatrix<T> &A, const T *x, const T &beta, T *y)
{
   // note :  the mkl version of sparse mv product does not seems any better than mine ...
   // At least on titan ...
#ifdef HAVE_MKL
   mkl_dcsrmv("N", const_cast<int *>(&A.nrows), const_cast<int *>(&A.ncolumns), const_cast<double *>(&alpha), "GXXC",
              //"SUNC",  //if matrix is Symmetric, accessing only Up, diagonal Non constant, Cindexing
              const_cast<double *>(A.values), const_cast<int *>(A.columns), const_cast<int *>(&A.rowindex[0]),
              const_cast<int *>(&A.rowindex[1]), const_cast<double *>(&x[0]), const_cast<double *>(&beta), &y[0]);

#else
   gemv(alpha, A, x, 1, beta, y, 1);
#endif
}

void gemv(const double &alpha, const xlinalg::xFixedCRSMatrix<double> &A, const xlinalg::xCSRVector &x, const double &beta,
          xlinalg::xCSRVector &y);

//! implementation of a blas style Matrix vector multiplication for xFixedCRSMatrix.
template <class T>
void gemv(const T &alpha, const xFixedCRSMatrix<T> &A, const T *x, const int &incx, const T &beta, T *y, const int &incy)
{
#define CHUNKSIZE 100
   int nrows = A.nrows;
   int i, j, jb, je, col;
   T tmp;
#pragma omp parallel num_threads(4) private(j, jb, je, col, tmp)
   {
#pragma omp for schedule(dynamic, CHUNKSIZE)
      for (i = 0; i < nrows; ++i)
      {
         tmp = 0;
         jb = A.rowindex[i];
         je = A.rowindex[i + 1];
         for (j = jb; j < je; ++j)
         {
            col = A.columns[j];
            tmp += A.values[j] * x[col * incx];
         }
         y[i * incy] = tmp * alpha + beta * y[i * incy];
      }
   }
}

template <class T>
void trsv(const char *UPLOW, const xFixedCRSMatrix<T> &A, T *x, const int &incx)
{
   if (UPLOW[0] == 'U')
      trsvUp(A, x, incx);
   else if (UPLOW[0] == 'L')
      trsvLow(A, x, incx);
   else
      throw;
}

//! Rransform a graph of the non zero term of a matrix into the array columns and rowindex using the convention of the CRS format.
/*! The graph on entry is suupposed to be an std::vector. each line of the vector, correspond to a line of the matrix.
    Each entry i of the vector contain on std::set of int, corresponding to the column of the non zero term of the line i of the
  matrix. the array columns is supposed to have at least size  nnz where nnz is the total number of non zero term in tha Matrix.
    the array rowindex is supposed to have a size n+1, where n is the number of line of the matrix.
  !*/
void graphToCRS(std::vector<std::set<size_t>> &graph, int *columns, int *rowindex);

/// inout values rows columns out rowindex;
template <class T>
void sortCoordinateInPlace(int nrows, int nnz, T *values, int *rows, int *columns, int *rowindex)
{
   int pos = 0;
   int nsorted = 0;
   int line = 0;
   rowindex[0] = 0;
   rowindex[nrows] = nnz;
   // sorting by row ...
   while (nsorted < nnz)
   {
      int *current = std::find(&rows[nsorted], &rows[nnz], line);
      std::cout << current << " " << &rows[nsorted] << " " << &rows[nnz] << std::endl;
      if (current == &rows[nnz])
      {
         line += 1;
         rowindex[line] = nsorted;
      }
      else
      {
         int k = std::distance(rows, current);
         std::swap(rows[nsorted], rows[k]);
         std::swap(columns[nsorted], columns[k]);
         std::swap(values[nsorted], values[k]);
         nsorted += 1;
      }
   }
   // sort each row
   for (int i = 0; i < nrows; ++i)
   {
      int starti = rowindex[i];
      int endi = rowindex[i + 1];
      int online = endi - starti;
      int sorted = 0;
      while (sorted < online)
      {
         int *jmin = std::min_element(&columns[starti + sorted], &columns[endi]);
         int k = std::distance(columns, jmin);
         std::swap(columns[k], columns[starti + sorted]);
         std::swap(values[k], values[starti + sorted]);
         ++sorted;
      }
   }
}

// void trsv(uplo,trans, diag, N, A, LDA, x incx)
void trsv(const char *UPLOW, const xFixedCRSMatrix<double> &A, xCSRVector &x);

template <class T>
void trsvLow(const xFixedCRSMatrix<T> &A, T *x, const int &incx)
{
   int nrows = A.nrows;
   for (int i = 0; i < nrows; ++i)
   {
      T tmp = 0;
      int j = A.rowindex[i];
      int je = A.rowindex[i + 1];
      int col = A.columns[j];
      while ((j < je) && (col < i))
      {
         tmp += A.values[j] * x[col * incx];
         ++j;
         col = A.columns[j];
      }
      assert(col == i);
      x[col * incx] = (x[col * incx] - tmp) / (A.values[j]);
   }
}

template <class T>
void trsvUp(const xFixedCRSMatrix<T> &A, T *x, const int &incx)
{
   int nrows = A.nrows;
   for (int i = nrows - 1; i > -1; --i)
   {
      T tmp = 0;
      int jb = A.rowindex[i];
      int j = A.rowindex[i + 1] - 1;
      int col = A.columns[j];
      while ((j >= jb) && (col > i))
      {
         tmp += A.values[j] * x[col * incx];
         --j;
         col = A.columns[j];
      }
      assert(col == i);
      x[col * incx] = (x[col * incx] - tmp) / A.values[j];
   }
}

xCSRVector operator*(const xFixedCRSMatrix<double> &A, xCSRVector &x);
}  // namespace xlinalg
#endif
