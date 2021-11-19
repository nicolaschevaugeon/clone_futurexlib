/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xDiagonalMatrix_
#define _xDiagonalMatrix_
#include <cmath>

#include "xBlasDef.h"
#include "xCSRVector.h"

namespace xlinalg
{
template <class T>
class xDiagonalMatrix
{
  public:
   xDiagonalMatrix(const int &_nrows, const int &_ncolumns)
       : nrows(_nrows), ncolumns(_ncolumns), ndiag(std::min(nrows, ncolumns)), values(new T[ndiag]), owner(true){};
   xDiagonalMatrix(const int &_nrows, const int &_ncolumns, T *_values)
       : nrows(_nrows), ncolumns(_ncolumns), ndiag(std::min(nrows, ncolumns)), values(_values), owner(false)
   {
   }

   xDiagonalMatrix(const xDiagonalMatrix<double> &M)
       : nrows(M.nrows), ncolumns(M.ncolumns), ndiag(std::min(nrows, ncolumns)), values(new T[ndiag]), owner(true)
   {
      for (int i = 0; i < ndiag; ++i) values[i] = M.values[i];
   };

   ~xDiagonalMatrix()
   {
      if (owner) delete[] values;
   };

   template <class M>
   explicit xDiagonalMatrix(const M &A)
       : nrows(A.nrows), ncolumns(A.ncolumns), ndiag(std::min(nrows, ncolumns)), values(new T[ndiag]), owner(true)
   {
      for (int i = 0; i < ndiag; ++i) values[i] = A(i, i);
   }

   int nrows;
   int ncolumns;
   int ndiag;
   T *values;
   bool owner;
};

void solve(const xDiagonalMatrix<double> &D, xCSRVector &x);
void solve(const xDiagonalMatrix<double> &D, const xCSRVector &b, xCSRVector &x);
void trsv(char *UPLOW, const xDiagonalMatrix<double> &D, xCSRVector &x);

template <class T>
void gemv(const T &alpha, const xDiagonalMatrix<T> &D, const T *x, const int &incx, const T &beta, T *y, const int &incy)
{
   int nrows = D.nrows;
   int ncolumns = D.ncolumns;
   int ndiag = std::min(nrows, ncolumns);
   for (int i = 0; i < ndiag; ++i)
   {
      y[i * incy] = beta * y[i * incy] + alpha * x[incx * i] * D.values[i];
   }
}

void gemv(const double &alpha, const xDiagonalMatrix<double> &D, const xCSRVector &x, const double &beta, xCSRVector &y);

}  // namespace xlinalg
#endif
