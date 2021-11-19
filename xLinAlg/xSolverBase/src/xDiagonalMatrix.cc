/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xDiagonalMatrix.h"
namespace xlinalg
{
void solve(const xDiagonalMatrix<double> &D, xCSRVector &x)
{
   int N = x.size();
   int K = 0;  // number of subdiagonals
   int LDA = 1;
   int INCX = 1;
   dtbsv_("U", "N", "N", &N, &K, const_cast<double *>(&(D.values[0])), &LDA, &x[0], &INCX);
}

void solve(const xDiagonalMatrix<double> &D, const xCSRVector &b, xCSRVector &x)
{
   x = b;
   solve(D, x);
}

void trsv(char *UPLOW, const xDiagonalMatrix<double> &D, xCSRVector &x) { solve(D, x); }

void gemv(const double &alpha, const xDiagonalMatrix<double> &D, const xCSRVector &x, const double &beta, xCSRVector &y)
{
   gemv<double>(alpha, D, &x[0], 1, beta, &y[0], 1);
}
}  // end namespace xlinalg
