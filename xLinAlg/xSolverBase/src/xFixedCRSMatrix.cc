/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xFixedCRSMatrix.h"
namespace xlinalg
{
void graphToCRS(std::vector<std::set<size_t>> &graph, int *columns, int *rowindex)
{
   int nlines = graph.size();
   rowindex[0] = 0;
   for (int i = 0; i < nlines; ++i)
   {
      rowindex[i + 1] = rowindex[i] + graph[i].size();
      std::copy(graph[i].begin(), graph[i].end(), &columns[rowindex[i]]);
   }
}

void gemv(const double &alpha, const xFixedCRSMatrix<double> &A, const xCSRVector &x, const double &beta, xCSRVector &y)
{
   gemv<double>(alpha, A, &x[0], beta, &y[0]);
   return;
}

void trsv(const char *UPLOW, const xFixedCRSMatrix<double> &A, xCSRVector &x) { xlinalg::trsv<double>(UPLOW, A, &x[0], 1); }

xCSRVector operator*(const xFixedCRSMatrix<double> &A, xCSRVector &x)
{
   xCSRVector y(A.nrows);
   gemv<double>(1, A, &x[0], 1, 0., &y[0], 0);
   return y;
}

}  // end namespace xlinalg
