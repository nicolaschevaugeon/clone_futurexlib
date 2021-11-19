/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include "xLinearSystemSolverDiagonal.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"

xLinearSystemSolverDiagonal::xLinearSystemSolverDiagonal() :
  mat(nullptr) {}

xLinearSystemSolverDiagonal::~xLinearSystemSolverDiagonal() = default;

void xLinearSystemSolverDiagonal::connectMatrix(xlinalg::xCSRMatrix& m)
{
  mat = &m;
}

int xLinearSystemSolverDiagonal::solve(xlinalg::xCSRVector& rhs, xlinalg::xCSRVector& sol)
{
  if (!mat) throw;
  int n=mat->GetNbUnknown();
  for (int i=0; i<n; ++i) {
    sol[i]=rhs[i]/mat->GetMatrix(i+1,i+1);
  }
  return 0;
}
