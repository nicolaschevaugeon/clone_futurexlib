/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef _SOLVER_DIAGONAL_H
#define _SOLVER_DIAGONAL_H

namespace xlinalg
{
class xCSRMatrix;
class xCSRVector;
}

class xLinearSystemSolverDiagonal
{
public:
  xLinearSystemSolverDiagonal();
  ~xLinearSystemSolverDiagonal();
  void connectMatrix(xlinalg::xCSRMatrix& m);
  int solve(xlinalg::xCSRVector& rhs, xlinalg::xCSRVector& sol);
private:
  xlinalg::xCSRMatrix* mat;
};

#endif





