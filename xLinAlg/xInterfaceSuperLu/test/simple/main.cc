/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#include <fstream>

#include "xCSRMatrix.h"
#include "xCSRVector.h"
#include "xLinearSystemSolverSuperLU.h"

using namespace xlinalg;
using namespace std;

int main()
{
   /*
   simple test : solve Ax = b
   where A is :
    1   0   2   0
    4   1   0   3
    0   0   1   0
    0   5   0   1
   and b is :
    1
    1
    1
    1
   solution is :
    -1.00000
    -0.14286
     1.00000
     1.71429
   */

   ofstream out("res.txt");
   xCSRMatrix A(4);
   A.AddMatrix(1, 1, 1.);
   A.AddMatrix(1, 3, 2.);
   A.AddMatrix(2, 1, 4.);
   A.AddMatrix(2, 2, 1.);
   A.AddMatrix(2, 4, 3.);
   A.AddMatrix(3, 3, 1.);
   A.AddMatrix(4, 2, 5.);
   A.AddMatrix(4, 4, 1.);

   xCSRVector b(4);
   b.AddVal(1, 1.);
   b.AddVal(2, 1.);
   b.AddVal(3, 1.);
   b.AddVal(4, 1.);

   xCSRVector x(4);
   xCSRVector ref(4);
   ref.AddVal(1, -1.);
   ref.AddVal(2, -0.142857);
   ref.AddVal(3, 1.);
   ref.AddVal(4, 1.71429);

   xLinearSystemSolverSuperLU<> solver;

   // Testing the solve operation
   solver.connectMatrix(A);
   solver.solve(b, x);

   out << "Numeric Solution Is :" << endl;
   x.OutputVector(out);
   out << "Ref Solution Is :" << endl;
   ref.OutputVector(out);

   // test a second solve to make sure that A and b are not modified.
   solver.solve(b, x);
   out << "Numeric Solution Is :" << endl;
   x.OutputVector(out);
   out << "Ref Solution Is :" << endl;
   ref.OutputVector(out);

   // test solve with multiples right-handsides;
   xDenseMatrix B(4, 3);
   B.AddMatrix(1, 1, 1.);
   B.AddMatrix(2, 1, 1.);
   B.AddMatrix(3, 1, 1.);
   B.AddMatrix(4, 1, 1.);

   B.AddMatrix(1, 2, 0.);
   B.AddMatrix(2, 2, 1.);
   B.AddMatrix(3, 2, 0.);
   B.AddMatrix(4, 2, 1.);

   B.AddMatrix(1, 3, 1.);
   B.AddMatrix(2, 3, 1.);
   B.AddMatrix(3, 3, 0.);
   B.AddMatrix(4, 3, 1.);

   xDenseMatrix REF(4, 3);
   REF.AddMatrix(1, 1, -1.);
   REF.AddMatrix(2, 1, -0.14286);
   REF.AddMatrix(3, 1, 1.);
   REF.AddMatrix(4, 1, 1.71429);

   REF.AddMatrix(1, 2, 0.);
   REF.AddMatrix(2, 2, 0.14286);
   REF.AddMatrix(3, 2, 0.);
   REF.AddMatrix(4, 2, 0.28571);

   REF.AddMatrix(1, 3, 1.);
   REF.AddMatrix(2, 3, 0.42857);
   REF.AddMatrix(3, 3, 0.);
   REF.AddMatrix(4, 3, -1.14286);

   xDenseMatrix X(4, 3);
   solver.solve(B, X);

   out << "Numeric Solution Is :" << endl;
   out << X << endl;
   out << "Ref Solution Is :" << endl;
   out << REF << endl;
}
