/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "xLapackInterface.h"
#include "xDenseMatrix.h"

using std::copy;
using std::cout;
using std::endl;
using std::ofstream;
using namespace xlinalg;




typedef std::ostream_iterator<double , char, 
                         std::char_traits<char> >     os_iter;

int main(int argc, char * argv[]){
  std::cout << "GGGGG"<< std::endl;
  xDenseMatrix a(3);
  xCSRVector v0(3);
  xCSRVector v1(3);
  xCSRVector v2(3);
  xCSRVector x(3);
  v0(0) =1.;
  v1(1) =1.;
  v2(2) =1.;

  a(0,0) = 1.;
  a(1,1) = 1.;
  a.AddMatrix(3,3, 1.);
  a(0,1) = -1.;
  a(1,0) = -1.;
  a(1,2) = -1.;
  a(2,1) = -1.;
  ofstream out("result_test");
  std::cout << "Before SOLVE " << std::endl;

  solve(a, v0, x );
  /*
  std::copy (x.begin(), x.end(), os_iter (out, " ") );
  out << endl;  
  solve(a, v1, x );
  std::copy (x.begin(), x.end(), os_iter (out, " ") );
  out << endl;  
  solve(a, v2, x );
  std::copy (x.begin(), x.end(), os_iter (out, " ") );
  out << endl << endl;
  

  xDenseLUFactor lua(a);
  solve(lua, v0, x);
  std::copy (x.begin(), x.end(), os_iter (out, " ") );
  out << endl;
  solve(lua, v1, x );
  std::copy (x.begin(), x.end(), os_iter (out, " ") ); 
  out << endl;
  solve(lua, v2, x );
  std::copy (x.begin(), x.end(), os_iter (out, " ") );
  out << endl;
*/

}
