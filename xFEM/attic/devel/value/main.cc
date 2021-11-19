/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "main.h" 
#include <vector>
#include <cassert>
#include <iostream>
#include "xValue.h"

using namespace xfem;
using namespace std;

int main(int argc, char *argv[])
{ 
  std::ofstream out("value.txt");
  out << " value class example " << endl;
  // v2 = v1
  xValueDouble v1;
  out << "initial value is ";
  v1.print(out);
  out << endl;
  xValueLinearCombination v2(1.0, &v1);
  v1.setVal(2.);
  v1.printVal(out);
  assert( v2.getVal() == v1.getVal() );
  assert( v2.getVal() == 2.);
  
  // v3 = 2 * v1 + 3
  xValueLinearCombination v3(2.0, &v1, 3.);
  assert(v3.getVal() == 2. * v1.getVal() + 3.);
  
  // v4 = v1 + 2 * v2 + 3 * v3 + 5
  std::vector<double> coeffs;
  coeffs.push_back(1.0);
  coeffs.push_back(2.0);
  coeffs.push_back(3.0);
  std::vector<xValue<double>*> vals;
  vals.push_back(&v1);
  vals.push_back(&v2);
  vals.push_back(&v3);
  
  
  xValueLinearCombination v4(coeffs, vals, 5.0);
  assert(v4.getVal() == v1.getVal() + 2. * v2.getVal() + 3. * v3.getVal() + 5.);
  return 0; 

}

