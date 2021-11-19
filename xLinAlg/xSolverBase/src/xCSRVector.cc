/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include <iomanip>
#include "xCSRVector.h"
#include "xBlasDef.h"

namespace xlinalg 
{

void xCSRVector::OutputOctaveFormat(const std::string &filename) const
{
  std::ofstream  out(filename.c_str());
  OutputOctaveFormat(out);
  out.close();
}

void xCSRVector::OutputOctaveFormat(std::ostream &out) const
{
    using namespace std;
    double val;
    int i;
    const int Asize = size();
    out << setiosflags(ios_base::scientific) ;
    int p=out.precision(20);
    for(i=0;i<Asize;i++)
    {
      val = Array[i];
      out  << val << " "; 
      out << std::endl;
    }
    out << resetiosflags(ios_base::scientific) ;
    out.precision(p);
    return;
}

  void axpy(const double &a, const xCSRVector & x,  xCSRVector & y){
   assert(x.size() == y.size());
   xCPPBlasDef<double>::axpy(x.size(), a, const_cast<double *> (&(x[0])), &(y[0]));
 }
 
  double nrm2( const xCSRVector & x){
   return xCPPBlasDef<double>::nrm2(x.size(), const_cast<double *>( &x[0]));
 }

  double dot( const xCSRVector & x, const xCSRVector &y){
   return xCPPBlasDef<double>::dot(x.size(), const_cast<double *>( &x[0]), const_cast<double *>( &y[0]));
 }
 
  void scal(const double& alpha,    xCSRVector & x){
      xCPPBlasDef<double>::scal(x.size(), alpha,  &x[0]);
 }



} // end of namespace

