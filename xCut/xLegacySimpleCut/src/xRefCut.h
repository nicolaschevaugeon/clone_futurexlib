/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/

#ifndef ___XREFCUT__H
#define ___XREFCUT__H


#include <vector>

namespace xcut
{

  class xRefMesh;

  // cutting function for 1D edge
  extern int  cutEdgeRefByLevelSet(const  std::vector<double>  &  , const std::vector < int > & , xRefMesh &  );
  // cutting function for 2D triangle
  extern int  cutTriRefByLevelSet(const  std::vector<double>  &  , const std::vector < int > & , xRefMesh &  );
  // cutting function for 3D tetraedron
  extern int cutTetRefByLevelSet( const std::vector < double > & , const std::vector < int > & , xRefMesh & );


} // end of namespace xfem

#endif
