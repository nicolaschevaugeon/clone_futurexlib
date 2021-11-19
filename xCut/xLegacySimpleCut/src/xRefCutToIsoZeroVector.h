/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/


#ifndef _REF_CUT_TO_ISO_ZERO_VECTOR_H
#define _REF_CUT_TO_ISO_ZERO_VECTOR_H

#include <vector>

namespace xfem
{
  class xLevelSet;
}


namespace xcut
{	
  class xRefCutToIsoZeroVector
  {
  public:
    int cutEdgeRefByLevelSet(const std::vector < double >  &, std::vector < double > &  ) const ;
    int cutTriRefByLevelSet(const std::vector < double >  &, std::vector < double > &  ) const ;
    int cutTetRefByLevelSet(const std::vector < double >  &, std::vector < double > &  ) const ;
  private:
    // for tet
    static unsigned char c1_node[6][2];
    static unsigned char c2_node[5][5];
    static double tet_node[4][3];
    static unsigned char size3d;

  };

}

#endif
