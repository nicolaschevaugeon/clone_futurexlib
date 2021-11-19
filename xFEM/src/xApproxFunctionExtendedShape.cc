/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xApproxFunctionExtendedShape.h"

namespace xfem
{

void xApproxFunctionExtendedShape::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const 
{
  res = asso[geo_integ->getCurrentIntegrationPointId()];
}



} // end of namespace


