/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <algorithm>
#include <cstdio>
#include <iostream>
// xinterface
#include "xAOMDEntityUtil.h"
// mapping
#include "xLagrangeMapping.h"
// xfem
#include "xLagrangeMappingBuilder.h"

using namespace xmapping;
namespace xfem
{
xMapping* xLagrangeMappingBuilder::buildMapping(AOMD::mEntity* e, int SystemOfCoordinates)
{
   const auto points = xinterface::aomd::getPoints(*e);
   const auto type = xinterface::aomd::getRefElementType(*e);
   if (SystemOfCoordinates != 0) return new xCylindricalCoordinatesLagrangeMapping(type, points);
   switch (type)
   {
      case xReferenceElementType::QUAD:
      case xReferenceElementType::HEX:
         break;
      case xReferenceElementType::VERTEX:
      case xReferenceElementType::EDGE:
      case xReferenceElementType::TRI:
      case xReferenceElementType::TET:
         return new xConstantLagrangeMapping(type, points);
         break;
      default:
         break;
   }
   return new xLagrangeMapping(type, points);
}

xMapping* xRegularCubeLagrangeMappingBuilder::buildMapping(AOMD::mEntity* e, int SystemOfCoordinates)
{
   if (SystemOfCoordinates != 0) throw;
   return new xRegularCubeLagrangeMapping(xinterface::aomd::getRefElementType(*e), xinterface::aomd::getPoints(*e));
}

xMapping* xCylindricalCoordinatesLagrangeMappingBuilder::buildMapping(AOMD::mEntity* e, int SystemOfCoordinates)
{
   if (SystemOfCoordinates != 0) throw;
   return new xCylindricalCoordinatesLagrangeMapping(xinterface::aomd::getRefElementType(*e), xinterface::aomd::getPoints(*e));
}

}  // namespace xfem
