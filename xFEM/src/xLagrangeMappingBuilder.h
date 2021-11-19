/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef XLAGRANGEMAPPINGBUILDER_
#define XLAGRANGEMAPPINGBUILDER_

#include "xMappingBuilder.h"

namespace xfem
{
class xLagrangeMappingBuilder : public xMappingBuilder
{
  public:
   xmapping::xMapping* buildMapping(AOMD::mEntity*, int SystemOfCoordinates = 0) override;
};

class xRegularCubeLagrangeMappingBuilder : public xMappingBuilder
{
  public:
   xmapping::xMapping* buildMapping(AOMD::mEntity* e, int SystemOfCoordinates = 0) override;
};

class xCylindricalCoordinatesLagrangeMappingBuilder : public xMappingBuilder
{
  public:
   xmapping::xMapping* buildMapping(AOMD::mEntity* e, int SystemOfCoordinates = 0) override;
};

}  // namespace xfem

#endif
