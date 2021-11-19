/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef XMAPPING_BUILDER_
#define XMAPPING_BUILDER_

// xinterface
#include "xAOMDEntityUtil.h"
// xmapping
#include "xMapping.h"

namespace xfem
{
class xMappingBuilder
{
  public:
   virtual xmapping::xMapping *buildMapping(AOMD::mEntity *, int SystemOfCoordinates = 0) = 0;
   // virtual std::string getMappingBuilderName() const = 0;
   virtual ~xMappingBuilder() = default;
};
}  // namespace xfem

#endif
