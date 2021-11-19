/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef XMAPPING_BUILDER_HOLDER_H
#define XMAPPING_BUILDER_HOLDER_H

#include <string>

// xtool
#include "xSingleton.h"
// xmapping
#include "xLagrangeMappingBuilder.h"
#include "xMapping.h"
#include "xMappingBuilder.h"

namespace xfem
{
class xMappingBuilderHolder;
using xMappingBuilderHolderSingleton = xtool::xSingleton<xMappingBuilderHolder>;

class xMappingBuilderHolder
{
  public:
   xMappingBuilderHolder();
   ~xMappingBuilderHolder();
   /// default builder is xLagrangeMappingBuilder.
   xMappingBuilder* getMappingBuilder() { return builder; }
   xmapping::xMapping* buildMapping(const AOMD::mEntity& e, int SystemOfCoordinates = 0) const;
   template <class T>
   void setMappingBuilder(T& other)
   {
      delete builder;
      builder = new T(other);
   }

  private:
   xMappingBuilder* builder = nullptr;
};

}  // namespace xfem

#endif
