#include "xMappingBuilderHolder.h"

#include "xLagrangeMapping.h"

using namespace xmapping;
namespace xfem
{
xMappingBuilderHolder::xMappingBuilderHolder() : builder(new xLagrangeMappingBuilder) {}
xMapping* xMappingBuilderHolder::buildMapping(const AOMD::mEntity& e, int SystemOfCoordinates) const
{
   return builder->buildMapping(const_cast<AOMD::mEntity*>(&e), SystemOfCoordinates);
}
xMappingBuilderHolder::~xMappingBuilderHolder() { delete builder; }
}  // namespace xfem
