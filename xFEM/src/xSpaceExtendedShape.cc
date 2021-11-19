/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xSpaceExtendedShape.h"

#include "xApproxFunctionExtendedShape.h"
using AOMD::mEntity;

namespace xfem
{
xSpaceExtendedShape::xSpaceExtendedShape(const std::string& space_name, const std::string& physical_name,
                                         xExtendShapeGeneratorBase* generator_)
    : xSpaceRegular(physical_name, xSpace::SCALAR), generator(generator_)
{
}

xSpaceExtendedShape::~xSpaceExtendedShape() = default;

void xSpaceExtendedShape::setCurentExtend(mEntity* e)
{
   curent_extend = generator->getExtendedShapeFcts(*e);
   if (!curent_extend)
   {
      curent_extend = generator->generateExtendShapeFcts(e);
   }
}

void xSpaceExtendedShape::getKeys(mEntity* e, femKeys* keys)
{
   setCurentExtend(e);
   keys->reserve(keys->size() + curent_extend->associated_keys_and_value.size());
   for (const auto& kv : curent_extend->associated_keys_and_value) keys->push_back(kv.first);
}

void xSpaceExtendedShape::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   setCurentExtend(e);
   const size_t nbfunc = curent_extend->associated_keys_and_value.size();
   keys->reserve(keys->size() + nbfunc);
   appro->reserve(appro->size() + nbfunc);
   for (const auto& kv : curent_extend->associated_keys_and_value)
   {
      keys->push_back(kv.first);
      appro->push_back(shapeFctPtr(new xApproxFunctionExtendedShape(kv.second)));
   }
}

}  // end of namespace xfem
