/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xForm.h"

#include <typeinfo>

#include "xGeomElem.h"
#include "xMapping.h"

namespace xfem
{
xForm::~xForm() = default;

bool integ_is_appro(xGeomElem* geo_integ, xGeomElem* geo_appro)
{
   auto& mapping_appro = *geo_appro->getMapping();
   auto& mapping_integ = *geo_integ->getMapping();
   return ((geo_appro->getEntity() == geo_integ->getEntity()) && (typeid(mapping_appro) == typeid(mapping_integ)));
}

}  // namespace xfem
