/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xCommandOnGeomElem.h"
#include "xGeomElem.h"
#include "xForm.h"

namespace xfem
{
void xCommandOnGeomElem::openApproxElem(xGeomElem* g_appro) { geom_appro = g_appro; }
void xCommandOnGeomElem::setIntegElem(xGeomElem* g_integ) { geom_integ = g_integ; }
void xCommandOnGeomElem::closeApproxElem(xGeomElem* g_appro) { }
xIntegrateFormCommand::xIntegrateFormCommand(xForm* f) : form(f) {}
void xIntegrateFormCommand::execute(xGeomElem* geo_integ) {form->accumulate(geo_integ);
  // std::cout <<  geo_integ << std::endl;
}
} // end of namespace
