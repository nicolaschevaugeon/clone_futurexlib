/*
   xfem : C++ Finite Element Library
   developed under the GNU Lesser General Public License
   See the NOTICE & LICENSE files for conditions.
*/

#include "OctreeUtilities.h"
#include "xLevelSet.h"
#include "xMesh.h"

using namespace xfem;
using namespace xoctree;

namespace xinterface
{
namespace xoctree
{
// Caution : this one is VERY inefficient...

xLevelSetTo_ls_xyz_t_Adaptator::xLevelSetTo_ls_xyz_t_Adaptator(xLevelSet &ls_) : ls(ls_), mesh(ls.getSupport().getMesh()) {}

double xLevelSetTo_ls_xyz_t_Adaptator::operator()(const double &x, const double &y, const double &z)
{
   Trellis_Util::mPoint p(x, y, z);
   std::set<AOMD::mEntity *> elts;
//    mesh->locateElement(p, elts);
   throw;//RECODER LA LIGNE CI-DESSUS AVEC LA NOUVELLE VERSION DE locateElement
   if (!elts.empty())
   {
      AOMD::mEntity *e = *(elts.begin());
      xfem::xGeomElem geo(e);
      geo.setUVWForXYZ(p);
      return ls.getVal(e, geo.getUVW());
   }
   else
   {
      throw;
   }
}

xPointToDoubleTo_ls_xyz_t_Adaptator::xPointToDoubleTo_ls_xyz_t_Adaptator(xPointToDouble &xpd_) : xpd(xpd_) {}

double xPointToDoubleTo_ls_xyz_t_Adaptator::operator()(const double &x, const double &y, const double &z)
{
   Trellis_Util::mPoint p(x, y, z);
   return xpd(p);
}

}  // end namespace xoctree
}  // end namespace xinterface
