/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xEvalStorage.h"

#include "xDebug.h"
#include "xGeomElem.h"

namespace xfem
{
xEvalStorage::xEvalStorage() : e_current_d(nullptr), e_current_v(nullptr), e_current_t(nullptr) {}

bool xEvalStorage::exist(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& v) const
{
   const bool debug = xdebug_flag;
   if (geo_integ->getEntity() == e_current_d && equal_uvw(geo_integ->getUVW(), uvw_current_d))
   {
      if (debug) std::cout << "the value double exist " << val_d << std::endl;
      v = val_d;
      return true;
   }
   return false;
}

void xEvalStorage::store(const xGeomElem* geo_appro, const xGeomElem* geo_integ, const double& v)
{
   const bool debug = xdebug_flag;
   e_current_d = geo_integ->getEntity();
   uvw_current_d = geo_integ->getUVW();
   val_d = v;
   if (debug) std::cout << "the value double is stored " << val_d << std::endl;
}

bool xEvalStorage::exist(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& v) const
{
   if (geo_integ->getEntity() == e_current_v && equal_uvw(geo_integ->getUVW(), uvw_current_v))
   {
      v = val_v;
      return true;
   }
   return false;
}

void xEvalStorage::store(const xGeomElem* geo_appro, const xGeomElem* geo_integ, const xtensor::xVector<>& v)
{
   e_current_v = geo_integ->getEntity();
   uvw_current_v = geo_integ->getUVW();
   val_v = v;
}

bool xEvalStorage::exist(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>& v) const
{
   if (geo_integ->getEntity() == e_current_t && equal_uvw(geo_integ->getUVW(), uvw_current_t))
   {
      v = val_t;
      return true;
   }
   return false;
}

void xEvalStorage::store(const xGeomElem* geo_appro, const xGeomElem* geo_integ, const xtensor::xTensor2<>& v)
{
   e_current_t = geo_integ->getEntity();
   uvw_current_t = geo_integ->getUVW();
   val_t = v;
}

bool xEvalStorage::equal_uvw(const xtensor::xPoint& a, const xtensor::xPoint& b) const
{
   return (a(0) == b(0) && a(1) == b(1) && a(2) == b(2));
}

}  // namespace xfem
