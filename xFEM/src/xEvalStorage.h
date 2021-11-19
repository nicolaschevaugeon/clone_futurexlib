/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef eval_storage__h
#define eval_storage__h

#include <iostream>

#include "xTensor2.h"
#include "xVector.h"

namespace AOMD
{
class mEntity;
}

namespace xfem
{
class xGeomElem;

class xEvalStorage
{
  public:
   xEvalStorage();
   bool exist(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& v) const;
   void store(const xGeomElem* geo_appro, const xGeomElem* geo_integ, const double& v);
   bool exist(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& v) const;
   void store(const xGeomElem* geo_appro, const xGeomElem* geo_integ, const xtensor::xVector<>& v);
   bool exist(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>& v) const;
   void store(const xGeomElem* geo_appro, const xGeomElem* geo_integ, const xtensor::xTensor2<>& v);
   bool equal_uvw(const xtensor::xPoint& a, const xtensor::xPoint& b) const;

   void reset()
   {
      e_current_d = nullptr;
      e_current_v = nullptr;
      e_current_t = nullptr;
      // uvw_current_d=xtensor::xPoint(0.,0.,0.);
   }

  private:
   AOMD::mEntity* e_current_d;
   xtensor::xPoint uvw_current_d;
   double val_d;
   ///
   AOMD::mEntity* e_current_v;
   xtensor::xPoint uvw_current_v;
   xtensor::xVector<> val_v;
   ///
   AOMD::mEntity* e_current_t;
   xtensor::xPoint uvw_current_t;
   xtensor::xTensor2<> val_t;
};  // end clas xEvalStorage
}  // namespace xfem

#endif
