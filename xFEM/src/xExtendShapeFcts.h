/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef X_EXTEND_SHAPE_FCT_H
#define X_EXTEND_SHAPE_FCT_H
#include <iostream>
#include <map>
#include <set>
// xfem
#include "xField.h"
#include "xFiniteElement.h"
#include "xGeomElem.h"
#include "xIntegrationRuleStored.h"
#include "xMesh.h"

namespace xfem
{
class xExtendShapeFcts
{
  public:
   typedef std::map<int, double> extGaussVal_t;
   typedef extGaussVal_t::iterator extGaussVal_it_t;
   typedef std::map<xValKey, extGaussVal_t> extShapeFct_t;
   typedef extShapeFct_t::iterator extShapeFct_it_t;
   extShapeFct_t associated_keys_and_value;
};

class xExtendShapeGeneratorBase
{
  public:
   virtual xExtendShapeFcts *generateExtendShapeFcts(AOMD::mEntity *e) = 0;
   const xExtendShapeFcts *getExtendedShapeFcts(const AOMD::mEntity &e) const { return data.getData(e); }
   virtual ~xExtendShapeGeneratorBase() = default;

  protected:
   xinterface::aomd::xAttachedDataManagerAOMD<xExtendShapeFcts> data;
};

template <typename DIST, typename VT = double>
class xExtendShapeGenerator : public xExtendShapeGeneratorBase
{
  public:
   template <typename ITERFRONT>
   xExtendShapeGenerator(ITERFRONT beg, ITERFRONT end, xEntityToEntity upper_, xfem::xField<VT> &vgama,
                         std::function<void(AOMD::mEntity *, const xtensor::xPoint &, xmapping::xMapping *,
                                            xfem::xIntegrationRuleStoredDataManager &)>
                             integration_order_function_object_,
                         xfem::xIntegrationRuleStoredDataManager &gp_storage_);
   xExtendShapeFcts *generateExtendShapeFcts(AOMD::mEntity *e) override;
   template <typename ITER>
   void cleanExtendShapeFcts(ITER begin, ITER end);

  private:
   DIST distance_kd_tree;
   xEntityToEntity upper;
   xfem::xField<VT> &vgama;
   xfem::xIntegrationRuleStoredDataManager &gp_storage;
   std::function<void(AOMD::mEntity *, const xtensor::xPoint &, xmapping::xMapping *,
                      xfem::xIntegrationRuleStoredDataManager &gp_manager)>
       integration_order_function_object;

   void cleanExtendShapeFctsElem(AOMD::mEntity *e);
};

}  // namespace xfem

#include "xExtendShapeFcts_imp.h"
#endif
