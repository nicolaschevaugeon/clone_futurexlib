/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef SHAPE_FUNCTIONS_EXTENDED_SHAPE_
#define SHAPE_FUNCTIONS_EXTENDED_SHAPE_

#include "xApproxFunction.h"
#include "xExtendShapeFcts.h"

namespace xfem
{
/// Specialization of xApproxFunction for Extended Shape function
class xApproxFunctionExtendedShape : public xApproxFunction
{
  public:
   xApproxFunctionExtendedShape(const xExtendShapeFcts::extGaussVal_t& asso_) : xApproxFunction(), asso(asso_) {}
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;

  private:
   mutable xExtendShapeFcts::extGaussVal_t asso;
};

}  // namespace xfem

#endif
