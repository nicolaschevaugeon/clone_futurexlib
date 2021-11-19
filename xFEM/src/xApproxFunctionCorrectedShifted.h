/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _CORSHIFT_hh_
#define _CORSHIFT_hh_

#include "xApproxFunction.h"
#include "xEval.h"

namespace xfem
{
/// \brief{Shifted and Corrected approx function}
/// \param[in] uvwInterpo_ : interpolant point for the classical shape function (needed to shift)
/// \param[in] ramp_ : ramp function
/// \param[in] gramp_ : grad of the ramp function
class xApproxFunctionEnrichedXFEMShiftedCorr : public xApproxFunction
{
  public:
   xApproxFunctionEnrichedXFEMShiftedCorr(shapeFctPtr b, shapeFctPtr enr, AOMD::mEntity* e_, xtensor::xPoint uvwInterpo_,
                                          xEval<double>* ramp_, xEval<xtensor::xVector<>>* gramp_)
       : base(b), enrichment(enr), enti(e_), ev_ramp(ramp_), ev_gramp(gramp_), uvwInterpo(uvwInterpo_)
   { /* std::cout<<"entity is"; enti->print();std::cout<<std::endl;*/
   }
   std::string name() override { return "xcApproxFunctionEnrichedXFEM SHIFTED CORRECTED"; }
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   void resetStorage() override { enrichment->resetStorage(); }

  protected:
   shapeFctPtr base;
   shapeFctPtr enrichment;
   AOMD::mEntity* enti;
   const xEval<double>* ev_ramp;
   const xEval<xtensor::xVector<>>* ev_gramp;
   const xtensor::xPoint uvwInterpo;
};

}  // namespace xfem
#endif
