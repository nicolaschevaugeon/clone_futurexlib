/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _complexPotential_
#define _complexPotential_

#include <complex>

// xtensor
#include "xTensor2.h"
#include "xVector.h"

// xfem
#include "xEval.h"
#include "xGeomElem.h"

using xfem::xEval;
using xfem::xGeomElem;

namespace xanalyticalsolution
{
class complexPotentialElasticity
{
  public:
   typedef std::complex<double> complex;
   enum hypothesis
   {
      PLANESTRAIN,
      PLANESTRESS
   };
   complexPotentialElasticity(double _E, double _nu, hypothesis _hyp);
   virtual ~complexPotentialElasticity();
   xtensor::xVector<> displacementCartesian(complex z) const;
   xtensor::xVector<> displacementAxisymetric(complex z) const;
   xtensor::xTensor2<> stressCartesian(complex z) const;
   xtensor::xTensor2<> strainCartesian(complex z) const;

  protected:
   virtual complex f(complex z) const = 0;
   virtual complex df(complex z) const = 0;
   virtual complex ddf(complex z) const = 0;
   virtual complex dg(complex z) const = 0;
   virtual complex ddg(complex z) const = 0;
   double E, nu;
   hypothesis hyp;
   double k2, mu;
};

class complexPotentialElasticityCircularHole : public complexPotentialElasticity
{
  public:
   complexPotentialElasticityCircularHole(double _E, double _nu, hypothesis _hyp, double _R, double _P);
   complex f(complex z) const override;
   complex df(complex z) const override;
   complex ddf(complex z) const override;
   complex dg(complex z) const override;
   complex ddg(complex z) const override;

  private:
   double R, P;
};

class xEvalDisplacementComplexPotentialElasticity : public xEval<xtensor::xVector<>>
{
  public:
   xEvalDisplacementComplexPotentialElasticity(const complexPotentialElasticity& _cPE) : cPE(_cPE) {}
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& val) const override
   {
      xtensor::xPoint xyz = geo_integ->getXYZ();
      complexPotentialElasticity::complex z(xyz(0), xyz(1));
      val = cPE.displacementCartesian(z);
   }

  private:
   const complexPotentialElasticity& cPE;
};

class xEvalStrainComplexPotentialElasticity : public xEval<xtensor::xTensor2<>>
{
  public:
   xEvalStrainComplexPotentialElasticity(const complexPotentialElasticity& _cPE) : cPE(_cPE) {}
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& val) const override
   {
      xtensor::xPoint xyz = geo_integ->getXYZ();
      complexPotentialElasticity::complex z(xyz(0), xyz(1));
      val = cPE.strainCartesian(z);
   }

  private:
   const complexPotentialElasticity& cPE;
};

class xEvalStressComplexPotentialElasticity : public xEval<xtensor::xTensor2<>>
{
  public:
   xEvalStressComplexPotentialElasticity(const complexPotentialElasticity& _cPE) : cPE(_cPE) {}
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& val) const override
   {
      xtensor::xPoint xyz = geo_integ->getXYZ();
      complexPotentialElasticity::complex z(xyz(0), xyz(1));
      val = cPE.stressCartesian(z);
   }

  private:
   const complexPotentialElasticity& cPE;
};

}  // namespace xanalyticalsolution
#endif
