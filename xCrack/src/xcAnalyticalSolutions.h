/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef _XCANALYTICALSOLUTIONS_H
#define _XCANALYTICALSOLUTIONS_H

#include <complex>
#include <string>

// xfem
#include "xEval.h"

namespace xfem
{
class xGeomElem;
}

namespace xcrack
{
class xcCrackBase;
}

using namespace xfem;
using namespace xcrack;

/// exact Solution evaluator for planar crack in infinite media.
/*!  Operator  return the stress */
class xcEvalExactCrackStressPlanarCrack : public xEval<xtensor::xTensor2<>>
{
  public:
   xcEvalExactCrackStressPlanarCrack(double K1p, double K2p, double K3p, const xcCrackBase& _fissure);
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& stress) const override;

  private:
   const double K1, K2, K3;
   const xcCrackBase& crack;
};

/// exact Solution for planar crack in infinite media
/*!  Operator  return the displacement */
class xcEvalExactDispPlanarCrack : public xEval<xtensor::xVector<>>
{
  public:
   xcEvalExactDispPlanarCrack(double K1p, double K2p, double K3p, const xcCrackBase& thecrack);
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& disp) const override;

  private:
   const double K1, K2, K3;
   void dispe(double r, double th, result_type& disp) const;
   void getval(const xtensor::xPoint&, result_type&) const;
   const xcCrackBase& crack;
};

/// Class to compute Exact solution for a circle arc crack
/*! The crack is a circle arc of radius 1, centered at (0. 0.), with angle from -theta to +theta
    The domain is an infinite plane, with uniaxial tension at infinity.
    plane strain is assumed.
*/
class xcExactCircularArcCrack
{
  public:
   /// constructor,
   /*! siginf is the uniaxial stress at infinity, alpha is the angle in the xy plane between
       x and the direction of the tension, theta describe the extend of the arc crack. the arc goes from -theta to +theta */
   xcExactCircularArcCrack(double siginf, double alpha, double theta);
   ///  given z, coordinate of a point in the complex plane, return the stress tensor with component express in the x y base;
   xtensor::xTensor2<> stressXY(std::complex<double> z) const;
   ///  given z, coordinate of a point in the complex plane, return the stress tensor with component express in the r alpha base;
   xtensor::xTensor2<> stressAXI(std::complex<double> z) const;

  private:
   std::complex<double> G(std::complex<double> z) const;
   std::complex<double> dGdz(std::complex<double> z) const;
   std::complex<double> d2Gd2z(std::complex<double> z) const;
   std::complex<double> dfdz(std::complex<double> z) const;
   std::complex<double> d2fd2z(std::complex<double> z) const;
   std::complex<double> d2gd2z(std::complex<double> z) const;
   std::complex<double> q(std::complex<double> z) const;

   const double siginf, alpha, theta, eps;
   std::vector<std::complex<double>> a;
   std::vector<std::complex<double>> b;
};

/// xEval that eval thet stress in XY axis state for an arc crack as defined in xcExactCircularArcCrack given on entry
class xcEvalExactStressXYCircularArcCrack : public xEval<xtensor::xTensor2<>>
{
  public:
   xcEvalExactStressXYCircularArcCrack(const xcExactCircularArcCrack& _exactsol);
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& sig) const override;

  private:
   const xcExactCircularArcCrack& exactsol;
};

/// xEval that eval thet stress in axisymetrical axis state for an arc crack as defined in xcExactCircularArcCrack given on entry
class xcEvalExactStressAXICircularArcCrack : public xEval<xtensor::xTensor2<>>
{
  public:
   xcEvalExactStressAXICircularArcCrack(const xcExactCircularArcCrack& _exactsol);
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& sig) const override;

  private:
   const xcExactCircularArcCrack& exactsol;
};

class xcExactPennyShapeCrack
{
  public:
   xcExactPennyShapeCrack(double young, double poisson, double a = 1., double p = 1., double tx = 0., double ty = 0.);
   void setParameters(const std::map<std::string, double>&);
   xtensor::xVector<> getDisplacement(xtensor::xPoint X) const;
   xtensor::xTensor2<> getStress(xtensor::xPoint X) const;

  private:
   double young, poisson, a, p, tx, ty;
};

class frameChange
{
  public:
   frameChange(xtensor::xPoint _center, xtensor::xVector<> _n);
   xtensor::xPoint pushForward(const xtensor::xPoint&) const;
   xtensor::xVector<> pushForward(const xtensor::xVector<>&) const;
   xtensor::xTensor2<> pushForward(const xtensor::xTensor2<>&) const;

   xtensor::xPoint pullBack(const xtensor::xPoint&) const;
   xtensor::xVector<> pullBack(const xtensor::xVector<>&) const;
   xtensor::xTensor2<> pullBack(const xtensor::xTensor2<>&) const;

  private:
   xtensor::xPoint center;
   xtensor::xVector<> n;
   xtensor::xTensor2<> Q;
};

class xcEvalExactPennyShapeCrackStress : public xEval<xtensor::xTensor2<>>
{
  public:
   xcEvalExactPennyShapeCrackStress(double young, double poisson, xtensor::xVector<> _normal, xtensor::xPoint _center,
                                    double _radius, double p = 1., double tx = 0., double ty = 0.);
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& sig) const override;
   void setParameters(const std::map<std::string, double>& para) { exactsol.setParameters(para); };

  private:
   const frameChange localframe;
   //  const xtensor::xVector<> normal;
   // const xtensor::xPoint center;
   const double radius;
   xcExactPennyShapeCrack exactsol;
};

class xcEvalExactPennyShapeCrackDisplacement : public xEval<xtensor::xVector<>>
{
  public:
   xcEvalExactPennyShapeCrackDisplacement(xtensor::xVector<> _normal = xtensor::xVector<>(0., 0., 1.),
                                          xtensor::xPoint _center = xtensor::xPoint(0., 0., 0.), double _radius = 1.,
                                          double p = 1., double tx = 0., double ty = 0.);
   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& sig) const override;
   void setParameters(const std::map<std::string, double>& para) { exactsol.setParameters(para); };

  private:
   const frameChange localframe;
   //  const xtensor::xVector<> normal;
   // const xtensor::xPoint center;
   const double radius;
   xcExactPennyShapeCrack exactsol;
};

/*class xcEvalPennyShapeCrackDisplacement :public xEval<xtensor::xVector<> >{

  };*/

#endif
