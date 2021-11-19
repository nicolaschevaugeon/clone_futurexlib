
/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/

  
#include "complexPotential.h"

namespace xanalyticalsolution
{

  complexPotentialElasticity::complexPotentialElasticity(double _E, double _nu, hypothesis _hyp):E(_E), nu(_nu), hyp(_hyp){
    mu  = E/(2.*(1.+nu));
    switch (hyp){
    case PLANESTRAIN:
      k2 = (1.-2.*nu) /(2.*(1.-nu));
      break;
    case PLANESTRESS:
      k2 = (1.-nu) /2.;
      break;
    }
  }

  complexPotentialElasticity::~complexPotentialElasticity()= default;;
  xtensor::xVector<> complexPotentialElasticity::displacementCartesian  (complex z) const {
    // 2*mu(u+iv)
    complex mu2upiv = (1+k2)/(1-k2) * f(z) - z * conj(df(z)) - conj(dg(z));
    return xtensor::xVector<>(real( mu2upiv)/2./mu, imag( mu2upiv)/2./mu );
  };
  
  xtensor::xVector<> complexPotentialElasticity:: displacementAxisymetric(complex z) const{
    // 2*mu(u+iv)
    complex mu2upiv = ((1+k2)/(1-k2) * f(z) - z * conj(df(z)) - conj(dg(z)))*exp( complex(0, -arg(z) ));
    return xtensor::xVector<>(real( mu2upiv)/2./mu, imag( mu2upiv)/2./mu );
  };

  xtensor::xTensor2<> complexPotentialElasticity::stressCartesian ( complex z) const {
    //sigmaxx + sigmayy
    double sxxpsyy = 4*real(df(z));
    //sigmayy - sigmaxx + 2isigmaxy
    complex syymsxxp2isxy = 2.*( conj(z)* ddf(z) + ddg(z) );
    double sxx = (sxxpsyy - real(syymsxxp2isxy))*0.5;
    double syy = (sxxpsyy + real(syymsxxp2isxy))*0.5;
    double sxy = 0.5*imag(syymsxxp2isxy);
    double szz = 0.;
    if (hyp == PLANESTRAIN)
      szz = nu * (sxx+syy);
    xtensor::xTensor2<> stress(0.);
    stress(0,0) = sxx;
    stress(0,1) = sxy;
    stress(1,0) = sxy;
    stress(1,1) = syy;
    stress(2,2) = szz;
    return stress;
  }; 
  
  xtensor::xTensor2<> complexPotentialElasticity::strainCartesian ( complex z) const{
    xtensor::xTensor2<> stress = stressCartesian(z);
    double epsxx=0., epsyy=0., epsxy=0., epszz=0.;
    switch (hyp){
    case PLANESTRESS :{
      epsxx = (1./E) * (stress(0,0) - nu*stress(1,1));
      epsyy = (1./E) * (stress(1,1) - nu*stress(0,0));
      epsxy = (1./E) * (1+nu) *stress(0,1);
      epszz = - nu / (1-nu) *(epsxx+epsyy);
      break;
    }
    case PLANESTRAIN:{
      epsxx = (1.+nu)/E * ((1-nu)* stress(0,0) - nu*stress(1,1));
      epsyy = (1.+nu)/E * ((1-nu)* stress(1,1) - nu*stress(0,0));
      epsxy = (1.+nu)/E * stress(0,1);
      epszz = 0.;
      break;
    }
    }
    xtensor::xTensor2<> strain(0.);
    strain(0,0) = epsxx;
    strain(1,1) = epsyy;
    strain(0,1) = epsxy;
    strain(1,0) = epsxy;
    strain(2,2) = epszz;
    return strain;
  }

  complexPotentialElasticityCircularHole::complexPotentialElasticityCircularHole(double _E, double _nu, hypothesis _hyp, double _R, double _P):complexPotentialElasticity(_E, _nu, _hyp), R(_R), P(_P){};

  complexPotentialElasticity::complex complexPotentialElasticityCircularHole::f(complex z) const{
    return P*R/4.*(z/R+ 2*R/z);
  }

  complexPotentialElasticity::complex complexPotentialElasticityCircularHole::df(complex z) const{
    return P*R/4.*(1/R  - 2*R/z/z);
  }

  complexPotentialElasticity::complex complexPotentialElasticityCircularHole::ddf(complex z) const{
    return P*R *R/z/z/z;
  }

  complexPotentialElasticity::complex complexPotentialElasticityCircularHole::dg(complex z) const{
    return -P*R/2.*(z/R+ R/z - pow(R, 3)*pow(z, -3) );
  }

  complexPotentialElasticity::complex complexPotentialElasticityCircularHole::ddg(complex z) const{
    return -P*R/2.*(1/R- R/z/z +3* pow(R, 3)*pow(z, -4) );
  }
  
}
