/*
     This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

// xfem
#include "xGeomElem.h"
#include "xMaterialManager.h"

// xcrack
#include "xcAnalyticalSolutions.h"
#include "xcCrackBase.h"

using std::string;

xcEvalExactCrackStressPlanarCrack::xcEvalExactCrackStressPlanarCrack(double K1p, double K2p, double K3p,
                                                                     const xcCrackBase& _fissure)
    : K1(K1p), K2(K2p), K3(K3p), crack(_fissure)
{
}

void xcEvalExactCrackStressPlanarCrack::operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                   result_type& stress) const
{
   // note c++20 standard we have std::numbers::pi
   const double PI = (3.1415926535897932384626433832795);
   xtensor::xPoint uvw = geo_appro->getUVW();
   AOMD::mEntity* e = geo_appro->getEntity();
   double th, r, side;
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);

   xMaterial* mat = xMaterialManagerSingleton::instance().getMaterial(geo_appro);
   const xTensors* properties = mat->getProperties();
   // double young   = properties->scalar("YOUNG_MODULUS");
   double poisson = properties->scalar("POISSON_RATIO");

   const double nu = poisson;
   // double c=cos(th);
   // double s=sin(th);
   double c2 = cos(th / 2.);
   double s2 = sin(th / 2.);
   double c32 = cos(3. * th / 2.);
   double s32 = sin(3. * th / 2.);
   double K1r = K1 / sqrt(2. * PI * r);
   double K2r = K2 / sqrt(2. * PI * r);
   double K3r = K3 / sqrt(2. * PI * r);

   stress(0, 0) = K1r * c2 * (1. - s2 * s32) - K2r * s2 * (2. + c2 * c32);
   stress(1, 0) = stress(0, 1) = K1r * c2 * s2 * c32 + K2r * c2 * (1. - s2 * s32);
   stress(1, 1) = K1r * c2 * (1. + s2 * s32) + K2r * s2 * c2 * c32;
   stress(0, 2) = stress(2, 0) = -K3r * s2;
   stress(1, 2) = stress(2, 1) = K3r * c2;
   stress(2, 2) = nu * (stress(0, 0) + stress(1, 1));

   crack.localToGlobal(e, uvw, stress);
}

xcEvalExactDispPlanarCrack::xcEvalExactDispPlanarCrack(double K1p, double K2p, double K3p, const xcCrackBase& thecrack)
    : K1(K1p), K2(K2p), K3(K3p), crack(thecrack)
{
}

void xcEvalExactDispPlanarCrack::dispe(double r, double th, result_type& disp) const
{
   // compute displacement for mu =1.;
   const double PI = (3.1415926535897932384626433832795);

   double mu = 1.;
   const double k = 3;
   double c = cos(th);
   // double s=sin(th);

   double c2 = cos(th / 2.);
   double s2 = sin(th / 2.);
   double K1r = K1 * sqrt(r / (2 * PI)) / (2 * mu);
   double K2r = K2 * sqrt(r / (2 * PI)) / (2 * mu);
   double K3r = 2. * K3 * sqrt(r / (2 * PI)) / (mu);

   disp(0) = K1r * c2 * (k - c) + K2r * s2 * (k + c + 2);
   disp(1) = K1r * s2 * (k - c) - K2r * c2 * (k + c - 2);
   disp(2) = K3r * s2;
}

void xcEvalExactDispPlanarCrack::getval(const xtensor::xPoint& xyz, result_type& disp) const { dispe(xyz(0), xyz(1), disp); }

void xcEvalExactDispPlanarCrack::operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                            result_type& disp) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   AOMD::mEntity* e = geo_appro->getEntity();
   double th, r, side;
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   result_type disp_coin, disp_coin2;
   dispe(r, th, disp);
   xMaterial* mat = xMaterialManagerSingleton::instance().getMaterial(geo_appro);
   const xTensors* properties = mat->getProperties();
   double young = properties->scalar("YOUNG_MODULUS");
   double poisson = properties->scalar("POISSON_RATIO");
   const double nu = poisson;
   const double E = young;
   // const double lam=nu*E/((1-2*nu)*(1+nu));
   const double mu = E / (2 * (1 + nu));
   disp /= mu;
   crack.localToGlobal(e, uvw, disp);
}

xcExactCircularArcCrack::xcExactCircularArcCrack(double _siginf, double _alpha, double _theta)
    : siginf(_siginf), alpha(_alpha), theta(_theta), eps(0.0001)
{
   a.reserve(3);
   b.reserve(4);
   std::complex<double> e2ialpha(cos(2 * alpha), sin(2 * alpha));

   b[0] = -siginf / 4. * e2ialpha;
   b[1] = -cos(theta) * b[0];
   b[3] = (3. * b[0] * pow(cos(theta), 2) - siginf + 2. * b[1] * cos(theta) - b[0]) / (2. * (cos(theta) - 3.));
   b[2] = -cos(theta) * b[3];

   a[0] = b[0];
   a[1] = 0.;
   a[2] = siginf / 4. - b[3];
}

xtensor::xTensor2<> xcExactCircularArcCrack::stressAXI(std::complex<double> z) const
{
   // $\Theta$ and $\Phi$ are the complex potential as defined in Elasticity, J.R.Barber
   // $\xi$ is $-f$ and $\theta$ is $-g'$

   std::complex<double> THETA = 2. * (dfdz(z) + conj(dfdz(z)));
   std::complex<double> PHIalpha = -2. * conj(z) / z * (z * conj(d2fd2z(z)) + conj(d2gd2z(z)));

   double sigrr = std::real(0.5 * (THETA + PHIalpha));
   double sigralpha = std::imag(0.5 * (THETA + PHIalpha));
   double sigalphaalpha = std::real(THETA) - sigrr;
   xtensor::xTensor2<> sig;
   sig(0, 0) = sigrr;
   sig(0, 1) = sigralpha;
   sig(1, 0) = sigralpha;
   sig(1, 1) = sigalphaalpha;
   sig(0, 2) = 0.;
   sig(1, 2) = 0.;
   sig(2, 2) = 0.;
   sig(2, 0) = 0.;
   sig(2, 1) = 0.;
   sig(2, 2) = 0.;
   return sig;
}

xtensor::xTensor2<> xcExactCircularArcCrack::stressXY(std::complex<double> z) const
{
   std::complex<double> THETA = 2. * (dfdz(z) + conj(dfdz(z)));
   std::complex<double> PHI = -2. * (z * conj(d2fd2z(z)) + conj(d2gd2z(z)));
   double sigxx = std::real(0.5 * (THETA + PHI));
   double sigxy = std::imag(0.5 * (THETA + PHI));
   double sigyy = std::real(THETA) - sigxx;
   xtensor::xTensor2<> sig;
   sig(0, 0) = sigxx;
   sig(0, 1) = sigxy;
   sig(1, 0) = sigxy;
   sig(1, 1) = sigyy;
   sig(0, 2) = 0.;
   sig(1, 2) = 0.;
   sig(2, 2) = 0.;
   sig(2, 0) = 0.;
   sig(2, 1) = 0.;
   sig(2, 2) = 0.;
   return sig;
};

std::complex<double> xcExactCircularArcCrack::G(std::complex<double> z) const
{
   const std::complex<double> i(0., 1.);
   const std::complex<double> one(1., 0.);
   std::complex<double> a = exp(i * theta);
   std::complex<double> b = exp(-i * theta);
   std::complex<double> g = sqrt((pow(z, 2) - 2. * z * cos(theta) + one));
   if (real(z - a) < 0 || abs(z) < 1.) g = -g;
   return g;
   std::complex<double> sqrtzmb =
       ((real(z - b) < 0 && imag(z - b) > 0) || (((pow(real(z), 2) + pow(imag(z), 2)) < 1.) && imag(z - b) > 0.)) ? sqrt(z - b)
                                                                                                                  : -sqrt(z - b);
   std::complex<double> sqrtzma = (real(z - a) < 0 && imag(z - a) > 0) ? sqrt(z - a) : -sqrt(z - a);
   return sqrtzmb * sqrtzma;
}

std::complex<double> xcExactCircularArcCrack::dGdz(std::complex<double> z) const
{
   std::complex<double> dGdz = -(cos(theta) - z) / G(z);
   return dGdz;
}

std::complex<double> xcExactCircularArcCrack::d2Gd2z(std::complex<double> z) const
{
   const std::complex<double> Gz = G(z);
   const std::complex<double> dGdzz = dGdz(z);
   std::complex<double> d2Gd2z = (Gz + cos(theta - z) * dGdzz) / Gz / Gz;
   return d2Gd2z;
}

std::complex<double> xcExactCircularArcCrack::dfdz(std::complex<double> z) const
{
   std::complex<double> dfdz(0., 0.);
   std::complex<double> g = G(z);
   if (abs(z) < eps)
   {
      dfdz += a[2] - b[2] - cos(theta) * b[1] - 0.5 * (3. * pow(cos(theta), 2) - 1.) * b[0];
      dfdz += (a[3] - b[3] - cos(theta) * b[2] - 0.5 * (3. * pow(cos(theta), 2) - 1.) * b[1]) * z;
      std::cout << "warning using dl df" << dfdz << std::endl;
      return dfdz;
   }

   std::vector<std::complex<double>> znm2(4);
   for (int i = 0; i < 4; ++i)
   {
      znm2[i] = pow(z, i - 2);
   }

   for (int i = 0; i < 4; ++i)
   {
      dfdz += b[i] * znm2[i];
   }
   dfdz /= g;
   for (int i = 0; i < 3; ++i)
   {
      dfdz += a[i] * znm2[i];
   }
   return dfdz;
}

std::complex<double> xcExactCircularArcCrack::q(std::complex<double> z) const
{
   std::complex<double> q(0., 0.);
   std::complex<double> g = G(z);
   if (abs(z) < eps)
   {
      // std::cout << "warning using dl q"  << std::endl;
      q += (b[2] + b[3] * z) / g;
      q += a[2];
      return q;
   }

   std::vector<std::complex<double>> znm2(4);
   for (int i = 0; i < 4; ++i)
   {
      znm2[i] = pow(z, i - 2);
   }

   for (int i = 0; i < 4; ++i)
   {
      q += b[i] * znm2[i];
   }
   q /= g;
   for (int i = 0; i < 3; ++i)
   {
      q -= a[i] * znm2[i];
   }
   return q;
}

std::complex<double> xcExactCircularArcCrack::d2fd2z(std::complex<double> z) const
{
   std::complex<double> d2fd2z(0., 0.);
   std::complex<double> Gz = G(z);
   std::complex<double> dGdzz = dGdz(z);
   if (abs(z) < eps)
   {
      // std::cout << "warning using dl d2fd2z"  << std::endl;
      d2fd2z = -(b[2] + b[3] * z) * dGdzz / pow(Gz, 2) + b[3] / Gz;  // this is wrong ...need to be corrected ..
      return d2fd2z;
   }
   std::vector<std::complex<double>> znm3(4);
   std::vector<std::complex<double>> znm2(4);
   for (int i = 0; i < 4; ++i)
   {
      znm3[i] = pow(z, i - 3);
      znm2[i] = pow(z, i - 2);
   }

   for (int i = 0; i < 4; ++i)
   {
      d2fd2z += b[i] * znm2[i];
   }
   d2fd2z *= (-dGdzz / pow(Gz, 2));

   std::complex<double> temp(0., 0.);
   for (int i = 0; i < 4; ++i)
   {
      temp += (i - 2.) * b[i] * znm3[i];
   }
   d2fd2z += (temp / Gz);

   for (int i = 0; i < 3; ++i)
   {
      d2fd2z += (i - 2.) * a[i] * znm3[i];
   }

   return d2fd2z;
}

std::complex<double> xcExactCircularArcCrack::d2gd2z(std::complex<double> z) const
{
   std::complex<double> d2gd2z(0., 0.);
   std::complex<double> qinvz = q(1. / conj(z));

   if (abs(z) < eps)
   {  // d2gd2z should be regular at 0, using 0 order term and first order of the dl
      double c = cos(theta);
      double s = sin(theta);
      std::complex<double> dfdivzz =
          -c * b[3] + 0.5 * b[2] * (c * c - 1.) + 0.5 * b[1] * (c - c * c * c) +
          b[0] * (-0.125 + 0.75 * c * c - 0.625 * c * c * c * c) +
          0.125 * s * s * ((7 * c * c * c - 3 * c) * b[0] + (-1. + 5. * c * c) * b[1] + 4. * c * b[2] + 4. * b[3]) * z;
      std::complex<double> d2fd2zdivz =
          (1.5 * c * c - 0.25 - 1.25 * c * c * c * c) * b[0] + (c - c * c * c) * b[1] + (1. - c * c) * b[2] - 2 * c * b[3] +
          (0.375 * s * s * ((7 * c * c * c - 3 * c) * b[0] + (5. * c * c - 1) * b[1] + 4. * c * b[2] + 4. * b[3])) * z;
      std::complex<double> qinvzz =
          (-c) * b[0] + (0.5 - 0.5 * c * c) * b[1] + (0.5 * c - 0.5 * c * c * c) * b[2] +
          (-0.125 + 0.75 * c * c - 0.625 * c * c * c * c) * b[3] - a[0] +
          0.125 * s * s * (4. * b[0] - b[2] + 4. * c * b[1] - 3. * c * b[3] + 5. * c * c * b[2] + 7. * b[3] * c * c * c) * z;
      d2gd2z = -qinvzz + dfdivzz - d2fd2zdivz;
      return d2gd2z;
   }
   d2gd2z = -conj(qinvz) / (z * z) + dfdz(z) / (z * z) - d2fd2z(z) / (z);
   return d2gd2z;
}

xcEvalExactStressXYCircularArcCrack::xcEvalExactStressXYCircularArcCrack(const xcExactCircularArcCrack& _exactsol)
    : exactsol(_exactsol){};
void xcEvalExactStressXYCircularArcCrack::operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                     result_type& sig) const
{
   xtensor::xPoint XYZ = geo_appro->getXYZ();
   double x = XYZ(0);
   double y = XYZ(1);
   std::complex<double> z(x, y);
   sig = exactsol.stressXY(z);
   return;
};

xcEvalExactStressAXICircularArcCrack::xcEvalExactStressAXICircularArcCrack(const xcExactCircularArcCrack& _exactsol)
    : exactsol(_exactsol){};
void xcEvalExactStressAXICircularArcCrack::operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                      result_type& sig) const
{
   xtensor::xPoint XYZ = geo_appro->getXYZ();
   double x = XYZ(0);
   double y = XYZ(1);
   std::complex<double> z(x, y);
   sig = exactsol.stressAXI(z);
   return;
};

xcExactPennyShapeCrack::xcExactPennyShapeCrack(double _young, double _poisson, double _a, double _p, double _tx, double _ty)
    : young(_young), poisson(_poisson), a(_a), p(_p), tx(_tx), ty(_ty)
{
}

xtensor::xVector<> xcExactPennyShapeCrack::getDisplacement(xtensor::xPoint X) const
{
   const double x = X(0);
   const double y = X(1);
   const double z = X(2);
   const double rho = sqrt((x) * (x) + (y) * (y));
   const double l1 = 0.5 * (sqrt((a + rho) * (a + rho) + z * z) - sqrt((a - rho) * (a - rho) + z * z));
   const double l2 = 0.5 * (sqrt((a + rho) * (a + rho) + z * z) + sqrt((a - rho) * (a - rho) + z * z));
   const double phi = atan2(y, x);
   double l2l2maa = l2 * l2 - a * a;
   const double pi = M_PI;
   /*
   xMaterial *mat = xMaterialManagerSingleton::instance().getMaterial(geo_appro);
   const xTensors* properties = mat->getProperties();
   const double young   = properties->scalar("YOUNG_MODULUS");
   const double nu = properties->scalar("POISSON_RATIO");
   */

   if (l2l2maa < 0.) l2l2maa = 0.;
   // double K1exa=2.*p*sqrt(a/pi);
   double GG = 0.5 * young / (1. + poisson);
   double u = p * rho / (2. * pi * GG) *
              ((1 - 2. * poisson) * (a * sqrt(l2 * l2 - a * a) / (l2 * l2) - asin(a / l2)) +
               2. * a * a * sqrt(z * z) * sqrt(a * a - l1 * l1) / (l2 * l2 * (l2 * l2 - l1 * l1)));
   double uz = p / (pi * GG) *
                   (2. * (1 - poisson) * (z / sqrt(z * z) * sqrt(a * a - l1 * l1) - z * asin(a / l2)) +
                    z * (asin(a / l2) - a * sqrt(l2 * l2 - a * a) / (l2 * l2 - l1 * l1))) +
               p * z / young;

   double ux = u * cos(phi) - poisson * p * x / young;
   double uy = u * sin(phi) - poisson * p * y / young;
   xtensor::xVector<> disp;
   disp(0) = ux;
   disp(1) = uy;
   disp(2) = uz;
   return disp;
}

void xcExactPennyShapeCrack::setParameters(const std::map<std::string, double>& para)
{
   std::map<string, double>::const_iterator it = para.begin();
   std::map<string, double>::const_iterator itend = para.end();
   while (it != itend)
   {
      std::string data = it->first;
      double val = it->second;
      if (data == "radius")
         a = val;
      else if (data == "pressure")
         p = val;
      else
      {
         std::cout << "Warning :  data " << data << " is not known by xcExactPennyShapeCrack" << std::endl;
      }
      ++it;
   }
};

xtensor::xTensor2<> xcExactPennyShapeCrack::getStress(xtensor::xPoint X) const
{
   const bool debug = xdebug_flag;

   const double x = X(0);
   const double y = X(1);
   const double z = X(2);
   const double rho = sqrt((x) * (x) + (y) * (y));
   const double l1 = 0.5 * (sqrt((a + rho) * (a + rho) + z * z) - sqrt((a - rho) * (a - rho) + z * z));
   const double l2 = 0.5 * (sqrt((a + rho) * (a + rho) + z * z) + sqrt((a - rho) * (a - rho) + z * z));
   const double phi = atan2(y, x);

   const double pi = M_PI;

   double l2l2maa = l2 * l2 - a * a;
   if (l2l2maa < 0.) l2l2maa = 0.;  // GEOMEtricaly It SHOULD >= 0. BUT ... sometime IT'S SLIGHTLY LESS THAN 0.
   double sqrtl2l2maa = sqrt(l2l2maa);
   double asinadl2 = ((a / l2) >= 1) ? pi / 2. : (((a / l2) <= -1) ? -pi / 2 : asin(a / l2));

   double zzdl2l2maa = 0.;      // z^2 / (l2^2 - a^2)
   double zzdsqrtl2l2maa = 0.;  // z^2 / sqrt(l2^2 - a^2)

   if ((fabs(z) < 1.e-6) && (rho < a))
   {
      zzdl2l2maa = 1. / (a * (1. / (a + rho) + 1. / ((a - rho))));
      // z^2 / (l2^2 - a^2) -> 1/a(1/(a+rho) + 1/(a-rho)) when z-> 0
      zzdsqrtl2l2maa = z * sqrt(zzdl2l2maa);
      // z^2 /sqrt (l2^2 - a^2) -> z / sqrt(a(1/(a+rho) + 1/(a-rho)))  when z-> 0
   }
   else
   {
      zzdl2l2maa = z * z / l2l2maa;
      zzdsqrtl2l2maa = z * z / sqrt(l2l2maa);
   }

   // double K1exa=2.*p*sqrt(a/pi);
   // double GG=0.5*young/(1.+nu);
   // if ((l2*l2 - a*a) < 0.) std::cout << l2 << " " << a << l2*l2 - a*a << std::endl;

   if (std::fabs(l1 - l2) < 1.e-6) std::cout << l1 << " " << l2 << " " << std::endl;
   if (std::fabs(l2) < 1.e-6) std::cout << l2 << std::endl;

   double sigma1 =
       2. * p / pi *
       ((1 + 2. * poisson) * (a * sqrtl2l2maa / (l2 * l2 - l1 * l1) - asinadl2) +
        a * zzdsqrtl2l2maa * (pow(l1, 4) + a * a * (2. * a * a + 2. * z * z - 3. * rho * rho)) / std::pow(l2 * l2 - l1 * l1, 3));

   double sigma2 = 2. * p / pi * a * l1 * l1 * sqrtl2l2maa / (l2 * l2 * (l2 * l2 - l1 * l1)) *
                   (1 - 2. * poisson +
                    (zzdl2l2maa * (a * a * (6. * l2 * l2 - 2. * l1 * l1 + rho * rho) - 5. * l2 * l2 * l2 * l2) /
                     ((l2 * l2 - l1 * l1) * (l2 * l2 - l1 * l1))));

   double sigmazz = 2. * p / pi *
                        (a * sqrtl2l2maa / (l2 * l2 - l1 * l1) - asinadl2 -
                         a * zzdsqrtl2l2maa * (l1 * l1 * l1 * l1 + a * a * (2. * a * a + 2 * z * z - 3. * rho * rho)) /
                             std::pow(l2 * l2 - l1 * l1, 3)) +
                    p;

   double tauz = -2. * p / pi * z * l1 * sqrtl2l2maa * (a * a * (4. * l2 * l2 - 5. * rho * rho) + pow(l1, 4)) /
                 (l2 * std::pow(l2 * l2 - l1 * l1, 3));
   if (debug)
   {
      if (std::isnan(sigma1)) std::cout << "sigma1 " << rho << " " << z << std::endl;
      if (std::isnan(sigma2)) std::cout << "sigma2 " << std::endl;
      if (std::isnan(tauz)) std::cout << "tauz " << std::endl;
      if (std::isnan(sigmazz)) std::cout << "sigmazz " << rho << " " << z << std::endl;
   }
   xtensor::xTensor2<> stress;
   stress(0, 0) = (sigma1 + sigma2 * cos(2. * phi)) * 0.5;
   stress(1, 1) = (sigma1 - sigma2 * cos(2. * phi)) * 0.5;
   stress(1, 0) = stress(0, 1) = sigma2 * sin(2. * phi) * 0.5;
   stress(0, 2) = stress(2, 0) = tauz * cos(phi);
   stress(1, 2) = stress(2, 1) = tauz * sin(phi);
   stress(2, 2) = sigmazz;
   // std::cout << stress(0,0) << " " <<  stress(1,1) << " " <<   stress(2,2) << " " <<  stress(1,0) << " "<<  stress(0,2)<< " "
   // << stress(1,2) << " "<< std::endl;
   return stress;
}

/*
xtensor::xTensor2<> xcExactPennyShapeCrack::getStress(xtensor::xPoint X) const {
  //xtensor::xTensor2<> xcExactPennyShapeCrack::getStressShearLoading(xtensor::xPoint X) const {
  const bool debug = 1;

  const double x = X(0);
  const double y = X(1);
  const double z = X(2);
  const double rho=sqrt((x)*(x)+(y)*(y));
  const double l1=0.5*(sqrt((a+rho)*(a+rho)+z*z)-sqrt((a-rho)*(a-rho)+z*z));
  const double l2=0.5*(sqrt((a+rho)*(a+rho)+z*z)+sqrt((a-rho)*(a-rho)+z*z));
  const double phi=atan2(y,x);

  const double pi=M_PI;

  double l2l2maa = l2*l2 - a*a;
  if (l2l2maa<0.) l2l2maa = 0.; // GEOMEtricaly It SHOULD >= 0. BUT ... sometime IT'S SLIGHTLY LESS THAN 0.
  double sqrtl2l2maa = sqrt(l2l2maa);
  double asinadl2 = ((a/l2) >= 1) ? pi/2. :(((a/l2) <= -1) ? -pi/2: asin(a/l2));



  double zzdl2l2maa =0.;     //z^2 / (l2^2 - a^2)
  double zzdsqrtl2l2maa = 0.; // z^2 / sqrt(l2^2 - a^2)


  if ((fabs(z)< 1.e-6)&&(rho < a))  {

      zzdl2l2maa =  1./(a* (1./(a+rho) + 1./((a-rho)) ) );
      // z^2 / (l2^2 - a^2) -> 1/a(1/(a+rho) + 1/(a-rho)) when z-> 0
      zzdsqrtl2l2maa =  z *sqrt ( zzdl2l2maa);
      //z^2 /sqrt (l2^2 - a^2) -> z / sqrt(a(1/(a+rho) + 1/(a-rho)))  when z-> 0
    }
  else{
    zzdl2l2maa = z*z/ l2l2maa;
    zzdsqrtl2l2maa = z*z / sqrt(l2l2maa);
  }

  //double K1exa=2.*p*sqrt(a/pi);
  //double GG=0.5*young/(1.+nu);
  // if ((l2*l2 - a*a) < 0.) std::cout << l2 << " " << a << l2*l2 - a*a << std::endl;

  if (std::fabs(l1-l2) < 1.e-6) std::cout << l1 << " " << l2 << " "  << std::endl;
  if  (std::fabs(l2)< 1.e-6) std::cout << l2  << std::endl;
  double nu = poisson;


  double sigma1= 2/pi/(2-nu)*(
                 - 2*(1+nu)*a*l1*sqrt(a*a-l1*l1)/(l2*(l2*l2-l1*l1))
                 + z*l1*sqrt(l2*l2-a*a)*(a*a*(4*l2*l2-5*rho*rho)+l1*l1*l1*l1)/(l2*pow(l2*l2-l1*l1,3))
                                          );

 if (debug ){
   if (isnan(sigma1))
     std::cout << "sigma1 " << rho << " " << z  << std::endl;
   // if (isnan(sigma2))
   //  std::cout << "sigma2 " << std::endl;
   //if (isnan(tauz))
   //std::cout << "tauz "   << std::endl;
   //if (isnan(sigmazz))
   //std::cout << "sigmazz " << rho << " " << z  << std::endl;
 }
  xtensor::xTensor2<> stress;
  stress(0,0) = sigma1;

  //std::cout << stress(0,0) << " " <<  stress(1,1) << " " <<   stress(2,2) << " " <<  stress(1,0) << " "<<  stress(0,2)<< " " <<
stress(1,2) << " "<< std::endl; return stress;
  }*/

xcEvalExactPennyShapeCrackStress::xcEvalExactPennyShapeCrackStress(double young, double poisson, xtensor::xVector<> _normal,
                                                                   xtensor::xPoint _center, double _radius, double _p, double _tx,
                                                                   double _ty)
    : localframe(_center, _normal), radius(_radius), exactsol(young, poisson, _radius, _p, _tx, _ty){};

// normal(_normal.norm()), center(_center), radius(_radius ), exactsol(_radius, _p, _tx, _ty){};

void xcEvalExactPennyShapeCrackStress::operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                  result_type& sig) const
{
   // std::cout << "Eval Stress Penny " << std::endl;

   xtensor::xPoint XYZ = geo_appro->getXYZ();
   xtensor::xPoint x = localframe.pushForward(XYZ);

   xtensor::xTensor2<> sigxy = exactsol.getStress(x);
   sig = localframe.pullBack(sigxy);
   return;
}

xcEvalExactPennyShapeCrackDisplacement::xcEvalExactPennyShapeCrackDisplacement(xtensor::xVector<> _normal,
                                                                               xtensor::xPoint _center, double _radius, double _p,
                                                                               double _tx, double _ty)
    : localframe(_center, _normal), radius(_radius), exactsol(_radius, _p, _tx, _ty){};

// normal(_normal.norm()), center(_center), radius(_radius ), exactsol(_radius, _p, _tx, _ty){};

void xcEvalExactPennyShapeCrackDisplacement::operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                        result_type& disp) const
{
   xtensor::xPoint XYZ = geo_appro->getXYZ();
   xtensor::xPoint x = localframe.pushForward(XYZ);
   xtensor::xVector<> displ = exactsol.getDisplacement(x);
   disp = localframe.pullBack(displ);
   return;
}

frameChange::frameChange(xtensor::xPoint _center, xtensor::xVector<> _n) : center(_center), n(_n.norm())
{
   xtensor::xVector<> eZ = n;
   xtensor::xVector<> ex(1., 0., 0.);
   xtensor::xVector<> ey(0., 1., 0.);
   xtensor::xVector<> ez(0., 0., 1.);
   xtensor::xVector<> eX = (ey % eZ);
   if (eX.mag() < 1.e-3) eX = (ex % eZ);
   eX.norm();
   xtensor::xVector<> eY = (eZ % eX);
   Q(0, 0) = ex * eX;
   Q(0, 1) = ey * eX;
   Q(0, 2) = ez * eX;
   Q(1, 0) = ex * eY;
   Q(1, 1) = ey * eY;
   Q(1, 2) = ez * eY;
   Q(2, 0) = ex * eZ;
   Q(2, 1) = ey * eZ;
   Q(2, 2) = ez * eZ;
   return;
}

xtensor::xPoint frameChange::pushForward(const xtensor::xPoint& X) const
{
   xtensor::xVector<> x(X(0) - center(0), X(1) - center(1), X(2) - center(2));
   xtensor::xVector<> r = Q * x;
   return xtensor::xPoint(r(0), r(1), r(2));
};

xtensor::xVector<> frameChange::pushForward(const xtensor::xVector<>& U) const { return Q * U; };

xtensor::xTensor2<> frameChange::pushForward(const xtensor::xTensor2<>& S) const { return Q * S * (!Q); };

xtensor::xPoint frameChange::pullBack(const xtensor::xPoint& x) const
{
   xtensor::xVector<> u(x(0), x(1), x(2));
   xtensor::xVector<> U = (!Q) * u;
   return xtensor::xPoint(U(0) + center(0), U(1) + center(1), U(2) + center(2));
};
xtensor::xVector<> frameChange::pullBack(const xtensor::xVector<>& u) const { return (!Q) * u; };

xtensor::xTensor2<> frameChange::pullBack(const xtensor::xTensor2<>& s) const { return (!Q) * s * Q; };
