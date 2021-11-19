/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xcFrontSpaces.h"

#include <complex>

// xfem
#include "xLevelSet.h"

using AOMD::mEntity;
using namespace xfem;

xcSpaceModalConstant::xcSpaceModalConstant(const std::string& a, TensorialType_t b, mEntity* e, const xLevelSet* _ls)
    : xSpaceRegular(a, b), keyedentity(e), ls(_ls)
{
}

void xcSpaceModalConstant::getKeys(mEntity* e, femKeys* keys)
{
   keys->push_back(xValKey(Phys, xKeyInfo::getGeomId("MODALCONSTANT"), keyedentity));
}

void xcSpaceModalConstant::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   getKeys(e, keys);
   appro->push_back(shapeFctPtr(new xcApproxFunctionModalConstant(ls)));
}

void xcApproxFunctionModalConstant::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   if (ls)
   {
      mEntity* e = geo_appro->getEntity();
      xtensor::xPoint uvw = geo_appro->getUVW();
      double temp;
      bool paramdef;
      paramdef = ls->getVal(e, uvw, temp);
      if (paramdef)
         res = 1.;
      else
         res = 0.;
   }
   else
   {
      res = 1.;
   }
}

void xcApproxFunctionModalConstant::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                            xtensor::xVector<>& res) const
{
   res = xtensor::xVector<>(0., 0., 0.);
   return;
}

xcSpaceModalFourier::xcSpaceModalFourier(const std::string& a, TensorialType_t b, int nb_mod, const xLevelSet& _ls_cos_s,
                                         const xLevelSet& _ls_sin_s)
    : xSpaceRegular(a, b),
      nb_modes(nb_mod),
      ls_cos_s(_ls_cos_s),
      ls_sin_s(_ls_sin_s),
      modalfourier(ls_cos_s, ls_sin_s, nb_modes),
      MODE_SIN_KEY(nb_mod - 1),
      MODE_COS_KEY(nb_mod)
{
   for (int i = 1; i < nb_mod; ++i)
   {
      char str[12];
      std::sprintf(str, "%d", i);
      std::string name("MODE_SIN_" + std::string(str));
      MODE_SIN_KEY[i - 1] = xKeyInfo::getGeomId(name);
   }
   for (int i = 0; i < nb_mod; ++i)
   {
      char str[12];
      std::sprintf(str, "%d", i);
      std::string name("MODE_COS_" + std::string(str));
      MODE_COS_KEY[i] = xKeyInfo::getGeomId(name);
   }
}

void xcSpaceModalFourier::getKeys(mEntity* e, femKeys* keys)
{
   for (int i = 1; i < nb_modes; ++i)
   {
      keys->push_back(xValKey(Phys, MODE_SIN_KEY[i - 1], *ls_cos_s.getSupport().begin(0)));
   }
   for (int i = 0; i < nb_modes; ++i)
   {
      keys->push_back(xValKey(Phys, MODE_COS_KEY[i], *ls_sin_s.getSupport().begin(0)));
   }
}

void xcSpaceModalFourier::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   // warnign if e is a vertex this function does not have to much sense
   getKeys(e, keys);

   for (int i = 1; i < nb_modes; ++i)
   {
      // appro->push_back(shapeFctPtr(new xcApproxFunctionModalFourier(ls_cos_s, ls_sin_s, i,
      // xcApproxFunctionModalFourier::SIN)));
      appro->push_back(shapeFctPtr(new xcApproxFunctionModalFourier(modalfourier, i, xcApproxFunctionModalFourier::SIN)));
   }
   for (int i = 0; i < nb_modes; ++i)
   {
      // appro->push_back(shapeFctPtr(new xcApproxFunctionModalFourier(ls_cos_s, ls_sin_s ,  i,
      // xcApproxFunctionModalFourier::COS)));
      appro->push_back(shapeFctPtr(new xcApproxFunctionModalFourier(modalfourier, i, xcApproxFunctionModalFourier::COS)));
   }
}

void ModalFourier::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ) const
{
   if ((geo_integ->getEntity() == e_current_v) && (equal_uvw(geo_integ->getUVW(), uvw_current_v))) return;
   e_current_v = geo_integ->getEntity();
   uvw_current_v = geo_integ->getUVW();

   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   double costheta;
   bool cosexist = ls_cos_s.getVal(e, uvw, costheta);
   double sintheta;
   bool sinexist = ls_sin_s.getVal(e, uvw, sintheta);
   if ((!sinexist) || (!cosexist))
   {
      std::fill(eintheta.begin(), eintheta.end(), std::complex<double>(0., 0.));
      return;
   }
   eintheta[0] = std::complex<double>(1., 0.);
   std::complex<double> eitheta(costheta, sintheta);
   // std::cout <<" eitheta "   << eitheta << " " <<eintheta[0]<<" " <<eitheta*eintheta[0] << std::endl;
   for (int i = 1; i < nb_modes; ++i)
   {
      eintheta[i] = eintheta[i - 1] * eitheta;
      //  std::cout << i << " " << eintheta[i] <<std::endl;
   }
}

void ModalFourier::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ) const
{
   if (geo_integ->getEntity() == e_current_g && equal_uvw(geo_integ->getUVW(), uvw_current_g)) return;
   e_current_g = geo_integ->getEntity();
   uvw_current_g = geo_integ->getUVW();
   getVal(geo_appro, geo_integ);
   xtensor::xVector<> gcostheta, gsintheta;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   bool havegcostheta = ls_cos_s.getGrad(e, uvw, gcostheta);
   bool havegsintheta = ls_sin_s.getGrad(e, uvw, gsintheta);
   if (!(havegcostheta && havegsintheta))
   {
      std::fill(gx.begin(), gx.end(), std::complex<double>(0., 0.));
      std::fill(gy.begin(), gy.end(), std::complex<double>(0., 0.));
      std::fill(gz.begin(), gz.end(), std::complex<double>(0., 0.));
      return;
   }
   gx[0] = std::complex<double>(0., 0.);
   gy[0] = std::complex<double>(0., 0.);
   gz[0] = std::complex<double>(0., 0.);
   std::complex<double> gthetax = std::complex<double>(gcostheta(0), gsintheta(0));
   std::complex<double> gthetay = std::complex<double>(gcostheta(1), gsintheta(1));
   std::complex<double> gthetaz = std::complex<double>(gcostheta(2), gsintheta(2));
   for (int i = 1; i < nb_modes; ++i)
   {
      gx[i] = (i * 1.) * eintheta[i - 1] * gthetax;
      gy[i] = (i * 1.) * eintheta[i - 1] * gthetay;
      gz[i] = (i * 1.) * eintheta[i - 1] * gthetaz;
   }
}

xcApproxFunctionModalFourier::xcApproxFunctionModalFourier(const ModalFourier& _modalfourier, int i, cos_sin_t t)
    : modalfourier(_modalfourier), ith(i), cos_or_sin(t)
{
}

void xcApproxFunctionModalFourier::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   modalfourier.getVal(geo_appro, geo_integ);
   switch (cos_or_sin)
   {
      case SIN:
         res = imag(modalfourier.eintheta[ith]);
         break;
      case COS:
         res = real(modalfourier.eintheta[ith]);
         break;
   }
   // std::cout << res << std::endl;
   return;
}

void xcApproxFunctionModalFourier::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                           xtensor::xVector<>& res) const
{
   modalfourier.getGrad(geo_appro, geo_integ);
   switch (cos_or_sin)
   {
      case SIN:
         res(0) = imag(modalfourier.gx[ith]);
         res(1) = imag(modalfourier.gy[ith]);
         res(2) = imag(modalfourier.gz[ith]);
         break;
      case COS:
         res(0) = real(modalfourier.gx[ith]);
         res(1) = real(modalfourier.gy[ith]);
         res(2) = real(modalfourier.gz[ith]);
         break;
   }
   return;
}

/*
xcApproxFunctionModalFourier::xcApproxFunctionModalFourier(const xLevelSet& cos_s_modal_,const xLevelSet& sin_s_modal_,   int i,
cos_sin_t t) : cos_s_modal(cos_s_modal_),  sin_s_modal( sin_s_modal_),  ith(i), cos_or_sin(t) {}


        void  xcApproxFunctionModalFourier::getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double& res)
const
{
        mEntity* e = geo_appro->getEntity();
        xtensor::xPoint uvw=geo_appro->getUVW();
        double costheta;
        bool cosexist = cos_s_modal.getVal(e, uvw, costheta);
        double sintheta;
        bool sinexist = sin_s_modal.getVal(e, uvw, sintheta);
        if ((!sinexist) || (!cosexist)) {
          res = 0.;
          return;
        }

        std::complex <double > eitheta (costheta, sintheta);
        std::complex <double > eintheta (std::pow(eitheta, ith));

        switch (cos_or_sin)
        {
                case SIN :
                        res = imag(eintheta);
                        break;
                case COS :
                        res = real(eintheta);
                        break;
        }
        return;
}

void  xcApproxFunctionModalFourier::getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&
res) const
{
        if (ith == 0) {
                res = xtensor::xVector<>(0., 0. , 0.);
                return;
        }

        mEntity* e = geo_appro->getEntity();
        xtensor::xPoint uvw=geo_appro->getUVW();
        double costheta;
        bool cosexist = cos_s_modal.getVal(e, uvw, costheta);
        double sintheta;
        bool sinexist = sin_s_modal.getVal(e, uvw, sintheta);
        if ((!sinexist) || (!cosexist)) {
          res = xtensor::xVector<>(0.);
          return;
        }
        xtensor::xVector<> gcostheta = cos_s_modal.getGrad(e, uvw);
        xtensor::xVector<> gsintheta = sin_s_modal.getGrad(e, uvw);
        std::complex<double > eitheta(costheta , sintheta);
        std::complex<double > einm1thetha(pow (eitheta, ith-1));
        std::complex<double > gx = (ith*1.) *  einm1thetha * std::complex<double >( gcostheta(0), gsintheta(0) );
        std::complex<double > gy = (ith *1.)*  einm1thetha * std::complex<double >( gcostheta(1), gsintheta(1) );
        std::complex<double > gz = (ith *1.)  *einm1thetha * std::complex<double >( gcostheta(2), gsintheta(2) );

        switch (cos_or_sin)
          {
          case SIN :
            res(0) = imag(gx);
            res(1) = imag(gy);
            res(2) = imag(gz);
            break;
          case COS :
            res(0) = real(gx);
            res(1) = real(gy);
            res(2) = real(gz);
            break;
          }
        return;
}*/

xcApproxFunctionModalLegendre::xcApproxFunctionModalLegendre(const xLevelSet& s_modal_, int i) : s_modal(s_modal_), ith(i) {}

xcSpaceModalLegendre::xcSpaceModalLegendre(const std::string& a, TensorialType_t b, int nb_mod, const xLevelSet& ls)
    : xSpaceRegular(a, b), nb_modes(nb_mod), curvi_coord_on_front_part(ls), MODE_LEGENDRE_KEY(nb_mod)
{
   for (int i = 0; i < nb_mod; ++i)
   {
      char str[12];
      std::sprintf(str, "%d", i);
      std::string name("MODE_LEGENDRE_" + std::string(str));
      MODE_LEGENDRE_KEY[i] = xKeyInfo::getGeomId(name);
   }
}

void xcSpaceModalLegendre::getKeys(mEntity* e, femKeys* keys)
{
   // test
   const xRegion& region = curvi_coord_on_front_part.getSupport();
   // std::cout << "in Keys Legendre" << std::endl;
   // std::cout << region.IsInRegion(e) << " "<< region.size(0)<< " "<< region.size(1)<< " "<<region.size(2)<< " "<<
   // region.size(3)  << std::endl;
   // if (region.IsInRegion(e))
   {
      for (int i = 0; i < nb_modes; ++i)
      {
         keys->push_back(xValKey(Phys, MODE_LEGENDRE_KEY[i], *(region.begin(0))));
      }
   }
}

void xcSpaceModalLegendre::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   // warnign if e is a vertex this function does not have to much sense
   getKeys(e, keys);
   for (int i = 0; i < nb_modes; ++i)
   {
      appro->push_back(shapeFctPtr(new xcApproxFunctionModalLegendre(curvi_coord_on_front_part, i)));
   }
}

void xcApproxFunctionModalLegendre::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   double s;
   if (s_modal.getVal(e, uvw, s))
   {
      switch (ith)
      {
         case 0:
            res = 1.;
            break;
         case 1:
            res = s;
            break;
         case 2:
            res = 0.5 * (3. * s * s - 1.);
            break;
         case 3:
            res = 0.5 * (5. * pow(s, 3) - 3. * s);
            break;
         case 4:
            res = 0.125 * (35. * pow(s, 4) - 30. * s * s + 3.);
            break;
         case 5:
            res = 0.125 * (63. * pow(s, 5) - 70. * pow(s, 3) + 15. * s);
            break;
         case 6:
            res = (1. / 16.) * (231. * pow(s, 6) - 315. * pow(s, 4) + 105. * s * s - 5.);
            break;
         default:
            //	std::cerr << " above order 6 not implemented yet " << endl;
            throw;
            break;
      }
   }
   else
      res = 0.;
   return;
}

void xcApproxFunctionModalLegendre::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                            xtensor::xVector<>& res) const
{
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   double s;
   if (s_modal.getVal(e, uvw, s))
   {
      res = s_modal.getGrad(e);

      switch (ith)
      {
         case 0:
            res *= 0.;
            break;
         case 1:
            res *= 1.;
            break;
         case 2:
            res *= 3. * s;
            break;
         case 3:
            res *= 0.5 * (15. * pow(s, 2) - 3.);
            break;
         case 4:
            res *= 0.125 * (35. * 4. * pow(s, 3) - 60. * s);
            break;
         case 5:
            res *= 0.125 * (63. * 5 * pow(s, 4) - 70. * 3 * pow(s, 2) + 15.);
            break;
         case 6:
            res *= (1. / 16.) * (231. * 6. * pow(s, 5) - 315. * 4. * pow(s, 3) + 210. * s);
            break;
         default:
            //	std::cerr << " above order 6 not implemented yet " << endl;
            throw;
            break;
      }
   }
   else
      res = xtensor::xVector<>(0.);
   return;
}

xcApproxFunctionModalPolynomeHermiteUniform::xcApproxFunctionModalPolynomeHermiteUniform(const xLevelSet& s_modal_, int i,
                                                                                         int nb_mod, double _mode_width)
    : s_modal(s_modal_), s_modal_sin(s_modal_), ith(i), total_number_of_modes(nb_mod), line_open(true), mode_width(_mode_width)
{
}

xcApproxFunctionModalPolynomeHermiteUniform::xcApproxFunctionModalPolynomeHermiteUniform(const xLevelSet& _s_modal_cos,
                                                                                         const xLevelSet& _s_modal_sin, int i,
                                                                                         int nb_mod, double _mode_width)
    : s_modal(_s_modal_cos),
      s_modal_sin(_s_modal_sin),
      ith(i),
      total_number_of_modes(nb_mod),
      line_open(false),
      mode_width(_mode_width)
{
}

xcSpaceModalPolynomeHermiteUniform::xcSpaceModalPolynomeHermiteUniform(const std::string& a, TensorialType_t b, int nb_mod,
                                                                       const xLevelSet& ls, double _mode_width)
    : xSpaceRegular(a, b),
      nb_modes(nb_mod),
      curvi_coord_on_front_part(ls),
      curvi_coord_on_front_part_sin(ls),
      MODE_POLYNOME_HERMITE_UNIFORM_KEY(nb_mod),
      line_open(true),
      mode_width(_mode_width)
{
   for (int i = 0; i < nb_mod; ++i)
   {
      char str[12];
      std::sprintf(str, "%d", i);
      std::string name("MODE_POLYNOME_HERMITE_UNIFORM_" + std::string(str));
      MODE_POLYNOME_HERMITE_UNIFORM_KEY[i] = xKeyInfo::getGeomId(name);
   }
}

xcSpaceModalPolynomeHermiteUniform::xcSpaceModalPolynomeHermiteUniform(const std::string& a, TensorialType_t b, int nb_mod,
                                                                       const xLevelSet& ls, const xLevelSet& ls_sin,
                                                                       double _mode_width)
    : xSpaceRegular(a, b),
      nb_modes(nb_mod),
      curvi_coord_on_front_part(ls),
      curvi_coord_on_front_part_sin(ls_sin),
      MODE_POLYNOME_HERMITE_UNIFORM_KEY(nb_mod),
      line_open(false),
      mode_width(_mode_width)
{
   for (int i = 0; i < nb_mod; ++i)
   {
      char str[12];
      std::sprintf(str, "%d", i);
      std::string name("MODE_POLYNOME_HERMITE_UNIFORM_" + std::string(str));
      MODE_POLYNOME_HERMITE_UNIFORM_KEY[i] = xKeyInfo::getGeomId(name);
   }
}

void xcSpaceModalPolynomeHermiteUniform::getKeys(mEntity* e, femKeys* keys)
{
   const xRegion& region = curvi_coord_on_front_part.getSupport();
   {
      for (int i = 0; i < nb_modes; ++i)
      {
         keys->push_back(xValKey(Phys, MODE_POLYNOME_HERMITE_UNIFORM_KEY[i], *(region.begin(0))));
      }
   }
}

void xcSpaceModalPolynomeHermiteUniform::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   // warnign if e is a vertex this function does not have to much sense
   getKeys(e, keys);
   for (int i = 0; i < nb_modes; ++i)
   {
      if (line_open)
         appro->push_back(
             shapeFctPtr(new xcApproxFunctionModalPolynomeHermiteUniform(curvi_coord_on_front_part, i, nb_modes, mode_width)));
      else
         appro->push_back(shapeFctPtr(new xcApproxFunctionModalPolynomeHermiteUniform(
             curvi_coord_on_front_part, curvi_coord_on_front_part_sin, i, nb_modes, mode_width)));
   }
}

void xcApproxFunctionModalPolynomeHermiteUniform::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                         double& res) const
{
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   double s;
   bool paramdef;
   if (line_open)
   {
      paramdef = s_modal.getVal(e, uvw, s);
   }
   else
   {
      double cos, sin;
      paramdef = s_modal.getVal(e, uvw, cos);
      paramdef = ((paramdef) && (s_modal_sin.getVal(e, uvw, sin)));
      s = atan2(sin, cos) / M_PI;
   }
   if (!paramdef)
   {
      res = 0.;
      return;
   }
   // s belongs to [-1,1]
   // modes made of two parts: x in [-1,0] and in [0,1]
   // ensuring partition of unity

   // repartition uniforme
   double l = mode_width;
   double xmid = -1 + ith * l / 2.;
   double x0 = xmid - l / 2;
   double x1 = xmid + l / 2;
   bool compute = false;
   if ((line_open) || (ith != 0))
   {
      if ((s > x0) && (s < x1))
      {  // domain of definition of mode ith
         compute = true;
      }
   }
   else
   {
      // should use modulo...
      double x0_periodic = x0 + 2.;
      double x1_periodic = x1 + 2.;
      if ((s > x0) && (s < x1))
      {  // domain of definition of mode ith for lineclose
         compute = true;
      }
      if ((s > x0_periodic) && (s < x1_periodic))
      {  // periodic condition
         compute = true;
         s = s - 2.;
      }
   }

   if (compute)
   {
      double coefb = -3.;
      double coefd = 1.;
      double coefa;
      if (s <= xmid)  // first part
         coefa = -2.;
      else  // second part
         coefa = 2.;
      res = coefa * pow(2 / l * (s - xmid), 3) + coefb * pow(2 / l * (s - xmid), 2) + coefd;
   }
   else
      res = 0.;

   //		cout << " ith = " << ith << " xmid " << xmid << " s = " << s << " line_open = " << line_open << " compute = " <<
   // compute << endl;
}

void xcApproxFunctionModalPolynomeHermiteUniform::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                          xtensor::xVector<>& res) const
{
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   double s;
   bool paramdef, paramdef2;
   double l = mode_width;
   double xmid = -1 + ith * l / 2.;
   double x0 = xmid - l / 2;
   double x1 = xmid + l / 2;

   if (line_open)
   {
      paramdef = s_modal.getVal(e, uvw, s);
      if (!paramdef)
      {
         res = xtensor::xVector<>(0.);
         return;
      }
      paramdef2 = s_modal.getGrad(e, uvw, res);
      if (!paramdef2)
      {
         res = xtensor::xVector<>(0.);
         return;
      }
      // res = s_modal.getGrad(e);
      if ((s > x0) && (s < x1))
      {
         double coefb = -3.;
         double coefa;
         if (s <= xmid)  // first part
            coefa = -2.;
         else  // second part
            coefa = 2.;
         res *= coefa * 3 * 2 / l * pow(2 / l * (s - xmid), 2) + 2 * coefb * 2 / l * 2 / l * (s - xmid);
      }
      else
         res *= 0.;
   }

   else
   {
      double cos, sin;
      paramdef = s_modal.getVal(e, uvw, cos);
      paramdef = ((paramdef) && (s_modal_sin.getVal(e, uvw, sin)));
      s = atan2(sin, cos) / M_PI;
      if (!paramdef)
      {
         res = xtensor::xVector<>(0.);
         return;
      }
      if (ith == 0)
      {  // premier mode (periodique)
         if ((s > (x0 + 2)) && (s < (x1 + 2)))
         {  // periodic condition
            s = s - 2.;
         }
      }
      if (!((s > x0) && (s < x1))) res = xtensor::xVector<>(0.);
      return;
      double coefb = -3.;
      double coefa = s <= xmid ? -2. : 2.;  // assemblage de 2 function de hermite

      double dHds = coefa * 3 * 2 / l * pow(2 / l * (s - xmid), 2) + 2 * coefb * 2 / l * 2 / l * (s - xmid);

      xtensor::xVector<> gcos = s_modal.getGrad(e);
      xtensor::xVector<> gsin = s_modal_sin.getGrad(e);
      double dHdcos = -sin / M_PI * dHds;
      double dHdsin = cos / M_PI * dHds;
      res = gcos * dHdcos + gsin * dHdsin;
      return;
   }

   /*

 res = s_modal.getGrad(e);

 // repartition uniforme
 double l = mode_width;
 double xmid = -1+ith*l/2.;
 double x0 = xmid-l/2;
 double x1 = xmid+l/2;
 bool compute=false;
 if ((line_open)||(ith!=0)){
         if ((s>x0)&&(s<x1)){// domain of definition of mode ith
                 compute=true;
         }
 }
 else{
         // should use modulo...
         double x0_periodic = x0+2.;
         double x1_periodic = x1+2.;
         if ((s>x0)&&(s<x1)){// domain of definition of mode ith for lineclose
                 compute=true;
         }
         if ((s>x0_periodic)&&(s<x1_periodic)){// periodic condition
                 compute=true;
                 s = s-2.;
         }
 }

 if (compute){
         double coefb = -3.;
         double coefa;
         if (s<=xmid)// first part
                 coefa = -2.;
         else// second part
                 coefa = 2.;
         res *= 	coefa*3*2/l*pow(2/l*(s-xmid),2) + 2*coefb*2/l*2/l*(s-xmid);
 }
 else
         res = 0.;
   */
}
