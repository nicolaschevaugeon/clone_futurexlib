/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "FunctionCrackXFEM.h"

#include <algorithm>
#include <cmath>

// xfem
#include "xElement.h"
#include "xField.h"
#include "xGeomElem.h"
#include "xValKey.h"

// xcrack
#include "lCrack.h"

using namespace xfem;
using namespace std;
using AOMD::mEntity;

namespace xcrack
{
#define debug (false)
#define PI 3.141592654
#define SCAL 1. / sqrt(2. * PI)

ScalarFunctionCrackXFEM_c::ScalarFunctionCrackXFEM_c(const xcCrackBase& crk) : crack(crk) {}

void ScalarFunctionCrackXFEM_c::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   res = (double)crack.sideOf(geo_appro, geo_integ);
   storage.store(geo_appro, geo_integ, res);
}

void ScalarFunctionCrackXFEM_c::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                        xtensor::xVector<>& res) const
{
   res[0] = 0.0;
   res[1] = 0.0;
   res[2] = 0.0;
}

bool xcVectorFunctionCrack::evaltrick = false;
xcVectorFunctionCrack::xcVectorFunctionCrack(const xcCrackBase& crk, double nu)
    : crack(crk),
      kappa(3. - 4. * nu),
      scalfac(1. / sqrt(2. * PI)){

          // std::cout << "bl " << scalfac;
      };

void xcVectorFunctionCrack::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                   xtensor::xVector<>& res) const
{
   // if (exist(geo_appro, geo_integ, res)) return;
   getValLocalAxis(geo_appro, geo_integ, res);
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   xtensor::xVector<> e1, e2, e3, u;
   xtensor::xTensor2<> ge1, ge2, ge3;
   crack.getLocalCurv(e, uvw, e1, e2, e3, ge1, ge2, ge3);
   axisChangeCurv T(e1, e2, e3, ge1, ge2, ge3);
   res = T.localToGlobal(res);
   res *= SCAL;
   // store(geo_appro, geo_integ, res);
   return;
};

void xcVectorFunctionCrack::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>& res,
                                   const xValKey& basekey) const
{
   mEntity* enrichedEntity = basekey.getEnti();
   // if (exist(geo_appro, geo_integ, res)) return;
   if (evaltrick)
   {
      // std::cout << "in eval trick " << std::endl;
      res(0) = 1;
      res(1) = 1;
      res(2) = 1;
      return;
   }
   getValLocalAxis(geo_appro, geo_integ, res);

   xtensor::xTensor2<> Q;
   crack.getCrackOrthoAxis(enrichedEntity, Q);
   Q.transpose();
   res = Q * res;

   /*
   xtensor::xVector<> e1, e2, e3;
   crack.getLocalSmoothOrthoAxis(enrichedEntity, e1, e2, e3);
   axisChange T(e1, e2, e3);
   res = T.localToGlobal(res);*/

   res *= SCAL;
   // store(geo_appro, geo_integ, res);
   return;
};

void xcVectorFunctionCrack::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                    xtensor::xTensor2<>& res) const
{
   // if (exist(geo_appro, geo_integ, res)) return;
   getGradLocalAxis(geo_appro, geo_integ, res);
   xtensor::xVector<> e1, e2, e3, u;
   xtensor::xTensor2<> ge1, ge2, ge3;
   getValLocalAxis(geo_appro, geo_integ, u);
   crack.getLocalCurv(geo_appro->getEntity(), geo_appro->getUVW(), e1, e2, e3, ge1, ge2, ge3);
   axisChangeCurv T(e1, e2, e3, ge1, ge2, ge3);
   res = T.localToGlobal(u, res);
   res *= SCAL;

   // store(geo_appro, geo_integ, res);
   return;
};

void xcVectorFunctionCrack::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xTensor2<>& res,
                                    const xValKey& basekey) const
{
   // if (exist(geo_appro, geo_integ, res)) return;
   mEntity* enrichedEntity = basekey.getEnti();
   getGradLocalAxis(geo_appro, geo_integ, res);
   //  xtensor::xVector<> e1, e2, e3;
   // crack.getLocalSmoothOrthoAxis(enrichedEntity, e1, e2, e3);
   // xtensor::xTensor2<> Q=xtensor::xTensor2<>(e1,e2,e3);
   xtensor::xTensor2<> Q;
   crack.getCrackOrthoAxis(enrichedEntity, Q);
   Q.transpose();
   xtensor::xTensor2<> G;
   crack.getCrackAxis(geo_appro->getEntity(), G);
   // xtensor::xVector<> g1, g2, g3;
   // g1=crack.getFieldt()->getGrad(geo_appro->getEntity());
   // g2=crack.getFieldn()->getGrad(geo_appro->getEntity());
   // g3=g1%g2;
   // xtensor::xTensor2<> G(g1,g2,g3);

   res = Q * res * G;

   res *= SCAL;

   // store(geo_appro, geo_integ, res);
   // getGradNum(geo_appro, geo_integ, res, enrichedEntity);
}

void xcVectorFunctionCrack::getGradNum(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                       xtensor::xTensor2<>& res, const xValKey& basekey) const
{
   xtensor::xVector<> u1, u2;
   double eps = 0.00001;
   xfem::xGeomElem geo_appro_copy(*geo_appro);
   xfem::xGeomElem geo_integ_copy(*geo_integ);
   xtensor::xPoint starting = geo_appro_copy.getXYZ();
   getVal(&geo_appro_copy, &geo_integ_copy, u1, basekey);
   xtensor::xPoint DX(starting);
   DX(0) += eps;
   geo_appro_copy.setUVWForXYZ(DX);
   geo_integ_copy.setUVWForXYZ(DX);
   getVal(&geo_appro_copy, &geo_integ_copy, u2, basekey);
   xtensor::xVector<> dudx = (u2 - u1) / eps;
   DX(0) -= eps;
   DX(1) += eps;
   geo_appro_copy.setUVWForXYZ(DX);
   geo_integ_copy.setUVWForXYZ(DX);
   getVal(&geo_appro_copy, &geo_integ_copy, u2, basekey);
   xtensor::xVector<> dudy = (u2 - u1) / eps;
   res(0, 0) = dudx(0);
   res(0, 1) = dudy(0);
   res(1, 0) = dudx(1);
   res(1, 1) = dudy(1);
   res(0, 2) = res(1, 2) = res(2, 2) = res(0, 2) = res(1, 2) = 0.;
   geo_appro_copy.setUVWForXYZ(starting);
   geo_integ_copy.setUVWForXYZ(starting);
   return;
};

xcVectorFunctionCrack1::xcVectorFunctionCrack1(const xcCrackBase& crk, double nu) : xcVectorFunctionCrack(crk, nu){};

void xcVectorFunctionCrack1::getValLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                             xtensor::xVector<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   double coshalfth = cos(0.5 * th);
   double sinhalfth = sin(0.5 * th);
   // double costh = cos(th);
   // double sinth = sin(th);
   double sqrtr = sqrt(r);
   res[0] = sqrtr * coshalfth * (kappa - 1. + 2. * sinhalfth * sinhalfth);
   res[1] = sqrtr * sinhalfth * (kappa + 1. - 2. * coshalfth * coshalfth);
   res[2] = 0.;
   return;
};

void xcVectorFunctionCrack1::getGradLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                              xtensor::xTensor2<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double dthbdth = side*(th>0?1:-1);
   double coshalfth = cos(0.5 * th);
   double sinhalfth = sin(0.5 * th);
   double sqrtr = sqrt(r);

   double drdx = cos(th);
   double drdy = sin(th);
   double dthdx = -sin(th) / r;
   double dthdy = cos(th) / r;

   double duxdr = 0.5 / sqrtr * coshalfth * (kappa - 1. + 2 * sinhalfth * sinhalfth);
   double duxdth = sqrtr * (-0.5 * sinhalfth * (kappa - 1. + 2 * sinhalfth * sinhalfth) + 2 * coshalfth * coshalfth * sinhalfth);
   double duydr = 0.5 / sqrtr * sinhalfth * (kappa + 1. - 2. * coshalfth * coshalfth);
   double duydth = sqrtr * (0.5 * coshalfth * (kappa + 1. - 2 * coshalfth * coshalfth) + 2 * sinhalfth * sinhalfth * coshalfth);

   double duxdx = duxdr * drdx + duxdth * dthdx;
   double duxdy = duxdr * drdy + duxdth * dthdy;
   double duydx = duydr * drdx + duydth * dthdx;
   double duydy = duydr * drdy + duydth * dthdy;

   res(0, 0) = duxdx;
   res(0, 1) = duxdy;
   res(0, 2) = 0.;
   res(1, 0) = duydx;
   res(1, 1) = duydy;
   res(1, 2) = 0.;
   res(2, 0) = 0.;
   res(2, 1) = 0.;
   res(2, 2) = 0.;
   return;
}

xcVectorFunctionCrack1r32::xcVectorFunctionCrack1r32(const xcCrackBase& crk, double nu) : xcVectorFunctionCrack(crk, nu) {}

void xcVectorFunctionCrack1r32::getValLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                xtensor::xVector<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   // double costh = cos(th);
   // double sinth = sin(th);
   double sqrtr = sqrt(r);
   res[0] = sqrtr * r * ((kappa + 3. / 2. - 1.) * cos(3. * th / 2.) - 3. / 2. * cos((3. / 2. - 2.) * th));
   res[1] = sqrtr * r * ((kappa - 3. / 2. + 1.) * sin(3. * th / 2.) + 3. / 2. * sin((3. / 2. - 2.) * th));
   res[2] = 0.;
   return;
};

void xcVectorFunctionCrack1r32::getGradLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                 xtensor::xTensor2<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double dthbdth = side*(th>0?1:-1);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   double sqrtr = sqrt(r);

   double drdx = cos(th);
   double drdy = sin(th);
   double dthdx = -sin(th) / r;
   double dthdy = cos(th) / r;

   double duxdr = 3. / 2. * sqrtr * ((kappa + 3. / 2. - 1.) * cos(3. * th / 2.) - 3. / 2. * cos((3. / 2. - 2) * th));
   double duxdth =
       sqrtr * r * (-3. / 2. * (kappa + 3. / 2. - 1.) * sin(3. * th / 2.) + 3. / 2. * (3. / 2. - 2.) * sin((3. / 2. - 2.) * th));
   double duydr = 3. / 2. * sqrtr * ((kappa - 3. / 2. + 1.) * sin(3. * th / 2.) + 3. / 2. * sin((3. / 2. - 2.) * th));
   double duydth =
       sqrtr * r * (3. / 2. * (kappa - 3. / 2. + 1.) * cos(3. * th / 2.) + 3. / 2. * (3. / 2. - 2.) * cos((3. / 2. - 2) * th));

   double duxdx = duxdr * drdx + duxdth * dthdx;
   double duxdy = duxdr * drdy + duxdth * dthdy;
   double duydx = duydr * drdx + duydth * dthdx;
   double duydy = duydr * drdy + duydth * dthdy;

   res(0, 0) = duxdx;
   res(0, 1) = duxdy;
   res(0, 2) = 0.;
   res(1, 0) = duydx;
   res(1, 1) = duydy;
   res(1, 2) = 0.;
   res(2, 0) = 0.;
   res(2, 1) = 0.;
   res(2, 2) = 0.;
   return;
}

xcVectorFunctionCrack1r::xcVectorFunctionCrack1r(const xcCrackBase& crk, double nu) : xcVectorFunctionCrack(crk, nu) {}
void xcVectorFunctionCrack1r::getValLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                              xtensor::xVector<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   // double costh = cos(th);
   // double sinth = sin(th);
   // double sqrtr = sqrt(r);
   res[0] = r * ((kappa + 2.) * cos(th) - cos(-1 * th));
   res[1] = r * ((kappa - 2.) * sin(th) + sin(-1. * th));
   res[2] = 0.;
   return;
};

void xcVectorFunctionCrack1r::getGradLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                               xtensor::xTensor2<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double dthbdth = side*(th>0?1:-1);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   // double sqrtr = sqrt(r);

   double drdx = cos(th);
   double drdy = sin(th);
   double dthdx = -sin(th) / r;
   double dthdy = cos(th) / r;

   double duxdr = ((kappa + 2.) * cos(th) - cos(-1 * th));
   double duxdth = r * (-(kappa + 2.) * sin(th) - sin(-1 * th));
   double duydr = ((kappa - 2.) * sin(th) + sin(-1. * th));
   double duydth = r * ((kappa - 2.) * cos(th) - cos(-1. * th));

   double duxdx = duxdr * drdx + duxdth * dthdx;
   double duxdy = duxdr * drdy + duxdth * dthdy;
   double duydx = duydr * drdx + duydth * dthdx;
   double duydy = duydr * drdy + duydth * dthdy;

   res(0, 0) = duxdx;
   res(0, 1) = duxdy;
   res(0, 2) = 0.;
   res(1, 0) = duydx;
   res(1, 1) = duydy;
   res(1, 2) = 0.;
   res(2, 0) = 0.;
   res(2, 1) = 0.;
   res(2, 2) = 0.;
   return;
}

xcVectorFunctionCrack2::xcVectorFunctionCrack2(const xcCrackBase& crk, double nu) : xcVectorFunctionCrack(crk, nu) {}

void xcVectorFunctionCrack2::getValLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                             xtensor::xVector<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   double coshalfth = cos(0.5 * th);
   double sinhalfth = sin(0.5 * th);
   // double costh = cos(th);
   // double sinth = sin(th);
   double sqrtr = sqrt(r);
   res[0] = sqrtr * sinhalfth * (kappa + 1. + 2. * coshalfth * coshalfth);
   res[1] = -sqrtr * coshalfth * (kappa - 1. - 2. * sinhalfth * sinhalfth);
   res[2] = 0.;
   return;
};

void xcVectorFunctionCrack2::getGradLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                              xtensor::xTensor2<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double dthbdth = side*(th>0?1:-1);
   double coshalfth = cos(0.5 * th);
   double sinhalfth = sin(0.5 * th);
   double sqrtr = sqrt(r);

   double drdx = cos(th);
   double drdy = sin(th);
   double dthdx = -sin(th) / r;
   double dthdy = cos(th) / r;

   double duxdr = 0.5 / sqrtr * sinhalfth * (kappa + 1. + 2. * coshalfth * coshalfth);
   double duxdth = sqrtr * (0.5 * coshalfth * (kappa + 1. + 2. * coshalfth * coshalfth) - 2 * sinhalfth * sinhalfth * coshalfth);
   double duydr = -0.5 / sqrtr * coshalfth * (kappa - 1. - 2. * sinhalfth * sinhalfth);
   double duydth = -sqrtr * (-0.5 * sinhalfth * (kappa - 1. - 2 * sinhalfth * sinhalfth) - 2 * coshalfth * coshalfth * sinhalfth);

   double duxdx = duxdr * drdx + duxdth * dthdx;
   double duxdy = duxdr * drdy + duxdth * dthdy;
   double duydx = duydr * drdx + duydth * dthdx;
   double duydy = duydr * drdy + duydth * dthdy;

   res(0, 0) = duxdx;
   res(0, 1) = duxdy;
   res(0, 2) = 0.;
   res(1, 0) = duydx;
   res(1, 1) = duydy;
   res(1, 2) = 0.;
   res(2, 0) = 0.;
   res(2, 1) = 0.;
   res(2, 2) = 0.;
   return;
}

xcVectorFunctionCrack2r32::xcVectorFunctionCrack2r32(const xcCrackBase& crk, double nu) : xcVectorFunctionCrack(crk, nu) {}

void xcVectorFunctionCrack2r32::getValLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                xtensor::xVector<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   // double costh = cos(th);
   // double sinth = sin(th);
   double sqrtr = sqrt(r);
   res[0] = sqrtr * r * (-(kappa + 3. / 2. + 1.) * sin(3. * th / 2.) + (3. / 2.) * sin((3. / 2. - 2.) * th));
   res[1] = sqrtr * r * ((kappa - 3. / 2. - 1.) * cos(3. * th / 2.) + (3. / 2.) * cos((3. / 2. - 2.) * th));
   res[2] = 0.;
   return;
};

void xcVectorFunctionCrack2r32::getGradLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                                 xtensor::xTensor2<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double dthbdth = side*(th>0?1:-1);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   double sqrtr = sqrt(r);

   double drdx = cos(th);
   double drdy = sin(th);
   double dthdx = -sin(th) / r;
   double dthdy = cos(th) / r;

   double duxdr = 3. / 2. * sqrtr * (-(kappa + 3. / 2. + 1.) * sin(3. * th / 2.) + (3. / 2.) * sin((3. / 2. - 2.) * th));
   double duxdth =
       sqrtr * r *
       (-(kappa + 3. / 2. + 1.) * (3. / 2.) * cos(3. * th / 2.) + (3. / 2.) * (3. / 2. - 2.) * cos((3. / 2. - 2.) * th));
   double duydr = 3. / 2. * sqrtr * ((kappa - 3. / 2. - 1.) * cos(3. * th / 2.) + 3. / 2. * cos((3. / 2. - 2.) * th));
   double duydth =
       sqrtr * r *
       (-(3. / 2.) * (kappa - 3. / 2. - 1.) * sin(3. * th / 2.) - (3. / 2.) * (3. / 2. - 2.) * sin((3. / 2. - 2.) * th));

   double duxdx = duxdr * drdx + duxdth * dthdx;
   double duxdy = duxdr * drdy + duxdth * dthdy;
   double duydx = duydr * drdx + duydth * dthdx;
   double duydy = duydr * drdy + duydth * dthdy;

   res(0, 0) = duxdx;
   res(0, 1) = duxdy;
   res(0, 2) = 0.;
   res(1, 0) = duydx;
   res(1, 1) = duydy;
   res(1, 2) = 0.;
   res(2, 0) = 0.;
   res(2, 1) = 0.;
   res(2, 2) = 0.;
   return;
}

xcVectorFunctionCrack2r::xcVectorFunctionCrack2r(const xcCrackBase& crk, double nu) : xcVectorFunctionCrack(crk, nu) {}

void xcVectorFunctionCrack2r::getValLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                              xtensor::xVector<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   // double costh = cos(th);
   // double sinth = sin(th);
   // double sqrtr = sqrt(r);
   res[0] = r * (-(kappa)*sin(th) + sin(-1. * th));
   res[1] = r * ((kappa)*cos(th) + cos(-1. * th));
   res[2] = 0.;
   return;
};

void xcVectorFunctionCrack2r::getGradLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                               xtensor::xTensor2<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double dthbdth = side*(th>0?1:-1);
   // double coshalfth = cos(0.5*th);
   // double sinhalfth = sin(0.5*th);
   // double sqrtr = sqrt(r);

   double drdx = cos(th);
   double drdy = sin(th);
   double dthdx = -sin(th) / r;
   double dthdy = cos(th) / r;

   double duxdr = (-(kappa)*sin(th) + sin(-1. * th));
   double duxdth = r * (-(kappa)*cos(th) - cos(-1. * th));
   double duydr = ((kappa)*cos(th) + cos(-1. * th));
   double duydth = r * (-(kappa)*sin(th) + sin(-1. * th));

   double duxdx = duxdr * drdx + duxdth * dthdx;
   double duxdy = duxdr * drdy + duxdth * dthdy;
   double duydx = duydr * drdx + duydth * dthdx;
   double duydy = duydr * drdy + duydth * dthdy;

   res(0, 0) = duxdx;
   res(0, 1) = duxdy;
   res(0, 2) = 0.;
   res(1, 0) = duydx;
   res(1, 1) = duydy;
   res(1, 2) = 0.;
   res(2, 0) = 0.;
   res(2, 1) = 0.;
   res(2, 2) = 0.;
   return;
}

xcVectorFunctionCrack3::xcVectorFunctionCrack3(const xcCrackBase& crk, double nu) : xcVectorFunctionCrack(crk, nu){};

void xcVectorFunctionCrack3::getValLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                             xtensor::xVector<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   double sinhalfth = sin(0.5 * th);
   double sqrtr = sqrt(r);
   res[0] = 0.;
   res[1] = 0.;
   res[2] = sqrtr * sinhalfth;
   return;
};

void xcVectorFunctionCrack3::getGradLocalAxis(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                              xtensor::xTensor2<>& res) const
{
   double r, th, side;
   xtensor::xPoint uvw = geo_appro->getUVW();
   mEntity* e = geo_appro->getEntity();
   crack.getLocalCoords(e, uvw, r, th);
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   // double dthbdth = side*(th>0?1:-1);
   double coshalfth = cos(0.5 * th);
   double sinhalfth = sin(0.5 * th);
   double sqrtr = sqrt(r);

   double drdx = cos(th);
   double drdy = sin(th);
   double dthdx = -sin(th) / r;
   double dthdy = cos(th) / r;

   double duzdr = 0.5 / sqrtr * sinhalfth;
   double duzdth = 0.5 * sqrtr * coshalfth;

   double duzdx = duzdr * drdx + duzdth * dthdx;
   double duzdy = duzdr * drdy + duzdth * dthdy;

   res(0, 0) = 0.;
   res(0, 1) = 0.;
   res(0, 2) = 0.;
   res(1, 0) = 0.;
   res(1, 1) = 0.;
   res(1, 2) = 0.;
   res(2, 0) = duzdx;
   res(2, 1) = duzdy;
   res(2, 2) = 0.;
   return;
}

bool ScalarFunctionCrackTip1XFEM_c::evaltrick = false;
ScalarFunctionCrackTip1XFEM_c::ScalarFunctionCrackTip1XFEM_c(const xcCrackBase& crk) : crack(crk) {}

void ScalarFunctionCrackTip1XFEM_c::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   if (evaltrick)
   {
      res = 1.;
      return;
   }
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip1 uvw " << uvw << " xyz " << xyz;
   }

   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   res = sin(0.5 * th) * sqrt(r);
   //  res = sin(0.5*th)*cos(th)*sqrt(r);

   if (debug)
   {
      cout << " fct=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

void ScalarFunctionCrackTip1XFEM_c::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                            xtensor::xVector<>& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip1 uvw " << uvw << " xyz " << xyz;
   }
   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = -side * fabs(th);
   res[0] = sin(0.5 * th) / (2 * sqrt(r));
   res[1] = cos(0.5 * th) / (2 * sqrt(r));
   //  res[0]=(sin(0.5*th)-sin(th)*cos(1.5*th))/(2*sqrt(r));
   //  res[1]=(cos(th)*cos(1.5*th))/(2*sqrt(r));

   res[2] = 0.0;
   crack.localToGlobal(e, uvw, res);
   if (debug)
   {
      cout << " grad=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

bool ScalarFunctionCrackTip2XFEM_c::evaltrick = false;
ScalarFunctionCrackTip2XFEM_c::ScalarFunctionCrackTip2XFEM_c(const xcCrackBase& crk) : crack(crk) {}

void ScalarFunctionCrackTip2XFEM_c::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   if (evaltrick)
   {
      res = 1.;
      return;
   }
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip2 uvw " << uvw << " xyz " << xyz;
   }
   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   res = (sqrt(r) * cos(0.5 * th));
   //  res = sqrt(r)*cos(0.5*th)*cos(th);

   if (debug)
   {
      cout << " fct=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

void ScalarFunctionCrackTip2XFEM_c::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                            xtensor::xVector<>& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip2 uvw " << uvw << " xyz " << xyz;
   }
   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   res[0] = cos(0.5 * th) / (2 * sqrt(r));
   res[1] = sin(0.5 * th) / (2 * sqrt(r));

   //  res[0]=(cos(0.5*th)+sin(th)*sin(1.5*th))/(2*sqrt(r));
   //  res[1]=-sin(1.5*th)*cos(th)/(2*sqrt(r));

   res[2] = 0.0;
   crack.localToGlobal(e, uvw, res);
   if (debug)
   {
      cout << " grad=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

bool ScalarFunctionCrackTip3XFEM_c::evaltrick = false;
ScalarFunctionCrackTip3XFEM_c::ScalarFunctionCrackTip3XFEM_c(const xcCrackBase& crk) : crack(crk) {}

void ScalarFunctionCrackTip3XFEM_c::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   if (evaltrick)
   {
      res = 1.;
      return;
   }
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip3 uvw " << uvw << " xyz " << xyz;
   }
   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   res = sin(0.5 * th) * sin(th) * sqrt(r);
   if (debug)
   {
      cout << " fct=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

void ScalarFunctionCrackTip3XFEM_c::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                            xtensor::xVector<>& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip3 uvw " << uvw << " xyz " << xyz;
   }
   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   res[0] = -(sin(th) * sin(1.5 * th)) / (2 * sqrt(r));
   res[1] = (sin(0.5 * th) + sin(1.5 * th) * cos(th)) / (2 * sqrt(r));
   res[2] = 0.0;
   crack.localToGlobal(e, uvw, res);
   if (debug)
   {
      cout << " grad=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

bool ScalarFunctionCrackTip4XFEM_c::evaltrick = false;
ScalarFunctionCrackTip4XFEM_c::ScalarFunctionCrackTip4XFEM_c(const xcCrackBase& crk) : crack(crk) {}

void ScalarFunctionCrackTip4XFEM_c::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   if (evaltrick)
   {
      res = 1.;
      return;
   }
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip4 uvw " << uvw << " xyz " << xyz;
   }
   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   res = sqrt(r) * sin(th) * cos(0.5 * th);
   //  res = (sqrt(r)*cos(0.5*th));
   if (debug)
   {
      cout << " fct=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

void ScalarFunctionCrackTip4XFEM_c::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                            xtensor::xVector<>& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   xtensor::xPoint xyz = geo_appro->getXYZ();
   if (debug)
   {
      e->print();
      cout << "tip4 uvw " << uvw << " xyz " << xyz;
   }
   double r, th, side;
   crack.getLocalCoords(e, uvw, r, th);
   if (debug)
   {
      cout << " r " << r << " th " << th;
   }
   side = crack.sideOf(geo_appro, geo_integ);
   th = side * fabs(th);

   res[0] = -sin(th) * cos(1.5 * th) / (2 * sqrt(r));
   res[1] = (cos(0.5 * th) + cos(1.5 * th) * cos(th)) / (2 * sqrt(r));

   //  res[0]=cos(0.5*th)/(2*sqrt(r));
   //  res[1]=sin(0.5*th)/(2*sqrt(r));

   res[2] = 0.0;
   crack.localToGlobal(e, uvw, res);
   if (debug)
   {
      cout << " grad=" << res << endl;
   }
   storage.store(geo_appro, geo_integ, res);
}

EnrichMod_c::EnrichMod_c(lCrack& crk, double r0p, double attp) : crack(&crk), r0(r0p), att(attp) {}

void EnrichMod_c::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   // xtensor::xPoint xyz=geo_appro->getXYZ();

   double r, th;
   crack->getLocalCoords(e, uvw, r, th);
   if ((r < r0) && (r > 0))
      res = 1 - att * exp(4 + (r0 * r0) / (r * (r - r0)));
   else
      res = 1.0;

   storage.store(geo_appro, geo_integ, res);
}

void EnrichMod_c::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();
   xtensor::xPoint uvw = geo_appro->getUVW();
   // xtensor::xPoint xyz=geo_appro->getXYZ();
   double r, th, side;
   crack->getLocalCoords(e, uvw, r, th);
   side = crack->sideOf(geo_appro, geo_integ);
   th = side * fabs(th);
   double gg;
   if ((r > 0) && (r < r0))
      gg = att * r0 * r0 * ((2 * r - r0) / (r * r * (r - r0) * (r - r0))) * exp(4 + (r0 * r0) / (r * (r - r0)));
   else
      gg = 0.;

   res[0] = cos(th) * gg;
   res[1] = sin(th) * gg;
   res[2] = 0.0;
   crack->localToGlobal(e, uvw, res);
   storage.store(geo_appro, geo_integ, res);
}

axisChangeCurv::axisChangeCurv(const xtensor::xVector<>& _e0, const xtensor::xVector<>& _e1, const xtensor::xVector<>& _e2,
                               const xtensor::xTensor2<>& _ge0, const xtensor::xTensor2<>& _ge1, const xtensor::xTensor2<>& _ge2)
    : e0(_e0), e1(_e1), e2(_e2), ge0(_ge0), ge1(_ge1), ge2(_ge2)
{
   Q = xtensor::xTensor2<>(e0, e1, e2);
   Q.transpose();
   Qinv = Q.invert();
}

xtensor::xVector<> axisChangeCurv::localToGlobal(const xtensor::xVector<>& v) const { return Q * v; }

xtensor::xTensor2<> axisChangeCurv::localToGlobal(const xtensor::xVector<>& v, const xtensor::xTensor2<>& gv) const
{
   return Q * gv * Qinv + ge0 * v(0) + ge1 * v(1) + ge2 * v(2);
}

axisChange::axisChange(const xtensor::xVector<>& _e0, const xtensor::xVector<>& _e1, const xtensor::xVector<>& _e2)
    : e0(_e0), e1(_e1), e2(_e2)
{
   Q = xtensor::xTensor2<>(e0, e1, e2);
   Q.transpose();
   Qinv = Q.invert();
   // std::cout << Q << std::endl;
}

xtensor::xVector<> axisChange::localToGlobal(const xtensor::xVector<>& v) const { return Q * v; }

xtensor::xTensor2<> axisChange::localToGlobal(const xtensor::xTensor2<>& gv) const { return Q * gv * Qinv; }

xcCrackHeavisideEnrichment::xcCrackHeavisideEnrichment(const xcCrackBase& _crack) : crack(_crack) {}

void xcCrackHeavisideEnrichment::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   res = crack.sideOf(geo_appro, geo_integ);
   return;
}

void xcCrackHeavisideEnrichment::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                         xtensor::xVector<>& res) const
{
   for (int i = 0; i < 3; ++i) res[i] = 0.;
}

#undef debug

}  // namespace xcrack
