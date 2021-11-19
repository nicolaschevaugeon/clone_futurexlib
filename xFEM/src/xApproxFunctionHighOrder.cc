/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xApproxFunctionHighOrder.h"

#include "xGeomElem.h"

using AOMD::mEntity;

namespace xfem
{
// double LagrangeScaled(const double coef, const std::vector<int > & n, const double *L);
double LagrangeScaled(const double coef, const int *n, const double *L, const int dim);

void dLagrangeScaled1(const int n, double &L, double &dL, double &fac);
void dLagrangeScaled2(const int *n, double *L, double *dL);
void dLagrangeScaled3(const int *n, double *L, double *dL);
void dLagrangeScaled4(const int *n, double *L, double *dL);

int nbShapeFunctionEdgeTotal(const int &order) { return order + 1; }

int nbShapeFunctionEdge(const int &order) { return std::max(order - 1, 0); }
int nbShapeFunctionTriTotal(const int &order) { return (order + 2) * (order + 1) / 2; }
int nbShapeFunctionTri(const int &order) { return (order > 2) ? (order - 2) * (order - 1) / 2 : 0; }

int nbShapeFunctionTetTotal(const int &order) { return (order + 3) * (order + 2) * (order + 1) / 6; }

int nbShapeFunctionTet(const int &order) { return (order > 3) ? (order - 3) * (order - 2) * (order - 1) / 6 : 0; }

long int factorial(long int n)
{
   if (n <= 1)
      return 1;
   else
      return n * factorial(n - 1);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// xSimplex Implementation                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////

xSimplexVertex::xSimplexVertex(int _node) : node(_node) {}

void xSimplexVertex::getL(const mEntity::mType &type, const xtensor::xPoint &uvw, double &L) const
{
   switch (type)
   {
      case mEntity::TET:
      {
         const double &u = uvw(0), v = uvw(1), w = uvw(2);
         switch (node)
         {
            case 0:
            {
               L = 1 - u - v - w;
               break;
            }
            case 1:
            {
               L = u;
               break;
            }
            case 2:
            {
               L = v;
               break;
            }
            case 3:
            {
               L = w;
               break;
            }
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::TRI:
      {
         const double &u = uvw(0), v = uvw(1);
         switch (node)
         {
            case 0:
            {
               L = 1 - u - v;
               break;
            }
            case 1:
            {
               L = u;
               break;
            }
            case 2:
            {
               L = v;
               break;
            }
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::EDGE:
      {
         const double &u = uvw(0);
         switch (node)
         {
            case 0:
            {
               L = (1 - u) * 0.5;
               break;
            }
            case 1:
            {
               L = (1 + u) * 0.5;
               break;
            }
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::VERTEX:
      {
         L = 1.;
         break;
      }
      default:
      {
         throw;
      }
   }
}

void xSimplexVertex::getdL(const mEntity::mType &type, double &dLdu, double &dLdv, double &dLdw) const
{
   switch (type)
   {
      case mEntity::TET:
      {
         switch (node)
         {
            case 0:
            {
               dLdu = -1;
               dLdv = -1;
               dLdw = -1;
               break;
            }
            case 1:
            {
               dLdu = 1;
               dLdv = 0;
               dLdw = 0;
               break;
            }
            case 2:
            {
               dLdu = 0;
               dLdv = 1;
               dLdw = 0;
               break;
            }
            case 3:
            {
               dLdu = 0;
               dLdv = 0;
               dLdw = 1;
               break;
            }
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::TRI:
      {
         switch (node)
         {
            case 0:
            {
               dLdu = -1;
               dLdv = -1;
               dLdw = 0;
               break;
            }
            case 1:
            {
               dLdu = 1;
               dLdv = 0;
               dLdw = 0;
               break;
            }
            case 2:
            {
               dLdu = 0;
               dLdv = 1;
               dLdw = 0;
               break;
            }
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::EDGE:
      {
         switch (node)
         {
            case 0:
            {
               dLdu = -0.5;
               dLdv = 0;
               dLdw = 0;
               break;
            }
            case 1:
            {
               dLdu = 0.5;
               dLdv = 0.;
               dLdw = 0;
               break;
            }
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::VERTEX:
      {
         dLdu = 0.;
         dLdv = 0;
         dLdw = 0;
         break;
      }
      default:
      {
         throw;
      }
   }
   return;
}

void xSimplexVertex::getLdL(const mEntity::mType &type, const xtensor::xPoint &uvw, double &L, double &dLdu, double &dLdv,
                            double &dLdw) const
{
   getL(type, uvw, L);
   getdL(type, dLdu, dLdv, dLdw);
   return;
}

xSimplexEdge::xSimplexEdge(int _edge) : edge(_edge) { lastpoint = nullptr; }
void xSimplexEdge::getL(const mEntity::mType &type, const xtensor::xPoint &uvw, double *L) const
{
   switch (type)
   {
      case mEntity::TET:
      {
         const double u = uvw(0), v = uvw(1), w = uvw(2);
         switch (edge)
         {
            case 0:
            {
               L[0] = 1. - u - v - w;
               L[1] = u;
               break;
            }  // 01
            case 1:
            {
               L[0] = u;
               L[1] = v;
               break;
            }  // 12
            case 2:
            {
               L[0] = v;
               L[1] = 1. - u - v - w;
               break;
            }  // 20
            case 3:
            {
               L[0] = 1. - u - v - w;
               L[1] = w;
               break;
            }  // 03
            case 4:
            {
               L[0] = u;
               L[1] = w;
               break;
            }  // 13
            case 5:
            {
               L[0] = v;
               L[1] = w;
               break;
            }  // 23
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::TRI:
      {
         const double u = uvw(0);
         const double v = uvw(1);
         switch (edge)
         {
            case 0:
            {
               L[0] = 1. - u - v;
               L[1] = u;
               break;
            }  // 01
            case 1:
            {
               L[0] = u;
               L[1] = v;
               break;
            }  // 12
            case 2:
            {
               L[0] = v;
               L[1] = 1. - u - v;
               break;
            }  // 21
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::EDGE:
      {
         const double u = uvw(0);
         L[0] = (1. - u) * 0.5;
         L[1] = (1. + u) * 0.5;
         break;
      }
      default:
      {
         throw;
      }
   }
   return;
}

void xSimplexEdge::getdL(const mEntity::mType &type, double *dLdu, double *dLdv, double *dLdw) const
{
   switch (type)
   {
      case mEntity::TET:
      {
         switch (edge)
         {
            case 0:
            {
               dLdu[0] = -1.;
               dLdu[1] = 1.;
               dLdv[0] = -1.;
               dLdv[1] = 0.;
               dLdw[0] = -1.;
               dLdw[1] = 0.;
               break;
            }  // 01
            case 1:
            {
               dLdu[0] = 1.;
               dLdu[1] = 0.;
               dLdv[0] = 0.;
               dLdv[1] = 1.;
               dLdw[0] = 0.;
               dLdw[1] = 0.;
               break;
            }  // 12
            case 2:
            {
               dLdu[0] = 0.;
               dLdu[1] = -1.;
               dLdv[0] = 1.;
               dLdv[1] = -1.;
               dLdw[0] = 0.;
               dLdw[1] = -1.;
               break;
            }  // 20
            case 3:
            {
               dLdu[0] = -1.;
               dLdu[1] = 0.;
               dLdv[0] = -1.;
               dLdv[1] = 0.;
               dLdw[0] = -1.;
               dLdw[1] = 1.;
               break;
            }  // 03
            case 4:
            {
               dLdu[0] = 1.;
               dLdu[1] = 0.;
               dLdv[0] = 0.;
               dLdv[1] = 0.;
               dLdw[0] = 0.;
               dLdw[1] = 1.;
               break;
            }  // 13
            case 5:
            {
               dLdu[0] = 0.;
               dLdu[1] = 0.;
               dLdv[0] = 1.;
               dLdv[1] = 0.;
               dLdw[0] = 0.;
               dLdw[1] = 1.;
               break;
            }  // 23
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::TRI:
      {
         switch (edge)
         {
            case 0:
            {
               dLdu[0] = -1.;
               dLdu[1] = 1.;
               dLdv[0] = -1.;
               dLdv[1] = 0.;
               dLdw[0] = 0.;
               dLdw[1] = 0.;
               break;
            }  // 01
            case 1:
            {
               dLdu[0] = 1.;
               dLdu[1] = 0.;
               dLdv[0] = 0.;
               dLdv[1] = 1.;
               dLdw[0] = 0.;
               dLdw[1] = 0.;
               break;
            }  // 12
            case 2:
            {
               dLdu[0] = 0.;
               dLdu[1] = -1.;
               dLdv[0] = 1.;
               dLdv[1] = -1.;
               dLdw[0] = 0.;
               dLdw[1] = 0.;
               break;
            }  // 20
            default:
            {
               throw;
            }
         }
         break;
      }
      case mEntity::EDGE:
      {
         dLdu[0] = -0.5;
         dLdu[1] = 0.5;
         dLdv[0] = 0.;
         dLdv[1] = 0.;
         dLdw[0] = 0.;
         dLdw[1] = 0.;
         break;
      }
      default:
      {
         throw;
      }
   }
   return;
}

void xSimplexEdge::getLdL(const mEntity::mType &type, const xtensor::xPoint &uvw, double *L, double *dLdu, double *dLdv,
                          double *dLdw) const
{
   getL(type, uvw, L);
   getdL(type, dLdu, dLdv, dLdw);
   return;
}

const int xSimplexTet::Tfv[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
const int xSimplexTet::Tev[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

xSimplexTri::xSimplexTri(int _face) : face(_face)
{
   if ((face < 0) || (face > 3)) throw;
}

void xSimplexTri::getL(const mEntity::mType &type, const xtensor::xPoint &uvw, double *L) const
{
   switch (type)
   {
      case mEntity::TET:
      {
         const double u = uvw(0), v = uvw(1), w = uvw(2);
         switch (face)
         {
            case 0:
            {
               L[0] = 1 - u - v - w;
               L[1] = u;
               L[2] = v;
               break;
            }  // 012
            case 1:
            {
               L[0] = 1 - u - v - w;
               L[1] = u;
               L[2] = w;
               break;
            }  // 013
            case 2:
            {
               L[0] = u;
               L[1] = v;
               L[2] = w;
               break;
            }  // 123
            case 3:
            {
               L[0] = 1 - u - v - w;
               L[1] = v;
               L[2] = w;
               break;
            }  // 023
            default:
               throw;
         }
         break;
      }
      case mEntity::TRI:
      {
         const double u = uvw(0);
         const double v = uvw(1);
         L[0] = 1 - u - v;
         L[1] = u;
         L[2] = v;
         break;
      }
      default:
      {
         throw;
      }
   }
}

void xSimplexTri::getLdL(const mEntity::mType &type, const xtensor::xPoint &uvw, double *L, double *dLdu, double *dLdv,
                         double *dLdw) const
{
   getL(type, uvw, L);
   getdL(type, dLdu, dLdv, dLdw);
}

void xSimplexTri::getdL(const mEntity::mType &type, double *dLdu, double *dLdv, double *dLdw) const
{
   switch (type)
   {
      case mEntity::TRI:
      {
         dLdu[0] = -1.;
         dLdu[1] = 1.;
         dLdu[2] = 0.;
         dLdv[0] = -1.;
         dLdv[1] = 0.;
         dLdv[2] = 1.;
         dLdw[0] = 0.;
         dLdw[1] = 0.;
         dLdw[2] = 0.;
         break;
      }
      case mEntity::TET:
      {
         switch (face)
         {
            case 0:
            {  // 012
               dLdu[0] = -1.;
               dLdu[1] = 1.;
               dLdu[2] = 0.;
               dLdv[0] = -1.;
               dLdv[1] = 0.;
               dLdv[2] = 1.;
               dLdw[0] = -1.;
               dLdw[1] = 0.;
               dLdw[2] = 0.;
               break;
            }
            case 1:
            {  // 013
               dLdu[0] = -1.;
               dLdu[1] = 1.;
               dLdu[2] = 0.;
               dLdv[0] = -1.;
               dLdv[1] = 0.;
               dLdv[2] = 0.;
               dLdw[0] = -1.;
               dLdw[1] = 0.;
               dLdw[2] = 1.;
               break;
            }
            case 2:
            {  // 123
               dLdu[0] = 1.;
               dLdu[1] = 0.;
               dLdu[2] = 0.;
               dLdv[0] = 0.;
               dLdv[1] = 1.;
               dLdv[2] = 0.;
               dLdw[0] = 0.;
               dLdw[1] = 0.;
               dLdw[2] = 1.;
               break;
            }
            case 3:
            {  // 023
               dLdu[0] = -1.;
               dLdu[1] = 0.;
               dLdu[2] = 0.;
               dLdv[0] = -1.;
               dLdv[1] = 1.;
               dLdv[2] = 0.;
               dLdw[0] = -1.;
               dLdw[1] = 0.;
               dLdw[2] = 1.;
               break;
            }
            default:
               throw;
         }

         break;
      }
      default:
      {
         throw;
      }
   }
}

xSimplexTet::xSimplexTet() = default;
void xSimplexTet::getL(const mEntity::mType &type, const xtensor::xPoint &uvw, double *L) const
{
   const double &u = uvw(0), v = uvw(1), w = uvw(2);
   L[0] = 1. - u - v - w;
   L[1] = u;
   L[2] = v;
   L[3] = w;
   return;
}

void xSimplexTet::getdL(const mEntity::mType &type, double *dLdu, double *dLdv, double *dLdw) const
{
   dLdu[0] = -1.;
   dLdu[1] = 1.;
   dLdu[2] = 0.;
   dLdu[3] = 0.;
   dLdv[0] = -1.;
   dLdv[1] = 0.;
   dLdv[2] = 1.;
   dLdv[3] = 0.;
   dLdw[0] = -1.;
   dLdw[1] = 0.;
   dLdw[2] = 0.;
   dLdw[3] = 1.;
   return;
}

void xSimplexTet::getLdL(const mEntity::mType &type, const xtensor::xPoint &uvw, double *L, double *dLdu, double *dLdv,
                         double *dLdw) const
{
   getL(type, uvw, L);
   getdL(type, dLdu, dLdv, dLdw);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// xApproxFunction Bernstein Implementation                                                  //
///////////////////////////////////////////////////////////////////////////////////////////////

xApproxFunctionScalarVertexBernstein::xApproxFunctionScalarVertexBernstein(int pL0, int _node) : xSimplexVertex(_node)
{
   powL = pL0;
}

void xApproxFunctionScalarVertexBernstein::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   double L;
   getL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L);
   res = pow(L, powL);
   storage.store(geo_appro, geo_integ, res);
}

void xApproxFunctionScalarVertexBernstein::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                   xtensor::xVector<> &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   storage.store(geo_appro, geo_integ, res);
   return;
}

void xApproxFunctionScalarVertexBernstein::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                        xtensor::xVector<> &res) const
{
   double L, dLdu, dLdv, dLdw;
   getLdL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L, dLdu, dLdv, dLdw);
   const double dfdL = powL * std::pow(L, std::max(powL - 1, 0));
   res = xtensor::xVector<>(dfdL * dLdu, dfdL * dLdv, dfdL * dLdw);
}

xApproxFunctionScalarEdgeBernstein::xApproxFunctionScalarEdgeBernstein(int pL0, int pL1, int _edge) : xSimplexEdge(_edge)
{
   powL[0] = pL0;
   powL[1] = pL1;
   int n = pL0 + pL1;
   factor = factorial(n) / factorial(pL0) / factorial(pL1);
   // factor =1;
}

void xApproxFunctionScalarEdgeBernstein::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   const xtensor::xPoint uvw = geo_appro->getUVW();
   const mEntity *e = geo_appro->getEntity();
   getVal_(e, uvw, res);
   storage.store(geo_appro, geo_integ, res);
}

void xApproxFunctionScalarEdgeBernstein::getVal_(const mEntity *e, const xtensor::xPoint &uvw, double &res) const
{
   double L[2];
   getL(e->getType(), uvw, L);
   res = 1.;
   for (int i = 0; i < 2; ++i)
   {
      res *= std::pow(L[i], powL[i]);
   }
   res *= factor;
}

void xApproxFunctionScalarEdgeBernstein::getVal_(const double &u, double &res) const
{
   res = 1.;
   double L[2] = {(1 - u) * 0.5, (1 + u) * 0.5};
   for (int i = 0; i < 2; ++i)
   {
      res *= std::pow(L[i], powL[i]);
   }
   res *= factor;
}

void xApproxFunctionScalarEdgeBernstein::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                      xtensor::xVector<> &res) const
{
   double L[2], dLdu[2], dLdv[2], dLdw[2];
   getLdL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L, dLdu, dLdv, dLdw);
   const double tmp = std::pow(L[0], powL[0] - 1) * std::pow(L[1], powL[1] - 1);
   double dfdL_0 = powL[0] * tmp * L[1];
   double dfdL_1 = powL[1] * tmp * L[0];

   res[0] = (dfdL_0 * dLdu[0] + dfdL_1 * dLdu[1]) * factor;  // dfdu
   res[1] = (dfdL_0 * dLdv[0] + dfdL_1 * dLdv[1]) * factor;  // dfdv
   res[2] = (dfdL_0 * dLdw[0] + dfdL_1 * dLdw[1]) * factor;  // dfdw
}

void xApproxFunctionScalarEdgeBernstein::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                 xtensor::xVector<> &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   storage.store(geo_appro, geo_integ, res);
}

double xApproxFunctionScalarEdgeBernstein::getIntegrale() const
{
   const int a = powL[0], b = powL[1];
   return 1. * factorial(a) * factorial(b) / (1. * factorial(a + b + 2));
}

xApproxFunctionScalarTriBernstein::xApproxFunctionScalarTriBernstein(int pL0, int pL1, int pL2, int _face) : xSimplexTri(_face)
{
   powL[0] = pL0;
   powL[1] = pL1;
   powL[2] = pL2;
   int n = pL0 + pL1 + pL2;
   factor = factorial(n) / factorial(pL0) / factorial(pL1) / factorial(pL2);
   // factor = 1.;
}

void xApproxFunctionScalarTriBernstein::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   double L[3];
   getL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L);
   res = std::pow(L[0], powL[0]) * std::pow(L[1], powL[1]) * std::pow(L[2], powL[2]);
   res *= factor;
   storage.store(geo_appro, geo_integ, res);
}

void xApproxFunctionScalarTriBernstein::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                xtensor::xVector<> &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   storage.store(geo_appro, geo_integ, res);
   return;
}

void xApproxFunctionScalarTriBernstein::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                     xtensor::xVector<> &res) const
{
   const xtensor::xPoint &uvw = geo_appro->getUVW();
   double L[3], dLdu[3], dLdv[3], dLdw[3];
   getLdL(geo_appro->getEntity()->getType(), uvw, L, dLdu, dLdv, dLdw);
   const double tmp = std::pow(L[0], powL[0] - 1) * std::pow(L[1], powL[1] - 1) * std::pow(L[2], powL[2] - 1);
   //    double dfdL[] = { powL[0]*tmp*L[1]*L[2], powL[1]*tmp*L[0]*L[2], dfdL[2] = powL[2]*tmp*L[0]*L[1]};
   double dfdL[3] = {powL[0] * tmp * L[1] * L[2], powL[1] * tmp * L[0] * L[2], powL[2] * tmp * L[0] * L[1]};
   res[0] = dfdL[0] * dLdu[0] + dfdL[1] * dLdu[1] + dfdL[2] * dLdu[2];
   res[1] = dfdL[0] * dLdv[0] + dfdL[1] * dLdv[1] + dfdL[2] * dLdv[2];
   res[2] = dfdL[0] * dLdw[0] + dfdL[1] * dLdw[1] + dfdL[2] * dLdw[2];
   res *= factor;
   return;
}

xApproxFunctionScalarTetBernstein::xApproxFunctionScalarTetBernstein(int pL0, int pL1, int pL2, int pL3) : xSimplexTet()
{
   powL[0] = pL0;
   powL[1] = pL1;
   powL[2] = pL2;
   powL[3] = pL3;
   int n = pL0 + pL1 + pL2 + pL3;
   factor = factorial(n) / factorial(pL0) / factorial(pL1) / factorial(pL2) / factorial(pL3);
}

void xApproxFunctionScalarTetBernstein::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   double L[4];
   getL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L);
   res = std::pow(L[0], powL[0]) * std::pow(L[1], powL[1]) * std::pow(L[2], powL[2]) * std::pow(L[3], powL[3]);
   res *= factor;
   storage.store(geo_appro, geo_integ, res);
}

void xApproxFunctionScalarTetBernstein::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                xtensor::xVector<> &res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   storage.store(geo_appro, geo_integ, res);
   return;
}

void xApproxFunctionScalarTetBernstein::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                     xtensor::xVector<> &res) const
{
   const xtensor::xPoint &uvw = geo_appro->getUVW();
   double L[4], dLdu[4], dLdv[4], dLdw[4];
   getLdL(geo_appro->getEntity()->getType(), uvw, L, dLdu, dLdv, dLdw);
   const double tmp =
       std::pow(L[0], powL[0] - 1) * std::pow(L[1], powL[1] - 1) * std::pow(L[2], powL[2] - 1) * std::pow(L[3], powL[3] - 1);
   double dfdL[] = {powL[0] * tmp * L[1] * L[2] * L[3], powL[1] * tmp * L[0] * L[2] * L[3], powL[2] * tmp * L[0] * L[1] * L[3],
                    powL[3] * tmp * L[0] * L[1] * L[2]};
   res[0] = dfdL[0] * dLdu[0] + dfdL[1] * dLdu[1] + dfdL[2] * dLdu[2] + dfdL[3] * dLdu[3];
   res[1] = dfdL[0] * dLdv[0] + dfdL[1] * dLdv[1] + dfdL[2] * dLdv[2] + dfdL[3] * dLdv[3];
   res[2] = dfdL[0] * dLdw[0] + dfdL[1] * dLdw[1] + dfdL[2] * dLdw[2] + dfdL[3] * dLdw[3];
   ;
   res *= factor;
   return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// xApproxFunction Lagrange Implementation                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////

std::map<int, std::vector<std::vector<int>>> xApproxFunctionScalarTetLagrange::noi;
std::map<int, std::vector<std::vector<int>>> xApproxFunctionScalarTriLagrange::noi;

xApproxFunctionScalarVertexLagrange::xApproxFunctionScalarVertexLagrange(int _order, int _node, double _coef)
    : xSimplexVertex(_node), order(_order), coef(_coef), coef_order(_coef * _order), n(order)
{
}

void xApproxFunctionScalarVertexLagrange::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   double L;
   getL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L);
   L *= order;
   res = LagrangeScaled(coef, &n, &L, 1);
   return;
}

void xApproxFunctionScalarVertexLagrange::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                       xtensor::xVector<> &res) const
{
   double L, dL, dLdu, dLdv, dLdw, fac;
   getLdL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L, dLdu, dLdv, dLdw);
   L *= order;
   dLagrangeScaled1(n, L, dL, fac);
   res[0] = dL * dLdu * coef_order;
   res[1] = dL * dLdv * coef_order;
   res[2] = dL * dLdw * coef_order;
   return;
}
void xApproxFunctionScalarVertexLagrange::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                  xtensor::xVector<> &res) const
{
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   return;
}

xApproxFunctionScalarEdgeLagrange::xApproxFunctionScalarEdgeLagrange(int _order, int _edge, int _ionedge, double _coef)
    : xSimplexEdge(_edge), order(_order), ionedge(_ionedge), coef(_coef), coef_order(_coef * _order)
{
   n[0] = order - 1 - ionedge;
   n[1] = 1 + ionedge;
}

void xApproxFunctionScalarEdgeLagrange::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   double L[2];
   getL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L);
   for (int i = 0; i < 2; ++i)
   {
      L[i] *= order;
   }
   res = LagrangeScaled(coef, n, L, 2);
   return;
}

void xApproxFunctionScalarEdgeLagrange::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                     xtensor::xVector<> &res) const
{
   double L[2], dL[2], dLdu[2], dLdv[2], dLdw[2];
   getLdL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L, dLdu, dLdv, dLdw);
   for (int i = 0; i < 2; ++i) L[i] *= order;
   dLagrangeScaled2(n, L, dL);
   res[0] = (dL[0] * dLdu[0] + dL[1] * dLdu[1]) * coef_order;
   res[1] = (dL[0] * dLdv[0] + dL[1] * dLdv[1]) * coef_order;
   res[2] = (dL[0] * dLdw[0] + dL[1] * dLdw[1]) * coef_order;
   return;
}
void xApproxFunctionScalarEdgeLagrange::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                xtensor::xVector<> &res) const
{
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   return;
}

void xApproxFunctionScalarEdgeLagrange::getNodePos(int order, int ionedge, int *nodepos)
{
   nodepos[0] = order - 1 - ionedge;
   nodepos[1] = 1 + ionedge;
}

void xApproxFunctionScalarEdgeLagrange::getNodeUvw(std::vector<double> &nodeuvw)
{
   nodeuvw.resize(3);
   nodeuvw[1] = 0.;
   nodeuvw[2] = 0.;
   double lambda = (1. / order) * this->n[0];

   // uvw = lambda * x1 + (1-lambda) * x2
   // ici x1=-1. et x2=1.
   nodeuvw[0] = 1. - 2 * lambda;
}

void xApproxFunctionScalarTriLagrange::getNodePos(int order, int ionface, int *nodepos)
{
   if (noi[order].size() == 0)
   {
      noi[order].resize(nbShapeFunctionTri(order));
      int index = 0;
      for (int i = 1; i <= order - 2; ++i)
      {
         for (int j = 1; j <= order - 2 - i + 1; ++j)
         {
            const int k = order - i - j;
            noi[order][index] = {i, j, k};
            ++index;
         }
      }
   }
   std::copy(noi[order][ionface].begin(), noi[order][ionface].end(), nodepos);
   return;
}

/*
  if (noi[order].size() == 0){
        noi[order].resize(  nbShapeFunctionTri(order) );
    }
    if (noi[order][ionface].size() == 0){
        noi[order][ionface].resize(3);
        std::vector <int > subfaces(order/3+1);
        int maxlevel = order/3;
        subfaces[0] = 0;
        for (int i = 1; i <= maxlevel; ++i){
            subfaces[i] = ((order-3*(i-1)) == 3 ) ? 1: (order-3*(i))*3;
        }
        for (int i = 1; i <= maxlevel; ++i){
            subfaces[i] += subfaces[i-1];
        }
        int level = 1;
        if (ionface > subfaces[order/3]) throw;
        while (ionface>= subfaces[level]) level++;
        int ionsubface = ionface-subfaces[level-1];
        if (ionsubface < 3){
            nodepos[ionsubface%3] = level + order-3*(level);
            nodepos[(ionsubface+1)%3] = level;
            nodepos[(ionsubface+2)%3] = level;
            for (int i =0; i <3; ++i) noi[order][ionface][i] = nodepos[i]  ;
            return;
        }
        int subedge    = (ionsubface -3)/(order-1-3*level);
        int ionsubedge = (ionsubface -3)%(order-1-3*level);
        int   nedge[2];
        xApproxFunctionScalarEdgeLagrange::getNodePos( order-3*level , ionsubedge, nedge);
        nodepos[  subedge %3]     = level + nedge[0];
        nodepos[ (subedge+1) %3]  = level + nedge[1];
        nodepos[ (subedge+2) %3 ] = level;
        for (int i =0; i <3; ++i) noi[order][ionface][i] = nodepos[i]  ;
        return;
    }
    else{
        for (int i =0; i <3; ++i) nodepos[i] = noi[order][ionface][i];
        for (int i =0; i <3; ++i){
          std::cout << nodepos[i]<< " ";
        }
        std::cout << std::endl;

        return;
    }
}
*/

void xApproxFunctionScalarTriLagrange::getNodeUvw(std::vector<double> &nodeuvw)
{
   // nodeuvw.resize(3);

   // uvw = lambda1 * x1 + lambda2 * x2 + (1-lambda1 - lambda2) * x3
   // ici x1=(0.,0.), x2=(1.,0.) et x3=(0.,1.)

   // const double lambda1 = (1./order)*this->n[0];
   const double lambda2 = (1. / order) * this->n[1];
   const double lambda3 = (1. / order) * this->n[2];
   nodeuvw = {lambda2, lambda3, 0.};
   // nodeuvw[0] = lambda2;
   // nodeuvw[1] = lambda3;
   // nodeuvw[2]=0.;
}

xApproxFunctionScalarTriLagrange::xApproxFunctionScalarTriLagrange(int _order, int _face, int _ionface, double _coef)
    : xSimplexTri(_face), order(_order), ionface(_ionface), coef(_coef), coef_order(_coef * _order)
{
   getNodePos(order, ionface, n);
}

void xApproxFunctionScalarTriLagrange::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   double L[3];
   getL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L);
   for (int i = 0; i < 3; ++i)
   {
      L[i] *= order;
   }
   res = LagrangeScaled(coef, n, L, 3);
   return;
}

void xApproxFunctionScalarTriLagrange::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                    xtensor::xVector<> &res) const
{
   double L[3], dL[3], dLdu[3], dLdv[3], dLdw[3];
   getLdL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L, dLdu, dLdv, dLdw);
   for (int i = 0; i < 3; ++i) L[i] *= order;
   dLagrangeScaled3(n, L, dL);
   res[0] = (dL[0] * dLdu[0] + dL[1] * dLdu[1] + dL[2] * dLdu[2]) * coef_order;
   res[1] = (dL[0] * dLdv[0] + dL[1] * dLdv[1] + dL[2] * dLdv[2]) * coef_order;
   res[2] = (dL[0] * dLdw[0] + dL[1] * dLdw[1] + dL[2] * dLdw[2]) * coef_order;
}
void xApproxFunctionScalarTriLagrange::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                               xtensor::xVector<> &res) const
{
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   return;
}

void xApproxFunctionScalarTetLagrange::getNodePos(int order, int iontet, int *nodepos)
{
   if (noi[order].size() == 0)
   {
      noi[order].resize(nbShapeFunctionTet(order));
   }
   if (noi[order][iontet].size() == 0)
   {
      noi[order][iontet].resize(4);
      int maxlevel = order / 4;
      std::vector<int> nnodeonleveli(maxlevel);
      for (int i = 0; i < maxlevel; ++i)
      {
         nnodeonleveli[i] = nbShapeFunctionTet(order - i * 4) - nbShapeFunctionTet(order - (i + 1) * 4);
      }
      int level = 0;
      int ionlevel = iontet;
      while (ionlevel >= nnodeonleveli[level])
      {
         ionlevel -= nnodeonleveli[level];
         ++level;
      }
      level += 1;
      nodepos[0] = nodepos[1] = nodepos[2] = nodepos[3] = level;

      if (ionlevel < 4)
      {
         nodepos[ionlevel] += order - 4 * level;
         for (int i = 0; i < 4; ++i) noi[order][iontet][i] = nodepos[i];
         return;
      }
      ionlevel -= 4;
      const int nnodeedge = nbShapeFunctionEdge(order - (level)*4);

      if (ionlevel < (nnodeedge * 6))
      {
         int subedge = ionlevel / (nnodeedge);
         int ionsubedge = ionlevel % (nnodeedge);
         int n0 = Tev[subedge][0];
         int n1 = Tev[subedge][1];
         int nedge[2];
         xApproxFunctionScalarEdgeLagrange::getNodePos(order - 4 * level, ionsubedge, nedge);
         nodepos[n0] += nedge[0];
         nodepos[n1] += nedge[1];
         for (int i = 0; i < 4; ++i) noi[order][iontet][i] = nodepos[i];
         return;
      }
      ionlevel -= nnodeedge * 6;
      int nnodeface = nbShapeFunctionTri(order - (level)*4);
      if (ionlevel < nnodeface * 4)
      {
         int subface = ionlevel / (nnodeface);
         int ionsubface = ionlevel % (nnodeface);
         int n0 = Tfv[subface][0];
         int n1 = Tfv[subface][1];
         int n2 = Tfv[subface][2];
         int nface[3];
         xApproxFunctionScalarTriLagrange::getNodePos(order - 4 * level, ionsubface, nface);
         nodepos[n0] += nface[0];
         nodepos[n1] += nface[1];
         nodepos[n2] += nface[2];
         for (int i = 0; i < 4; ++i) noi[order][iontet][i] = nodepos[i];
         return;
      }
      else
         throw;
   }
   else
   {
      std::copy(noi[order][iontet].begin(), noi[order][iontet].end(), nodepos);
   }
}

xApproxFunctionScalarTetLagrange::xApproxFunctionScalarTetLagrange(int _order, int iontet, double _coef)
    : order(_order), coef(_coef), coef_order(_coef * _order)
{
   getNodePos(order, iontet, n);
}

void xApproxFunctionScalarTetLagrange::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   double L[4];
   getL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L);
   for (int i = 0; i < 4; ++i)
   {
      L[i] *= order;
   }
   res = LagrangeScaled(coef, n, L, 4);
   return;
}

void xApproxFunctionScalarTetLagrange::getGradLocal(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                                    xtensor::xVector<> &res) const
{
   double L[4], dL[4], dLdu[4], dLdv[4], dLdw[4];
   getLdL(geo_appro->getEntity()->getType(), geo_appro->getUVW(), L, dLdu, dLdv, dLdw);
   for (int i = 0; i < 4; ++i) L[i] *= order;
   dLagrangeScaled4(n, L, dL);
   res[0] = (dL[0] * dLdu[0] + dL[1] * dLdu[1] + dL[2] * dLdu[2] + dL[3] * dLdu[3]) * coef_order;
   res[1] = (dL[0] * dLdv[0] + dL[1] * dLdv[1] + dL[2] * dLdv[2] + dL[3] * dLdv[3]) * coef_order;
   res[2] = (dL[0] * dLdw[0] + dL[1] * dLdw[1] + dL[2] * dLdw[2] + dL[3] * dLdw[3]) * coef_order;
   return;
}

void xApproxFunctionScalarTetLagrange::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ,
                                               xtensor::xVector<> &res) const
{
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   return;
}

double LagrangeScaled(const double coef, const int *n, const double *L, const int dim)
{
   double res = coef;
   for (int i = 0; i < dim; ++i)
   {
      const double Li = L[i];
      const int kmax = n[i];
      for (int k = 0; k < kmax; ++k)
      {
         res *= (Li - k) / (kmax - k);
      }
   }
   return res;
}

void dLagrangeScaled1(const int kmax, double &L, double &dL, double &fac)
{
   int k, l;
   const double one = 1.0;
   switch (kmax)
   {
      case 1:
      {
         dL = one;
         fac = L;
         break;
      }
      case 2:
      {
         const double Li = L;
         const double oneontwo = 0.5;
         dL = Li - oneontwo;
         fac = oneontwo * Li * (Li - one);
         break;
      }
      case 3:
      {
         const double Li = L;
         const double oneonsix = 1. / 6.;
         const double v1 = Li * oneonsix;
         const double v2 = Li - one;
         const double v3 = Li - 2.;
         const double v4 = v2 * v3 * oneonsix;
         dL = v1 * (v2 + v3) + v4;
         fac = v4 * Li;
         break;
      }
      case 4:
      {
         const double Li = L;
         const double oneontwentyfor = 1. / 24.;
         const double v1 = Li * oneontwentyfor;
         const double v2 = Li - 3.;
         const double v3 = Li - one;
         const double v4 = Li - 2.;
         const double v5 = v3 * v4;
         const double v6 = v1 * v2;
         dL = (v1 + v2 * oneontwentyfor) * v5 + (v3 + v4) * v6;
         fac = v5 * v6;
         break;
      }
      default:
      {
         double tmp[XAPPROXFUNCTIONLAGRANGEMAXORDER];
         const double Li = L;
         double facl = 1.;
         double tmpl;
         for (k = 0; k < kmax; ++k) facl *= (Li - k) / (kmax - k);
         tmp[0] = one / kmax;
         tmpl = Li * tmp[0];
         for (l = 1; l < kmax; ++l) tmp[l] = tmpl;
         for (k = 1; k < kmax; ++k)
         {
            tmpl = (Li - k) / (kmax - k);
            for (l = 0; l < kmax; ++l)
            {
               if (l != k)
                  tmp[l] *= tmpl;
               else
                  tmp[l] /= (kmax - k);
            }
         }
         tmpl = tmp[0];
         for (k = 1; k < kmax; ++k)
         {
            tmpl += tmp[k];
         }
         dL = tmpl;
         fac = facl;
      }
   }
}
void dLagrangeScaled2(const int *n, double *L, double *dL)
{
   double fac[2];
   dLagrangeScaled1(n[0], L[0], dL[0], fac[0]);
   dLagrangeScaled1(n[1], L[1], dL[1], fac[1]);
   dL[0] *= fac[1];
   dL[1] *= fac[0];
}

void dLagrangeScaled3(const int *n, double *L, double *dL)
{
   double fac[3];
   dLagrangeScaled1(n[0], L[0], dL[0], fac[0]);
   dLagrangeScaled1(n[1], L[1], dL[1], fac[1]);
   dLagrangeScaled1(n[2], L[2], dL[2], fac[2]);
   dL[0] *= fac[1] * fac[2];
   dL[1] *= fac[2] * fac[0];
   dL[2] *= fac[0] * fac[1];
}

void dLagrangeScaled4(const int *n, double *L, double *dL)
{
   double fac[4];
   dLagrangeScaled1(n[0], L[0], dL[0], fac[0]);
   dLagrangeScaled1(n[1], L[1], dL[1], fac[1]);
   dLagrangeScaled1(n[2], L[2], dL[2], fac[2]);
   dLagrangeScaled1(n[3], L[3], dL[3], fac[3]);
   double v1 = fac[2] * fac[3];
   dL[0] *= fac[1] * v1;
   dL[1] *= v1 * fac[0];
   v1 = fac[0] * fac[1];
   dL[2] *= fac[3] * v1;
   dL[3] *= v1 * fac[2];
}

/// regular spaced lagrange. u  between 0 and 1.
/*!
    n: order of the lagrange basis,
    i: node on which lag =1. all other node, lag = 0.
    0---1---2---  ---i--- ---n-1---n
    u:  0/n-1/n-2/n-    -i/n-  -i/(n-1)-n/n
  */
double lag(int n, int i, double u)
{
   double res = 1.;
   if (n == 0) return res;
   double ui = (1. * i) / n;
   for (int k = 0; k < (n + 1); ++k)
   {
      if (k != i)
      {
         double rk = (1. * k) / n;
         res *= (u - rk) / (ui - rk);
      }
   }
   return res;
}

}  // end namespace xfem
