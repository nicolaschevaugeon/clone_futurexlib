/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

// std
#include <iostream>

// xtensor
#include "xTensor2.h"
#include "xTensorOperations.h"
#include "xVector.h"
// xmapping
#include "xMapping.h"
// xfem
#include "xApproxFunction.h"
#include "xDebug.h"
#include "xElement.h"
#include "xGeomElem.h"
#include "xLevelSet.h"
#include "xMesh.h"

namespace xfem
{
using AOMD::mEntity;
using std::cerr;
using std::cout;
using std::endl;

xApproxFunction::xApproxFunction() = default;
void xApproxFunction::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const
{
   cerr << "getVal not coded yet\n";
   assert(0);
   throw;
}
void xApproxFunction::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const
{
   cerr << "getGrad not coded yet\n";
   assert(0);
   throw;
}
void xApproxFunction::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const
{
   cerr << "getGradLocal not coded yet\n";
   assert(0);
}
void xApproxFunction::getHessian(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const
{
   cerr << "getGrad not coded yet\n";
   assert(0);
   throw;
}

void xApproxFunction::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const
{
   cerr << "getVal not coded yet\n";
   assert(0);
   throw;
}
void xApproxFunction::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const
{
   cerr << "getGrad not coded yet\n";
   assert(0);
   throw;
}

void xApproxFunction::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&, const xValKey&) const
{
   cerr << "getVal not coded yet\n";
   assert(0);
   throw;
}
void xApproxFunction::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&, const xValKey&) const
{
   cerr << "getGrad not coded yet\n";
   assert(0);
   throw;
}

void xApproxFunction::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const
{
   cerr << "getGradLocal not coded yet\n";
   assert(0);
}
void xApproxFunction::getGradLocalAxpy(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>& y,
                                       const double& a) const
{
   cerr << "getGradLocalAxpy not coded yet\n";
   assert(0);
}
void xApproxFunction::getHessian(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor4<>&) const
{
   cerr << "getGrad not coded yet\n";
   assert(0);
}

void xApproxFunctionConstant::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const { res = 1.; }

void xApproxFunctionConstant::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   res(0) = 0.;
   res(1) = 0.;
   res(2) = 0.;
}

xApproxFunctionVector::xApproxFunctionVector(shapeFctPtr sf, unsigned int comp) : scalarFunction(sf), component(comp) {}

void xApproxFunctionVector::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   res *= 0.;
   scalarFunction->getVal(geo_appro, geo_integ, res(component));
}

void xApproxFunctionVector::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>& res) const
{
   xtensor::xVector<> res1;
   scalarFunction->getGrad(geo_appro, geo_integ, res1);

   // Les xtensor::xTensor2<> sont deja initialises a zero: inutile de recommencer.
   // Matrice 3x3: pas la peine de faire deux petites boucles imbriquees : gros overhead...
   // On deroule la boucle
   res *= 0.;
   res(component, 0) = res1(0);
   res(component, 1) = res1(1);
   res(component, 2) = res1(2);
}

xApproxFunctionSummed::xApproxFunctionSummed(std::vector<shapeFctPtr>& approx) : approx_functions(approx) {}
void xApproxFunctionSummed::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   std::vector<shapeFctPtr>::const_iterator ita = approx_functions.begin();
   std::vector<shapeFctPtr>::const_iterator itae = approx_functions.end();

   double tmp;
   (*ita)->getVal(geo_appro, geo_integ, res);

   for (++ita; ita != itae; ++ita)
   {
      (*ita)->getVal(geo_appro, geo_integ, tmp);
      res += tmp;
   }

   return;
}
void xApproxFunctionSummed::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   std::vector<shapeFctPtr>::const_iterator ita = approx_functions.begin();
   std::vector<shapeFctPtr>::const_iterator itae = approx_functions.end();

   (*ita)->getGradLocal(geo_appro, geo_integ, res);
   xtensor::xVector<> tmp_vect;

   for (++ita; ita != itae; ++ita)
   {
      (*ita)->getGradLocal(geo_appro, geo_integ, tmp_vect);
      res += tmp_vect;
   }

   return;
}
void xApproxFunctionSummed::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   xtensor::xVector<> gl;
   getGradLocal(geo_appro, geo_integ, gl);
   res = (geo_appro)->PushBack(gl);
   return;
}

void xApproxFunctionEnrichedXFEM::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   const bool debug = xdebug_flag;
   double b, e;
   base->getVal(geo_appro, geo_integ, b);
   enrichment->getVal(geo_appro, geo_integ, e);
   res = b * e;

   if (debug)
   {
      cout << " enrichment at point " << geo_appro->getXYZ() << " is " << e;
      cout << " ScalarxApproxFunctionEnrichedXFEM::getVal " << res << endl;
   }
}

void xApproxFunctionEnrichedXFEM::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   const bool debug = xdebug_flag;
   xtensor::xVector<> bg, eg;
   double es, bs;

   base->getGrad(geo_appro, geo_integ, bg);
   enrichment->getVal(geo_appro, geo_integ, es);
   enrichment->getGrad(geo_appro, geo_integ, eg);
   base->getVal(geo_appro, geo_integ, bs);
   res = bg * es + eg * bs;

   if (debug)
   {
      cout << " enrichment  at point " << geo_appro->getXYZ() << " is " << es << endl;
      cout << " grad enrichment " << eg << endl;
      cout << " ScalarxApproxFunctionEnrichedXFEM::getGrad " << res << endl;
   }

   return;
}

void xApproxFunctionEnrichedXFEM::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   xtensor::xVector<> bg;
   double es;
   base->getVal(geo_appro, geo_integ, bg);
   enrichment->getVal(geo_appro, geo_integ, es);
   res = bg * es;
   return;
}

void xApproxFunctionEnrichedXFEM::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>& res) const
{
   const bool debug = xdebug_flag;
   xtensor::xTensor2<> bt;
   xtensor::xVector<> eg, bv;
   double es;
   base->getGrad(geo_appro, geo_integ, bt);
   enrichment->getVal(geo_appro, geo_integ, es);
   enrichment->getGrad(geo_appro, geo_integ, eg);
   base->getVal(geo_appro, geo_integ, bv);
   res = bt * es + tensor_product(bv, eg);
   if (debug)
   {
      cout << " enrichment  at point " << geo_appro->getXYZ() << " is " << es << endl;
      cout << " grad enrichment " << eg << endl;
      cout << "val classical " << endl;
      cout << bv << endl;
      cout << "grad classical " << endl;
      cout << bt << endl;

      cout << "first part " << endl;
      cout << bt * es << endl;
      cout << "second part " << endl;
      cout << xtensor::tensor_product(eg, bv) << endl;
   }
}

void xApproxFunctionEnrichedVectorXFEM::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                               xtensor::xVector<>& res) const
{
   xtensor::xVector<> ev;
   double bs;
   base->getVal(geo_appro, geo_integ, bs);
   enrichment->getVal(geo_appro, geo_integ, ev);
   res = ev * bs;
}

void xApproxFunctionEnrichedVectorXFEM::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                xtensor::xTensor2<>& res) const
{
   xtensor::xTensor2<> egt;
   xtensor::xVector<> ev, bgv;
   double bs;
   base->getVal(geo_appro, geo_integ, bs);
   base->getGrad(geo_appro, geo_integ, bgv);

   enrichment->getVal(geo_appro, geo_integ, ev);
   enrichment->getGrad(geo_appro, geo_integ, egt);

   res = egt * bs + tensor_product(ev, bgv);
   return;
}

void xApproxFunctionEnrichedVector2XFEM::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                xtensor::xVector<>& res) const
{
   xtensor::xVector<> ev;
   double bs;
   base->getVal(geo_appro, geo_integ, bs);
   enrichment->getVal(geo_appro, geo_integ, ev, basekey);
   res = ev * bs;
}

void xApproxFunctionEnrichedVector2XFEM::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                 xtensor::xTensor2<>& res) const
{
   xtensor::xTensor2<> egt;
   xtensor::xVector<> ev, bgv;
   double bs;
   base->getVal(geo_appro, geo_integ, bs);
   base->getGrad(geo_appro, geo_integ, bgv);

   enrichment->getVal(geo_appro, geo_integ, ev, basekey);
   enrichment->getGrad(geo_appro, geo_integ, egt, basekey);
   res = egt * bs + tensor_product(ev, bgv);
   return;
}

void xApproxFunctionScalarHierarchical::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   const bool debug = xdebug_flag;
   res = xfem::HierarchicalApproxFunction(ith, geo_appro);
   if (debug) cout << " res of xApproxFunctionScalarHierarchical " << res << endl;
}

void xApproxFunctionScalarHierarchical::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                xtensor::xVector<>& res) const
{
   xtensor::xVector<> in = xfem::GradHierarchicalApproxFunction(ith, geo_appro);
   res = (geo_appro)->PushBack(in);
}

void xApproxFunctionScalarPieceWiseHierarchical::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   // care must be taken in evaluating the function because it is discontinuous
   const bool debug = xdebug_flag;
   auto cdg = geo_integ->getCDGxyz();
   const xmapping::xMapping* mapping = geo_appro->getMapping();
   double u, v, w;
   mapping->invert(cdg(0), cdg(1), cdg(2), u, v, w);

   if (ith == 1 && u <= 0.0)
      res = 1.;
   else if (ith == 2 && u >= 0.0)
      res = 1.;
   else
      res = 0.;

   if (debug)
   {
      cout << " us is " << u << endl;
      cout << " res of xApproxFunctionScalarPieceWiseHierarchical " << res << endl;
   }
}
void xApproxFunctionScalarPieceWiseHierarchical::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                         xtensor::xVector<>& res) const
{
   res(0) = res(1) = res(2) = 0.0;
}

void xApproxFunctionVectorHierarchical::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                               xtensor::xVector<>& res) const
{
   res(0) = 0.;
   res(1) = 0.;
   res(2) = 0.;
   res(comp) = xfem::HierarchicalApproxFunction(ith, geo_appro);
}

void xApproxFunctionVectorHierarchical::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                xtensor::xTensor2<>& res) const
{
   xtensor::xVector<> s = xfem::GradHierarchicalApproxFunction(ith, geo_appro);
   geo_appro->PushBack(s);
   for (int j = 0; j < 3; j++)
   {
      for (int i = 0; i < 3; i++)
      {
         if (i == comp)
            res(i, j) = s(j);
         else
            res(i, j) = 0.;
      }
   }
}

void xApproxFunctionVectorHierarchical::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                     xtensor::xTensor2<>& res) const
{
   xtensor::xVector<> s = xfem::GradHierarchicalApproxFunction(ith, geo_appro);
   for (int j = 0; j < 3; j++)
   {
      for (int i = 0; i < 3; i++)
      {
         if (i == comp)
            res(i, j) = s(j);
         else
            res(i, j) = 0.;
      }
   }
}

void xApproxFunctionVectorHierarchical::getGradLocalAxpy(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                         xtensor::xTensor2<>& y, const double& a) const
{
   xtensor::xVector<> x = xfem::GradHierarchicalApproxFunction(ith, geo_appro);
   for (int j = 0; j < 3; j++) y(comp, j) += x(j) * a;
}

double HierarchicalApproxFunction(int NumNode, const xGeomElem* geo)
{
   auto uvw = geo->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);
   mEntity* e = geo->getEntity();
   switch (e->getType())
   {
      case mEntity::VERTEX:
         return 1.;
      case mEntity::EDGE:
         return HierarchicalApproxFunctionEdge(NumNode, u);
         break;

      case mEntity::TRI:
      {
         double temp = HierarchicalApproxFunctionTriangle(NumNode, u, v);
         switch (NumNode)
         {
            case 7:
            {
               mEntity* edge = e->get(1, 0);
               mEntity* node1 = e->get(0, 0);
               if (node1 == edge->get(0, 0))
                  return temp;
               else
                  return -temp;
            }
            case 8:
            {
               mEntity* edge = e->get(1, 1);
               mEntity* node1 = e->get(0, 1);
               if (node1 == edge->get(0, 0))
                  return temp;
               else
                  return -temp;
            }
            case 9:
            {
               mEntity* edge = e->get(1, 2);
               mEntity* node1 = e->get(0, 2);
               if (node1 == edge->get(0, 0))
                  return temp;
               else
                  return -temp;
            }
            default:
               return temp;
         }
      }
      case mEntity::QUAD:
         return HierarchicalApproxFunctionQuadrilateral(NumNode, u, v);
      case mEntity::TET:
         return HierarchicalApproxFunctionTetrahedron(NumNode, u, v, w);
      case mEntity::HEX:
         return HierarchicalApproxFunctionHexahedron(NumNode, u, v, w);
      case mEntity::PRISM:
         return HierarchicalApproxFunctionPrism(NumNode, u, v, w);
      default:
         throw;
         return 0.;
   }
   throw;
   return 0.;
}

xtensor::xVector<> GradHierarchicalApproxFunction(int NumNode, const xGeomElem* geo)
{
   double s[3];
   auto uvw = geo->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);
   mEntity* e = geo->getEntity();
   switch (e->getType())
   {
      case mEntity::EDGE:
         GradHierarchicalApproxFunctionEdge(NumNode, u, s);
         break;
      case mEntity::TRI:
      {
         GradHierarchicalApproxFunctionTriangle(NumNode, u, v, s);
         switch (NumNode)
         {
            case 7:
            {
               mEntity* edge = e->get(1, 0);
               mEntity* node1 = e->get(0, 0);
               if (node1 != edge->get(0, 0))
                  for (int i = 0; i < 3; i++) s[i] = -s[i];
               break;
            }
            case 8:
            {
               mEntity* edge = e->get(1, 1);
               mEntity* node1 = e->get(0, 1);
               if (node1 != edge->get(0, 0))
                  for (int i = 0; i < 3; i++) s[i] = -s[i];
               ;
               break;
            }
            case 9:
            {
               mEntity* edge = e->get(1, 2);
               mEntity* node1 = e->get(0, 2);
               if (node1 != edge->get(0, 0))
                  for (int i = 0; i < 3; i++) s[i] = -s[i];
               break;
            }
            default:
               break;
         }
         break;
      }
      case mEntity::QUAD:
         GradHierarchicalApproxFunctionQuadrilateral(NumNode, u, v, s);
         break;
      case mEntity::TET:
         GradHierarchicalApproxFunctionTetrahedron(NumNode, u, v, w, s);
         break;
      case mEntity::HEX:
         GradHierarchicalApproxFunctionHexahedron(NumNode, u, v, w, s);
         break;
      case mEntity::PRISM:
         GradHierarchicalApproxFunctionPrism(NumNode, u, v, w, s);
         break;
      default:
         assert(0);
         break;
   }
   return xtensor::xVector<>(s[0], s[1], s[2]);
}

double HierarchicalApproxFunctionEdge(int ith, double u)
{
   // attention il faut que les fcts ci-dessous soit les traces
   // des fonctions sur elt 2D
   // u in [-1 1]
   switch (ith)
   {
      case 1:
         return 0.5 * (1. - u);
      case 2:
         return 0.5 * (1. + u);
      case 3:
         return 0.25 * Psi(u, 2) * (1. - u * u);
      case 4:
         return 0.25 * Psi(u, 3) * (1. - u * u);
      default:
         std::cout
             << "Warning , only four shape functions are implemented for HierarchicalApproxFunctionEdge (order 3) returning 0."
             << std::endl;
         assert(0);
         return 0.;
   }
}

double HierarchicalApproxFunctionTriangle(int ith, double u, double v)
{
   switch (ith)
   {
      case 1:
         return 1. - u - v;
      case 2:
         return u;
      case 3:
         return v;
      case 4:
         return (1. - u - v) * u * Psi(u - (1. - u - v), 2);
      case 5:
         return u * v * Psi(v - u, 2);
      case 6:
         return v * (1. - u - v) * Psi((1. - u - v) - v, 2);
      case 7:
         return (1. - u - v) * u * Psi(u - (1. - u - v), 3);
      case 8:
         return u * v * Psi(v - u, 3);
      case 9:
         return v * (1. - u - v) * Psi((1. - u - v) - v, 3);
      case 10:
         return (1. - u - v) * u * v;
      default:
         assert(0);
         return 0.;
         break;
   }
}

double HierarchicalApproxFunctionQuadrilateral(int ith, double u, double v)
{
   switch (ith)
   {
      case 1:
         return 0.25 * (1. - u) * (1. - v);
      case 2:
         return 0.25 * (1. + u) * (1. - v);
      case 3:
         return 0.25 * (1. + u) * (1. + v);
      case 4:
         return 0.25 * (1. - u) * (1. + v);
      case 5:
         return 0.5 * (1. - v) * Phi(u, 2);
      case 6:
         return 0.5 * (1. + u) * Phi(v, 2);
      case 7:
         return 0.5 * (1. + v) * Phi(u, 2);
      case 8:
         return 0.5 * (1. - u) * Phi(v, 2);
      case 9:
         return Phi(u, 2) * Phi(v, 2);
      case 10:
         return 0.5 * (1. - v) * Phi(u, 3);
      case 11:
         return 0.5 * (1. + u) * Phi(v, 3);
      case 12:
         return 0.5 * (1. + v) * Phi(u, 3);
      case 13:
         return 0.5 * (1. - u) * Phi(v, 3);
      case 14:
         return Phi(u, 2) * Phi(v, 3);
      case 15:
         return Phi(u, 3) * Phi(v, 2);
      case 16:
         return Phi(u, 3) * Phi(v, 3);
      default:
         cerr << "ith " << ith << " not valid\n";
         assert(0);
         return 0.;
         break;
   }
}
double HierarchicalApproxFunctionTetrahedron(int ith, double u, double v, double w)
{
   switch (ith)
   {
      case 1:
         return 1. - u - v - w;
      case 2:
         return u;
      case 3:
         return v;
      case 4:
         return w;
      case 5:
         return (1. - u - v - w) * u;
      case 6:
         return u * v;
      case 7:
         return (1. - u - v - w) * v;
      case 8:
         return (1. - u - v - w) * w;
      case 9:
         return u * w;
      case 10:
         return v * w;
      default:
         assert(0);
         return 0.;
         break;
   }
}

double HierarchicalApproxFunctionHexahedron(int ith, double u, double v, double w)
{
   switch (ith)
   {
      case 1:
         return (1. - u) * (1. - v) * (1. - w) * 0.125;
      case 2:
         return (1. + u) * (1. - v) * (1. - w) * 0.125;
      case 3:
         return (1. + u) * (1. + v) * (1. - w) * 0.125;
      case 4:
         return (1. - u) * (1. + v) * (1. - w) * 0.125;
      case 5:
         return (1. - u) * (1. - v) * (1. + w) * 0.125;
      case 6:
         return (1. + u) * (1. - v) * (1. + w) * 0.125;
      case 7:
         return (1. + u) * (1. + v) * (1. + w) * 0.125;
      case 8:
         return (1. - u) * (1. + v) * (1. + w) * 0.125;
      default:
         assert(0);
         return 0.;
         break;
   }
}

double HierarchicalApproxFunctionPrism(int ith, double u, double v, double w)
{
   switch (ith)
   {
      case 1:
         return (1. - u - v) * (1. - w) * 0.5;
      case 2:
         return u * (1. - w) * 0.5;
      case 3:
         return v * (1. - w) * 0.5;
      case 4:
         return (1. - u - v) * (1. + w) * 0.5;
      case 5:
         return u * (1. + w) * 0.5;
      case 6:
         return v * (1. + w) * 0.5;
      default:
         return 0.;
   }
}

void GradHierarchicalApproxFunctionEdge(int ith, double u, double s[])
{
   switch (ith)
   {
      case 1:
         s[0] = -0.5;
         s[1] = 0.;
         s[2] = 0.;
         break;
      case 2:
         s[0] = 0.5;
         s[1] = 0.;
         s[2] = 0.;
         break;
      case 3:
         s[0] = 0.25 * Psi(u, 2) * (1. - 2 * u);
         s[1] = 0.;
         s[2] = 0.;
         break;
      case 4:
         s[0] = 0.25 * Psi(u, 3) * (1. - 2 * u) + 0.25 * (1. - u * u) * Dpsi(u, 3);
         s[1] = 0.;
         s[2] = 0.;
         break;
      default:
         assert(0);
         return;
         break;
   }
}

void GradHierarchicalApproxFunctionTriangle(int ith, double u, double v, double s[])
{
   switch (ith)
   {
         // Degree 1
      case 1:
         s[0] = -1.;
         s[1] = -1.;
         s[2] = 0.;
         break;
      case 2:
         s[0] = 1.;
         s[1] = 0.;
         s[2] = 0.;
         break;
      case 3:
         s[0] = 0.;
         s[1] = 1.;
         s[2] = 0.;
         break;
         // Degre 2
      case 4:
      {
         const double psi = Psi(u - (1. - u - v), 2);
         const double dpsi = Dpsi(u - (1. - u - v), 2);
         s[0] = (1. - u - v) * u * (2.) * dpsi + (1. - 2. * u - v) * psi;
         s[1] = (1. - u - v) * u * dpsi - u * psi;
         s[2] = 0.;
      }
      break;
      case 5:
      {
         const double psi = Psi(v - u, 2);
         const double dpsi = Dpsi(v - u, 2);
         s[0] = u * v * (-1) * dpsi + v * psi;
         s[1] = u * v * dpsi + u * psi;
         s[2] = 0.;
      }
      break;
      case 6:
      {
         const double psi = Psi((1. - u - v) - v, 2);
         const double dpsi = Dpsi((1. - u - v) - v, 2);
         s[0] = v * (1. - u - v) * (-1) * dpsi - v * psi;
         s[1] = v * (1. - u - v) * (-2) * dpsi + (1. - u - 2. * v) * psi;
         s[2] = 0.;
      }
      break;
         // Degre 3
      case 7:
      {
         const double psi = Psi(u - (1. - u - v), 3);
         const double dpsi = Dpsi(u - (1. - u - v), 3);
         s[0] = (1. - u - v) * u * (2) * dpsi + (1. - 2. * u - v) * psi;
         s[1] = (1. - u - v) * u * dpsi - u * psi;
         s[2] = 0.;
      }
      break;
      case 8:
      {
         const double psi = Psi(v - u, 3);
         const double dpsi = Dpsi(v - u, 3);
         s[0] = u * v * (-1) * dpsi + v * psi;
         s[1] = u * v * dpsi + u * psi;
         s[2] = 0.;
      }
      break;
      case 9:
      {
         const double psi = Psi((1. - u - v) - v, 3);
         const double dpsi = Dpsi((1. - u - v) - v, 3);
         s[0] = v * (1. - u - v) * (-1) * dpsi - v * psi;
         s[1] = v * (1. - u - v) * (-2) * dpsi + (1. - u - 2. * v) * psi;
         s[2] = 0.;
      }
      break;
      case 10:
         s[0] = (1. - 2. * u - v) * v;
         s[1] = (1. - u - 2. * v) * u;
         s[2] = 0.;
         break;
      default:
         assert(0);
         return;
         break;
   }
}

void GradHierarchicalApproxFunctionQuadrilateral(int ith, double u, double v, double s[])
{
   switch (ith)
   {
         // Degree 1
      case 1:
         s[0] = -0.25 * (1. - v);
         s[1] = -0.25 * (1. - u);
         s[2] = 0.;
         break;
      case 2:
         s[0] = 0.25 * (1. - v);
         s[1] = -0.25 * (1. + u);
         s[2] = 0.;
         break;
      case 3:
         s[0] = 0.25 * (1. + v);
         s[1] = 0.25 * (1. + u);
         s[2] = 0.;
         break;
      case 4:
         s[0] = -0.25 * (1. + v);
         s[1] = 0.25 * (1. - u);
         s[2] = 0.;
         break;
         // Degre 2
      case 5:
         s[0] = 0.5 * (1 - v) * Dphi(u, 2);
         s[1] = -0.5 * Phi(u, 2);
         s[2] = 0.;
         break;
      case 6:
         s[0] = 0.5 * Phi(v, 2);
         s[1] = 0.5 * (1 + u) * Dphi(v, 2);
         s[2] = 0.;
         break;
      case 7:
         s[0] = 0.5 * (1 + v) * Dphi(u, 2);
         s[1] = 0.5 * Phi(u, 2);
         s[2] = 0.;
         break;
      case 8:
         s[0] = -0.5 * Phi(v, 2);
         s[1] = 0.5 * (1 - u) * Dphi(v, 2);
         s[2] = 0.;
         break;
      case 9:
         s[0] = Dphi(u, 2) * Phi(v, 2);
         s[1] = Phi(u, 2) * Dphi(v, 2);
         s[2] = 0.0;
         break;
         // what follows is close to dangerous
         // never tested
         //     case 9  :	Phi(u, 2, phi); Dphi(u, 2, dphi);
         //       s[0] = 0.5*(1-v)*dphi; s[1] = -0.5*phi; s[2] = 0.; break;
         //       // Degre 3
         //     case 10 :   Phi(u, 2, phiu); Phi(v, 2, phiv);
         //       Dphi(u, 2, dphiu); Dphi(v, 2, dphiv);
         //       s[0] = dphiu*phiv; s[1] = phiu*dphiv; s[2] = 0.0; break;
         //     case 11 :	Phi(v, 2, phi); Dphi(v, 2, dphi);
         //       s[0] = 0.5*phi; s[1] = 0.5*(1+u)*dphi; s[2] = 0.; break;
         //     case 12 :	Phi(u, 2, phi); Dphi(u, 2, dphi);
         //       s[0] = 0.5*(1+v)*dphi;    s[1] = 0.5*phi; s[2] = 0.; break;
         //     case 13 :	Phi(v, 2, phi); Dphi(v, 2, dphi);
         //       s[0] = -0.5*phi; s[1] = 0.5*(1-u)*dphi; s[2] = 0.; break;
         //     case 14 :   Phi(u, 2, phiu); Phi(v, 3, phiv);
         //       Dphi(u, 2, dphiu); Dphi(v, 3, dphiv);
         //       s[0] = dphiu*phiv; s[1] = phiu*dphiv; s[2] = 0.0; break;
         //     case 15 :   Phi(u, 3, phiu); Phi(v, 2, phiv);
         //       Dphi(u, 3, dphiu); Dphi(v, 2, dphiv);
         //       s[0] = dphiu*phiv; s[1] = phiu*dphiv; s[2] = 0.0; break;
         //     case 16 :   Phi(u, 3, phiu); Phi(v, 3, phiv);
         //       Dphi(u, 3, dphiu); Dphi(v, 3, dphiv);
         //       s[0] = dphiu*phiv; s[1] = phiu*dphiv; s[2] = 0.0; break;
      default:
         assert(0);
         break;
   }
}

void GradHierarchicalApproxFunctionTetrahedron(int ith, double u, double v, double w, double s[])
{
   switch (ith)
   {
      case 1:
         s[0] = -1.;
         s[1] = -1.;
         s[2] = -1.;
         break;
      case 2:
         s[0] = 1.;
         s[1] = 0.;
         s[2] = 0.;
         break;
      case 3:
         s[0] = 0.;
         s[1] = 1.;
         s[2] = 0.;
         break;
      case 4:
         s[0] = 0.;
         s[1] = 0.;
         s[2] = 1.;
         break;
      case 5:
         s[0] = (1. - 2. * u - v - w);
         s[1] = -u;
         s[2] = -u;
         break;
      case 6:
         s[0] = v;
         s[1] = u;
         s[2] = 0.;
         break;
      case 7:
         s[0] = -v;
         s[1] = (1. - u - 2. * v - w);
         s[2] = -v;
         break;
      case 8:
         s[0] = -w;
         s[1] = -w;
         s[2] = (1. - u - v - 2. * w);
         break;
      case 9:
         s[0] = w;
         s[1] = 0.;
         s[2] = u;
         break;
      case 10:
         s[0] = 0.;
         s[1] = w;
         s[2] = v;
         break;
      default:
         assert(0);
         break;
   }
}

void GradHierarchicalApproxFunctionHexahedron(int ith, double u, double v, double w, double s[])
{
   switch (ith)
   {
      case 1:
         s[0] = -0.125 * (1. - v) * (1. - w);
         s[1] = -0.125 * (1. - u) * (1. - w);
         s[2] = -0.125 * (1. - u) * (1. - v);
         break;
      case 2:
         s[0] = 0.125 * (1. - v) * (1. - w);
         s[1] = -0.125 * (1. + u) * (1. - w);
         s[2] = -0.125 * (1. + u) * (1. - v);
         break;
      case 3:
         s[0] = 0.125 * (1. + v) * (1. - w);
         s[1] = 0.125 * (1. + u) * (1. - w);
         s[2] = -0.125 * (1. + u) * (1. + v);
         break;
      case 4:
         s[0] = -0.125 * (1. + v) * (1. - w);
         s[1] = 0.125 * (1. - u) * (1. - w);
         s[2] = -0.125 * (1. - u) * (1. + v);
         break;
      case 5:
         s[0] = -0.125 * (1. - v) * (1. + w);
         s[1] = -0.125 * (1. - u) * (1. + w);
         s[2] = 0.125 * (1. - u) * (1. - v);
         break;
      case 6:
         s[0] = 0.125 * (1. - v) * (1. + w);
         s[1] = -0.125 * (1. + u) * (1. + w);
         s[2] = 0.125 * (1. + u) * (1. - v);
         break;
      case 7:
         s[0] = 0.125 * (1. + v) * (1. + w);
         s[1] = 0.125 * (1. + u) * (1. + w);
         s[2] = 0.125 * (1. + u) * (1. + v);
         break;
      case 8:
         s[0] = -0.125 * (1. + v) * (1. + w);
         s[1] = 0.125 * (1. - u) * (1. + w);
         s[2] = 0.125 * (1. - u) * (1. + v);
         break;
      default:
         assert(0);
         break;
   }
}

void GradHierarchicalApproxFunctionPrism(int ith, double u, double v, double w, double s[])
{
   switch (ith)
   {
      case 1:
         s[0] = -0.5 * (1. - w);
         s[1] = -0.5 * (1. - w);
         s[2] = -0.5 * (1. - u - v);
         break;
      case 2:
         s[0] = 0.5 * (1. - w);
         s[1] = 0.;
         s[2] = -0.5 * u;
         break;
      case 3:
         s[0] = 0.;
         s[1] = 0.5 * (1. - w);
         s[2] = -0.5 * v;
         break;
      case 4:
         s[0] = -0.5 * (1. + w);
         s[1] = -0.5 * (1. + w);
         s[2] = 0.5 * (1. - u - v);
         break;
      case 5:
         s[0] = 0.5 * (1. + w);
         s[1] = 0.;
         s[2] = 0.5 * u;
         break;
      case 6:
         s[0] = 0.;
         s[1] = 0.5 * (1. + w);
         s[2] = 0.5 * v;
         break;
      default:
         assert(0);
         break;
   }
}

double Psi(double x, int n)
{
   switch (n)
   {
      case 2:
         return -sqrt(6.);
         break;
      case 3:
         return -sqrt(10.) * x;
         break;
      case 4:
         return -sqrt(7. / 8.) * (5. * x * x - 1.);
         break;
      default:
         assert((n >= 2) && (n <= 4));
         return 0.;
         break;
   }
}
double Dpsi(double x, int n)
{
   switch (n)
   {
      case 2:
         return 0.;
         break;
      case 3:
         return -sqrt(10.);
         break;
      case 4:
         return -sqrt(7. / 8.) * 10. * x;
         break;
      default:
         assert((n >= 2) && (n <= 4));
         return 0.;
         break;
   }
}
double Phi(double x, int n) { return (LegendrePol(x, n) - LegendrePol(x, n - 2)) / sqrt(2. * (2. * ((double)n) - 1)); }
double Dphi(double x, int n) { return sqrt((2. * ((double)n) - 1) / 2.) * LegendrePol(x, n - 1); }
double LegendrePol(double x, int n)
{
   switch (n)
   {
      case 0:
         return 1.;
         break;
      case 1:
         return x;
         break;
      case 2:
         return (3. * x * x - 1.) / 2.;
         break;
      case 3:
         return (5. * x * x * x - 3. * x) / 2.;
         break;
      case 4:
         return (35. * x * x * x * x - 30. * x * x + 3.) / 8.;
         break;
      case 5:
         return (63. * x * x * x * x * x - 70. * x * x * x + 15. * x) / 8.;
         break;
      case 6:
         return (231. * x * x * x * x * x * x - 315. * x * x * x * x + 105. * x * x - 5.) / 16.;
         break;
      case 7:
         return (429. * x * x * x * x * x * x * x - 693. * x * x * x * x * x + 315. * x * x * x - 35. * x) / 16.;
         break;
      default:
         assert((n >= 0) && (n <= 8));
         return 0.;
         break;
   }
}

double SimplexApproxFunctionQuadrilateral(int ith, double u, double v)
{
   double s0 = 0.0;
   int count = 0;

   if (v >= u)
   {
      switch (ith)
      {
         case 3:
            s0 += (u + 1.0) / 2.0;
            break;
         case 4:
            s0 += 1.0 - (u + 1.0) / 2.0 + (v - 1.0) / 2.0;
            break;
         case 1:
            s0 += -(v - 1.0) / 2.0;
            break;
         default:
            s0 += 0.;
            break;
      }
      count++;
   }

   if (u >= v)
   {
      switch (ith)
      {
         case 3:
            s0 += (v + 1.0) / 2.0;
            break;
         case 1:
            s0 += -(u - 1.0) / 2.0;
            break;
         case 2:
            s0 += 1.0 + (u - 1.0) / 2.0 - (v + 1.0) / 2.0;
            break;
         default:
            s0 += 0.;
            break;
      }
      count++;
   }

   if (count == 0)
   {
      assert(0);
   }
   return s0 / (double)count;
}

inline void vecMvec(double vec1[3], double vec2[3], double res[3])
{
   res[0] = vec1[0] - vec2[0];
   res[1] = vec1[1] - vec2[1];
   res[2] = vec1[2] - vec2[2];
   return;
}
inline void determinant(double a[3], double b[3], double c[3], double& vol)
{
   double res[3];

   res[0] = a[1] * b[2] - a[2] * b[1];
   res[1] = a[2] * b[0] - a[0] * b[2];
   res[2] = a[0] * b[1] - a[1] * b[0];

   vol = (res[0] * c[0] + res[1] * c[1] + res[2] * c[2]);
}

double SimplexApproxFunctionHexahedron(int ith, double u, double v, double w)
{
   double V145 = 2.0 / 3.0 * (u + 1.0);
   double V548 = V145;

   double V673 = 2.0 / 3.0 * (-u + 1.0);
   double V632 = V673;

   double V615 = 2.0 / 3.0 * (v + 1.0);
   double V621 = V615;

   double V743 = 2.0 / 3.0 * (-v + 1.0);
   double V847 = V743;

   double V241 = 2.0 / 3.0 * (w + 1.0);
   double V342 = V241;

   double V658 = 2.0 / 3.0 * (-w + 1.0);
   double V687 = V658;

   double vol;
   double a[3], b[3], c[3];

   // codinate of the nodes
   double N[8][3] = {{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0}, {1.0, 1.0, -1.0}, {-1.0, 1.0, -1.0},
                     {-1.0, -1.0, 1.0},  {1.0, -1.0, 1.0},  {1.0, 1.0, 1.0},  {-1.0, 1.0, 1.0}};

   double uvw[3] = {u, v, w};

   vecMvec(N[4 - 1], N[6 - 1], a);
   vecMvec(uvw, N[6 - 1], c);

   vecMvec(N[1 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V416 = vol / 6.0;
   double V461 = -V416;

   vecMvec(N[2 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V426 = vol / 6.0;
   double V462 = -V426;

   vecMvec(N[3 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V436 = vol / 6.0;
   double V463 = -V436;

   vecMvec(N[5 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V456 = vol / 6.0;
   double V465 = -V456;

   vecMvec(N[7 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V476 = vol / 6.0;
   double V467 = -V476;

   vecMvec(N[8 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V486 = vol / 6.0;
   double V468 = -V486;

   double vol_tetra = 8.0 / 6.0;

   // int con_2;
   int count = 0;
   double s0 = 0.0;

   double tolerance = 0.0;

   if (V486 >= tolerance && V467 >= tolerance)
   {  // && V687 >0. && V847 >0. ){
      count++;
      // 4-8-6-7   4-8-6  4-6-7  6-8-7  8-4-7
      switch (ith)
      {
         case 7:
            s0 += V486 / vol_tetra;
            break;
         case 8:
            s0 += V467 / vol_tetra;
            break;
         case 4:
            s0 += V687 / vol_tetra;
            break;
         case 6:
            s0 += V847 / vol_tetra;
            break;
         default:
            s0 += 0.;
            break;
      }
      //     con_2=1;
   }

   if (V476 >= tolerance && V463 >= tolerance)
   {  //&& V673 > 0. && V743 > 0. ){
      count++;
      // 4-7-6-3   4-7-6  4-6-3  6-7-3  7-4-3
      switch (ith)
      {
         case 3:
            s0 += V476 / vol_tetra;
            break;
         case 7:
            s0 += V463 / vol_tetra;
            break;
         case 4:
            s0 += V673 / vol_tetra;
            break;
         case 6:
            s0 += V743 / vol_tetra;
            break;
         default:
            s0 += 0.;
            break;
      }
      //     con_2=2;
   }

   if (V456 >= tolerance && V468 >= tolerance)
   {  //&& V658 > 0. && V548 > 0. ){
      count++;
      // 4-5-6-8   4-5-6  4-6-8  6-5-8  5-4-8
      switch (ith)
      {
         case 8:
            s0 += V456 / vol_tetra;
            break;
         case 5:
            s0 += V468 / vol_tetra;
            break;
         case 4:
            s0 += V658 / vol_tetra;
            break;
         case 6:
            s0 += V548 / vol_tetra;
            break;
         default:
            s0 += 0.;
            break;
      }
      //     con_2=3;
   }

   if (V436 >= tolerance && V462 >= tolerance)
   {  //&& V632 > 0. && V342 > 0. ){
      count++;
      // 4-3-6-2   4-3-6  4-6-2  6-3-2  3-4-2
      switch (ith)
      {
         case 2:
            s0 += V436 / vol_tetra;
            break;
         case 3:
            s0 += V462 / vol_tetra;
            break;
         case 4:
            s0 += V632 / vol_tetra;
            break;
         case 6:
            s0 += V342 / vol_tetra;
            break;
         default:
            s0 += 0.;
            break;
      }
      //     con_2=4;
   }

   if (V426 >= tolerance && V461 >= tolerance)
   {  //&& V621 > 0. && V241 > 0. ){
      count++;
      // 4-2-6-1   4-2-6  4-6-1  6-2-1  2-4-1
      switch (ith)
      {
         case 1:
            s0 += V426 / vol_tetra;
            break;
         case 2:
            s0 += V461 / vol_tetra;
            break;
         case 4:
            s0 += V621 / vol_tetra;
            break;
         case 6:
            s0 += V241 / vol_tetra;
            break;
         default:
            s0 += 0.;
            break;
      }
      // con_2=5;
   }

   if (V416 >= tolerance && V465 >= tolerance)
   {  //&& V615 > 0. && V145 > 0. ){
      count++;
      // 4-1-6-5   4-1-6  4-6-5  6-1-5  1-4-5
      switch (ith)
      {
         case 1:
            s0 += V465 / vol_tetra;
            break;
         case 4:
            s0 += V615 / vol_tetra;
            break;
         case 5:
            s0 += V416 / vol_tetra;
            break;
         case 6:
            s0 += V145 / vol_tetra;
            break;
         default:
            s0 += 0.;
            break;
      }
      //     con_2=6;
   }
   if (count == 0)
   {
      assert(0);
   }
   return s0 / (double)count;
}

void GradSimplexApproxFunctionQuadrilateral(int ith, double u, double v, double s[])
{
   // double s0 = 0.0;
   // int count = 0;

   if (u < v)
   {
      switch (ith)
      {
         case 3:
            s[0] = 0.5;
            s[1] = 0.0;
            s[2] = 0.0;
            break;
         case 4:
            s[0] = -0.5;
            s[1] = 0.5;
            s[2] = 0.0;
            break;
         case 1:
            s[0] = 0.0;
            s[1] = -0.5;
            s[2] = 0.0;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else if (u > v)
   {
      switch (ith)
      {
         case 3:
            s[0] = 0.0;
            s[1] = 0.5;
            s[2] = 0.0;
            break;
         case 1:
            s[0] = -0.5;
            s[1] = 0.0;
            s[2] = 0.0;
            break;
         case 2:
            s[0] = 0.5;
            s[1] = -0.5;
            s[2] = 0.0;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else
   {
      assert(0);
   }
}

void GradSimplexApproxFunctionHexahedron(int ith, double u, double v, double w, double s[])
{
   double vol;
   double a[3], b[3], c[3];

   // volumes of tetrahedron inside Hexadedron.
   // codinate of the nodes
   double N[8][3] = {{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0}, {1.0, 1.0, -1.0}, {-1.0, 1.0, -1.0},
                     {-1.0, -1.0, 1.0},  {1.0, -1.0, 1.0},  {1.0, 1.0, 1.0},  {-1.0, 1.0, 1.0}};

   double uvw[3] = {u, v, w};

   vecMvec(N[4 - 1], N[6 - 1], a);
   vecMvec(uvw, N[6 - 1], c);

   vecMvec(N[1 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V416 = vol / 6.0;
   double V461 = -V416;

   vecMvec(N[2 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V426 = vol / 6.0;
   double V462 = -V426;

   vecMvec(N[3 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V436 = vol / 6.0;
   double V463 = -V436;

   vecMvec(N[5 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V456 = vol / 6.0;
   double V465 = -V456;

   vecMvec(N[7 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V476 = vol / 6.0;
   double V467 = -V476;

   vecMvec(N[8 - 1], N[6 - 1], b);
   determinant(a, b, c, vol);
   double V486 = vol / 6.0;
   double V468 = -V486;

   if (V486 > 0. && V467 > 0.)
   {  // && V687 >=0. && V847 >=0. ){
      // 4-8-6-7   4-8-6  4-6-7  6-8-7  8-4-7
      switch (ith)
      {
         case 4:
            s[0] = 0.0;
            s[1] = 0.0;
            s[2] = -0.5;
            break;
         case 6:
            s[0] = 0.0;
            s[1] = -0.5;
            s[2] = 0.0;
            break;
         case 7:
            s[0] = 0.5;
            s[1] = 0.5;
            s[2] = 0.0;
            break;
         case 8:
            s[0] = -0.5;
            s[1] = 0.0;
            s[2] = 0.5;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else if (V476 > 0. && V463 > 0.)
   {  //&& V673 >= 0. && V743 >= 0. ){
      // 4-7-6-3   4-7-6  4-6-3  6-7-3  7-4-3
      switch (ith)
      {
         case 3:
            s[0] = 0.5;
            s[1] = 0.0;
            s[2] = -0.5;
            break;
         case 4:
            s[0] = -0.5;
            s[1] = 0.0;
            s[2] = 0.0;
            break;
         case 6:
            s[0] = 0.0;
            s[1] = -0.5;
            s[2] = 0.0;
            break;
         case 7:
            s[0] = 0.0;
            s[1] = 0.5;
            s[2] = 0.5;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else if (V456 > 0. && V468 > 0.)
   {  //&& V658 >= 0. && V548 >= 0. ){
      // 4-5-6-8   4-5-6  4-6-8  6-5-8  5-4-8
      switch (ith)
      {
         case 4:
            s[0] = 0.0;
            s[1] = 0.0;
            s[2] = -0.5;
            break;
         case 5:
            s[0] = -0.5;
            s[1] = -0.5;
            s[2] = 0.0;
            break;
         case 6:
            s[0] = 0.5;
            s[1] = 0.0;
            s[2] = 0.0;
            break;
         case 8:
            s[0] = 0.0;
            s[1] = 0.5;
            s[2] = 0.5;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else if (V436 > 0. && V462 > 0.)
   {  //&& V632 >= 0. && V342 >= 0. ){
      // 4-3-6-2   4-3-6  4-6-2  6-3-2  3-4-2
      switch (ith)
      {
         case 2:
            s[0] = 0.0;
            s[1] = -0.5;
            s[2] = -0.5;
            break;
         case 3:
            s[0] = 0.5;
            s[1] = 0.5;
            s[2] = 0.0;
            break;
         case 4:
            s[0] = -0.5;
            s[1] = 0.0;
            s[2] = 0.0;
            break;
         case 6:
            s[0] = 0.0;
            s[1] = 0.0;
            s[2] = 0.5;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else if (V426 > 0. && V461 > 0.)
   {  //&& V621 >= 0. && V241 >= 0. ){
      // 4-2-6-1   4-2-6  4-6-1  6-2-1  2-4-1
      switch (ith)
      {
         case 1:
            s[0] = -0.5;
            s[1] = -0.5;
            s[2] = 0.0;
            break;
         case 2:
            s[0] = 0.5;
            s[1] = 0.0;
            s[2] = -0.5;
            break;
         case 4:
            s[0] = 0.0;
            s[1] = 0.5;
            s[2] = 0.0;
            break;
         case 6:
            s[0] = 0.0;
            s[1] = 0.0;
            s[2] = 0.5;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else if (V416 > 0. && V465 > 0.)
   {  //&& V621 >= 0. && V241 >= 0. ){
      // 4-1-6-5   4-1-6  4-6-5  6-1-5  1-4-5
      switch (ith)
      {
         case 1:
            s[0] = 0.0;
            s[1] = -0.5;
            s[2] = -0.5;
            break;
         case 4:
            s[0] = 0.0;
            s[1] = 0.5;
            s[2] = 0.0;
            break;
         case 5:
            s[0] = -0.5;
            s[1] = 0.0;
            s[2] = 0.5;
            break;
         case 6:
            s[0] = 0.5;
            s[1] = 0.0;
            s[2] = 0.0;
            break;
         default:
            s[0] = s[1] = s[2] = 0.0;
            break;
      }
   }
   else
   {
      printf(" u v w = %12.8e , %12.8e, %12.8e is on an subelement boudary.\n", u, v, w);
      assert(0);
   }
}

xScalarFunctionDerivDiscXFEM::xScalarFunctionDerivDiscXFEM(const xLevelSet& ls) : field(ls) {}

void xScalarFunctionDerivDiscXFEM::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();

   if (!field.isDefinedOnElement(e))
   {
      res = 0.;
      storage.store(geo_appro, geo_integ, res);
      return;
   }

   xElement elem(e);
   elem.setUvw(geo_appro->getUVW());

   std::vector<double> vals = field.getVals(e);

   if (*std::min_element(vals.begin(), vals.end()) * *std::max_element(vals.begin(), vals.end()) >= 0.)
   {
      res = 0.;
      storage.store(geo_appro, geo_integ, res);
      return;
   }

   double f = elem.getInterpoSca(vals);
   std::transform(vals.begin(), vals.end(), vals.begin(), xtool::xAbs<double>());
   double fa = elem.getInterpoSca(vals);

   res = fa - fabs(f);
   storage.store(geo_appro, geo_integ, res);
}

void xScalarFunctionDerivDiscXFEM::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   const bool debug = xdebug_flag;
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_appro->getEntity();

   if (debug) cout << " xScalarFunctionDerivDiscXFEM::getGrad for element " << endl;
   if (debug) e->print();

   if (!field.isDefinedOnElement(e))
   {
      res = xtensor::xVector<>(0., 0., 0.);
      storage.store(geo_appro, geo_integ, res);
      return;
   }

   if (debug)
   {
      cout << " field is defined on the element and level set values are " << endl;
      std::vector<double> ls_v = field.getVals(e);
      std::copy(ls_v.begin(), ls_v.end(), std::ostream_iterator<double>(cout, " "));
   }

   xElement elem(e);
   elem.setUvw(geo_appro->getUVW());

   std::vector<double> vals = field.getVals(e);

   if (*std::min_element(vals.begin(), vals.end()) * *std::max_element(vals.begin(), vals.end()) >= 0.)
   {
      res = xtensor::xVector<>(0., 0., 0.);
      storage.store(geo_appro, geo_integ, res);
      return;
   }

   double valLS = elem.getInterpoSca(vals);
   double s = elem.getInterpoSca(vals) > 0.0 ? 1. : -1.;
   std::transform(vals.begin(), vals.end(), vals.begin(), xtool::xAbs<double>());
   xtensor::xVector<> fa = elem.getGradInterpoSca(vals);
   xtensor::xVector<> f = field.getGrad(e);

   if (fabs(valLS) < 1e-3)
   {
      if (geo_integ != nullptr)
      {
         elem.xyz2uvw(geo_integ->getCDGxyz());
      }
      else
      {
         elem.setUvw(geo_appro->getCDGuvw());
      }

      std::vector<double> vals = field.getVals(e);
      s = elem.getInterpoSca(vals) > 0.0 ? 1. : -1.;
      std::transform(vals.begin(), vals.end(), vals.begin(), xtool::xAbs<double>());
      fa = elem.getGradInterpoSca(vals);
      f = field.getGrad(e);
   }

   xtensor::xVector<> resbis(fa - f * s);
   res = resbis;

   storage.store(geo_appro, geo_integ, res);
   return;
}

xScalarFunctionRegularizedHeaviside::xScalarFunctionRegularizedHeaviside(const xEntityFilter& _f_in_isozero,
                                                                         const xEntityFilter& _f_subelement_regularized,
                                                                         const xEntityFilter& _f_value_one)
    : filter_in_iso_zero(_f_in_isozero), filter_subelement_regularized(_f_subelement_regularized), filter_value_one(_f_value_one)
{
}

void xScalarFunctionRegularizedHeaviside::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_integ->getEntity();

   if (filter_subelement_regularized(e))
   {  // requires enrichment, linear part
      std::vector<double> vals;
      for (int inode = 0; inode < e->size(0); inode++)
      {
         mEntity* v = e->get(0, inode);
         vals.push_back(filter_in_iso_zero(v));
      }
      xElement elem(e);
      elem.setUvw(geo_integ->getUVW());
      res = elem.getInterpoSca(vals);
   }
   else if (filter_value_one(e))
      res = 1.;
   else
      res = 0.;
   storage.store(geo_appro, geo_integ, res);
   return;
}

void xScalarFunctionRegularizedHeaviside::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                  xtensor::xVector<>& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   mEntity* e = geo_integ->getEntity();

   if (filter_subelement_regularized(e))
   {  // requires enrichment, linear part
      std::vector<double> vals;
      for (int inode = 0; inode < e->size(0); inode++)
      {
         mEntity* v = e->get(0, inode);
         vals.push_back(filter_in_iso_zero(v));
      }
      xElement elem(e);
      elem.setUvw(geo_integ->getUVW());
      res = elem.getGradInterpoSca(vals);
   }
   else
      res = xtensor::xVector<>(0., 0., 0.);
   storage.store(geo_appro, geo_integ, res);
   return;
}

xScalarFunctionDiscXFEM::xScalarFunctionDiscXFEM(const xLevelSet& ls_) : ls(ls_) {}

void xScalarFunctionDiscXFEM::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   if (storage.exist(geo_appro, geo_integ, res)) return;
   res = (double)ls.side_of(geo_appro, geo_integ);
   storage.store(geo_appro, geo_integ, res);
}

void xScalarFunctionDiscXFEM::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const
{
   res[0] = 0.0;
   res[1] = 0.0;
   res[2] = 0.0;
}

// <<<<<<< bizarre, apres le update du 09/01/09
// //   std::vector<double> vals = field->getVals(e);

// //   if (*std::min_element(vals.begin(), vals.end()) * *std::max_element(vals.begin(), vals.end()) >= 0.)
// //     {
// //       res = xtensor::xVector<>(0.,0.,0.);
// //       //store(geo_appro, geo_integ, res);
// //       return;
// //     }

//   xtensor::xVector<> f =  elem.getGradInterpoSca(vals);
//   std::transform(vals.begin(), vals.end(), vals.begin(), xAbs());
//   xtensor::xVector<> fa = elem.getGradInterpoSca(vals);

//   int s = surf->side_of(geo_appro, geo_integ);
//   std::cout << geo_appro << " " << geo_integ << " " << s << std::endl;
//   xtensor::xPoint p = geo_appro->getUVW();
//   std::cout<<  p(0)  << " " <<  p(1)<<" "  << p(2) << std::endl;
//   // std::cout << vals[0] << " " << vals[1] <<  " " << vals[2] << std::endl;
//   //std::cout << fa << " " << f << std::endl;
//   xtensor::xVector<> resbis(f - f * s);
//   res=resbis;

//   store(geo_appro, geo_integ, res);
//   return;
// }
// <<<<

/* needed data :
 *  xMesh::get_const_was_created_by().
 *  xMesh::get_levelset_value_tag()
 *  xMesh::get_cut_edge_tag()
 *  xMesh::get_refined_elements_tag()
 */

void xApproxFunctionEnrichedXFEMShifted::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   const bool debug = xdebug_flag;
   double b, e;
   base->getVal(geo_appro, geo_integ, b);
   cout << geo_integ->getXYZ() << endl;
   enrichment->getVal(geo_appro, geo_integ, e);

   xGeomElem geo_entity(enti);
   double ee;
   enrichment->getVal(geo_appro, &geo_entity, ee);
   if (enti) enti->print();

   res = b * (e - ee);

   if (debug)
   {
      cout << " enrichment at point " << geo_appro->getXYZ() << " is " << e;
      cout << " ScalarxApproxFunctionEnrichedXFEM::getVal " << res << endl;
   }
}

void xApproxFunctionEnrichedXFEMShifted::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                 xtensor::xVector<>& res) const
{
   const bool debug = xdebug_flag;
   xtensor::xVector<> bg, eg;
   double es, bs;

   base->getGrad(geo_appro, geo_integ, bg);
   enrichment->getVal(geo_appro, geo_integ, es);
   enrichment->getGrad(geo_appro, geo_integ, eg);
   base->getVal(geo_appro, geo_integ, bs);

   xGeomElem geo_entity(enti);
   double ee;
   enrichment->getVal(geo_appro, &geo_entity, ee);

   res = bg * (es - ee) + eg * bs;

   if (debug)
   {
      cout << " enrichment  at point " << geo_appro->getXYZ() << " is " << es << endl;
      cout << " grad enrichment " << eg << endl;
      cout << " ScalarxApproxFunctionEnrichedXFEM::getGrad " << res << endl;
   }

   return;
}

void xApproxFunctionEnrichedXFEMShifted::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                xtensor::xVector<>& res) const
{
   xtensor::xVector<> bg;
   double es;
   base->getVal(geo_appro, geo_integ, bg);
   enrichment->getVal(geo_appro, geo_integ, es);

   xGeomElem geo_entity(enti);
   xGeomElem geo_appro2(geo_appro->getEntity());
   geo_appro2.setUVWForXYZ(geo_entity.getXYZ());
   double ee;
   enrichment->getVal(&geo_appro2, &geo_entity, ee);

   res = bg * (es - ee);
   return;
}

void xApproxFunctionEnrichedXFEMShifted::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                 xtensor::xTensor2<>& res) const
{
   const bool debug = xdebug_flag;
   xtensor::xTensor2<> bt;
   xtensor::xVector<> eg, bv;
   double es;
   base->getGrad(geo_appro, geo_integ, bt);
   enrichment->getVal(geo_appro, geo_integ, es);
   enrichment->getGrad(geo_appro, geo_integ, eg);
   base->getVal(geo_appro, geo_integ, bv);

   xGeomElem geo_entity(enti);
   xGeomElem geo_appro2(geo_appro->getEntity());
   geo_appro2.setUVWForXYZ(geo_entity.getXYZ());
   double ee;
   enrichment->getVal(&geo_appro2, &geo_entity, ee);

   res = bt * (es - ee) + tensor_product(bv, eg);
   if (debug)
   {
      cout << " enrichment  at point " << geo_appro->getXYZ() << " is " << es << endl;
      cout << " grad enrichment " << eg << endl;
      cout << "val classical " << endl;
      cout << bv << endl;
      cout << "grad classical " << endl;
      cout << bt << endl;

      cout << "first part " << endl;
      cout << bt * es << endl;
      cout << "second part " << endl;
      cout << xtensor::tensor_product(eg, bv) << endl;
   }
}

}  // namespace xfem
