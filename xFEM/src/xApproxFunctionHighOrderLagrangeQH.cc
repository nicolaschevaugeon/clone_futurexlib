/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <cmath>

#include "xApproxFunctionHighOrder.h"
#include "xApproxFunctionHighOrderQH.h"
#include "xGeomElem.h"

using AOMD::mEntity;
using std::cout;
using std::endl;

// Max number of tabulated Lagrange functions
#define XAPPROXFUNCTIONLAGRANGEQHMAXORDER 15

namespace xfem
{
///////////////////////////////////////////////////////////////////////////////////////////////
// xApproxFunction Lagrange Implementation     QUADS AND HEXAHEDRA                           //
///////////////////////////////////////////////////////////////////////////////////////////////

xApproxFunctionScalarVertexLagrangeQH::xApproxFunctionScalarVertexLagrangeQH(int _order, int _node) : order(_order), inode(_node)
{
   if (order > XAPPROXFUNCTIONLAGRANGEQHMAXORDER)
   {
      cout << "Tabulated functions up to order XAPPROXFUNCTIONLAGRANGEQHMAXORDER !\n";
      throw;
   }
}

const int xApproxFunctionScalarVertexLagrangeQH::Vindices[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                                                   {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

void xApproxFunctionScalarVertexLagrangeQH::correctUvwForCurrentEntity(const mEntity::mType& type, const xtensor::xPoint& uvw,
                                                                       xtensor::xPoint& uvwCorr)
{
   switch (type)
   {
      case mEntity::HEX:
      {
         uvwCorr = uvw;
         break;
      }
      case mEntity::QUAD:
      {
         uvwCorr = {uvw(0), uvw(1), -1.};
         break;
      }
      case mEntity::EDGE:
      {
         uvwCorr = {uvw(0), -1., -1.};
         break;
      }
      case mEntity::VERTEX:
      {
         uvwCorr = {-1., -1., -1.};
         break;
      }
      default:
      {
         throw;
      }
   }
}

void xApproxFunctionScalarVertexLagrangeQH::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvwLoc;
   correctUvwForCurrentEntity(geo_appro->getEntity()->getType(), geo_appro->getUVW(), uvwLoc);
   res = lagPoly(order, Vindices[inode][0] * order, uvwLoc(0)) * lagPoly(order, Vindices[inode][1] * order, uvwLoc(1)) *
         lagPoly(order, Vindices[inode][2] * order, uvwLoc(2));
   return;
}

void xApproxFunctionScalarVertexLagrangeQH::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                    xtensor::xVector<>& res) const
{
   xtensor::xPoint uvwLoc;
   correctUvwForCurrentEntity(geo_appro->getEntity()->getType(), geo_appro->getUVW(), uvwLoc);
   res[0] = dLagPoly(order, Vindices[inode][0] * order, uvwLoc(0)) * lagPoly(order, Vindices[inode][1] * order, uvwLoc(1)) *
            lagPoly(order, Vindices[inode][2] * order, uvwLoc(2));
   res[1] = dLagPoly(order, Vindices[inode][1] * order, uvwLoc(1)) * lagPoly(order, Vindices[inode][0] * order, uvwLoc(0)) *
            lagPoly(order, Vindices[inode][2] * order, uvwLoc(2));
   if (geo_appro->getEntity()->getType() == mEntity::HEX)
   {
      res[2] = dLagPoly(order, Vindices[inode][2] * order, uvwLoc(2)) * lagPoly(order, Vindices[inode][0] * order, uvwLoc(0)) *
               lagPoly(order, Vindices[inode][1] * order, uvwLoc(1));
   }

   res = (geo_appro)->PushBack(res);
   return;
}

xApproxFunctionScalarEdgeLagrangeQH::xApproxFunctionScalarEdgeLagrangeQH(int _order, int _edge, int _ionedge)
    : order(_order), ionedge(_ionedge), edge(_edge)
{
   if (order > XAPPROXFUNCTIONLAGRANGEQHMAXORDER)
   {
      cout << "Tabulated functions up to order XAPPROXFUNCTIONLAGRANGEQHMAXORDER !\n";
      throw;
   }
}

/// The shape function is the 2D,3D shape function associated to the edge
/// and NOT the 2D,3D shape function evaluated on the edge (it as an evolution IN the element)
///
void xApproxFunctionScalarEdgeLagrangeQH::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);

   // On corrige le cas de l'edge (sinon, on a v et w nuls et donc la fonction de forme aussi)
   if (geo_appro->getEntity()->getType() == mEntity::EDGE)
   {
      v = -1.;
      w = -1.;
   }

   //   -1      0      1
   //    *------*------*
   //    a      b      c
   // lagPoly(order, 0  , ) associe au noeud a
   // lagPoly(order, order  , ) associe au noeud c

   switch (edge)
   {
      case 0:  // 2D et 3D ici...
         res = lagPoly(order, 1 + ionedge, u) * lagPoly(order, 0, v);
         break;
      case 1:
         res = lagPoly(order, 1 + ionedge, v) * lagPoly(order, order, u);
         break;
      case 2:
         res = lagPoly(order, 1 + ionedge, u) * lagPoly(order, order, v);
         break;
      case 3:
         res = lagPoly(order, 1 + ionedge, v) * lagPoly(order, 0, u);
         break;
      case 4:  // Only 3D here...
         res = lagPoly(order, 1 + ionedge, w) * lagPoly(order, 0, u) * lagPoly(order, 0, v);
         break;
      case 5:
         res = lagPoly(order, 1 + ionedge, w) * lagPoly(order, order, u) * lagPoly(order, 0, v);
         break;
      case 6:
         res = lagPoly(order, 1 + ionedge, w) * lagPoly(order, order, u) * lagPoly(order, order, v);
         break;
      case 7:
         res = lagPoly(order, 1 + ionedge, w) * lagPoly(order, 0, u) * lagPoly(order, order, v);
         break;
      case 8:
         res = lagPoly(order, 1 + ionedge, u) * lagPoly(order, 0, v) * lagPoly(order, order, w);
         break;
      case 9:
         res = lagPoly(order, 1 + ionedge, v) * lagPoly(order, order, u) * lagPoly(order, order, w);
         break;
      case 10:
         res = lagPoly(order, 1 + ionedge, u) * lagPoly(order, order, v) * lagPoly(order, order, w);
         break;
      case 11:
         res = lagPoly(order, 1 + ionedge, v) * lagPoly(order, 0, u) * lagPoly(order, order, w);
         break;
      default:
         res = 0.;
         break;
   }

   // Correction pour les Hex:
   if (geo_appro->getEntity()->getType() == mEntity::HEX && edge < 4)
   {
      res *= lagPoly(order, 0, w);
   }

   //    cout<<"xApproxFunctionScalarEdgeLagrangeQH "<<res<<endl;
   return;
}

void xApproxFunctionScalarEdgeLagrangeQH::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                       xtensor::xVector<>& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);
   bool tridim{false};
   if (geo_appro->getEntity()->getLevel() == 3) tridim = true;

   // On corrige le cas de l'edge (sinon, on a v et w nuls et donc la fonction de forme aussi)
   if (geo_appro->getEntity()->getType() == mEntity::EDGE)
   {
      v = -1.;
      w = -1.;
   }

   //   -1      0      1
   //    *------*------*
   //    a      b      c
   // lagPoly(order, 0  , ) associe au noeud a
   // lagPoly(order, order  , ) associe au noeud c

   switch (edge)
   {
      case 0:  // 2D et 3D ici...
         res[0] = dLagPoly(order, 1 + ionedge, u) * lagPoly(order, 0, v);
         res[1] = lagPoly(order, 1 + ionedge, u) * dLagPoly(order, 0, v);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= lagPoly(order, 0, w);
            res[1] *= lagPoly(order, 0, w);
            res[2] = lagPoly(order, 1 + ionedge, u) * lagPoly(order, 0, v) * dLagPoly(order, 0, w);
         }

         break;
      case 1:
         res[0] = lagPoly(order, 1 + ionedge, v) * dLagPoly(order, order, u);
         res[1] = dLagPoly(order, 1 + ionedge, v) * lagPoly(order, order, u);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= lagPoly(order, 0, w);
            res[1] *= lagPoly(order, 0, w);
            res[2] = lagPoly(order, 1 + ionedge, v) * lagPoly(order, order, u) * dLagPoly(order, 0, w);
         }

         break;
      case 2:
         res[0] = dLagPoly(order, 1 + ionedge, u) * lagPoly(order, order, v);
         res[1] = lagPoly(order, 1 + ionedge, u) * dLagPoly(order, order, v);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= lagPoly(order, 0, w);
            res[1] *= lagPoly(order, 0, w);
            res[2] = lagPoly(order, 1 + ionedge, u) * lagPoly(order, order, v) * dLagPoly(order, 0, w);
         }

         break;
      case 3:
         res[0] = lagPoly(order, 1 + ionedge, v) * dLagPoly(order, 0, u);
         res[1] = dLagPoly(order, 1 + ionedge, v) * lagPoly(order, 0, u);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= lagPoly(order, 0, w);
            res[1] *= lagPoly(order, 0, w);
            res[2] = lagPoly(order, 1 + ionedge, v) * lagPoly(order, 0, u) * dLagPoly(order, 0, w);
         }

         break;
      case 4:  // Only 3D here...
         res[0] = lagPoly(order, 1 + ionedge, w) * dLagPoly(order, 0, u) * lagPoly(order, 0, v);
         res[1] = lagPoly(order, 1 + ionedge, w) * lagPoly(order, 0, u) * dLagPoly(order, 0, v);
         res[2] = dLagPoly(order, 1 + ionedge, w) * lagPoly(order, 0, u) * lagPoly(order, 0, v);
         break;
      case 5:
         res[0] = lagPoly(order, 1 + ionedge, w) * dLagPoly(order, order, u) * lagPoly(order, 0, v);
         res[1] = lagPoly(order, 1 + ionedge, w) * lagPoly(order, order, u) * dLagPoly(order, 0, v);
         res[2] = dLagPoly(order, 1 + ionedge, w) * lagPoly(order, order, u) * lagPoly(order, 0, v);
         break;
      case 6:
         res[0] = lagPoly(order, 1 + ionedge, w) * dLagPoly(order, order, u) * lagPoly(order, order, v);
         res[1] = lagPoly(order, 1 + ionedge, w) * lagPoly(order, order, u) * dLagPoly(order, order, v);
         res[2] = dLagPoly(order, 1 + ionedge, w) * lagPoly(order, order, u) * lagPoly(order, order, v);
         break;
      case 7:
         res[0] = lagPoly(order, 1 + ionedge, w) * dLagPoly(order, 0, u) * lagPoly(order, order, v);
         res[1] = lagPoly(order, 1 + ionedge, w) * lagPoly(order, 0, u) * dLagPoly(order, order, v);
         res[2] = dLagPoly(order, 1 + ionedge, w) * lagPoly(order, 0, u) * lagPoly(order, order, v);
         break;
      case 8:
         res[0] = dLagPoly(order, 1 + ionedge, u) * lagPoly(order, 0, v) * lagPoly(order, order, w);
         res[1] = lagPoly(order, 1 + ionedge, u) * dLagPoly(order, 0, v) * lagPoly(order, order, w);
         res[2] = lagPoly(order, 1 + ionedge, u) * lagPoly(order, 0, v) * dLagPoly(order, order, w);
         break;
      case 9:
         res[0] = lagPoly(order, 1 + ionedge, v) * dLagPoly(order, order, u) * lagPoly(order, order, w);
         res[1] = dLagPoly(order, 1 + ionedge, v) * lagPoly(order, order, u) * lagPoly(order, order, w);
         res[2] = lagPoly(order, 1 + ionedge, v) * lagPoly(order, order, u) * dLagPoly(order, order, w);
         break;
      case 10:
         res[0] = dLagPoly(order, 1 + ionedge, u) * lagPoly(order, order, v) * lagPoly(order, order, w);
         res[1] = lagPoly(order, 1 + ionedge, u) * dLagPoly(order, order, v) * lagPoly(order, order, w);
         res[2] = lagPoly(order, 1 + ionedge, u) * lagPoly(order, order, v) * dLagPoly(order, order, w);
         break;
      case 11:
         res[0] = lagPoly(order, 1 + ionedge, v) * dLagPoly(order, 0, u) * lagPoly(order, order, w);
         res[1] = dLagPoly(order, 1 + ionedge, v) * lagPoly(order, 0, u) * lagPoly(order, order, w);
         res[2] = lagPoly(order, 1 + ionedge, v) * lagPoly(order, 0, u) * dLagPoly(order, order, w);
         break;
      default:
         res[0] = 0.;
         res[1] = 0.;
         res[2] = 0.;
         break;
   }
}

void xApproxFunctionScalarEdgeLagrangeQH::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                  xtensor::xVector<>& res) const
{
   getGradLocal(geo_appro, geo_integ, res);
   res = (geo_appro)->PushBack(res);
}

void xApproxFunctionScalarQuadLagrange::setNodeIndices(int order, int ionface, std::vector<std::vector<int>>& indices)
{
   if (indices.empty())
   {
      indices.resize((order - 1) * (order - 1));
      std::vector<int> idx(3);
      for (int j = 1; j < order; ++j)
      {
         for (int i = 1; i < order; ++i)
         {
            idx[0] = i;
            idx[1] = j;
            idx[2] = 0;
            indices[(i - 1) + (order - 1) * (j - 1)] = idx;
         }
      }
   }
}

xApproxFunctionScalarQuadLagrange::xApproxFunctionScalarQuadLagrange(int _order, int _face, int _ionface)
    : order(_order), ionface(_ionface), face(_face)
{
   if (order > XAPPROXFUNCTIONLAGRANGEQHMAXORDER)
   {
      cout << "Tabulated functions up to order XAPPROXFUNCTIONLAGRANGEQHMAXORDER !\n";
      throw;
   }

   setNodeIndices(order, ionface, indices);
}

void xApproxFunctionScalarQuadLagrange::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);

   switch (face)
   {
      case 0:  // 2D et 3D ici...
         res = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], v);

         break;
      case 1:  // Only 3D here...
         res = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, 0, v);
         break;
      case 2:
         res = lagPoly(order, indices[ionface][0], v) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, order, u);
         break;
      case 3:
         res = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, order, v);
         break;
      case 4:
         res = lagPoly(order, indices[ionface][0], v) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, 0, u);
         break;
      case 5:
         res = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], v) * lagPoly(order, order, w);
         break;
      default:
         res = 0.;
         break;
   }

   // Correction pour les Hex:
   if (geo_appro->getEntity()->getType() == mEntity::HEX && face < 1)
   {
      res *= lagPoly(order, 0, w);
   }

   return;
}

void xApproxFunctionScalarQuadLagrange::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                     xtensor::xVector<>& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);

   bool tridim{false};
   if (geo_appro->getEntity()->getLevel() == 3) tridim = true;

   switch (face)
   {
      case 0:  // 2D et 3D ici...
         res[0] = dLagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], v);
         res[1] = lagPoly(order, indices[ionface][0], u) * dLagPoly(order, indices[ionface][1], v);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= lagPoly(order, 0, w);
            res[1] *= lagPoly(order, 0, w);
            res[2] = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], v) * dLagPoly(order, 0, w);
         }
         break;
      case 1:  // Only 3D here...
         res[0] = dLagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, 0, v);
         res[1] = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], w) * dLagPoly(order, 0, v);
         res[2] = lagPoly(order, indices[ionface][0], u) * dLagPoly(order, indices[ionface][1], w) * lagPoly(order, 0, v);
         break;
      case 2:
         res[0] = lagPoly(order, indices[ionface][0], v) * lagPoly(order, indices[ionface][1], w) * dLagPoly(order, order, u);
         res[1] = dLagPoly(order, indices[ionface][0], v) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, order, u);
         res[2] = lagPoly(order, indices[ionface][0], v) * dLagPoly(order, indices[ionface][1], w) * lagPoly(order, order, u);
         break;
      case 3:
         res[0] = dLagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, order, v);
         res[1] = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], w) * dLagPoly(order, order, v);
         res[2] = lagPoly(order, indices[ionface][0], u) * dLagPoly(order, indices[ionface][1], w) * lagPoly(order, order, v);
         break;
      case 4:
         res[0] = lagPoly(order, indices[ionface][0], v) * lagPoly(order, indices[ionface][1], w) * dLagPoly(order, 0, u);
         res[1] = dLagPoly(order, indices[ionface][0], v) * lagPoly(order, indices[ionface][1], w) * lagPoly(order, 0, u);
         res[2] = lagPoly(order, indices[ionface][0], v) * dLagPoly(order, indices[ionface][1], w) * lagPoly(order, 0, u);
         break;
      case 5:
         res[0] = dLagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], v) * lagPoly(order, order, w);
         res[1] = lagPoly(order, indices[ionface][0], u) * dLagPoly(order, indices[ionface][1], v) * lagPoly(order, order, w);
         res[2] = lagPoly(order, indices[ionface][0], u) * lagPoly(order, indices[ionface][1], v) * dLagPoly(order, order, w);
         break;
      default:
         res = 0.;
         break;
   }

   return;
}

void xApproxFunctionScalarQuadLagrange::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                xtensor::xVector<>& res) const
{
   getGradLocal(geo_appro, geo_integ, res);
   res = (geo_appro)->PushBack(res);
}

void xApproxFunctionScalarHexLagrange::setNodeIndices(int order, int ionhex, std::vector<std::vector<int>>& indices)
{
   if (indices.empty())
   {
      indices.resize((order - 1) * (order - 1) * (order - 1));
      std::vector<int> idx(3);
      for (int k = 1; k < order; ++k)
      {
         for (int j = 1; j < order; ++j)
         {
            for (int i = 1; i < order; ++i)
            {
               idx[0] = i;
               idx[1] = j;
               idx[2] = k;
               indices[(i - 1) + (order - 1) * (j - 1) + (order - 1) * (order - 1) * (k - 1)] = idx;
            }
         }
      }
   }
}

xApproxFunctionScalarHexLagrange::xApproxFunctionScalarHexLagrange(int _order, int _ionhex) : order(_order), ionhex(_ionhex)
{
   setNodeIndices(order, ionhex, indices);
}

void xApproxFunctionScalarHexLagrange::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvwLoc = geo_appro->getUVW();
   res = lagPoly(order, indices[ionhex][0], uvwLoc(0)) * lagPoly(order, indices[ionhex][1], uvwLoc(1)) *
         lagPoly(order, indices[ionhex][2], uvwLoc(2));
   return;
}

void xApproxFunctionScalarHexLagrange::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                               xtensor::xVector<>& res) const
{
   xtensor::xPoint uvwLoc = geo_appro->getUVW();

   res[0] = dLagPoly(order, indices[ionhex][0], uvwLoc(0)) * lagPoly(order, indices[ionhex][1], uvwLoc(1)) *
            lagPoly(order, indices[ionhex][2], uvwLoc(2));
   res[1] = dLagPoly(order, indices[ionhex][1], uvwLoc(1)) * lagPoly(order, indices[ionhex][0], uvwLoc(0)) *
            lagPoly(order, indices[ionhex][2], uvwLoc(2));
   res[2] = dLagPoly(order, indices[ionhex][2], uvwLoc(2)) * lagPoly(order, indices[ionhex][0], uvwLoc(0)) *
            lagPoly(order, indices[ionhex][1], uvwLoc(1));

   res = (geo_appro)->PushBack(res);
   return;
}

// Tabulated Lagrange polynomial functions (for efficiency)
// Thank you sympy !
std::function<double(double)> lagPolys[136] = {
    [](double u) -> double { return 1.; },
    [](double u) -> double { return -0.5 * u + 0.5; },
    [](double u) -> double { return 0.5 * u + 0.5; },
    [](double u) -> double { return 0.5 * u * (u - 1.0); },
    [](double u) -> double { return -1.0 * u * u + 1.0; },
    [](double u) -> double { return 0.5 * u * (u + 1.0); },
    [](double u) -> double { return -0.5625 * u * u * u + 0.5625 * u * u + 0.0625 * u - 0.0625; },
    [](double u) -> double { return 1.124999999999999778 * (u - 1.0) * (u + 1.0) * (1.5 * u - 0.49999999999999988898); },
    [](double u) -> double { return -1.124999999999999778 * (u - 1.0) * (u + 1.0) * (1.5 * u + 0.5); },
    [](double u) -> double { return 0.5 * (0.75 * u + 0.25) * (u + 1.0) * (1.5 * u - 0.49999999999999983347); },
    [](double u) -> double { return 0.5 * u * (0.666666666666667 * u - 0.33333333333333331483) * (u - 1.0) * (2.0 * u + 1.0); },
    [](double u) -> double { return -2.6666666666666665186 * u * (u - 1.0) * (1.0 * u - 0.5) * (u + 1.0); },
    [](double u) -> double { return 4.0 * u * u * u * u - 5.0 * u * u + 1.0; },
    [](double u) -> double { return -2.6666666666666665186 * u * (u - 1.0) * (1.0 * u + 0.5) * (u + 1.0); },
    [](double u) -> double { return 0.5 * u * (0.666666666666667 * u + 0.33333333333333331483) * (u + 1.0) * (2.0 * u - 1.0); },
    [](double u) -> double {
       return -0.813802083333333 * u * u * u * u * u + 0.813802083333333 * u * u * u * u + 0.325520833333333 * u * u * u -
              0.325520833333333 * u * u - 0.01171875 * u + 0.01171875;
    },
    [](double u) -> double {
       return 1.5625 * (0.833333333333333 * u - 0.5) * (u - 1.0) * (u + 1.0) * (1.25 * u - 0.24999999999999994449) *
              (2.5 * u + 0.49999999999999988898);
    },
    [](double u) -> double {
       return -1.0416666666666667407 * (u - 1.0) * (u + 1.0) * (1.25 * u - 0.75000000000000011102) * (2.5 * u + 1.5) *
              (2.5 * u - 0.5);
    },
    [](double u) -> double {
       return 1.0416666666666667407 * (u - 1.0) * (u + 1.0) * (1.25 * u + 0.75) * (2.5 * u - 1.499999999999999778) *
              (2.5 * u + 0.5);
    },
    [](double u) -> double {
       return -1.562500000000000222 * (0.833333333333333 * u + 0.49999999999999994449) * (u - 1.0) * (u + 1.0) *
              (1.25 * u + 0.24999999999999994449) * (2.5 * u - 0.49999999999999972244);
    },
    [](double u) -> double {
       return 0.5 * (0.625 * u + 0.375) * (0.833333333333333 * u + 0.16666666666666662966) * (u + 1.0) *
              (1.25 * u - 0.24999999999999994449) * (2.5 * u - 1.5000000000000004441);
    },
    [](double u) -> double {
       return 0.5 * u * (0.6 * u - 0.4000000000000000222) * (0.75 * u - 0.24999999999999994449) * (u - 1.0) * (1.5 * u + 0.5) *
              (3.0 * u + 2.0000000000000008882);
    },
    [](double u) -> double {
       return -2.7000000000000006217 * u * (0.75 * u - 0.5) * (u - 1.0) * (1.0 * u - 0.33333333333333325932) * (u + 1.0) *
              (3.0 * u + 1.0);
    },
    [](double u) -> double {
       return 3.3749999999999986677 * u * (u - 1.0) * (1.0 * u - 0.66666666666666674068) * (u + 1.0) *
              (1.5 * u - 0.49999999999999988898) * (3.0 * u + 2.0);
    },
    [](double u) -> double {
       return -1.0 * (u - 1.0) * (u + 1.0) * (1.5 * u - 1.0) * (1.5 * u + 1.0) * (3.0 * u + 1.0) * (3.0 * u - 1.0);
    },
    [](double u) -> double {
       return 3.3750000000000008882 * u * (u - 1.0) * (1.0 * u + 0.66666666666666674068) * (u + 1.0) * (1.5 * u + 0.5) *
              (3.0 * u - 1.9999999999999993339);
    },
    [](double u) -> double {
       return -2.7000000000000001776 * u * (0.75 * u + 0.5) * (u - 1.0) * (1.0 * u + 0.33333333333333337034) * (u + 1.0) *
              (3.0 * u - 0.99999999999999933387);
    },
    [](double u) -> double {
       return 0.5 * u * (0.6 * u + 0.4000000000000000222) * (0.75 * u + 0.25) * (u + 1.0) * (1.5 * u - 0.49999999999999983347) *
              (3.0 * u - 2.0000000000000008882);
    },
    [](double u) -> double {
       return -1.27657335069444 * u * u * u * u * u * u * u + 1.27657335069444 * u * u * u * u * u * u +
              0.911838107638889 * u * u * u * u * u - 0.911838107638889 * u * u * u * u - 0.137706163194445 * u * u * u +
              0.137706163194444 * u * u + 0.00244140625000001 * u - 0.00244140625;
    },
    [](double u) -> double {
       return 2.0416666666666665186 * (0.7 * u - 0.5) * (0.875 * u - 0.375) * (u - 1.0) * (u + 1.0) *
              (1.16666666666667 * u - 0.1666666666666666019) * (1.75 * u + 0.25000000000000011102) * (3.5 * u + 1.5);
    },
    [](double u) -> double {
       return -1.2249999999999998668 * (0.875 * u - 0.62499999999999988898) * (u - 1.0) * (u + 1.0) *
              (1.16666666666667 * u - 0.5) * (1.75 * u - 0.24999999999999988898) * (3.5 * u + 0.50000000000000022204) *
              (3.5 * u + 2.5);
    },
    [](double u) -> double {
       return 1.0208333333333334814 * (u - 1.0) * (u + 1.0) * (1.16666666666667 * u - 0.83333333333333325932) *
              (1.75 * u - 0.75) * (1.75 * u + 1.25) * (3.5 * u - 0.49999999999999977796) * (3.5 * u + 1.5);
    },
    [](double u) -> double {
       return -1.0208333333333332593 * (u - 1.0) * (u + 1.0) * (1.16666666666667 * u + 0.83333333333333337034) *
              (1.75 * u - 1.249999999999999778) * (1.75 * u + 0.75) * (3.5 * u - 1.4999999999999995559) *
              (3.5 * u + 0.50000000000000022204);
    },
    [](double u) -> double {
       return 1.2249999999999998668 * (0.875 * u + 0.625) * (u - 1.0) * (u + 1.0) * (1.16666666666667 * u + 0.5) *
              (1.75 * u + 0.25000000000000005551) * (3.5 * u - 0.49999999999999961142) * (3.5 * u - 2.5000000000000008882);
    },
    [](double u) -> double {
       return -2.0416666666666660745 * (0.7 * u + 0.50000000000000011102) * (0.875 * u + 0.375) * (u - 1.0) * (u + 1.0) *
              (1.16666666666667 * u + 0.16666666666666674068) * (1.75 * u - 0.24999999999999988898) *
              (3.5 * u - 1.5000000000000008882);
    },
    [](double u) -> double {
       return 0.5 * (0.583333333333333 * u + 0.41666666666666662966) * (0.7 * u + 0.2999999999999999889) *
              (0.875 * u + 0.12500000000000005551) * (u + 1.0) * (1.16666666666667 * u - 0.16666666666666657415) *
              (1.75 * u - 0.75) * (3.5 * u - 2.4999999999999986677);
    },
    [](double u) -> double {
       return 0.5 * u * (0.571428571428571 * u - 0.42857142857142854764) * (0.666666666666667 * u - 0.33333333333333331483) *
              (0.8 * u - 0.2000000000000000111) * (u - 1.0) * (1.33333333333333 * u + 0.33333333333333331483) * (2.0 * u + 1.0) *
              (4.0 * u + 3.0);
    },
    [](double u) -> double {
       return -3.0476190476190474499 * u * (0.666666666666667 * u - 0.5) * (0.8 * u - 0.4000000000000000222) * (u - 1.0) *
              (1.0 * u - 0.25) * (u + 1.0) * (2.0 * u + 0.5) * (4.0 * u + 2.0);
    },
    [](double u) -> double {
       return 2.6666666666666665186 * u * (0.8 * u - 0.60000000000000008882) * (u - 1.0) * (1.0 * u - 0.5) * (u + 1.0) *
              (1.33333333333333 * u - 0.33333333333333331483) * (4.0 * u + 1.0) * (4.0 * u + 3.0);
    },
    [](double u) -> double {
       return -4.2666666666666666075 * u * (u - 1.0) * (1.0 * u - 0.75) * (u + 1.0) *
              (1.33333333333333 * u - 0.66666666666666662966) * (2.0 * u - 0.5) * (2.0 * u + 1.5) * (4.0 * u + 2.0);
    },
    [](double u) -> double {
       return 113.777777777778 * u * u * u * u * u * u * u * u - 213.333333333333 * u * u * u * u * u * u -
              1.4210854715202e-14 * u * u * u * u * u + 121.333333333333 * u * u * u * u + 7.105427357601e-15 * u * u * u -
              22.7777777777778 * u * u + 1.0;
    },
    [](double u) -> double {
       return -4.2666666666666666075 * u * (u - 1.0) * (1.0 * u + 0.75) * (u + 1.0) *
              (1.33333333333333 * u + 0.66666666666666662966) * (2.0 * u - 1.5) * (2.0 * u + 0.5) * (4.0 * u - 2.0);
    },
    [](double u) -> double {
       return 2.6666666666666665186 * u * (0.8 * u + 0.60000000000000008882) * (u - 1.0) * (1.0 * u + 0.5) * (u + 1.0) *
              (1.33333333333333 * u + 0.33333333333333331483) * (4.0 * u - 3.0) * (4.0 * u - 1.0);
    },
    [](double u) -> double {
       return -3.0476190476190474499 * u * (0.666666666666667 * u + 0.5) * (0.8 * u + 0.4000000000000000222) * (u - 1.0) *
              (1.0 * u + 0.25) * (u + 1.0) * (2.0 * u - 0.5) * (4.0 * u - 2.0);
    },
    [](double u) -> double {
       return 0.5 * u * (0.571428571428571 * u + 0.42857142857142854764) * (0.666666666666667 * u + 0.33333333333333331483) *
              (0.8 * u + 0.2000000000000000111) * (u + 1.0) * (1.33333333333333 * u - 0.33333333333333331483) * (2.0 * u - 1.0) *
              (4.0 * u - 3.0);
    },
    [](double u) -> double {
       return -2.08520900181362 * u * u * u * u * u * u * u * u * u + 2.08520900181362 * u * u * u * u * u * u * u * u +
              2.16243896484375 * u * u * u * u * u * u * u - 2.16243896484375 * u * u * u * u * u * u -
              0.627374267578125 * u * u * u * u * u + 0.627374267578125 * u * u * u * u + 0.0506783621651786 * u * u * u -
              0.0506783621651786 * u * u - 0.0005340576171875 * u + 0.00053405761718750032526;
    },
    [](double u) -> double {
       return 2.53125 * (0.642857142857143 * u - 0.5) * (0.75 * u - 0.41666666666666662966) * (0.9 * u - 0.29999999999999987788) *
              (u - 1.0) * (u + 1.0) * (1.125 * u - 0.12500000000000005551) * (1.5 * u + 0.16666666666666674068) *
              (2.25 * u + 0.75000000000000011102) * (4.5 * u + 2.5);
    },
    [](double u) -> double {
       return -1.4464285714285711748 * (0.75 * u - 0.58333333333333325932) * (0.9 * u - 0.5) * (u - 1.0) * (u + 1.0) *
              (1.125 * u - 0.37499999999999988898) * (1.5 * u - 0.16666666666666671293) * (2.25 * u + 0.25000000000000011102) *
              (4.5 * u + 1.500000000000000222) * (4.5 * u + 3.5);
    },
    [](double u) -> double {
       return 1.124999999999999778 * (0.9 * u - 0.69999999999999984457) * (u - 1.0) * (u + 1.0) * (1.125 * u - 0.625) *
              (1.5 * u - 0.49999999999999988898) * (2.25 * u - 0.25000000000000005551) * (2.25 * u + 1.75) *
              (4.5 * u + 0.50000000000000022204) * (4.5 * u + 2.5);
    },
    [](double u) -> double {
       return -1.0124999999999999556 * (u - 1.0) * (u + 1.0) * (1.125 * u - 0.87499999999999988898) *
              (1.5 * u - 0.83333333333333325932) * (1.5 * u + 1.1666666666666667407) * (2.25 * u - 0.74999999999999977796) *
              (2.25 * u + 1.25) * (4.5 * u - 0.5) * (4.5 * u + 1.500000000000000222);
    },
    [](double u) -> double {
       return 1.0124999999999999556 * (u - 1.0) * (u + 1.0) * (1.125 * u + 0.875) * (1.5 * u + 0.83333333333333325932) *
              (1.5 * u - 1.1666666666666669627) * (2.25 * u + 0.74999999999999988898) * (2.25 * u - 1.25) * (4.5 * u + 0.5) *
              (4.5 * u - 1.5000000000000004441);
    },
    [](double u) -> double {
       return -1.124999999999999778 * (0.9 * u + 0.69999999999999995559) * (u - 1.0) * (u + 1.0) * (1.125 * u + 0.625) *
              (1.5 * u + 0.5) * (2.25 * u - 1.749999999999999778) * (2.25 * u + 0.25000000000000011102) *
              (4.5 * u - 2.4999999999999991118) * (4.5 * u - 0.50000000000000055511);
    },
    [](double u) -> double {
       return 1.4464285714285711748 * (0.75 * u + 0.58333333333333325932) * (0.9 * u + 0.5) * (u - 1.0) * (u + 1.0) *
              (1.125 * u + 0.37500000000000005551) * (1.5 * u + 0.16666666666666671293) * (2.25 * u - 0.25000000000000011102) *
              (4.5 * u - 1.4999999999999991118) * (4.5 * u - 3.5000000000000017764);
    },
    [](double u) -> double {
       return -2.5312499999999991118 * (0.642857142857143 * u + 0.5) * (0.75 * u + 0.41666666666666668517) *
              (0.9 * u + 0.2999999999999999889) * (u - 1.0) * (u + 1.0) * (1.125 * u + 0.12500000000000005551) *
              (1.5 * u - 0.16666666666666679619) * (2.25 * u - 0.74999999999999977796) * (4.5 * u - 2.5000000000000017764);
    },
    [](double u) -> double {
       return 0.5 * (0.5625 * u + 0.4375) * (0.642857142857143 * u + 0.35714285714285715079) * (0.75 * u + 0.25) *
              (0.9 * u + 0.10000000000000003331) * (u + 1.0) * (1.125 * u - 0.12500000000000005551) *
              (1.5 * u - 0.49999999999999983347) * (2.25 * u - 1.25) * (4.5 * u - 3.4999999999999982236);
    },
    [](double u) -> double {
       return 0.5 * u * (0.555555555555556 * u - 0.44444444444444447528) * (0.625 * u - 0.37500000000000005551) *
              (0.714285714285714 * u - 0.28571428571428564291) * (0.833333333333333 * u - 0.16666666666666662966) * (u - 1.0) *
              (1.25 * u + 0.24999999999999994449) * (1.66666666666667 * u + 0.66666666666666674068) * (2.5 * u + 1.5) *
              (5.0 * u + 4.0000000000000008882);
    },
    [](double u) -> double {
       return -3.472222222222222765 * u * (0.625 * u - 0.5) * (0.714285714285714 * u - 0.42857142857142860315) *
              (0.833333333333333 * u - 0.33333333333333325932) * (u - 1.0) * (1.0 * u - 0.19999999999999995559) * (u + 1.0) *
              (1.66666666666667 * u + 0.33333333333333320381) * (2.5 * u + 1.0) * (5.0 * u + 2.9999999999999986677);
    },
    [](double u) -> double {
       return 2.6041666666666669627 * u * (0.714285714285714 * u - 0.57142857142857150787) * (0.833333333333333 * u - 0.5) *
              (u - 1.0) * (u + 1.0) * (1.0 * u - 0.4000000000000000222) * (1.25 * u - 0.24999999999999994449) *
              (2.5 * u + 0.49999999999999988898) * (5.0 * u + 3.9999999999999986677) * (5.0 * u + 2.0000000000000004441);
    },
    [](double u) -> double {
       return -2.9761904761904762751 * u * (0.833333333333333 * u - 0.66666666666666662966) * (u - 1.0) *
              (1.0 * u - 0.60000000000000008882) * (u + 1.0) * (1.25 * u - 0.49999999999999988898) *
              (1.66666666666667 * u - 0.33333333333333325932) * (2.5 * u + 2.0) * (5.0 * u + 0.99999999999999944489) *
              (5.0 * u + 3.0000000000000004441);
    },
    [](double u) -> double {
       return 5.2083333333333339255 * u * (u - 1.0) * (1.0 * u - 0.80000000000000004441) * (u + 1.0) *
              (1.25 * u - 0.75000000000000011102) * (1.66666666666667 * u + 1.3333333333333332593) *
              (1.66666666666667 * u - 0.66666666666666662966) * (2.5 * u + 1.5) * (2.5 * u - 0.5) *
              (5.0 * u + 1.9999999999999993339);
    },
    [](double u) -> double {
       return -678.168402777778 * u * u * u * u * u * u * u * u * u * u -
              1.13686837721616e-13 * u * u * u * u * u * u * u * u * u + 1491.97048611111 * u * u * u * u * u * u * u * u -
              1110.02604166667 * u * u * u * u * u * u + 2.8421709430404e-14 * u * u * u * u * u +
              331.814236111111 * u * u * u * u - 36.5902777777778 * u * u - 1.77635683940025e-15 * u + 1.0;
    },
    [](double u) -> double {
       return 5.2083333333333348136 * u * (u - 1.0) * (1.0 * u + 0.80000000000000004441) * (u + 1.0) * (1.25 * u + 0.75) *
              (1.66666666666667 * u - 1.3333333333333332593) * (1.66666666666667 * u + 0.66666666666666674068) *
              (2.5 * u - 1.499999999999999778) * (2.5 * u + 0.5) * (5.0 * u - 2.0);
    },
    [](double u) -> double {
       return -2.9761904761904762751 * u * (0.833333333333333 * u + 0.66666666666666674068) * (u - 1.0) * (u + 1.0) *
              (1.0 * u + 0.60000000000000008882) * (1.25 * u + 0.5) * (1.66666666666667 * u + 0.33333333333333331483) *
              (2.5 * u - 1.9999999999999993339) * (5.0 * u - 2.9999999999999977796) * (5.0 * u - 1.0);
    },
    [](double u) -> double {
       return 2.6041666666666669627 * u * (0.714285714285714 * u + 0.57142857142857139685) *
              (0.833333333333333 * u + 0.49999999999999994449) * (u - 1.0) * (1.0 * u + 0.4000000000000000222) * (u + 1.0) *
              (1.25 * u + 0.24999999999999994449) * (2.5 * u - 0.49999999999999972244) * (5.0 * u - 1.9999999999999977796) *
              (5.0 * u - 4.0000000000000008882);
    },
    [](double u) -> double {
       return -3.472222222222222765 * u * (0.625 * u + 0.5) * (0.714285714285714 * u + 0.42857142857142854764) *
              (0.833333333333333 * u + 0.33333333333333331483) * (u - 1.0) * (1.0 * u + 0.19999999999999995559) * (u + 1.0) *
              (1.66666666666667 * u - 0.33333333333333320381) * (2.5 * u - 0.99999999999999944489) *
              (5.0 * u - 3.0000000000000008882);
    },
    [](double u) -> double {
       return 0.5 * u * (0.555555555555556 * u + 0.44444444444444447528) * (0.625 * u + 0.375) *
              (0.714285714285714 * u + 0.28571428571428575394) * (0.833333333333333 * u + 0.16666666666666662966) * (u + 1.0) *
              (1.25 * u - 0.24999999999999994449) * (1.66666666666667 * u - 0.66666666666666640761) *
              (2.5 * u - 1.5000000000000004441) * (5.0 * u - 4.0000000000000008882);
    },
    [](double u) -> double {
       return -3.49006782020421 * u * u * u * u * u * u * u * u * u * u * u +
              3.49006782020421 * u * u * u * u * u * u * u * u * u * u + 4.75918339118756 * u * u * u * u * u * u * u * u * u -
              4.75918339118756 * u * u * u * u * u * u * u * u - 2.09246740835684 * u * u * u * u * u * u * u +
              2.09246740835684 * u * u * u * u * u * u + 0.340444737725367 * u * u * u * u * u -
              0.340444737725367 * u * u * u * u - 0.0172130633157397 * u * u * u + 0.0172130633157397 * u * u +
              0.000120162963867189 * u - 0.00012016296386718741868;
    },
    [](double u) -> double {
       return 3.0249999999999994671 * (0.611111111111111 * u - 0.5) * (0.6875 * u - 0.43750000000000005551) *
              (0.785714285714286 * u - 0.35714285714285715079) * (0.916666666666667 * u - 0.25) * (u - 1.0) * (u + 1.0) *
              (1.1 * u - 0.099999999999999922284) * (1.375 * u + 0.12500000000000005551) * (1.83333333333333 * u + 0.5) *
              (2.75 * u + 1.2500000000000004441) * (5.5 * u + 3.5000000000000008882);
    },
    [](double u) -> double {
       return -1.6805555555555558023 * (0.6875 * u - 0.5625) * (0.785714285714286 * u - 0.5) *
              (0.916666666666667 * u - 0.41666666666666674068) * (u - 1.0) * (u + 1.0) * (1.1 * u - 0.2999999999999999889) *
              (1.375 * u - 0.12499999999999990286) * (1.83333333333333 * u + 0.16666666666666674068) * (2.75 * u + 0.75) *
              (5.5 * u + 2.5000000000000008882) * (5.5 * u + 4.5000000000000008882);
    },
    [](double u) -> double {
       return 1.2604166666666667407 * (0.785714285714286 * u - 0.6428571428571427937) *
              (0.916666666666667 * u - 0.58333333333333337034) * (u - 1.0) * (u + 1.0) * (1.1 * u - 0.5) * (1.375 * u - 0.375) *
              (1.83333333333333 * u - 0.16666666666666651864) * (2.75 * u + 0.25000000000000011102) *
              (2.75 * u + 2.2500000000000004441) * (5.5 * u + 1.4999999999999993339) * (5.5 * u + 3.5000000000000008882);
    },
    [](double u) -> double {
       return -1.0803571428571427937 * (0.916666666666667 * u - 0.75) * (u - 1.0) * (u + 1.0) *
              (1.1 * u - 0.70000000000000006661) * (1.375 * u - 0.625) * (1.83333333333333 * u - 0.5) *
              (1.83333333333333 * u + 1.5) * (2.75 * u + 1.75) * (2.75 * u - 0.24999999999999986122) *
              (5.5 * u + 2.4999999999999995559) * (5.5 * u + 0.50000000000000033307);
    },
    [](double u) -> double {
       return 1.0083333333333335258 * (u - 1.0) * (u + 1.0) * (1.1 * u - 0.9000000000000000222) * (1.375 * u - 0.875) *
              (1.375 * u + 1.125) * (1.83333333333333 * u - 0.83333333333333325932) *
              (1.83333333333333 * u + 1.1666666666666667407) * (2.75 * u - 0.75) * (2.75 * u + 1.25) *
              (5.5 * u - 0.49999999999999972244) * (5.5 * u + 1.5000000000000004441);
    },
    [](double u) -> double {
       return -1.0083333333333333037 * (u - 1.0) * (u + 1.0) * (1.1 * u + 0.9000000000000000222) * (1.375 * u - 1.125) *
              (1.375 * u + 0.87500000000000011102) * (1.83333333333333 * u - 1.1666666666666665186) *
              (1.83333333333333 * u + 0.83333333333333348136) * (2.75 * u - 1.249999999999999778) *
              (2.75 * u + 0.75000000000000022204) * (5.5 * u - 1.4999999999999993339) * (5.5 * u + 0.50000000000000033307);
    },
    [](double u) -> double {
       return 1.0803571428571427937 * (0.916666666666667 * u + 0.75) * (u - 1.0) * (u + 1.0) *
              (1.1 * u + 0.70000000000000006661) * (1.375 * u + 0.625) * (1.83333333333333 * u - 1.5) *
              (1.83333333333333 * u + 0.5) * (2.75 * u - 1.749999999999999778) * (2.75 * u + 0.25000000000000011102) *
              (5.5 * u - 2.4999999999999995559) * (5.5 * u - 0.49999999999999938938);
    },
    [](double u) -> double {
       return -1.2604166666666667407 * (0.785714285714286 * u + 0.6428571428571427937) *
              (0.916666666666667 * u + 0.58333333333333337034) * (u - 1.0) * (u + 1.0) * (1.1 * u + 0.5) * (1.375 * u + 0.375) *
              (1.83333333333333 * u + 0.16666666666666668517) * (2.75 * u - 0.24999999999999969469) *
              (2.75 * u - 2.2500000000000004441) * (5.5 * u - 3.4999999999999995559) * (5.5 * u - 1.4999999999999993339);
    },
    [](double u) -> double {
       return 1.6805555555555560243 * (0.6875 * u + 0.5625) * (0.785714285714286 * u + 0.49999999999999988898) *
              (0.916666666666667 * u + 0.41666666666666662966) * (u - 1.0) * (u + 1.0) * (1.1 * u + 0.29999999999999993339) *
              (1.375 * u + 0.12500000000000002776) * (1.83333333333333 * u - 0.16666666666666646313) *
              (2.75 * u - 0.74999999999999966693) * (5.5 * u - 2.4999999999999995559) *
              (5.50000000000001 * u - 4.5000000000000044409);
    },
    [](double u) -> double {
       return -3.0249999999999994671 * (0.611111111111111 * u + 0.5) * (0.6875 * u + 0.4375) *
              (0.785714285714286 * u + 0.35714285714285715079) * (0.916666666666667 * u + 0.25) * (u - 1.0) * (u + 1.0) *
              (1.1 * u + 0.10000000000000004718) * (1.375 * u - 0.12499999999999988898) * (1.83333333333333 * u - 0.5) *
              (2.75 * u - 1.2500000000000004441) * (5.50000000000001 * u - 3.5000000000000039968);
    },
    [](double u) -> double {
       return 0.5 * (0.55 * u + 0.4500000000000000111) * (0.611111111111111 * u + 0.38888888888888889506) *
              (0.6875 * u + 0.3125) * (0.785714285714286 * u + 0.21428571428571427382) *
              (0.916666666666667 * u + 0.083333333333333370341) * (u + 1.0) * (1.1 * u - 0.099999999999999894529) *
              (1.375 * u - 0.375) * (1.83333333333333 * u - 0.83333333333333348136) * (2.75 * u - 1.7500000000000008882) *
              (5.5 * u - 4.4999999999999982236);
    },
    [](double u) -> double {
       return 0.5 * u * (0.545454545454546 * u - 0.45454545454545458583) * (0.6 * u - 0.4000000000000000222) *
              (0.666666666666667 * u - 0.33333333333333331483) * (0.75 * u - 0.24999999999999994449) *
              (0.857142857142857 * u - 0.14285714285714290472) * (u - 1.0) * (1.2 * u + 0.19999999999999995559) *
              (1.5 * u + 0.5) * (2.0 * u + 1.0) * (3.0 * u + 2.0000000000000008882) * (6.0 * u + 5.0000000000000017764);
    },
    [](double u) -> double {
       return -3.9272727272727281367 * u * (0.6 * u - 0.5) * (0.666666666666667 * u - 0.44444444444444447528) *
              (0.75 * u - 0.37499999999999994449) * (0.857142857142857 * u - 0.28571428571428569843) * (u - 1.0) *
              (1.0 * u - 0.16666666666666674068) * (u + 1.0) * (1.5 * u + 0.24999999999999991673) *
              (2.0 * u + 0.66666666666666674068) * (3.0 * u + 1.499999999999999778) * (6.0 * u + 4.0000000000000017764);
    },
    [](double u) -> double {
       return 2.7000000000000006217 * u * (0.666666666666667 * u - 0.5555555555555554692) * (0.75 * u - 0.5) *
              (0.857142857142857 * u - 0.42857142857142854764) * (u - 1.0) * (1.0 * u - 0.33333333333333325932) * (u + 1.0) *
              (1.2 * u - 0.20000000000000003886) * (2.0 * u + 0.33333333333333320381) * (3.0 * u + 1.0) *
              (6.0 * u + 2.9999999999999986677) * (6.0 * u + 5.0000000000000017764);
    },
    [](double u) -> double {
       return -2.6666666666666665186 * u * (0.75 * u - 0.625) * (0.857142857142857 * u - 0.57142857142857150787) * (u - 1.0) *
              (1.0 * u - 0.5) * (u + 1.0) * (1.2 * u - 0.39999999999999996669) * (1.5 * u - 0.25000000000000005551) *
              (3.0 * u + 0.49999999999999983347) * (3.0 * u + 2.4999999999999995559) * (6.0 * u + 3.9999999999999986677) *
              (6.0 * u + 2.0000000000000008882);
    },
    [](double u) -> double {
       return 3.3749999999999986677 * u * (0.857142857142857 * u - 0.71428571428571430157) * (u - 1.0) *
              (1.0 * u - 0.66666666666666674068) * (u + 1.0) * (1.2 * u - 0.5999999999999999778) *
              (1.5 * u - 0.49999999999999988898) * (2.0 * u - 0.33333333333333342585) * (2.0 * u + 1.6666666666666667407) *
              (3.0 * u + 2.0) * (6.0 * u + 0.99999999999999933387) * (6.0 * u + 3.0000000000000008882);
    },
    [](double u) -> double {
       return -6.1714285714285734841 * u * (u - 1.0) * (u + 1.0) * (1.0 * u - 0.83333333333333348136) *
              (1.2 * u - 0.80000000000000004441) * (1.5 * u + 1.249999999999999778) * (1.5 * u - 0.75) *
              (2.0 * u + 1.3333333333333332593) * (2.0 * u - 0.66666666666666662966) * (3.0 * u - 0.50000000000000011102) *
              (3.0 * u + 1.499999999999999778) * (6.0 * u + 1.9999999999999993339);
    },
    [](double u) -> double {
       return 1.0 * (u - 1.0) * (u + 1.0) * (1.2 * u + 1.0) * (1.2 * u - 1.0) * (1.5 * u - 1.0) * (1.5 * u + 1.0) *
              (2.0 * u - 1.0) * (2.0 * u + 1.0) * (3.0 * u + 1.0) * (3.0 * u - 1.0) * (6.0 * u - 1.0) * (6.0 * u + 1.0);
    },
    [](double u) -> double {
       return -6.1714285714285699314 * u * (u - 1.0) * (1.0 * u + 0.83333333333333337034) * (u + 1.0) *
              (1.2 * u + 0.79999999999999993339) * (1.5 * u + 0.74999999999999988898) * (1.5 * u - 1.250000000000000222) *
              (2.0 * u + 0.66666666666666662966) * (2.0 * u - 1.3333333333333334814) * (3.0 * u + 0.49999999999999983347) *
              (3.0 * u - 1.5000000000000004441) * (6.00000000000001 * u - 2.0000000000000013323);
    },
    [](double u) -> double {
       return 3.3750000000000008882 * u * (0.857142857142857 * u + 0.7142857142857144126) * (u - 1.0) *
              (1.0 * u + 0.66666666666666674068) * (u + 1.0) * (1.2 * u + 0.60000000000000008882) * (1.5 * u + 0.5) *
              (2.0 * u - 1.6666666666666665186) * (2.0 * u + 0.33333333333333331483) * (3.0 * u - 1.9999999999999993339) *
              (6.0 * u - 2.9999999999999986677) * (6.00000000000001 * u - 1.0000000000000013323);
    },
    [](double u) -> double {
       return -2.6666666666666665186 * u * (0.75 * u + 0.62499999999999988898) *
              (0.857142857142857 * u + 0.57142857142857150787) * (u - 1.0) * (1.0 * u + 0.5) * (u + 1.0) *
              (1.2 * u + 0.4000000000000000222) * (1.5 * u + 0.24999999999999994449) * (3.0 * u - 2.5000000000000004441) *
              (3.0 * u - 0.50000000000000033307) * (6.0 * u - 3.9999999999999986677) * (6.0 * u - 1.9999999999999986677);
    },
    [](double u) -> double {
       return 2.7000000000000001776 * u * (0.666666666666667 * u + 0.55555555555555558023) * (0.75 * u + 0.5) *
              (0.857142857142857 * u + 0.42857142857142854764) * (u - 1.0) * (1.0 * u + 0.33333333333333337034) * (u + 1.0) *
              (1.2 * u + 0.19999999999999995559) * (2.0 * u - 0.33333333333333348136) * (3.0 * u - 0.99999999999999933387) *
              (6.0 * u - 2.9999999999999986677) * (6.00000000000001 * u - 5.0000000000000035527);
    },
    [](double u) -> double {
       return -3.9272727272727268044 * u * (0.6 * u + 0.50000000000000011102) * (0.666666666666667 * u + 0.44444444444444447528) *
              (0.75 * u + 0.375) * (0.857142857142857 * u + 0.28571428571428575394) * (u - 1.0) * (u + 1.0) *
              (1.0 * u + 0.16666666666666665741) * (1.5 * u - 0.25000000000000016653) * (2.0 * u - 0.66666666666666651864) *
              (3.0 * u - 1.5000000000000004441) * (6.00000000000001 * u - 4.0000000000000044409);
    },
    [](double u) -> double {
       return 0.5 * u * (0.545454545454545 * u + 0.45454545454545453032) * (0.6 * u + 0.4000000000000000222) *
              (0.666666666666667 * u + 0.33333333333333331483) * (0.75 * u + 0.25) *
              (0.857142857142857 * u + 0.14285714285714284921) * (u + 1.0) * (1.2 * u - 0.20000000000000012212) *
              (1.5 * u - 0.49999999999999983347) * (2.0 * u - 1.0) * (3.0 * u - 2.0000000000000008882) *
              (6.0 * u - 4.9999999999999973355);
    },
    [](double u) -> double {
       return -5.93735914776246 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              5.93735914776246 * u * u * u * u * u * u * u * u * u * u * u * u +
              10.0478385577519 * u * u * u * u * u * u * u * u * u * u * u -
              10.0478385577519 * u * u * u * u * u * u * u * u * u * u - 5.97519393523114 * u * u * u * u * u * u * u * u * u +
              5.97519393523114 * u * u * u * u * u * u * u * u + 1.51908430700509 * u * u * u * u * u * u * u -
              1.5190843070051 * u * u * u * u * u * u - 0.159890189950148 * u * u * u * u * u +
              0.159890189950148 * u * u * u * u + 0.00554794553267719 * u * u * u - 0.00554794553267721 * u * u -
              2.75373458862304e-5 * u + 0.000027537345886230424704;
    },
    [](double u) -> double {
       return 3.5208333333333330373 * (0.590909090909091 * u - 0.5) * (0.65 * u - 0.4500000000000000111) *
              (0.722222222222222 * u - 0.38888888888888895057) * (0.8125 * u - 0.3125) *
              (0.928571428571428 * u - 0.21428571428571430157) * (u - 1.0) * (u + 1.0) *
              (1.08333333333333 * u - 0.083333333333333287074) * (1.3 * u + 0.099999999999999922284) *
              (1.625 * u + 0.37499999999999994449) * (2.16666666666667 * u + 0.8333333333333331483) * (3.25 * u + 1.75) *
              (6.5 * u + 4.5);
    },
    [](double u) -> double {
       return -1.9204545454545456362 * (0.65 * u - 0.55000000000000004441) * (0.722222222222222 * u - 0.5) *
              (0.8125 * u - 0.43750000000000005551) * (0.928571428571429 * u - 0.35714285714285715079) * (u - 1.0) * (u + 1.0) *
              (1.08333333333333 * u - 0.25000000000000005551) * (1.3 * u - 0.099999999999999936162) *
              (1.625 * u + 0.12499999999999991673) * (2.16666666666667 * u + 0.49999999999999988898) * (3.25 * u + 1.25) *
              (6.5 * u + 3.5) * (6.5 * u + 5.5);
    },
    [](double u) -> double {
       return 1.4083333333333332149 * (0.722222222222222 * u - 0.61111111111111116045) * (0.8125 * u - 0.5625) *
              (0.928571428571428 * u - 0.5) * (u - 1.0) * (u + 1.0) * (1.08333333333333 * u - 0.41666666666666668517) *
              (1.3 * u - 0.30000000000000004441) * (1.625 * u - 0.12499999999999993061) *
              (2.16666666666667 * u + 0.16666666666666654639) * (3.25 * u + 0.74999999999999988898) * (3.25 * u + 2.75) *
              (6.5 * u + 2.5) * (6.5 * u + 4.5);
    },
    [](double u) -> double {
       return -1.1736111111111111605 * (0.8125 * u - 0.68750000000000011102) * (0.928571428571429 * u - 0.64285714285714290472) *
              (u - 1.0) * (u + 1.0) * (1.08333333333333 * u - 0.58333333333333337034) * (1.3 * u - 0.5) *
              (1.625 * u - 0.37500000000000011102) * (2.16666666666667 * u + 1.8333333333333332593) *
              (2.16666666666667 * u - 0.16666666666666657415) * (3.25 * u + 0.24999999999999983347) * (3.25 * u + 2.25) *
              (6.5 * u + 1.499999999999999778) * (6.5 * u + 3.5);
    },
    [](double u) -> double {
       return 1.0562499999999999112 * (0.928571428571428 * u - 0.78571428571428569843) * (u - 1.0) * (u + 1.0) *
              (1.08333333333333 * u - 0.75000000000000011102) * (1.3 * u - 0.70000000000000006661) * (1.625 * u + 1.375) *
              (1.625 * u - 0.625) * (2.16666666666667 * u - 0.50000000000000011102) *
              (2.16666666666667 * u + 1.499999999999999778) * (3.25 * u + 1.75) * (3.25 * u - 0.24999999999999988898) *
              (6.5 * u + 0.49999999999999966693) * (6.5 * u + 2.5);
    },
    [](double u) -> double {
       return -1.0059523809523809312 * (u - 1.0) * (u + 1.0) * (1.08333333333333 * u - 0.91666666666666674068) *
              (1.3 * u + 1.0999999999999998668) * (1.3 * u - 0.9000000000000000222) * (1.625 * u - 0.87500000000000011102) *
              (1.625 * u + 1.125) * (2.16666666666667 * u + 1.1666666666666665186) *
              (2.16666666666667 * u - 0.83333333333333337034) * (3.25 * u - 0.75000000000000022204) * (3.25 * u + 1.25) *
              (6.5 * u + 1.499999999999999778) * (6.5 * u - 0.5);
    },
    [](double u) -> double {
       return 1.0059523809523809312 * (u - 1.0) * (u + 1.0) * (1.08333333333333 * u + 0.91666666666666674068) *
              (1.3 * u - 1.1000000000000000888) * (1.3 * u + 0.9000000000000000222) * (1.625 * u - 1.125) *
              (1.625 * u + 0.87500000000000011102) * (2.16666666666667 * u - 1.1666666666666665186) *
              (2.16666666666667 * u + 0.83333333333333337034) * (3.25 * u - 1.25) * (3.25 * u + 0.75000000000000011102) *
              (6.49999999999999 * u - 1.4999999999999991118) * (6.5 * u + 0.5);
    },
    [](double u) -> double {
       return -1.0562500000000001332 * (0.928571428571428 * u + 0.7857142857142855874) * (u - 1.0) * (u + 1.0) *
              (1.08333333333333 * u + 0.74999999999999988898) * (1.3 * u + 0.69999999999999984457) *
              (1.625 * u - 1.375000000000000222) * (1.625 * u + 0.625) * (2.16666666666667 * u + 0.49999999999999988898) *
              (2.16666666666667 * u - 1.500000000000000222) * (3.25 * u - 1.750000000000000222) *
              (3.25 * u + 0.24999999999999983347) * (6.49999999999999 * u - 0.49999999999999927836) *
              (6.5 * u - 2.5000000000000013323);
    },
    [](double u) -> double {
       return 1.1736111111111111605 * (0.8125 * u + 0.6875) * (0.928571428571429 * u + 0.64285714285714290472) * (u - 1.0) *
              (u + 1.0) * (1.08333333333333 * u + 0.58333333333333337034) * (1.3 * u + 0.5) * (1.625 * u + 0.375) *
              (2.16666666666667 * u - 1.8333333333333330373) * (2.16666666666667 * u + 0.16666666666666657415) *
              (3.25 * u - 2.25) * (3.25 * u - 0.24999999999999983347) * (6.49999999999999 * u - 3.4999999999999977796) *
              (6.5 * u - 1.5000000000000015543);
    },
    [](double u) -> double {
       return -1.408333333333333437 * (0.722222222222222 * u + 0.61111111111111104943) * (0.8125 * u + 0.5625) *
              (0.928571428571428 * u + 0.49999999999999988898) * (u - 1.0) * (u + 1.0) *
              (1.08333333333333 * u + 0.41666666666666657415) * (1.3 * u + 0.29999999999999987788) *
              (1.625 * u + 0.12499999999999991673) * (2.16666666666667 * u - 0.16666666666666651864) *
              (3.25 * u - 2.7500000000000004441) * (3.25 * u - 0.75000000000000022204) *
              (6.49999999999999 * u - 2.4999999999999977796) * (6.5 * u - 4.5000000000000026645);
    },
    [](double u) -> double {
       return 1.9204545454545456362 * (0.65 * u + 0.55000000000000004441) * (0.722222222222222 * u + 0.5) *
              (0.8125 * u + 0.4375) * (0.928571428571429 * u + 0.35714285714285715079) * (u - 1.0) * (u + 1.0) *
              (1.08333333333333 * u + 0.25) * (1.3 * u + 0.099999999999999936162) * (1.625 * u - 0.12499999999999991673) *
              (2.16666666666667 * u - 0.50000000000000022204) * (3.25 * u - 1.25) *
              (6.49999999999999 * u - 5.4999999999999964473) * (6.5 * u - 3.5000000000000031086);
    },
    [](double u) -> double {
       return -3.5208333333333357018 * (0.590909090909091 * u + 0.49999999999999994449) * (0.65 * u + 0.44999999999999990008) *
              (0.722222222222222 * u + 0.38888888888888883955) * (0.8125 * u + 0.3125) *
              (0.928571428571428 * u + 0.21428571428571419055) * (u - 1.0) * (u + 1.0) *
              (1.08333333333333 * u + 0.083333333333333273196) * (1.3 * u - 0.099999999999999922284) *
              (1.625 * u - 0.37500000000000011102) * (2.16666666666667 * u - 0.83333333333333303727) *
              (3.25 * u - 1.750000000000000222) * (6.49999999999999 * u - 4.4999999999999964473);
    },
    [](double u) -> double {
       return 0.5 * (0.541666666666667 * u + 0.45833333333333331483) * (0.590909090909091 * u + 0.40909090909090911614) *
              (0.65 * u + 0.3499999999999999778) * (0.722222222222222 * u + 0.2777777777777777346) *
              (0.8125 * u + 0.18749999999999997224) * (0.928571428571429 * u + 0.071428571428571382973) * (u + 1.0) *
              (1.08333333333333 * u - 0.083333333333333273196) * (1.3 * u - 0.30000000000000009992) * (1.625 * u - 0.625) *
              (2.16666666666667 * u - 1.1666666666666669627) * (3.25 * u - 2.25) * (6.5 * u - 5.5000000000000044409);
    },
    [](double u) -> double {
       return 0.5 * u * (0.538461538461538 * u - 0.46153846153846156408) * (0.583333333333333 * u - 0.41666666666666662966) *
              (0.636363636363636 * u - 0.36363636363636359095) * (0.7 * u - 0.2999999999999999889) *
              (0.777777777777778 * u - 0.2222222222222222654) * (0.875 * u - 0.12499999999999994449) * (u - 1.0) *
              (1.16666666666667 * u + 0.16666666666666674068) * (1.4 * u + 0.39999999999999996669) * (1.75 * u + 0.75) *
              (2.33333333333333 * u + 1.3333333333333330373) * (3.5 * u + 2.5) * (7.0 * u + 6.0000000000000035527);
    },
    [](double u) -> double {
       return -4.3974358974358986885 * u * (0.583333333333333 * u - 0.5) * (0.636363636363636 * u - 0.4545454545454544748) *
              (0.7 * u - 0.39999999999999996669) * (0.777777777777778 * u - 0.33333333333333331483) *
              (0.875 * u - 0.25000000000000005551) * (u - 1.0) * (1.0 * u - 0.1428571428571427937) * (u + 1.0) *
              (1.4 * u + 0.20000000000000006661) * (1.75 * u + 0.49999999999999988898) * (2.33333333333333 * u + 1.0) *
              (3.5 * u + 1.9999999999999991118) * (7.0 * u + 4.9999999999999982236);
    },
    [](double u) -> double {
       return 2.8583333333333325044 * u * (0.636363636363636 * u - 0.54545454545454541417) * (0.7 * u - 0.5) *
              (0.777777777777778 * u - 0.44444444444444447528) * (0.875 * u - 0.375) * (u - 1.0) *
              (1.0 * u - 0.28571428571428580945) * (u + 1.0) * (1.16666666666667 * u - 0.1666666666666666019) *
              (1.75 * u + 0.25000000000000011102) * (2.33333333333333 * u + 0.66666666666666651864) * (3.5 * u + 1.5) *
              (7.0 * u + 3.9999999999999982236) * (7.0 * u + 5.9999999999999982236);
    },
    [](double u) -> double {
       return -2.5984848484848481753 * u * (0.7 * u - 0.5999999999999999778) * (0.777777777777778 * u - 0.55555555555555558023) *
              (0.875 * u - 0.5) * (u - 1.0) * (1.0 * u - 0.42857142857142860315) * (u + 1.0) *
              (1.16666666666667 * u - 0.33333333333333342585) * (1.4 * u - 0.19999999999999992784) *
              (2.33333333333333 * u + 0.33333333333333353687) * (3.5 * u + 2.9999999999999991118) * (3.5 * u + 1.0) *
              (7.0 * u + 4.9999999999999982236) * (7.0 * u + 3.0000000000000017764);
    },
    [](double u) -> double {
       return 2.8583333333333329485 * u * (0.777777777777778 * u - 0.66666666666666662966) *
              (0.875 * u - 0.62499999999999988898) * (u - 1.0) * (1.0 * u - 0.57142857142857139685) * (u + 1.0) *
              (1.16666666666667 * u - 0.5) * (1.4 * u - 0.4000000000000000222) * (1.75 * u - 0.24999999999999988898) *
              (2.33333333333333 * u + 2.0) * (3.5 * u + 0.50000000000000022204) * (3.5 * u + 2.5) *
              (7.0 * u + 1.9999999999999991118) * (7.0 * u + 4.0000000000000017764);
    },
    [](double u) -> double {
       return -3.8111111111111113381 * u * (0.875 * u - 0.75) * (u - 1.0) * (u + 1.0) * (1.0 * u - 0.71428571428571430157) *
              (1.16666666666667 * u - 0.66666666666666662966) * (1.4 * u - 0.5999999999999999778) *
              (1.75 * u - 0.50000000000000011102) * (1.75 * u + 1.5) * (2.33333333333333 * u + 1.6666666666666665186) *
              (2.33333333333333 * u - 0.33333333333333325932) * (3.5 * u + 2.0) * (7.0 * u + 2.9999999999999991118) *
              (7.0 * u + 1.0000000000000008882);
    },
    [](double u) -> double {
       return 7.1458333333333312609 * u * (u - 1.0) * (1.0 * u - 0.8571428571428572063) * (u + 1.0) *
              (1.16666666666667 * u - 0.83333333333333325932) * (1.4 * u - 0.79999999999999993339) *
              (1.4 * u + 1.1999999999999999556) * (1.75 * u - 0.75) * (1.75 * u + 1.25) *
              (2.33333333333333 * u - 0.66666666666666662966) * (2.33333333333333 * u + 1.3333333333333337034) *
              (3.5 * u - 0.49999999999999977796) * (3.5 * u + 1.5) * (7.0 * u + 2.0000000000000008882);
    },
    [](double u) -> double {
       return -1.0 * (u - 1.0) * (u + 1.0) * (1.16666666666667 * u - 1.0) * (1.16666666666667 * u + 1.0) * (1.4 * u + 1.0) *
              (1.4 * u - 1.0) * (1.75 * u - 1.0) * (1.75 * u + 1.0) * (2.33333333333333 * u - 1.0) *
              (2.33333333333333 * u + 1.0) * (3.5 * u - 1.0) * (3.5 * u + 1.0) * (7.0 * u + 1.0) * (7.0 * u - 1.0);
    },
    [](double u) -> double {
       return 7.14583333333333659 * u * (u - 1.0) * (1.0 * u + 0.8571428571428572063) * (u + 1.0) *
              (1.16666666666667 * u + 0.83333333333333337034) * (1.4 * u - 1.1999999999999997335) *
              (1.4 * u + 0.80000000000000004441) * (1.75 * u - 1.249999999999999778) * (1.75 * u + 0.75) *
              (2.33333333333333 * u - 1.3333333333333330373) * (2.33333333333333 * u + 0.6666666666666668517) *
              (3.5 * u - 1.4999999999999995559) * (3.5 * u + 0.50000000000000022204) *
              (6.99999999999999 * u - 1.9999999999999984457);
    },
    [](double u) -> double {
       return -3.8111111111111095617 * u * (0.875 * u + 0.75) * (u - 1.0) * (1.0 * u + 0.71428571428571430157) * (u + 1.0) *
              (1.16666666666667 * u + 0.66666666666666651864) * (1.4 * u + 0.59999999999999986677) *
              (1.75 * u + 0.49999999999999988898) * (1.75 * u - 1.5) * (2.33333333333333 * u + 0.33333333333333331483) *
              (2.33333333333333 * u - 1.6666666666666671848) * (3.5 * u - 2.0000000000000008882) *
              (6.99999999999999 * u - 0.99999999999999844569) * (7.0 * u - 3.0000000000000017764);
    },
    [](double u) -> double {
       return 2.8583333333333329485 * u * (0.777777777777778 * u + 0.66666666666666662966) * (0.875 * u + 0.625) * (u - 1.0) *
              (1.0 * u + 0.57142857142857139685) * (u + 1.0) * (1.16666666666667 * u + 0.5) * (1.4 * u + 0.39999999999999996669) *
              (1.75 * u + 0.25000000000000005551) * (2.33333333333333 * u - 2.0) * (3.5 * u - 0.49999999999999961142) *
              (3.5 * u - 2.5000000000000008882) * (7.0 * u - 4.0000000000000017764) * (7.0 * u - 2.0000000000000017764);
    },
    [](double u) -> double {
       return -2.5984848484848477312 * u * (0.7 * u + 0.5999999999999999778) * (0.777777777777778 * u + 0.55555555555555569125) *
              (0.875 * u + 0.5) * (u - 1.0) * (1.0 * u + 0.42857142857142860315) * (u + 1.0) *
              (1.16666666666667 * u + 0.33333333333333331483) * (1.4 * u + 0.20000000000000006661) *
              (2.33333333333333 * u - 0.3333333333333331483) * (3.5 * u - 2.9999999999999991118) *
              (3.5 * u - 1.0000000000000008882) * (7.0 * u - 5.0000000000000017764) * (7.0 * u - 3.0000000000000017764);
    },
    [](double u) -> double {
       return 2.8583333333333325044 * u * (0.636363636363636 * u + 0.5454545454545455252) * (0.7 * u + 0.50000000000000011102) *
              (0.777777777777778 * u + 0.44444444444444447528) * (0.875 * u + 0.375) * (u - 1.0) * (u + 1.0) *
              (1.0 * u + 0.28571428571428575394) * (1.16666666666667 * u + 0.16666666666666674068) *
              (1.75 * u - 0.24999999999999988898) * (2.33333333333333 * u - 0.66666666666666718477) *
              (3.5 * u - 1.5000000000000008882) * (6.99999999999999 * u - 5.9999999999999937828) *
              (7.0 * u - 4.0000000000000017764);
    },
    [](double u) -> double {
       return -4.3974358974358986885 * u * (0.583333333333333 * u + 0.5) * (0.636363636363636 * u + 0.4545454545454544748) *
              (0.7 * u + 0.39999999999999996669) * (0.777777777777778 * u + 0.33333333333333331483) * (0.875 * u + 0.25) *
              (u - 1.0) * (1.0 * u + 0.14285714285714290472) * (u + 1.0) * (1.4 * u - 0.19999999999999987232) *
              (1.75 * u - 0.50000000000000022204) * (2.33333333333333 * u - 1.0) * (3.5 * u - 1.9999999999999991118) *
              (6.99999999999999 * u - 4.9999999999999937828);
    },
    [](double u) -> double {
       return 0.5 * u * (0.538461538461538 * u + 0.46153846153846156408) * (0.583333333333333 * u + 0.41666666666666662966) *
              (0.636363636363636 * u + 0.36363636363636359095) * (0.7 * u + 0.2999999999999999889) *
              (0.777777777777778 * u + 0.22222222222222223764) * (0.875 * u + 0.12500000000000005551) * (u + 1.0) *
              (1.16666666666667 * u - 0.16666666666666657415) * (1.4 * u - 0.40000000000000018874) * (1.75 * u - 0.75) *
              (2.33333333333333 * u - 1.3333333333333330373) * (3.5 * u - 2.4999999999999986677) *
              (7.0 * u - 6.0000000000000035527);
    },
    [](double u) -> double {
       return -10.2192574368469 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              10.2192574368469 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              20.6656094834015 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              20.6656094834015 * u * u * u * u * u * u * u * u * u * u * u * u -
              15.5589077621699 * u * u * u * u * u * u * u * u * u * u * u +
              15.5589077621698 * u * u * u * u * u * u * u * u * u * u + 5.46598445852192 * u * u * u * u * u * u * u * u * u -
              5.46598445852192 * u * u * u * u * u * u * u * u - 0.919786148837637 * u * u * u * u * u * u * u +
              0.919786148837635 * u * u * u * u * u * u + 0.0680742051565283 * u * u * u * u * u -
              0.0680742051565286 * u * u * u * u - 0.00172319182373183 * u * u * u + 0.00172319182373183 * u * u +
              6.39259815216064e-6 * u - 6.3925981521606402961e-6;
    },
    [](double u) -> double {
       return 4.017857142857144126 * (0.576923076923077 * u - 0.5) * (0.625 * u - 0.45833333333333337034) *
              (0.681818181818182 * u - 0.40909090909090911614) * (0.75 * u - 0.34999999999999992228) *
              (0.833333333333333 * u - 0.2777777777777777346) * (0.9375 * u - 0.18749999999999994449) * (u - 1.0) * (u + 1.0) *
              (1.07142857142857 * u - 0.071428571428571410729) * (1.25 * u + 0.08333333333333331483) *
              (1.5 * u + 0.29999999999999987788) * (1.875 * u + 0.62500000000000011102) * (2.5 * u + 1.1666666666666667407) *
              (3.75 * u + 2.2499999999999995559) * (7.5 * u + 5.5000000000000017764);
    },
    [](double u) -> double {
       return -2.163461538461538769 * (0.625 * u - 0.54166666666666674068) * (0.681818181818182 * u - 0.5) *
              (0.75 * u - 0.4500000000000000111) * (0.833333333333333 * u - 0.38888888888888883955) *
              (0.9375 * u - 0.31249999999999994449) * (u - 1.0) * (u + 1.0) * (1.07142857142857 * u - 0.21428571428571421831) *
              (1.25 * u - 0.08333333333333331483) * (1.5 * u + 0.099999999999999963918) * (1.875 * u + 0.37499999999999983347) *
              (2.5 * u + 0.83333333333333348136) * (3.75 * u + 1.7499999999999995559) * (7.5 * u + 4.4999999999999973355) *
              (7.5 * u + 6.5000000000000017764);
    },
    [](double u) -> double {
       return 1.5625 * (0.681818181818182 * u - 0.59090909090909093937) * (0.75 * u - 0.54999999999999993339) *
              (0.833333333333333 * u - 0.5) * (0.9375 * u - 0.4375) * (u - 1.0) * (u + 1.0) *
              (1.07142857142857 * u - 0.35714285714285715079) * (1.25 * u - 0.24999999999999994449) *
              (1.5 * u - 0.099999999999999977796) * (1.875 * u + 0.12499999999999997224) * (2.5 * u + 0.49999999999999988898) *
              (3.75 * u + 3.2499999999999995559) * (3.75 * u + 1.2500000000000004441) * (7.5 * u + 5.4999999999999973355) *
              (7.5 * u + 3.5000000000000008882);
    },
    [](double u) -> double {
       return -1.2784090909090908283 * (0.75 * u - 0.64999999999999991118) * (0.833333333333333 * u - 0.61111111111111116045) *
              (0.9375 * u - 0.5625) * (u - 1.0) * (u + 1.0) * (1.07142857142857 * u - 0.5) * (1.25 * u - 0.41666666666666657415) *
              (1.5 * u - 0.29999999999999993339) * (1.875 * u - 0.12499999999999997224) * (2.5 * u + 0.16666666666666662966) *
              (2.5 * u + 2.1666666666666669627) * (3.75 * u + 0.74999999999999966693) * (3.75 * u + 2.7499999999999995559) *
              (7.5 * u + 2.5000000000000008882) * (7.5 * u + 4.5000000000000008882);
    },
    [](double u) -> double {
       return 1.124999999999999778 * (0.833333333333333 * u - 0.72222222222222220989) * (0.9375 * u - 0.68749999999999988898) *
              (u - 1.0) * (u + 1.0) * (1.07142857142857 * u - 0.6428571428571427937) * (1.25 * u - 0.58333333333333325932) *
              (1.5 * u - 0.49999999999999988898) * (1.875 * u - 0.37499999999999988898) * (1.875 * u + 1.625) *
              (2.5 * u - 0.16666666666666662966) * (2.5 * u + 1.8333333333333334814) * (3.75 * u + 0.24999999999999988898) *
              (3.75 * u + 2.2500000000000004441) * (7.5 * u + 1.4999999999999988898) * (7.5 * u + 3.5000000000000008882);
    },
    [](double u) -> double {
       return -1.0416666666666667407 * (0.9375 * u - 0.8125) * (u - 1.0) * (u + 1.0) *
              (1.07142857142857 * u - 0.78571428571428580945) * (1.25 * u - 0.75000000000000011102) *
              (1.5 * u + 1.2999999999999998224) * (1.5 * u - 0.70000000000000006661) * (1.875 * u + 1.374999999999999778) *
              (1.875 * u - 0.625) * (2.5 * u + 1.5) * (2.5 * u - 0.5) * (3.75 * u + 1.7499999999999995559) * (3.75 * u - 0.25) *
              (7.5 * u + 2.4999999999999986677) * (7.5 * u + 0.5);
    },
    [](double u) -> double {
       return 1.0044642857142855874 * (u - 1.0) * (u + 1.0) * (1.07142857142857 * u - 0.92857142857142860315) *
              (1.25 * u - 0.91666666666666674068) * (1.25 * u + 1.0833333333333334814) * (1.5 * u - 0.9000000000000000222) *
              (1.5 * u + 1.0999999999999998668) * (1.875 * u + 1.125) * (1.875 * u - 0.875) * (2.5 * u + 1.1666666666666667407) *
              (2.5 * u - 0.83333333333333325932) * (3.75 * u + 1.249999999999999778) * (3.75 * u - 0.75) * (7.5 * u - 0.5) *
              (7.5 * u + 1.5);
    },
    [](double u) -> double {
       return -1.0044642857142855874 * (u - 1.0) * (u + 1.0) * (1.07142857142857 * u + 0.92857142857142860315) *
              (1.25 * u - 1.0833333333333334814) * (1.25 * u + 0.91666666666666674068) * (1.5 * u - 1.0999999999999998668) *
              (1.5 * u + 0.89999999999999991118) * (1.875 * u - 1.125) * (1.875 * u + 0.875) *
              (2.5 * u + 0.83333333333333348136) * (2.5 * u - 1.1666666666666665186) * (3.75 * u - 1.25) * (3.75 * u + 0.75) *
              (7.5 * u - 1.5) * (7.5 * u + 0.5);
    },
    [](double u) -> double {
       return 1.0416666666666667407 * (0.9375 * u + 0.8125) * (u - 1.0) * (u + 1.0) *
              (1.07142857142857 * u + 0.78571428571428580945) * (1.25 * u + 0.75) * (1.5 * u - 1.2999999999999998224) *
              (1.5 * u + 0.69999999999999995559) * (1.875 * u - 1.374999999999999778) * (1.875 * u + 0.62500000000000011102) *
              (2.5 * u - 1.499999999999999778) * (2.5 * u + 0.5) * (3.75 * u - 1.75) * (3.75 * u + 0.25) * (7.5 * u - 2.5) *
              (7.5 * u - 0.5);
    },
    [](double u) -> double {
       return -1.124999999999999778 * (0.833333333333333 * u + 0.72222222222222232091) * (0.9375 * u + 0.6875) * (u - 1.0) *
              (u + 1.0) * (1.07142857142857 * u + 0.64285714285714290472) * (1.25 * u + 0.58333333333333337034) *
              (1.5 * u + 0.5) * (1.875 * u - 1.624999999999999778) * (1.875 * u + 0.375) * (2.5 * u - 1.8333333333333328152) *
              (2.5 * u + 0.16666666666666665741) * (3.75 * u - 2.2499999999999991118) * (3.75 * u - 0.25) * (7.5 * u - 3.5) *
              (7.5 * u - 1.5);
    },
    [](double u) -> double {
       return 1.2784090909090908283 * (0.75 * u + 0.6500000000000000222) * (0.833333333333333 * u + 0.61111111111111116045) *
              (0.9375 * u + 0.56250000000000011102) * (u - 1.0) * (u + 1.0) * (1.07142857142857 * u + 0.50000000000000011102) *
              (1.25 * u + 0.41666666666666674068) * (1.5 * u + 0.30000000000000004441) * (1.875 * u + 0.125) *
              (2.5 * u - 2.1666666666666660745) * (2.5 * u - 0.16666666666666665741) * (3.75 * u - 2.7499999999999986677) *
              (3.75 * u - 0.75) * (7.49999999999999 * u - 4.4999999999999946709) * (7.5 * u - 2.5);
    },
    [](double u) -> double {
       return -1.562500000000000222 * (0.681818181818182 * u + 0.59090909090909093937) * (0.75 * u + 0.54999999999999993339) *
              (0.833333333333333 * u + 0.49999999999999994449) * (0.9375 * u + 0.43749999999999988898) * (u - 1.0) * (u + 1.0) *
              (1.07142857142857 * u + 0.35714285714285709528) * (1.25 * u + 0.24999999999999994449) *
              (1.5 * u + 0.099999999999999963918) * (1.875 * u - 0.12499999999999994449) * (2.5 * u - 0.49999999999999972244) *
              (3.75 * u - 1.2499999999999988898) * (3.75 * u - 3.2500000000000008882) *
              (7.49999999999999 * u - 3.4999999999999942268) * (7.5 * u - 5.5000000000000017764);
    },
    [](double u) -> double {
       return 2.163461538461538769 * (0.625 * u + 0.54166666666666674068) * (0.681818181818182 * u + 0.5) *
              (0.75 * u + 0.44999999999999990008) * (0.833333333333333 * u + 0.38888888888888883955) *
              (0.9375 * u + 0.31249999999999994449) * (u - 1.0) * (u + 1.0) * (1.07142857142857 * u + 0.21428571428571421831) *
              (1.25 * u + 0.08333333333333331483) * (1.5 * u - 0.099999999999999963918) * (1.875 * u - 0.37499999999999983347) *
              (2.5 * u - 0.83333333333333281523) * (3.75 * u - 1.7499999999999986677) * (7.5 * u - 6.5000000000000017764) *
              (7.5 * u - 4.5000000000000017764);
    },
    [](double u) -> double {
       return -4.017857142857144126 * (0.576923076923077 * u + 0.5) * (0.625 * u + 0.45833333333333337034) *
              (0.681818181818182 * u + 0.40909090909090906063) * (0.75 * u + 0.3499999999999999778) *
              (0.833333333333333 * u + 0.27777777777777779011) * (0.9375 * u + 0.18749999999999994449) * (u - 1.0) * (u + 1.0) *
              (1.07142857142857 * u + 0.071428571428571410729) * (1.25 * u - 0.08333333333333331483) *
              (1.5 * u - 0.29999999999999987788) * (1.875 * u - 0.62499999999999966693) * (2.5 * u - 1.1666666666666660745) *
              (3.75 * u - 2.2500000000000008882) * (7.5 * u - 5.5000000000000017764);
    },
    [](double u) -> double {
       return 0.5 * (0.535714285714286 * u + 0.46428571428571430157) * (0.576923076923077 * u + 0.42307692307692307265) *
              (0.625 * u + 0.375) * (0.681818181818182 * u + 0.31818181818181817677) * (0.75 * u + 0.25) *
              (0.833333333333333 * u + 0.16666666666666662966) * (0.9375 * u + 0.062499999999999986122) * (u + 1.0) *
              (1.07142857142857 * u - 0.071428571428571410729) * (1.25 * u - 0.24999999999999994449) *
              (1.5 * u - 0.49999999999999983347) * (1.875 * u - 0.87499999999999955591) * (2.5 * u - 1.5000000000000004441) *
              (3.75 * u - 2.7500000000000008882) * (7.5 * u - 6.5000000000000017764);
    }};

/// regular spaced lagrange. u  between -1 and 1.
/*!
 *  n: order of the lagrange basis,
 *  i: node on which lag =1. all other node, lag = 0.
 *  -1---1---2---  ---i--- ---n-1---n
 *  u:  -1+0*2/n_ -1+2/n_ -1+2*2/n-    _-1+i*2/n-  _-1+2n/n
 */
double lagPoly(int n, int i, double u)
{
   //     double res = 1.;
   //     if (n==0) return res;
   //     double ui = -1.+(2.*i)/n;
   //     for (int k = 0; k < (n+1); ++k) {
   //         if (k!=i) {
   //             double rk  =  -1.+(2.*k)/n;
   //             res *= (u - rk)/(ui-rk);
   //           }
   //       }
   //     return res;

   //    cout<<"lagpoly "<<n*(n+1)/2 + i<<endl;

   return lagPolys[n * (n + 1) / 2 + i](u);
}

// Tabulated Lagrange polynomial derivative functions (for efficiency)
// Thank you sympy !
std::function<double(double)> dLagPolys[136] = {
    [](double u) -> double { return 0.; },
    [](double u) -> double { return -0.50000000000000000000; },
    [](double u) -> double { return 0.50000000000000000000; },
    [](double u) -> double { return 1.0 * u - 0.5; },
    [](double u) -> double { return -2.0 * u; },
    [](double u) -> double { return 1.0 * u + 0.5; },
    [](double u) -> double { return -1.6875 * u * u + 1.125 * u + 0.0625; },
    [](double u) -> double { return 5.0625 * u * u - 1.125 * u - 1.6875; },
    [](double u) -> double { return -5.0625 * u * u - 1.125 * u + 1.6875; },
    [](double u) -> double { return 1.6875 * u * u + 1.125 * u - 0.0625; },
    [](double u) -> double {
       return 2.66666666666667 * u * u * u - 2.0 * u * u - 0.333333333333333 * u + 0.16666666666666665741;
    },
    [](double u) -> double { return -10.6666666666667 * u * u * u + 4.0 * u * u + 5.33333333333333 * u - 1.3333333333333332593; },
    [](double u) -> double { return u * (16.0 * u * u - 10.0); },
    [](double u) -> double { return -10.6666666666667 * u * u * u - 4.0 * u * u + 5.33333333333333 * u + 1.3333333333333332593; },
    [](double u) -> double {
       return 2.66666666666667 * u * u * u + 2.0 * u * u - 0.333333333333333 * u - 0.16666666666666665741;
    },
    [](double u) -> double {
       return -4.06901041666667 * u * u * u * u + 3.25520833333333 * u * u * u + 0.9765625 * u * u - 0.651041666666667 * u -
              0.01171875;
    },
    [](double u) -> double {
       return 20.3450520833333 * u * u * u * u - 9.765625 * u * u * u - 12.6953125 * u * u + 5.078125 * u +
              0.16276041666666654639;
    },
    [](double u) -> double {
       return -40.6901041666667 * u * u * u * u + 6.51041666666667 * u * u * u + 33.203125 * u * u - 4.42708333333334 * u -
              2.9296875000000008882;
    },
    [](double u) -> double {
       return 40.6901041666667 * u * u * u * u + 6.51041666666666 * u * u * u - 33.203125 * u * u - 4.42708333333333 * u +
              2.9296875000000008882;
    },
    [](double u) -> double {
       return -20.3450520833333 * u * u * u * u - 9.765625 * u * u * u + 12.6953125 * u * u + 5.078125 * u -
              0.16276041666666651864;
    },
    [](double u) -> double {
       return 4.06901041666667 * u * u * u * u + 3.25520833333333 * u * u * u - 0.9765625 * u * u - 0.651041666666667 * u +
              0.01171875;
    },
    [](double u) -> double {
       return 6.075 * u * u * u * u * u - 5.0625 * u * u * u * u - 2.25 * u * u * u + 1.6875 * u * u + 0.1 * u -
              0.050000000000000002776;
    },
    [](double u) -> double {
       return -36.45 * u * u * u * u * u + 20.25 * u * u * u * u + 27.0 * u * u * u - 13.5 * u * u - 1.35 * u +
              0.4500000000000000111;
    },
    [](double u) -> double {
       return 91.1249999999999 * u * u * u * u * u - 25.3125 * u * u * u * u - 87.75 * u * u * u + 21.9375 * u * u + 13.5 * u -
              2.25;
    },
    [](double u) -> double {
       return -121.5 * u * u * u * u * u - 1.24344978758018e-14 * u * u * u * u + 126.0 * u * u * u +
              1.24344978758018e-14 * u * u - 24.5 * u - 1.3322676295501878485e-15;
    },
    [](double u) -> double {
       return 91.125 * u * u * u * u * u + 25.3125 * u * u * u * u - 87.75 * u * u * u - 21.9375 * u * u + 13.5 * u + 2.25;
    },
    [](double u) -> double {
       return -36.45 * u * u * u * u * u - 20.25 * u * u * u * u + 27.0 * u * u * u + 13.5 * u * u - 1.35 * u -
              0.4500000000000000111;
    },
    [](double u) -> double {
       return 6.075 * u * u * u * u * u + 5.0625 * u * u * u * u - 2.25 * u * u * u - 1.6875 * u * u + 0.1 * u +
              0.050000000000000002776;
    },
    [](double u) -> double {
       return -8.93601345486111 * u * u * u * u * u * u + 7.65944010416667 * u * u * u * u * u +
              4.55919053819445 * u * u * u * u - 3.64735243055555 * u * u * u - 0.413118489583334 * u * u +
              0.275412326388889 * u + 0.0024414062500000173472;
    },
    [](double u) -> double {
       return 62.5520941840278 * u * u * u * u * u * u - 38.2972005208333 * u * u * u * u * u - 53.7984483506944 * u * u * u * u +
              30.7419704861111 * u * u * u + 5.57151692708334 * u * u - 2.65310329861111 * u - 0.033496093750000191513;
    },
    [](double u) -> double {
       return -187.656282552083 * u * u * u * u * u * u + 68.9349609374999 * u * u * u * u * u + 205.16357421875 * u * u * u * u -
              70.3417968749999 * u * u * u - 43.51142578125 * u * u + 12.4318359375 * u + 0.27913411458333386994;
    },
    [](double u) -> double {
       return 312.760470920139 * u * u * u * u * u * u - 38.2972005208333 * u * u * u * u * u - 378.412814670139 * u * u * u * u +
              43.2471788194444 * u * u * u + 105.568522135417 * u * u - 10.0541449652778 * u - 4.18701171875;
    },
    [](double u) -> double {
       return -312.760470920139 * u * u * u * u * u * u - 38.2972005208334 * u * u * u * u * u +
              378.412814670139 * u * u * u * u + 43.2471788194445 * u * u * u - 105.568522135417 * u * u - 10.0541449652778 * u +
              4.1870117187499982236;
    },
    [](double u) -> double {
       return 187.656282552083 * u * u * u * u * u * u + 68.9349609375 * u * u * u * u * u - 205.16357421875 * u * u * u * u -
              70.341796875 * u * u * u + 43.51142578125 * u * u + 12.4318359375 * u - 0.27913411458333292625;
    },
    [](double u) -> double {
       return -62.5520941840278 * u * u * u * u * u * u - 38.2972005208334 * u * u * u * u * u +
              53.7984483506945 * u * u * u * u + 30.7419704861111 * u * u * u - 5.57151692708333 * u * u - 2.65310329861111 * u +
              0.033496093749999893141;
    },
    [](double u) -> double {
       return 8.9360134548611 * u * u * u * u * u * u + 7.65944010416666 * u * u * u * u * u - 4.55919053819444 * u * u * u * u -
              3.64735243055555 * u * u * u + 0.413118489583333 * u * u + 0.275412326388889 * u - 0.0024414062499999817854;
    },
    [](double u) -> double {
       return 13.0031746031746 * u * u * u * u * u * u * u - 11.3777777777778 * u * u * u * u * u * u -
              8.53333333333333 * u * u * u * u * u + 7.11111111111111 * u * u * u * u + 1.24444444444444 * u * u * u -
              0.933333333333333 * u * u - 0.0285714285714286 * u + 0.014285714285714285268;
    },
    [](double u) -> double {
       return -104.025396825397 * u * u * u * u * u * u * u + 68.2666666666667 * u * u * u * u * u * u +
              102.4 * u * u * u * u * u - 64.0 * u * u * u * u - 17.0666666666667 * u * u * u + 9.6 * u * u +
              0.406349206349206 * u - 0.1523809523809523947;
    },
    [](double u) -> double {
       return 364.088888888889 * u * u * u * u * u * u * u - 159.288888888889 * u * u * u * u * u * u -
              443.733333333333 * u * u * u * u * u + 184.888888888889 * u * u * u * u + 120.177777777778 * u * u * u -
              45.0666666666667 * u * u - 3.2 * u + 0.80000000000000004441;
    },
    [](double u) -> double {
       return -728.177777777778 * u * u * u * u * u * u * u + 159.288888888889 * u * u * u * u * u * u +
              989.866666666667 * u * u * u * u * u - 206.222222222222 * u * u * u * u - 347.022222222222 * u * u * u +
              65.0666666666667 * u * u + 25.6 * u - 3.2000000000000001776;
    },
    [](double u) -> double {
       return u * (910.222222222222 * u * u * u * u * u * u - 1280.0 * u * u * u * u - 1.4210854715202e-14 * u * u * u +
                   485.333333333333 * u * u - 7.105427357601e-15 * u - 45.5555555555556);
    },
    [](double u) -> double {
       return -728.177777777778 * u * u * u * u * u * u * u - 159.288888888889 * u * u * u * u * u * u +
              989.866666666667 * u * u * u * u * u + 206.222222222222 * u * u * u * u - 347.022222222222 * u * u * u -
              65.0666666666666 * u * u + 25.6 * u + 3.2000000000000001776;
    },
    [](double u) -> double {
       return 364.088888888889 * u * u * u * u * u * u * u + 159.288888888889 * u * u * u * u * u * u -
              443.733333333333 * u * u * u * u * u - 184.888888888889 * u * u * u * u + 120.177777777778 * u * u * u +
              45.0666666666667 * u * u - 3.2 * u - 0.80000000000000004441;
    },
    [](double u) -> double {
       return -104.025396825397 * u * u * u * u * u * u * u - 68.2666666666667 * u * u * u * u * u * u +
              102.4 * u * u * u * u * u + 64.0 * u * u * u * u - 17.0666666666667 * u * u * u - 9.6 * u * u +
              0.406349206349206 * u + 0.1523809523809523947;
    },
    [](double u) -> double {
       return 13.0031746031746 * u * u * u * u * u * u * u + 11.3777777777778 * u * u * u * u * u * u -
              8.53333333333333 * u * u * u * u * u - 7.11111111111111 * u * u * u * u + 1.24444444444444 * u * u * u +
              0.933333333333333 * u * u - 0.0285714285714286 * u - 0.014285714285714285268;
    },
    [](double u) -> double {
       return -18.7668810163225 * u * u * u * u * u * u * u * u + 16.6816720145089 * u * u * u * u * u * u * u +
              15.1370727539062 * u * u * u * u * u * u - 12.9746337890625 * u * u * u * u * u - 3.13687133789062 * u * u * u * u +
              2.5094970703125 * u * u * u + 0.152035086495536 * u * u - 0.101356724330357 * u - 0.00053405761718750086736;
    },
    [](double u) -> double {
       return 168.901929146903 * u * u * u * u * u * u * u * u - 116.771704101562 * u * u * u * u * u * u * u -
              188.132189941406 * u * u * u * u * u * u + 125.421459960937 * u * u * u * u * u + 44.2499084472656 * u * u * u * u -
              27.5332763671875 * u * u * u - 2.24634312220982 * u * u + 1.1647705078125 * u + 0.007945469447544652522;
    },
    [](double u) -> double {
       return -675.607716587611 * u * u * u * u * u * u * u * u + 333.633440290178 * u * u * u * u * u * u * u +
              908.224365234375 * u * u * u * u * u * u - 432.48779296875 * u * u * u * u * u - 301.940551757812 * u * u * u * u +
              134.19580078125 * u * u * u + 17.3147670200893 * u * u - 6.41287667410714 * u - 0.062292480468750091593;
    },
    [](double u) -> double {
       return 1576.41800537109 * u * u * u * u * u * u * u * u - 467.08681640625 * u * u * u * u * u * u * u -
              2361.38334960937 * u * u * u * u * u * u + 674.68095703125 * u * u * u * u * u + 984.310180664062 * u * u * u * u -
              262.48271484375 * u * u * u - 105.248583984375 * u * u + 23.38857421875 * u + 0.40374755859375;
    },
    [](double u) -> double {
       return -2364.62700805664 * u * u * u * u * u * u * u * u + 233.543408203124 * u * u * u * u * u * u * u +
              3723.71989746093 * u * u * u * u * u * u - 354.639990234374 * u * u * u * u * u - 1724.74530029297 * u * u * u * u +
              153.310693359375 * u * u * u + 243.528002929687 * u * u - 18.039111328125 * u - 5.4505920410156214473;
    },
    [](double u) -> double {
       return 2364.62700805664 * u * u * u * u * u * u * u * u + 233.543408203126 * u * u * u * u * u * u * u -
              3723.71989746094 * u * u * u * u * u * u - 354.639990234376 * u * u * u * u * u + 1724.74530029297 * u * u * u * u +
              153.310693359375 * u * u * u - 243.528002929687 * u * u - 18.039111328125 * u + 5.4505920410156232236;
    },
    [](double u) -> double {
       return -1576.41800537109 * u * u * u * u * u * u * u * u - 467.08681640625 * u * u * u * u * u * u * u +
              2361.38334960937 * u * u * u * u * u * u + 674.68095703125 * u * u * u * u * u - 984.310180664062 * u * u * u * u -
              262.48271484375 * u * u * u + 105.248583984375 * u * u + 23.38857421875 * u - 0.40374755859375011102;
    },
    [](double u) -> double {
       return 675.607716587611 * u * u * u * u * u * u * u * u + 333.633440290179 * u * u * u * u * u * u * u -
              908.224365234375 * u * u * u * u * u * u - 432.48779296875 * u * u * u * u * u + 301.940551757812 * u * u * u * u +
              134.19580078125 * u * u * u - 17.3147670200893 * u * u - 6.41287667410714 * u + 0.062292480468750077716;
    },
    [](double u) -> double {
       return -168.901929146903 * u * u * u * u * u * u * u * u - 116.771704101563 * u * u * u * u * u * u * u +
              188.132189941406 * u * u * u * u * u * u + 125.421459960938 * u * u * u * u * u - 44.2499084472656 * u * u * u * u -
              27.5332763671875 * u * u * u + 2.24634312220982 * u * u + 1.1647705078125 * u - 0.0079454694475446369095;
    },
    [](double u) -> double {
       return 18.7668810163225 * u * u * u * u * u * u * u * u + 16.6816720145089 * u * u * u * u * u * u * u -
              15.1370727539062 * u * u * u * u * u * u - 12.9746337890625 * u * u * u * u * u + 3.13687133789062 * u * u * u * u +
              2.5094970703125 * u * u * u - 0.152035086495535 * u * u - 0.101356724330357 * u + 0.00053405761718750086736;
    },
    [](double u) -> double {
       return 26.9114445546737 * u * u * u * u * u * u * u * u * u - 24.2203000992064 * u * u * u * u * u * u * u * u -
              25.8349867724868 * u * u * u * u * u * u * u + 22.6056134259259 * u * u * u * u * u * u +
              7.05295138888889 * u * u * u * u * u - 5.87745949074074 * u * u * u * u - 0.564925044091711 * u * u * u +
              0.423693783068783 * u * u + 0.00793650793650794 * u - 0.0039682539682539680337;
    },
    [](double u) -> double {
       return -269.114445546737 * u * u * u * u * u * u * u * u * u + 193.762400793651 * u * u * u * u * u * u * u * u +
              335.854828042328 * u * u * u * u * u * u * u - 235.09837962963 * u * u * u * u * u * u -
              103.081597222222 * u * u * u * u * u + 68.7210648148148 * u * u * u * u + 8.68744488536155 * u * u * u -
              5.21246693121693 * u * u - 0.124007936507936 * u + 0.049603174603174558788;
    },
    [](double u) -> double {
       return 1211.01500496032 * u * u * u * u * u * u * u * u * u - 653.948102678572 * u * u * u * u * u * u * u * u -
              1782.61408730159 * u * u * u * u * u * u * u + 935.872395833334 * u * u * u * u * u * u +
              708.0078125 * u * u * u * u * u - 354.00390625 * u * u * u * u - 67.0882936507937 * u * u * u +
              30.1897321428571 * u * u + 0.992063492063492 * u - 0.2976190476190476164;
    },
    [](double u) -> double {
       return -3229.37334656085 * u * u * u * u * u * u * u * u * u + 1162.5744047619 * u * u * u * u * u * u * u * u +
              5270.3373015873 * u * u * u * u * u * u * u - 1844.61805555556 * u * u * u * u * u * u -
              2539.0625 * u * u * u * u * u + 846.354166666667 * u * u * u * u + 361.193783068783 * u * u * u -
              108.358134920635 * u * u - 5.95238095238095 * u + 1.1904761904761893554;
    },
    [](double u) -> double {
       return 5651.40335648148 * u * u * u * u * u * u * u * u * u - 1017.25260416667 * u * u * u * u * u * u * u * u -
              9765.625 * u * u * u * u * u * u * u + 1708.984375 * u * u * u * u * u * u + 5257.16145833334 * u * u * u * u * u -
              876.193576388889 * u * u * u * u - 965.856481481482 * u * u * u + 144.878472222222 * u * u + 41.6666666666667 * u -
              4.1666666666666669627;
    },
    [](double u) -> double {
       return -6781.68402777778 * u * u * u * u * u * u * u * u * u + 5.6843418860808e-14 * u * u * u * u * u * u * u * u +
              11935.7638888889 * u * u * u * u * u * u * u + 1.13686837721616e-13 * u * u * u * u * u * u -
              6660.15625 * u * u * u * u * u - 2.8421709430404e-13 * u * u * u * u + 1327.25694444445 * u * u * u -
              7.105427357601e-15 * u * u - 73.1805555555556 * u - 2.2204460492503130808e-16;
    },
    [](double u) -> double {
       return 5651.40335648148 * u * u * u * u * u * u * u * u * u + 1017.25260416667 * u * u * u * u * u * u * u * u -
              9765.625 * u * u * u * u * u * u * u - 1708.984375 * u * u * u * u * u * u + 5257.16145833334 * u * u * u * u * u +
              876.193576388889 * u * u * u * u - 965.856481481482 * u * u * u - 144.878472222222 * u * u + 41.6666666666667 * u +
              4.1666666666666669627;
    },
    [](double u) -> double {
       return -3229.37334656085 * u * u * u * u * u * u * u * u * u - 1162.5744047619 * u * u * u * u * u * u * u * u +
              5270.3373015873 * u * u * u * u * u * u * u + 1844.61805555555 * u * u * u * u * u * u -
              2539.0625 * u * u * u * u * u - 846.354166666666 * u * u * u * u + 361.193783068783 * u * u * u +
              108.358134920635 * u * u - 5.95238095238095 * u - 1.1904761904761893554;
    },
    [](double u) -> double {
       return 1211.01500496032 * u * u * u * u * u * u * u * u * u + 653.948102678571 * u * u * u * u * u * u * u * u -
              1782.61408730158 * u * u * u * u * u * u * u - 935.872395833332 * u * u * u * u * u * u +
              708.007812499999 * u * u * u * u * u + 354.00390625 * u * u * u * u - 67.0882936507935 * u * u * u -
              30.1897321428571 * u * u + 0.99206349206349 * u + 0.2976190476190471168;
    },
    [](double u) -> double {
       return -269.114445546737 * u * u * u * u * u * u * u * u * u - 193.762400793651 * u * u * u * u * u * u * u * u +
              335.854828042328 * u * u * u * u * u * u * u + 235.09837962963 * u * u * u * u * u * u -
              103.081597222222 * u * u * u * u * u - 68.7210648148148 * u * u * u * u + 8.68744488536155 * u * u * u +
              5.21246693121693 * u * u - 0.124007936507936 * u - 0.049603174603174558788;
    },
    [](double u) -> double {
       return 26.9114445546737 * u * u * u * u * u * u * u * u * u + 24.2203000992064 * u * u * u * u * u * u * u * u -
              25.8349867724868 * u * u * u * u * u * u * u - 22.6056134259259 * u * u * u * u * u * u +
              7.05295138888889 * u * u * u * u * u + 5.87745949074074 * u * u * u * u - 0.56492504409171 * u * u * u -
              0.423693783068783 * u * u + 0.00793650793650794 * u + 0.0039682539682539680337;
    },
    [](double u) -> double {
       return -38.3907460222463 * u * u * u * u * u * u * u * u * u * u + 34.9006782020421 * u * u * u * u * u * u * u * u * u +
              42.832650520688 * u * u * u * u * u * u * u * u - 38.0734671295004 * u * u * u * u * u * u * u -
              14.6472718584979 * u * u * u * u * u * u + 12.5548044501411 * u * u * u * u * u + 1.70222368862684 * u * u * u * u -
              1.36177895090147 * u * u * u - 0.0516391899472191 * u * u + 0.0344261266314794 * u + 0.00012016296386718869262;
    },
    [](double u) -> double {
       return 422.29820624471 * u * u * u * u * u * u * u * u * u * u - 314.106103818379 * u * u * u * u * u * u * u * u * u -
              585.37955711607 * u * u * u * u * u * u * u * u + 425.730586993506 * u * u * u * u * u * u * u +
              222.792714058205 * u * u * u * u * u * u - 156.244241027832 * u * u * u * u * u - 27.2800036973634 * u * u * u * u +
              17.8560024200924 * u * u * u + 0.845615040208931 * u * u - 0.461244567386688 * u - 0.0019745296902126951988;
    },
    [](double u) -> double {
       return -2111.49103122355 * u * u * u * u * u * u * u * u * u * u + 1221.52373707147 * u * u * u * u * u * u * u * u * u +
              3383.77939113436 * u * u * u * u * u * u * u * u - 1914.05702932852 * u * u * u * u * u * u * u -
              1572.10380285758 * u * u * u * u * u * u + 857.511165195042 * u * u * u * u * u + 214.310070750995 * u * u * u * u -
              109.103308745961 * u * u * u - 6.9415034521194 * u * u + 2.94488025241428 * u + 0.01632009233747226673;
    },
    [](double u) -> double {
       return 6334.47309367064 * u * u * u * u * u * u * u * u * u * u - 2617.55086515316 * u * u * u * u * u * u * u * u * u -
              11179.3217858996 * u * u * u * u * u * u * u * u + 4516.89769127256 * u * u * u * u * u * u * u +
              6117.16327353584 * u * u * u * u * u * u - 2383.31036631266 * u * u * u * u * u - 1104.38344236909 * u * u * u * u +
              401.593979043305 * u * u * u + 40.1335732051305 * u * u - 12.1616888500395 * u - 0.095962142944336636941;
    },
    [](double u) -> double {
       return -12668.9461873413 * u * u * u * u * u * u * u * u * u * u + 3141.06103818379 * u * u * u * u * u * u * u * u * u +
              23729.2883884612 * u * u * u * u * u * u * u * u - 5752.55476083907 * u * u * u * u * u * u * u -
              14454.5445972019 * u * u * u * u * u * u + 3378.98445129395 * u * u * u * u * u + 3246.35402937793 * u * u * u * u -
              708.295424591549 * u * u * u - 209.202613455152 * u * u + 38.0368388100276 * u + 0.53312301635742409545;
    },
    [](double u) -> double {
       return 17736.5246622778 * u * u * u * u * u * u * u * u * u * u - 1465.82848448577 * u * u * u * u * u * u * u * u * u -
              34180.4551155091 * u * u * u * u * u * u * u * u + 2762.05697903104 * u * u * u * u * u * u * u +
              21938.5296078491 * u * u * u * u * u * u - 1709.49581359864 * u * u * u * u * u - 5490.51979884395 * u * u * u * u +
              399.310530825015 * u * u * u + 468.487994232178 * u * u - 28.3932117716471 * u - 6.7173500061035182895;
    },
    [](double u) -> double {
       return -17736.5246622778 * u * u * u * u * u * u * u * u * u * u - 1465.82848448577 * u * u * u * u * u * u * u * u * u +
              34180.455115509 * u * u * u * u * u * u * u * u + 2762.05697903103 * u * u * u * u * u * u * u -
              21938.5296078491 * u * u * u * u * u * u - 1709.49581359863 * u * u * u * u * u + 5490.51979884395 * u * u * u * u +
              399.310530825014 * u * u * u - 468.487994232178 * u * u - 28.3932117716471 * u + 6.7173500061035174014;
    },
    [](double u) -> double {
       return 12668.9461873413 * u * u * u * u * u * u * u * u * u * u + 3141.06103818379 * u * u * u * u * u * u * u * u * u -
              23729.2883884611 * u * u * u * u * u * u * u * u - 5752.55476083907 * u * u * u * u * u * u * u +
              14454.5445972019 * u * u * u * u * u * u + 3378.98445129394 * u * u * u * u * u - 3246.35402937793 * u * u * u * u -
              708.295424591549 * u * u * u + 209.202613455152 * u * u + 38.0368388100276 * u - 0.53312301635741832229;
    },
    [](double u) -> double {
       return -6334.47309367063 * u * u * u * u * u * u * u * u * u * u - 2617.55086515315 * u * u * u * u * u * u * u * u * u +
              11179.3217858996 * u * u * u * u * u * u * u * u + 4516.89769127255 * u * u * u * u * u * u * u -
              6117.16327353583 * u * u * u * u * u * u - 2383.31036631266 * u * u * u * u * u + 1104.38344236909 * u * u * u * u +
              401.593979043305 * u * u * u - 40.1335732051304 * u * u - 12.1616888500395 * u + 0.095962142944335221406;
    },
    [](double u) -> double {
       return 2111.49103122355 * u * u * u * u * u * u * u * u * u * u + 1221.52373707147 * u * u * u * u * u * u * u * u * u -
              3383.77939113435 * u * u * u * u * u * u * u * u - 1914.05702932852 * u * u * u * u * u * u * u +
              1572.10380285757 * u * u * u * u * u * u + 857.511165195041 * u * u * u * u * u - 214.310070750995 * u * u * u * u -
              109.103308745961 * u * u * u + 6.94150345211937 * u * u + 2.94488025241428 * u - 0.016320092337471902438;
    },
    [](double u) -> double {
       return -422.29820624471 * u * u * u * u * u * u * u * u * u * u - 314.106103818379 * u * u * u * u * u * u * u * u * u +
              585.37955711607 * u * u * u * u * u * u * u * u + 425.730586993506 * u * u * u * u * u * u * u -
              222.792714058205 * u * u * u * u * u * u - 156.244241027832 * u * u * u * u * u + 27.2800036973634 * u * u * u * u +
              17.8560024200924 * u * u * u - 0.845615040208928 * u * u - 0.461244567386688 * u + 0.0019745296902126548665;
    },
    [](double u) -> double {
       return 38.3907460222463 * u * u * u * u * u * u * u * u * u * u + 34.9006782020421 * u * u * u * u * u * u * u * u * u -
              42.832650520688 * u * u * u * u * u * u * u * u - 38.0734671295005 * u * u * u * u * u * u * u +
              14.6472718584979 * u * u * u * u * u * u + 12.5548044501411 * u * u * u * u * u - 1.70222368862683 * u * u * u * u -
              1.36177895090147 * u * u * u + 0.051639189947219 * u * u + 0.0344261266314794 * u - 0.00012016296386718560265;
    },
    [](double u) -> double {
       return 54.5329870129871 * u * u * u * u * u * u * u * u * u * u * u -
              49.9885714285714 * u * u * u * u * u * u * u * u * u * u - 69.4285714285715 * u * u * u * u * u * u * u * u * u +
              62.4857142857143 * u * u * u * u * u * u * u * u + 28.6971428571429 * u * u * u * u * u * u * u -
              25.11 * u * u * u * u * u * u - 4.46785714285714 * u * u * u * u * u + 3.72321428571429 * u * u * u * u +
              0.228095238095239 * u * u * u - 0.171071428571429 * u * u - 0.00216450216450216 * u + 0.0010822510822510822501;
    },
    [](double u) -> double {
       return -654.395844155844 * u * u * u * u * u * u * u * u * u * u * u +
              499.885714285714 * u * u * u * u * u * u * u * u * u * u + 999.771428571429 * u * u * u * u * u * u * u * u * u -
              749.828571428572 * u * u * u * u * u * u * u * u - 455.451428571429 * u * u * u * u * u * u * u +
              332.1 * u * u * u * u * u * u + 74.6742857142857 * u * u * u * u * u - 51.8571428571429 * u * u * u * u -
              3.90857142857142 * u * u * u + 2.44285714285714 * u * u + 0.0374025974025974 * u - 0.015584415584415584402;
    },
    [](double u) -> double {
       return 3599.17714285714 * u * u * u * u * u * u * u * u * u * u * u -
              2199.49714285714 * u * u * u * u * u * u * u * u * u * u - 6248.57142857142 * u * u * u * u * u * u * u * u * u +
              3749.14285714285 * u * u * u * u * u * u * u * u + 3338.12571428571 * u * u * u * u * u * u * u -
              1947.24 * u * u * u * u * u * u - 602.678571428571 * u * u * u * u * u + 334.821428571428 * u * u * u * u +
              33.0685714285714 * u * u * u - 16.5342857142857 * u * u - 0.321428571428571 * u + 0.10714285714285713691;
    },
    [](double u) -> double {
       return -11997.2571428571 * u * u * u * u * u * u * u * u * u * u * u +
              5498.74285714285 * u * u * u * u * u * u * u * u * u * u + 22772.5714285714 * u * u * u * u * u * u * u * u * u -
              10247.6571428571 * u * u * u * u * u * u * u * u - 13978.2857142857 * u * u * u * u * u * u * u +
              6115.5 * u * u * u * u * u * u + 3097.02857142857 * u * u * u * u * u - 1290.42857142857 * u * u * u * u -
              189.295238095238 * u * u * u + 70.9857142857143 * u * u + 1.9047619047619 * u - 0.47619047619047616404;
    },
    [](double u) -> double {
       return 26993.8285714285 * u * u * u * u * u * u * u * u * u * u * u -
              8248.11428571428 * u * u * u * u * u * u * u * u * u * u - 54362.5714285714 * u * u * u * u * u * u * u * u * u +
              16308.7714285714 * u * u * u * u * u * u * u * u + 36866.5714285714 * u * u * u * u * u * u * u -
              10752.75 * u * u * u * u * u * u - 9793.18928571428 * u * u * u * u * u + 2720.33035714286 * u * u * u * u +
              861.87857142857 * u * u * u - 215.469642857143 * u * u - 9.64285714285713 * u + 1.6071428571428563181;
    },
    [](double u) -> double {
       return -43190.1257142857 * u * u * u * u * u * u * u * u * u * u * u +
              6598.49142857143 * u * u * u * u * u * u * u * u * u * u + 89979.4285714286 * u * u * u * u * u * u * u * u * u -
              13496.9142857143 * u * u * u * u * u * u * u * u - 64718.5371428571 * u * u * u * u * u * u * u +
              9438.11999999999 * u * u * u * u * u * u + 19236.3428571429 * u * u * u * u * u - 2671.71428571429 * u * u * u * u -
              2183.45142857143 * u * u * u + 272.931428571429 * u * u + 61.7142857142857 * u - 5.1428571428571432378;
    },
    [](double u) -> double {
       return 50388.48 * u * u * u * u * u * u * u * u * u * u * u +
              7.50333128962666e-12 * u * u * u * u * u * u * u * u * u * u - 106142.4 * u * u * u * u * u * u * u * u * u -
              1.45519152283669e-11 * u * u * u * u * u * u * u * u + 77837.76 * u * u * u * u * u * u * u +
              5.91171556152403e-12 * u * u * u * u * u * u - 24015.42 * u * u * u * u * u - 1.25055521493778e-12 * u * u * u * u +
              2962.96 * u * u * u - 2.8421709430404e-14 * u * u - 107.38 * u + 3.1086244689504383132e-15;
    },
    [](double u) -> double {
       return -43190.1257142857 * u * u * u * u * u * u * u * u * u * u * u -
              6598.49142857143 * u * u * u * u * u * u * u * u * u * u + 89979.4285714286 * u * u * u * u * u * u * u * u * u +
              13496.9142857143 * u * u * u * u * u * u * u * u - 64718.5371428572 * u * u * u * u * u * u * u -
              9438.12000000002 * u * u * u * u * u * u + 19236.3428571429 * u * u * u * u * u + 2671.71428571429 * u * u * u * u -
              2183.45142857143 * u * u * u - 272.931428571429 * u * u + 61.7142857142857 * u + 5.1428571428571459023;
    },
    [](double u) -> double {
       return 26993.8285714286 * u * u * u * u * u * u * u * u * u * u * u +
              8248.11428571429 * u * u * u * u * u * u * u * u * u * u - 54362.5714285715 * u * u * u * u * u * u * u * u * u -
              16308.7714285714 * u * u * u * u * u * u * u * u + 36866.5714285715 * u * u * u * u * u * u * u +
              10752.75 * u * u * u * u * u * u - 9793.1892857143 * u * u * u * u * u - 2720.33035714286 * u * u * u * u +
              861.878571428572 * u * u * u + 215.469642857143 * u * u - 9.64285714285717 * u - 1.6071428571428592047;
    },
    [](double u) -> double {
       return -11997.2571428571 * u * u * u * u * u * u * u * u * u * u * u -
              5498.74285714285 * u * u * u * u * u * u * u * u * u * u + 22772.5714285714 * u * u * u * u * u * u * u * u * u +
              10247.6571428571 * u * u * u * u * u * u * u * u - 13978.2857142857 * u * u * u * u * u * u * u -
              6115.5 * u * u * u * u * u * u + 3097.02857142857 * u * u * u * u * u + 1290.42857142857 * u * u * u * u -
              189.295238095238 * u * u * u - 70.9857142857143 * u * u + 1.90476190476191 * u + 0.47619047619047616404;
    },
    [](double u) -> double {
       return 3599.17714285714 * u * u * u * u * u * u * u * u * u * u * u +
              2199.49714285714 * u * u * u * u * u * u * u * u * u * u - 6248.57142857143 * u * u * u * u * u * u * u * u * u -
              3749.14285714286 * u * u * u * u * u * u * u * u + 3338.12571428571 * u * u * u * u * u * u * u +
              1947.24 * u * u * u * u * u * u - 602.678571428571 * u * u * u * u * u - 334.821428571428 * u * u * u * u +
              33.0685714285714 * u * u * u + 16.5342857142857 * u * u - 0.321428571428571 * u - 0.10714285714285713691;
    },
    [](double u) -> double {
       return -654.395844155845 * u * u * u * u * u * u * u * u * u * u * u -
              499.885714285715 * u * u * u * u * u * u * u * u * u * u + 999.771428571431 * u * u * u * u * u * u * u * u * u +
              749.828571428573 * u * u * u * u * u * u * u * u - 455.45142857143 * u * u * u * u * u * u * u -
              332.100000000001 * u * u * u * u * u * u + 74.6742857142859 * u * u * u * u * u + 51.857142857143 * u * u * u * u -
              3.90857142857144 * u * u * u - 2.44285714285715 * u * u + 0.0374025974025975 * u + 0.015584415584415641648;
    },
    [](double u) -> double {
       return 54.532987012987 * u * u * u * u * u * u * u * u * u * u * u +
              49.9885714285714 * u * u * u * u * u * u * u * u * u * u - 69.4285714285714 * u * u * u * u * u * u * u * u * u -
              62.4857142857143 * u * u * u * u * u * u * u * u + 28.6971428571429 * u * u * u * u * u * u * u +
              25.11 * u * u * u * u * u * u - 4.46785714285714 * u * u * u * u * u - 3.72321428571429 * u * u * u * u +
              0.228095238095238 * u * u * u + 0.171071428571429 * u * u - 0.00216450216450216 * u - 0.0010822510822510822501;
    },
    [](double u) -> double {
       return -77.185668920912 * u * u * u * u * u * u * u * u * u * u * u * u +
              71.2483097731496 * u * u * u * u * u * u * u * u * u * u * u +
              110.52622413527 * u * u * u * u * u * u * u * u * u * u - 100.478385577519 * u * u * u * u * u * u * u * u * u -
              53.7767454170802 * u * u * u * u * u * u * u * u + 47.8015514818491 * u * u * u * u * u * u * u +
              10.6335901490356 * u * u * u * u * u * u - 9.11450584203057 * u * u * u * u * u -
              0.799450949750739 * u * u * u * u + 0.639560759800592 * u * u * u + 0.0166438365980316 * u * u -
              0.0110958910653544 * u - 0.000027537345886230353554;
    },
    [](double u) -> double {
       return 1003.41369597186 * u * u * u * u * u * u * u * u * u * u * u * u -
              783.731407504645 * u * u * u * u * u * u * u * u * u * u * u -
              1677.98903914456 * u * u * u * u * u * u * u * u * u * u + 1290.76079934197 * u * u * u * u * u * u * u * u * u +
              891.730808333971 * u * u * u * u * u * u * u * u - 670.703513960594 * u * u * u * u * u * u * u -
              185.400756822672 * u * u * u * u * u * u + 134.46648297029 * u * u * u * u * u + 14.317232515348 * u * u * u * u -
              9.69166508731248 * u * u * u - 0.301607516881878 * u * u + 0.170137573625676 * u + 0.00049999627200039903959;
    },
    [](double u) -> double {
       return -6020.48217583114 * u * u * u * u * u * u * u * u * u * u * u * u +
              3847.40872775008 * u * u * u * u * u * u * u * u * u * u * u +
              11273.6748617976 * u * u * u * u * u * u * u * u * u * u - 7095.3198430894 * u * u * u * u * u * u * u * u * u -
              6827.23875399529 * u * u * u * u * u * u * u * u + 4201.37769476633 * u * u * u * u * u * u * u +
              1546.81058593138 * u * u * u * u * u * u - 917.887600442795 * u * u * u * u * u - 125.249515095548 * u * u * u * u +
              69.3689622067649 * u * u * u + 2.69402409562243 * u * u - 1.24339573644112 * u - 0.0044814480675591256836;
    },
    [](double u) -> double {
       return 22075.1013113808 * u * u * u * u * u * u * u * u * u * u * u * u -
              10972.239705065 * u * u * u * u * u * u * u * u * u * u * u -
              44873.6469989198 * u * u * u * u * u * u * u * u * u * u + 21966.1209085622 * u * u * u * u * u * u * u * u * u +
              30598.1655043273 * u * u * u * u * u * u * u * u - 14645.275796943 * u * u * u * u * u * u * u -
              8103.3670212434 * u * u * u * u * u * u + 3740.01554826619 * u * u * u * u * u + 722.721425654761 * u * u * u * u -
              311.326152589743 * u * u * u - 16.2180509291754 * u * u + 5.82186443611428 * u + 0.027163062776837942902;
    },
    [](double u) -> double {
       return -55187.7532784521 * u * u * u * u * u * u * u * u * u * u * u * u +
              19593.2851876162 * u * u * u * u * u * u * u * u * u * u * u +
              118815.690945416 * u * u * u * u * u * u * u * u * u * u - 41543.9478830126 * u * u * u * u * u * u * u * u * u -
              88727.6167482259 * u * u * u * u * u * u * u * u + 30334.2279481114 * u * u * u * u * u * u * u +
              27300.2281967779 * u * u * u * u * u * u - 9000.07522970701 * u * u * u * u * u - 3117.69432991404 * u * u * u * u +
              959.290563050472 * u * u * u + 78.1463689614856 * u * u - 20.037530502945 * u - 0.13309900760650617002;
    },
    [](double u) -> double {
       return 99337.9559012138 * u * u * u * u * u * u * u * u * u * u * u * u -
              21160.7480026254 * u * u * u * u * u * u * u * u * u * u * u -
              221826.131839488 * u * u * u * u * u * u * u * u * u * u + 46536.950735557 * u * u * u * u * u * u * u * u * u +
              175929.418674991 * u * u * u * u * u * u * u * u - 36088.0858820495 * u * u * u * u * u * u * u -
              60228.9353427386 * u * u * u * u * u * u + 11913.4157820802 * u * u * u * u * u + 8513.15939135312 * u * u * u * u -
              1571.66019532673 * u * u * u - 366.738530369562 * u * u + 56.421312364548 * u + 0.66549503803252973988;
    },
    [](double u) -> double {
       return -132450.607868285 * u * u * u * u * u * u * u * u * u * u * u * u +
              9404.77689005578 * u * u * u * u * u * u * u * u * u * u * u +
              301073.434544477 * u * u * u * u * u * u * u * u * u * u - 21054.0863317818 * u * u * u * u * u * u * u * u * u -
              246002.123229429 * u * u * u * u * u * u * u * u + 16820.6579985935 * u * u * u * u * u * u * u +
              88889.1105727594 * u * u * u * u * u * u - 5860.82047732476 * u * u * u * u * u - 14029.9075635347 * u * u * u * u +
              863.378926986752 * u * u * u + 801.865198754811 * u * u - 41.1212922438364 * u - 7.9859404563903844121;
    },
    [](double u) -> double {
       return 132450.607868285 * u * u * u * u * u * u * u * u * u * u * u * u +
              9404.77689005572 * u * u * u * u * u * u * u * u * u * u * u -
              301073.434544477 * u * u * u * u * u * u * u * u * u * u - 21054.0863317815 * u * u * u * u * u * u * u * u * u +
              246002.123229429 * u * u * u * u * u * u * u * u + 16820.6579985934 * u * u * u * u * u * u * u -
              88889.1105727594 * u * u * u * u * u * u - 5860.82047732478 * u * u * u * u * u + 14029.9075635347 * u * u * u * u +
              863.378926986745 * u * u * u - 801.865198754811 * u * u - 41.1212922438364 * u + 7.9859404563903853003;
    },
    [](double u) -> double {
       return -99337.9559012138 * u * u * u * u * u * u * u * u * u * u * u * u -
              21160.7480026254 * u * u * u * u * u * u * u * u * u * u * u +
              221826.131839488 * u * u * u * u * u * u * u * u * u * u + 46536.9507355568 * u * u * u * u * u * u * u * u * u -
              175929.418674991 * u * u * u * u * u * u * u * u - 36088.0858820494 * u * u * u * u * u * u * u +
              60228.9353427385 * u * u * u * u * u * u + 11913.4157820801 * u * u * u * u * u - 8513.15939135311 * u * u * u * u -
              1571.66019532673 * u * u * u + 366.738530369562 * u * u + 56.421312364548 * u - 0.66549503803253051704;
    },
    [](double u) -> double {
       return 55187.7532784521 * u * u * u * u * u * u * u * u * u * u * u * u +
              19593.2851876161 * u * u * u * u * u * u * u * u * u * u * u -
              118815.690945416 * u * u * u * u * u * u * u * u * u * u - 41543.9478830125 * u * u * u * u * u * u * u * u * u +
              88727.6167482259 * u * u * u * u * u * u * u * u + 30334.2279481114 * u * u * u * u * u * u * u -
              27300.2281967779 * u * u * u * u * u * u - 9000.075229707 * u * u * u * u * u + 3117.69432991403 * u * u * u * u +
              959.290563050471 * u * u * u - 78.1463689614857 * u * u - 20.037530502945 * u + 0.13309900760650608675;
    },
    [](double u) -> double {
       return -22075.1013113808 * u * u * u * u * u * u * u * u * u * u * u * u -
              10972.239705065 * u * u * u * u * u * u * u * u * u * u * u +
              44873.6469989198 * u * u * u * u * u * u * u * u * u * u + 21966.1209085621 * u * u * u * u * u * u * u * u * u -
              30598.1655043273 * u * u * u * u * u * u * u * u - 14645.275796943 * u * u * u * u * u * u * u +
              8103.3670212434 * u * u * u * u * u * u + 3740.01554826618 * u * u * u * u * u - 722.72142565476 * u * u * u * u -
              311.326152589743 * u * u * u + 16.2180509291755 * u * u + 5.82186443611427 * u - 0.027163062776838074741;
    },
    [](double u) -> double {
       return 6020.48217583114 * u * u * u * u * u * u * u * u * u * u * u * u +
              3847.40872775008 * u * u * u * u * u * u * u * u * u * u * u -
              11273.6748617976 * u * u * u * u * u * u * u * u * u * u - 7095.31984308939 * u * u * u * u * u * u * u * u * u +
              6827.23875399529 * u * u * u * u * u * u * u * u + 4201.37769476633 * u * u * u * u * u * u * u -
              1546.81058593138 * u * u * u * u * u * u - 917.887600442796 * u * u * u * u * u + 125.249515095548 * u * u * u * u +
              69.3689622067649 * u * u * u - 2.69402409562244 * u * u - 1.24339573644112 * u + 0.0044814480675591621128;
    },
    [](double u) -> double {
       return -1003.41369597185 * u * u * u * u * u * u * u * u * u * u * u * u -
              783.731407504644 * u * u * u * u * u * u * u * u * u * u * u +
              1677.98903914456 * u * u * u * u * u * u * u * u * u * u + 1290.76079934197 * u * u * u * u * u * u * u * u * u -
              891.730808333971 * u * u * u * u * u * u * u * u - 670.703513960594 * u * u * u * u * u * u * u +
              185.400756822672 * u * u * u * u * u * u + 134.46648297029 * u * u * u * u * u - 14.317232515348 * u * u * u * u -
              9.69166508731246 * u * u * u + 0.301607516881881 * u * u + 0.170137573625675 * u - 0.00049999627200040131642;
    },
    [](double u) -> double {
       return 77.1856689209121 * u * u * u * u * u * u * u * u * u * u * u * u +
              71.2483097731496 * u * u * u * u * u * u * u * u * u * u * u -
              110.526224135271 * u * u * u * u * u * u * u * u * u * u - 100.478385577519 * u * u * u * u * u * u * u * u * u +
              53.7767454170804 * u * u * u * u * u * u * u * u + 47.8015514818491 * u * u * u * u * u * u * u -
              10.6335901490357 * u * u * u * u * u * u - 9.11450584203058 * u * u * u * u * u +
              0.799450949750738 * u * u * u * u + 0.639560759800592 * u * u * u - 0.0166438365980317 * u * u -
              0.0110958910653544 * u + 0.000027537345886230516184;
    },
    [](double u) -> double {
       return 108.91614058026 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              101.136416253098 * u * u * u * u * u * u * u * u * u * u * u * u -
              173.37671357674 * u * u * u * u * u * u * u * u * u * u * u +
              158.928654112011 * u * u * u * u * u * u * u * u * u * u + 97.3032576195988 * u * u * u * u * u * u * u * u * u -
              87.5729318576389 * u * u * u * u * u * u * u * u - 23.5267740483539 * u * u * u * u * u * u * u +
              20.5859272923097 * u * u * u * u * u * u + 2.39914737654321 * u * u * u * u * u - 1.99928948045268 * u * u * u * u -
              0.0851725589225593 * u * u * u + 0.0638794191919192 * u * u + 0.000582750582750583 * u - 0.00029137529137529137504;
    },
    [](double u) -> double {
       return -1524.82596812363 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              1213.63699503718 * u * u * u * u * u * u * u * u * u * u * u * u +
              2774.02741722783 * u * u * u * u * u * u * u * u * u * u * u -
              2179.59297067901 * u * u * u * u * u * u * u * u * u * u - 1686.58979873971 * u * u * u * u * u * u * u * u * u +
              1301.08355902778 * u * u * u * u * u * u * u * u + 427.869562757201 * u * u * u * u * u * u * u -
              320.902172067901 * u * u * u * u * u * u - 44.8543296682099 * u * u * u * u * u + 32.0388069058642 * u * u * u * u +
              1.61499041339319 * u * u * u - 1.03820812289562 * u * u - 0.0111046361046361 * u + 0.0047591297591297546804;
    },
    [](double u) -> double {
       return 9911.3687928036 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              6573.86705645137 * u * u * u * u * u * u * u * u * u * u * u * u -
              19938.322061325 * u * u * u * u * u * u * u * u * u * u * u +
              13054.8537306295 * u * u * u * u * u * u * u * u * u * u + 13525.1528091242 * u * u * u * u * u * u * u * u * u -
              8694.74109157985 * u * u * u * u * u * u * u * u - 3704.14305877057 * u * u * u * u * u * u * u +
              2315.08941173161 * u * u * u * u * u * u + 406.62305941358 * u * u * u * u * u - 242.037535365226 * u * u * u * u -
              14.9918139730639 * u * u * u + 8.0313289141414 * u * u + 0.103939393939394 * u - 0.03712121212121208963;
    },
    [](double u) -> double {
       return -39645.4751712145 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              21036.3745806444 * u * u * u * u * u * u * u * u * u * u * u * u +
              85994.8499340629 * u * u * u * u * u * u * u * u * u * u * u -
              45044.9213940329 * u * u * u * u * u * u * u * u * u * u - 64609.3630594136 * u * u * u * u * u * u * u * u * u +
              33227.6724305556 * u * u * u * u * u * u * u * u + 20020.9064711934 * u * u * u * u * u * u * u -
              10010.4532355967 * u * u * u * u * u * u - 2393.01982445988 * u * u * u * u * u + 1139.5332497428 * u * u * u * u +
              92.2664225589226 * u * u * u - 39.5427525252525 * u * u - 0.649621212121215 * u + 0.18560606060606069101;
    },
    [](double u) -> double {
       return 109025.05672084 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              43387.522572579 * u * u * u * u * u * u * u * u * u * u * u * u -
              249835.844264082 * u * u * u * u * u * u * u * u * u * u * u +
              98149.7959608891 * u * u * u * u * u * u * u * u * u * u + 203785.455874646 * u * u * u * u * u * u * u * u * u -
              78602.9615516492 * u * u * u * u * u * u * u * u - 71526.7641711677 * u * u * u * u * u * u * u +
              26822.5365641879 * u * u * u * u * u * u + 10211.0954436728 * u * u * u * u * u - 3646.81980131172 * u * u * u * u -
              435.950529835391 * u * u * u + 140.126956018518 * u * u + 3.17592592592593 * u - 0.68055555555555524716;
    },
    [](double u) -> double {
       return -218050.11344168 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              57850.030096772 * u * u * u * u * u * u * u * u * u * u * u * u +
              518743.127021605 * u * u * u * u * u * u * u * u * u * u * u -
              135861.295172325 * u * u * u * u * u * u * u * u * u * u - 448762.62414159 * u * u * u * u * u * u * u * u * u +
              115396.103350694 * u * u * u * u * u * u * u * u + 173634.552165638 * u * u * u * u * u * u * u -
              43408.6380414094 * u * u * u * u * u * u - 29497.857414159 * u * u * u * u * u + 7023.29938432356 * u * u * u * u +
              1767.2519212963 * u * u * u - 378.696840277778 * u * u - 14.2916666666667 * u + 2.0416666666666674068;
    },
    [](double u) -> double {
       return 327075.170162519 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              43387.5225725792 * u * u * u * u * u * u * u * u * u * u * u * u -
              795278.985176505 * u * u * u * u * u * u * u * u * u * u * u +
              104143.676630257 * u * u * u * u * u * u * u * u * u * u + 712551.755548322 * u * u * u * u * u * u * u * u * u -
              91613.7971419273 * u * u * u * u * u * u * u * u - 292778.750522762 * u * u * u * u * u * u * u +
              36597.3438153453 * u * u * u * u * u * u + 55711.5090162037 * u * u * u * u * u - 6632.32250192902 * u * u * u * u -
              4300.88652777778 * u * u * u + 460.809270833333 * u * u + 85.75 * u - 6.125;
    },
    [](double u) -> double {
       return -373800.19447145 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              7.95807864051312e-11 * u * u * u * u * u * u * u * u * u * u * u * u +
              915429.047685184 * u * u * u * u * u * u * u * u * u * u * u +
              1.29148247651756e-10 * u * u * u * u * u * u * u * u * u * u -
              829802.180979938 * u * u * u * u * u * u * u * u * u - 6.00266503170133e-11 * u * u * u * u * u * u * u * u +
              347899.712654321 * u * u * u * u * u * u * u + 5.0931703299284e-11 * u * u * u * u * u * u -
              68791.7901967592 * u * u * u * u * u - 7.27595761418343e-12 * u * u * u * u + 5781.56141975309 * u * u * u -
              1.70530256582424e-13 * u * u - 148.156111111111 * u - 4.4408920985006261617e-15;
    },
    [](double u) -> double {
       return 327075.170162519 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              43387.5225725791 * u * u * u * u * u * u * u * u * u * u * u * u -
              795278.985176504 * u * u * u * u * u * u * u * u * u * u * u -
              104143.676630257 * u * u * u * u * u * u * u * u * u * u + 712551.755548321 * u * u * u * u * u * u * u * u * u +
              91613.797141927 * u * u * u * u * u * u * u * u - 292778.750522762 * u * u * u * u * u * u * u -
              36597.3438153453 * u * u * u * u * u * u + 55711.5090162037 * u * u * u * u * u + 6632.32250192901 * u * u * u * u -
              4300.88652777777 * u * u * u - 460.809270833333 * u * u + 85.75 * u + 6.1249999999999973355;
    },
    [](double u) -> double {
       return -218050.113441679 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              57850.0300967722 * u * u * u * u * u * u * u * u * u * u * u * u +
              518743.127021604 * u * u * u * u * u * u * u * u * u * u * u +
              135861.295172325 * u * u * u * u * u * u * u * u * u * u - 448762.624141589 * u * u * u * u * u * u * u * u * u -
              115396.103350694 * u * u * u * u * u * u * u * u + 173634.552165638 * u * u * u * u * u * u * u +
              43408.6380414095 * u * u * u * u * u * u - 29497.8574141589 * u * u * u * u * u - 7023.29938432356 * u * u * u * u +
              1767.25192129629 * u * u * u + 378.696840277777 * u * u - 14.2916666666666 * u - 2.0416666666666642982;
    },
    [](double u) -> double {
       return 109025.05672084 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              43387.5225725791 * u * u * u * u * u * u * u * u * u * u * u * u -
              249835.844264082 * u * u * u * u * u * u * u * u * u * u * u -
              98149.7959608893 * u * u * u * u * u * u * u * u * u * u + 203785.455874646 * u * u * u * u * u * u * u * u * u +
              78602.9615516494 * u * u * u * u * u * u * u * u - 71526.7641711677 * u * u * u * u * u * u * u -
              26822.5365641879 * u * u * u * u * u * u + 10211.0954436728 * u * u * u * u * u + 3646.81980131173 * u * u * u * u -
              435.950529835391 * u * u * u - 140.126956018519 * u * u + 3.17592592592592 * u + 0.68055555555555602432;
    },
    [](double u) -> double {
       return -39645.4751712145 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              21036.3745806444 * u * u * u * u * u * u * u * u * u * u * u * u +
              85994.8499340629 * u * u * u * u * u * u * u * u * u * u * u +
              45044.921394033 * u * u * u * u * u * u * u * u * u * u - 64609.3630594136 * u * u * u * u * u * u * u * u * u -
              33227.6724305556 * u * u * u * u * u * u * u * u + 20020.9064711934 * u * u * u * u * u * u * u +
              10010.4532355967 * u * u * u * u * u * u - 2393.01982445988 * u * u * u * u * u - 1139.5332497428 * u * u * u * u +
              92.2664225589227 * u * u * u + 39.5427525252526 * u * u - 0.64962121212121 * u - 0.18560606060606082979;
    },
    [](double u) -> double {
       return 9911.36879280363 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              6573.86705645138 * u * u * u * u * u * u * u * u * u * u * u * u -
              19938.3220613251 * u * u * u * u * u * u * u * u * u * u * u -
              13054.8537306295 * u * u * u * u * u * u * u * u * u * u + 13525.1528091242 * u * u * u * u * u * u * u * u * u +
              8694.74109157987 * u * u * u * u * u * u * u * u - 3704.14305877058 * u * u * u * u * u * u * u -
              2315.08941173161 * u * u * u * u * u * u + 406.62305941358 * u * u * u * u * u + 242.037535365226 * u * u * u * u -
              14.991813973064 * u * u * u - 8.03132891414142 * u * u + 0.103939393939394 * u + 0.037121212121212159019;
    },
    [](double u) -> double {
       return -1524.82596812363 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              1213.63699503717 * u * u * u * u * u * u * u * u * u * u * u * u +
              2774.02741722783 * u * u * u * u * u * u * u * u * u * u * u +
              2179.59297067901 * u * u * u * u * u * u * u * u * u * u - 1686.58979873971 * u * u * u * u * u * u * u * u * u -
              1301.08355902778 * u * u * u * u * u * u * u * u + 427.869562757201 * u * u * u * u * u * u * u +
              320.902172067901 * u * u * u * u * u * u - 44.8543296682098 * u * u * u * u * u - 32.0388069058641 * u * u * u * u +
              1.61499041339319 * u * u * u + 1.03820812289562 * u * u - 0.0111046361046361 * u - 0.004759129759129751211;
    },
    [](double u) -> double {
       return 108.916140580259 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              101.136416253098 * u * u * u * u * u * u * u * u * u * u * u * u -
              173.37671357674 * u * u * u * u * u * u * u * u * u * u * u -
              158.928654112011 * u * u * u * u * u * u * u * u * u * u + 97.3032576195987 * u * u * u * u * u * u * u * u * u +
              87.5729318576387 * u * u * u * u * u * u * u * u - 23.5267740483539 * u * u * u * u * u * u * u -
              20.5859272923097 * u * u * u * u * u * u + 2.3991473765432 * u * u * u * u * u + 1.99928948045267 * u * u * u * u -
              0.085172558922559 * u * u * u - 0.0638794191919192 * u * u + 0.000582750582750583 * u + 0.00029137529137529137504;
    },
    [](double u) -> double {
       return -153.288861552704 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              143.069604115857 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              268.652923284219 * u * u * u * u * u * u * u * u * u * u * u * u -
              247.987313800818 * u * u * u * u * u * u * u * u * u * u * u -
              171.147985383868 * u * u * u * u * u * u * u * u * u * u + 155.589077621699 * u * u * u * u * u * u * u * u * u +
              49.1938601266974 * u * u * u * u * u * u * u * u - 43.7278756681754 * u * u * u * u * u * u * u -
              6.43850304186347 * u * u * u * u * u * u + 5.51871689302581 * u * u * u * u * u +
              0.340371025782643 * u * u * u * u - 0.272296820626114 * u * u * u - 0.0051695754711955 * u * u +
              0.00344638364746366 * u + 6.3925981521606242025e-6;
    },
    [](double u) -> double {
       return 2299.33292329055 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              1859.90485350614 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              4525.76847686493 * u * u * u * u * u * u * u * u * u * u * u * u +
              3620.61478149194 * u * u * u * u * u * u * u * u * u * u * u +
              3100.66804688956 * u * u * u * u * u * u * u * u * u * u - 2442.95058239783 * u * u * u * u * u * u * u * u * u -
              932.858995523077 * u * u * u * u * u * u * u * u + 718.646929884444 * u * u * u * u * u * u * u +
              125.531991943717 * u * u * u * u * u * u - 93.2523368724755 * u * u * u * u * u - 6.74061658610773 * u * u * u * u +
              4.67349416636801 * u * u * u + 0.10311165597251 * u * u - 0.0595756234507835 * u - 0.0001276628329203682681;
    },
    [](double u) -> double {
       return -16095.3304630338 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              11016.3595169209 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              34656.2271036643 * u * u * u * u * u * u * u * u * u * u * u * u -
              23459.5998855574 * u * u * u * u * u * u * u * u * u * u * u -
              26069.2530511213 * u * u * u * u * u * u * u * u * u * u + 17379.5020340809 * u * u * u * u * u * u * u * u * u +
              8398.05180606032 * u * u * u * u * u * u * u * u - 5474.28562172821 * u * u * u * u * u * u * u -
              1180.87478913367 * u * u * u * u * u * u + 742.264153169734 * u * u * u * u * u + 65.0208144654316 * u * u * u * u -
              38.1455444863863 * u * u * u - 1.00613058126949 * u * u + 0.491886061953974 * u + 0.0012481415813619463513;
    },
    [](double u) -> double {
       return 69746.43200648 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              39058.0019236288 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              160923.101047248 * u * u * u * u * u * u * u * u * u * u * u * u +
              89126.640580014 * u * u * u * u * u * u * u * u * u * u * u +
              132284.056495077 * u * u * u * u * u * u * u * u * u * u - 72154.9399064056 * u * u * u * u * u * u * u * u * u -
              47133.2660942738 * u * u * u * u * u * u * u * u + 25137.741916946 * u * u * u * u * u * u * u +
              7141.36487804354 * u * u * u * u * u * u - 3672.70193727953 * u * u * u * u * u - 411.031336696956 * u * u * u * u +
              197.295041614538 * u * u * u + 6.49067815151308 * u * u - 2.59627126060523 * u - 0.0080795337756474558399;
    },
    [](double u) -> double {
       return -209239.29601944 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              91135.3378218006 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              508559.983777027 * u * u * u * u * u * u * u * u * u * u * u * u -
              219071.993011642 * u * u * u * u * u * u * u * u * u * u * u -
              450196.996098384 * u * u * u * u * u * u * u * u * u * u + 190992.665011435 * u * u * u * u * u * u * u * u * u +
              178073.881626395 * u * u * u * u * u * u * u * u - 73867.6842302085 * u * u * u * u * u * u * u -
              30853.4004259854 * u * u * u * u * u * u + 12341.3601703942 * u * u * u * u * u + 1942.72907196798 * u * u * u * u -
              725.285520201379 * u * u * u - 31.9704085636836 * u * u + 9.94634933092377 * u + 0.040067891989435597266;
    },
    [](double u) -> double {
       return 460326.451242768 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              143212.673719972 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              1161386.58735768 * u * u * u * u * u * u * u * u * u * u * u * u +
              357349.719186978 * u * u * u * u * u * u * u * u * u * u * u +
              1087414.28619916 * u * u * u * u * u * u * u * u * u * u - 329519.480666411 * u * u * u * u * u * u * u * u * u -
              469383.917848712 * u * u * u * u * u * u * u * u + 139076.716399618 * u * u * u * u * u * u * u +
              93671.241654083 * u * u * u * u * u * u - 26763.2119011666 * u * u * u * u * u - 7398.49768712051 * u * u * u * u +
              1972.93271656547 * u * u * u + 135.571570253798 * u * u - 30.1270156119551 * u - 0.17277275025844554546;
    },
    [](double u) -> double {
       return -767210.75207128 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              143212.673719972 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              1982927.22676082 * u * u * u * u * u * u * u * u * u * u * u * u -
              366078.872632767 * u * u * u * u * u * u * u * u * u * u * u -
              1927226.33697226 * u * u * u * u * u * u * u * u * u * u + 350404.78854041 * u * u * u * u * u * u * u * u * u +
              883993.408143787 * u * u * u * u * u * u * u * u - 157154.383670006 * u * u * u * u * u * u * u -
              196223.456086591 * u * u * u * u * u * u + 33638.3067577012 * u * u * u * u * u + 19114.820899087 * u * u * u * u -
              3058.37134385393 * u * u * u - 589.252214080521 * u * u + 78.5669618774028 * u + 0.79987384378909909266;
    },
    [](double u) -> double {
       return 986413.824091647 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              61376.8601657025 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              2579874.02229836 * u * u * u * u * u * u * u * u * u * u * u * u +
              158761.478295283 * u * u * u * u * u * u * u * u * u * u * u +
              2554450.3628875 * u * u * u * u * u * u * u * u * u * u - 154815.173508333 * u * u * u * u * u * u * u * u * u -
              1208367.72255086 * u * u * u * u * u * u * u * u + 71606.9761511626 * u * u * u * u * u * u * u +
              283469.963399694 * u * u * u * u * u * u - 16198.2836228396 * u * u * u * u * u - 30884.5022440489 * u * u * u * u +
              1647.17345301594 * u * u * u + 1265.08007605313 * u * u - 56.2257811579169 * u - 9.2556830495595967534;
    },
    [](double u) -> double {
       return -986413.824091647 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              61376.8601657021 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              2579874.02229836 * u * u * u * u * u * u * u * u * u * u * u * u +
              158761.478295284 * u * u * u * u * u * u * u * u * u * u * u -
              2554450.3628875 * u * u * u * u * u * u * u * u * u * u - 154815.173508334 * u * u * u * u * u * u * u * u * u +
              1208367.72255086 * u * u * u * u * u * u * u * u + 71606.9761511626 * u * u * u * u * u * u * u -
              283469.963399694 * u * u * u * u * u * u - 16198.2836228396 * u * u * u * u * u + 30884.5022440489 * u * u * u * u +
              1647.17345301595 * u * u * u - 1265.08007605313 * u * u - 56.2257811579169 * u + 9.2556830495595932007;
    },
    [](double u) -> double {
       return 767210.752071281 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              143212.673719973 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              1982927.22676082 * u * u * u * u * u * u * u * u * u * u * u * u -
              366078.872632768 * u * u * u * u * u * u * u * u * u * u * u +
              1927226.33697226 * u * u * u * u * u * u * u * u * u * u + 350404.788540411 * u * u * u * u * u * u * u * u * u -
              883993.408143787 * u * u * u * u * u * u * u * u - 157154.383670007 * u * u * u * u * u * u * u +
              196223.456086591 * u * u * u * u * u * u + 33638.3067577013 * u * u * u * u * u - 19114.8208990871 * u * u * u * u -
              3058.37134385393 * u * u * u + 589.252214080521 * u * u + 78.5669618774028 * u - 0.79987384378910042493;
    },
    [](double u) -> double {
       return -460326.451242768 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              143212.673719973 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              1161386.58735768 * u * u * u * u * u * u * u * u * u * u * u * u +
              357349.719186978 * u * u * u * u * u * u * u * u * u * u * u -
              1087414.28619916 * u * u * u * u * u * u * u * u * u * u - 329519.480666412 * u * u * u * u * u * u * u * u * u +
              469383.917848713 * u * u * u * u * u * u * u * u + 139076.716399619 * u * u * u * u * u * u * u -
              93671.2416540831 * u * u * u * u * u * u - 26763.2119011666 * u * u * u * u * u + 7398.49768712052 * u * u * u * u +
              1972.93271656547 * u * u * u - 135.571570253798 * u * u - 30.1270156119551 * u + 0.17277275025844568423;
    },
    [](double u) -> double {
       return 209239.29601944 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              91135.3378218005 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              508559.983777027 * u * u * u * u * u * u * u * u * u * u * u * u -
              219071.993011642 * u * u * u * u * u * u * u * u * u * u * u +
              450196.996098384 * u * u * u * u * u * u * u * u * u * u + 190992.665011436 * u * u * u * u * u * u * u * u * u -
              178073.881626395 * u * u * u * u * u * u * u * u - 73867.6842302083 * u * u * u * u * u * u * u +
              30853.4004259854 * u * u * u * u * u * u + 12341.3601703942 * u * u * u * u * u - 1942.72907196798 * u * u * u * u -
              725.285520201379 * u * u * u + 31.9704085636835 * u * u + 9.94634933092376 * u - 0.040067891989435333588;
    },
    [](double u) -> double {
       return -69746.4320064799 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              39058.0019236287 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              160923.101047247 * u * u * u * u * u * u * u * u * u * u * u * u +
              89126.6405800138 * u * u * u * u * u * u * u * u * u * u * u -
              132284.056495076 * u * u * u * u * u * u * u * u * u * u - 72154.9399064054 * u * u * u * u * u * u * u * u * u +
              47133.2660942736 * u * u * u * u * u * u * u * u + 25137.7419169459 * u * u * u * u * u * u * u -
              7141.36487804351 * u * u * u * u * u * u - 3672.70193727952 * u * u * u * u * u + 411.031336696954 * u * u * u * u +
              197.295041614538 * u * u * u - 6.49067815151304 * u * u - 2.59627126060522 * u + 0.0080795337756474107371;
    },
    [](double u) -> double {
       return 16095.3304630338 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              11016.3595169209 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              34656.2271036643 * u * u * u * u * u * u * u * u * u * u * u * u -
              23459.5998855574 * u * u * u * u * u * u * u * u * u * u * u +
              26069.2530511213 * u * u * u * u * u * u * u * u * u * u + 17379.5020340809 * u * u * u * u * u * u * u * u * u -
              8398.05180606032 * u * u * u * u * u * u * u * u - 5474.28562172821 * u * u * u * u * u * u * u +
              1180.87478913367 * u * u * u * u * u * u + 742.264153169733 * u * u * u * u * u - 65.0208144654313 * u * u * u * u -
              38.1455444863864 * u * u * u + 1.00613058126949 * u * u + 0.491886061953974 * u - 0.001248141581361937244;
    },
    [](double u) -> double {
       return -2299.33292329055 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              1859.90485350613 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              4525.76847686493 * u * u * u * u * u * u * u * u * u * u * u * u +
              3620.61478149194 * u * u * u * u * u * u * u * u * u * u * u -
              3100.66804688956 * u * u * u * u * u * u * u * u * u * u - 2442.95058239783 * u * u * u * u * u * u * u * u * u +
              932.858995523075 * u * u * u * u * u * u * u * u + 718.646929884443 * u * u * u * u * u * u * u -
              125.531991943717 * u * u * u * u * u * u - 93.2523368724754 * u * u * u * u * u + 6.7406165861077 * u * u * u * u +
              4.67349416636801 * u * u * u - 0.10311165597251 * u * u - 0.0595756234507835 * u + 0.00012766283292036813258;
    },
    [](double u) -> double {
       return 153.288861552703 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              143.069604115856 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              268.652923284219 * u * u * u * u * u * u * u * u * u * u * u * u -
              247.987313800818 * u * u * u * u * u * u * u * u * u * u * u +
              171.147985383868 * u * u * u * u * u * u * u * u * u * u + 155.589077621698 * u * u * u * u * u * u * u * u * u -
              49.1938601266973 * u * u * u * u * u * u * u * u - 43.7278756681754 * u * u * u * u * u * u * u +
              6.43850304186343 * u * u * u * u * u * u + 5.51871689302581 * u * u * u * u * u -
              0.340371025782643 * u * u * u * u - 0.272296820626114 * u * u * u + 0.0051695754711955 * u * u +
              0.00344638364746366 * u - 6.3925981521605699924e-6;
    }};

double dLagPoly(int n, int i, double u)
{
   //     double res = 0.;
   //     if (n==0) return 0.;
   //     double ui = -1.+(2.*i)/n;
   //
   //     for(int l =0; l< (n+1); ++l) {
   //         if(l!=i) {
   //             double rl  =  -1.+(2.*l)/n;
   //             double resLoc=1./(ui-rl);
   //
   //             for (int k = 0; k < (n+1); ++k) {
   //                 if (k!=i && k!=l) {
   //                     double rk  =  -1.+(2.*k)/n;
   //                     resLoc *= (u - rk)/(ui-rk);
   //                   }
   //               }
   //             res+=resLoc;
   //           }
   //       }
   //
   //     return res;

   return dLagPolys[n * (n + 1) / 2 + i](u);
}

}  // end namespace xfem
