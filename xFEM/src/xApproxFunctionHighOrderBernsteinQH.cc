/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <boost/math/special_functions/binomial.hpp>
#include <cmath>

#include "xApproxFunctionHighOrder.h"
#include "xApproxFunctionHighOrderQH.h"
#include "xGeomElem.h"

using AOMD::mEntity;
using std::cout;
using std::endl;
using xtensor::xPoint;

namespace xfem
{
///////////////////////////////////////////////////////////////////////////////////////////////
// xApproxFunction Bernstein Implementation     QUADS AND HEXAHEDRA                           //
///////////////////////////////////////////////////////////////////////////////////////////////

xApproxFunctionScalarVertexBernsteinQH::xApproxFunctionScalarVertexBernsteinQH(int _order, int _node)
    : order(_order), inode(_node)
{
}

const int xApproxFunctionScalarVertexBernsteinQH::Vindices[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                                                    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

void xApproxFunctionScalarVertexBernsteinQH::correctUvwForCurrentEntity(const mEntity::mType& type, const xtensor::xPoint& uvw,
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

/// Nothing to change here.
void xApproxFunctionScalarVertexBernsteinQH::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvwCorr;
   correctUvwForCurrentEntity(geo_appro->getEntity()->getType(), geo_appro->getUVW(), uvwCorr);

   res = bernPoly(Vindices[inode][0] * order, order, uvwCorr(0)) * bernPoly(Vindices[inode][1] * order, order, uvwCorr(1)) *
         bernPoly(Vindices[inode][2] * order, order, uvwCorr(2));

   return;
}

void xApproxFunctionScalarVertexBernsteinQH::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                          xtensor::xVector<>& res) const
{
   xtensor::xPoint uvwCorr;
   correctUvwForCurrentEntity(geo_appro->getEntity()->getType(), geo_appro->getUVW(), uvwCorr);

   res[0] = dBernPoly(Vindices[inode][0] * order, order, uvwCorr(0)) * bernPoly(Vindices[inode][1] * order, order, uvwCorr(1)) *
            bernPoly(Vindices[inode][2] * order, order, uvwCorr(2));
   res[1] = dBernPoly(Vindices[inode][1] * order, order, uvwCorr(1)) * bernPoly(Vindices[inode][0] * order, order, uvwCorr(0)) *
            bernPoly(Vindices[inode][2] * order, order, uvwCorr(2));
   if (geo_appro->getEntity()->getType() == mEntity::HEX)
   {
      res[2] = dBernPoly(Vindices[inode][2] * order, order, uvwCorr(2)) *
               bernPoly(Vindices[inode][0] * order, order, uvwCorr(0)) * bernPoly(Vindices[inode][1] * order, order, uvwCorr(1));
   }

   return;
}

void xApproxFunctionScalarVertexBernsteinQH::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                     xtensor::xVector<>& res) const
{
   getGradLocal(geo_appro, geo_integ, res);
   res = (geo_appro)->PushBack(res);
   return;
}

xApproxFunctionScalarEdgeBernsteinQH::xApproxFunctionScalarEdgeBernsteinQH(int _order, int _edge, int _ionedge)
    : order(_order), ionedge(_ionedge), edge(_edge)
{
}

/// For an edge, v and w are meaningless (for the value, NOT for the gradient)
void xApproxFunctionScalarEdgeBernsteinQH::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);

   //   -1      0      1
   //    *------*------*
   //    a      b      c
   // bernPoly(0, order  , ) associe au noeud a
   // bernPoly(order, order  , ) associe au noeud c

   // On corrige le cas de l'edge (sinon, on a v et w nuls et donc la fonction de forme aussi)
   if (geo_appro->getEntity()->getType() == mEntity::EDGE)
   {
      v = -1.;
      w = -1.;
   }

   switch (edge)
   {
      case 0:  // 2D et 3D ici...
         res = bernPoly(1 + ionedge, order, u) * bernPoly(0, order, v);
         break;
      case 1:
         res = bernPoly(1 + ionedge, order, v) * bernPoly(order, order, u);
         break;
      case 2:
         res = bernPoly(1 + ionedge, order, u) * bernPoly(order, order, v);
         break;
      case 3:
         res = bernPoly(1 + ionedge, order, v) * bernPoly(0, order, u);
         break;
      case 4:  // Only 3D here...
         res = bernPoly(1 + ionedge, order, w) * bernPoly(0, order, u) * bernPoly(0, order, v);
         break;
      case 5:
         res = bernPoly(1 + ionedge, order, w) * bernPoly(order, order, u) * bernPoly(0, order, v);
         break;
      case 6:
         res = bernPoly(1 + ionedge, order, w) * bernPoly(order, order, u) * bernPoly(order, order, v);
         break;
      case 7:
         res = bernPoly(1 + ionedge, order, w) * bernPoly(0, order, u) * bernPoly(order, order, v);
         break;
      case 8:
         res = bernPoly(1 + ionedge, order, u) * bernPoly(0, order, v) * bernPoly(order, order, w);
         break;
      case 9:
         res = bernPoly(1 + ionedge, order, v) * bernPoly(order, order, u) * bernPoly(order, order, w);
         break;
      case 10:
         res = bernPoly(1 + ionedge, order, u) * bernPoly(order, order, v) * bernPoly(order, order, w);
         break;
      case 11:
         res = bernPoly(1 + ionedge, order, v) * bernPoly(0, order, u) * bernPoly(order, order, w);
         break;
      default:
         res = 0.;
         break;
   }

   // Correction pour les Hex:
   if (geo_appro->getEntity()->getType() == mEntity::HEX && edge < 4)
   {
      res *= bernPoly(0, order, w);
   }

   return;
}

void xApproxFunctionScalarEdgeBernsteinQH::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                        xtensor::xVector<>& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);
   bool tridim{false};
   if (geo_appro->getEntity()->getLevel() == 3) tridim = true;

   //   -1      0      1
   //    *------*------*
   //    a      b      c
   // bernPoly(0, order  , ) associe au noeud a
   // bernPoly(order, order  , ) associe au noeud c

   // On corrige le cas de l'edge (sinon, on a v et w nuls et donc la fonction de forme aussi)
   if (geo_appro->getEntity()->getType() == mEntity::EDGE)
   {
      v = -1.;
      w = -1.;
   }

   switch (edge)
   {
      case 0:  // 2D et 3D ici...
         res[0] = dBernPoly(1 + ionedge, order, u) * bernPoly(0, order, v);
         res[1] = bernPoly(1 + ionedge, order, u) * dBernPoly(0, order, v);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= bernPoly(0, order, w);
            res[1] *= bernPoly(0, order, w);
            res[2] = bernPoly(1 + ionedge, order, u) * bernPoly(0, order, v) * dBernPoly(0, order, w);
         }

         break;
      case 1:
         res[0] = bernPoly(1 + ionedge, order, v) * dBernPoly(order, order, u);
         res[1] = dBernPoly(1 + ionedge, order, v) * bernPoly(order, order, u);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= bernPoly(0, order, w);
            res[1] *= bernPoly(0, order, w);
            res[2] = bernPoly(1 + ionedge, order, v) * bernPoly(order, order, u) * dBernPoly(0, order, w);
         }

         break;
      case 2:
         res[0] = dBernPoly(1 + ionedge, order, u) * bernPoly(order, order, v);
         res[1] = bernPoly(1 + ionedge, order, u) * dBernPoly(order, order, v);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= bernPoly(0, order, w);
            res[1] *= bernPoly(0, order, w);
            res[2] = bernPoly(1 + ionedge, order, u) * bernPoly(order, order, v) * dBernPoly(0, order, w);
         }

         break;
      case 3:
         res[0] = bernPoly(1 + ionedge, order, v) * dBernPoly(0, order, u);
         res[1] = dBernPoly(1 + ionedge, order, v) * bernPoly(0, order, u);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= bernPoly(0, order, w);
            res[1] *= bernPoly(0, order, w);
            res[2] = bernPoly(1 + ionedge, order, v) * bernPoly(0, order, u) * dBernPoly(0, order, w);
         }

         break;
      case 4:  // Only 3D here...
         res[0] = bernPoly(1 + ionedge, order, w) * dBernPoly(0, order, u) * bernPoly(0, order, v);
         res[1] = bernPoly(1 + ionedge, order, w) * bernPoly(0, order, u) * dBernPoly(0, order, v);
         res[2] = dBernPoly(1 + ionedge, order, w) * bernPoly(0, order, u) * bernPoly(0, order, v);
         break;
      case 5:
         res[0] = bernPoly(1 + ionedge, order, w) * dBernPoly(order, order, u) * bernPoly(0, order, v);
         res[1] = bernPoly(1 + ionedge, order, w) * bernPoly(order, order, u) * dBernPoly(0, order, v);
         res[2] = dBernPoly(1 + ionedge, order, w) * bernPoly(order, order, u) * bernPoly(0, order, v);
         break;
      case 6:
         res[0] = bernPoly(1 + ionedge, order, w) * dBernPoly(order, order, u) * bernPoly(order, order, v);
         res[1] = bernPoly(1 + ionedge, order, w) * bernPoly(order, order, u) * dBernPoly(order, order, v);
         res[2] = dBernPoly(1 + ionedge, order, w) * bernPoly(order, order, u) * bernPoly(order, order, v);
         break;
      case 7:
         res[0] = bernPoly(1 + ionedge, order, w) * dBernPoly(0, order, u) * bernPoly(order, order, v);
         res[1] = bernPoly(1 + ionedge, order, w) * bernPoly(0, order, u) * dBernPoly(order, order, v);
         res[2] = dBernPoly(1 + ionedge, order, w) * bernPoly(0, order, u) * bernPoly(order, order, v);
         break;
      case 8:
         res[0] = dBernPoly(1 + ionedge, order, u) * bernPoly(0, order, v) * bernPoly(order, order, w);
         res[1] = bernPoly(1 + ionedge, order, u) * dBernPoly(0, order, v) * bernPoly(order, order, w);
         res[2] = bernPoly(1 + ionedge, order, u) * bernPoly(0, order, v) * dBernPoly(order, order, w);
         break;
      case 9:
         res[0] = bernPoly(1 + ionedge, order, v) * dBernPoly(order, order, u) * bernPoly(order, order, w);
         res[1] = dBernPoly(1 + ionedge, order, v) * bernPoly(order, order, u) * bernPoly(order, order, w);
         res[2] = bernPoly(1 + ionedge, order, v) * bernPoly(order, order, u) * dBernPoly(order, order, w);
         break;
      case 10:
         res[0] = dBernPoly(1 + ionedge, order, u) * bernPoly(order, order, v) * bernPoly(order, order, w);
         res[1] = bernPoly(1 + ionedge, order, u) * dBernPoly(order, order, v) * bernPoly(order, order, w);
         res[2] = bernPoly(1 + ionedge, order, u) * bernPoly(order, order, v) * dBernPoly(order, order, w);
         break;
      case 11:
         res[0] = bernPoly(1 + ionedge, order, v) * dBernPoly(0, order, u) * bernPoly(order, order, w);
         res[1] = dBernPoly(1 + ionedge, order, v) * bernPoly(0, order, u) * bernPoly(order, order, w);
         res[2] = bernPoly(1 + ionedge, order, v) * bernPoly(0, order, u) * dBernPoly(order, order, w);
         break;
      default:
         res[0] = 0.;
         res[1] = 0.;
         res[2] = 0.;
         break;
   }
}

void xApproxFunctionScalarEdgeBernsteinQH::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                   xtensor::xVector<>& res) const
{
   getGradLocal(geo_appro, geo_integ, res);
   res = (geo_appro)->PushBack(res);
}

void xApproxFunctionScalarQuadBernstein::setNodeIndices(int order, int ionface, std::vector<std::vector<int>>& indices)
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

xApproxFunctionScalarQuadBernstein::xApproxFunctionScalarQuadBernstein(int _order, int _face, int _ionface)
    : order(_order), ionface(_ionface), face(_face)
{
   setNodeIndices(order, ionface, indices);
}

void xApproxFunctionScalarQuadBernstein::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);

   switch (face)
   {
      case 0:  // 2D et 3D ici...
         res = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, v);

         break;
      case 1:  // Only 3D here...
         res = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, w) * bernPoly(0, order, v);
         break;
      case 2:
         res = bernPoly(indices[ionface][0], order, v) * bernPoly(indices[ionface][1], order, w) * bernPoly(order, order, u);
         break;
      case 3:
         res = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, w) * bernPoly(order, order, v);
         break;
      case 4:
         res = bernPoly(indices[ionface][0], order, v) * bernPoly(indices[ionface][1], order, w) * bernPoly(0, order, u);
         break;
      case 5:
         res = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, v) * bernPoly(order, order, w);
         break;
      default:
         res = 0.;
         break;
   }

   // Correction pour les Hex:
   if (geo_appro->getEntity()->getType() == mEntity::HEX && face < 1)
   {
      res *= bernPoly(0, order, w);
   }

   return;
}

void xApproxFunctionScalarQuadBernstein::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
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
         res[0] = dBernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, v);
         res[1] = bernPoly(indices[ionface][0], order, u) * dBernPoly(indices[ionface][1], order, v);
         res[2] = 0.;

         if (tridim)
         {
            res[0] *= bernPoly(0, order, w);
            res[1] *= bernPoly(0, order, w);
            res[2] = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, v) * dBernPoly(0, order, w);
         }
         break;
      case 1:  // Only 3D here...
         res[0] = dBernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, w) * bernPoly(0, order, v);
         res[1] = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, w) * dBernPoly(0, order, v);
         res[2] = bernPoly(indices[ionface][0], order, u) * dBernPoly(indices[ionface][1], order, w) * bernPoly(0, order, v);
         break;
      case 2:
         res[0] = bernPoly(indices[ionface][0], order, v) * bernPoly(indices[ionface][1], order, w) * dBernPoly(order, order, u);
         res[1] = dBernPoly(indices[ionface][0], order, v) * bernPoly(indices[ionface][1], order, w) * bernPoly(order, order, u);
         res[2] = bernPoly(indices[ionface][0], order, v) * dBernPoly(indices[ionface][1], order, w) * bernPoly(order, order, u);
         break;
      case 3:
         res[0] = dBernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, w) * bernPoly(order, order, v);
         res[1] = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, w) * dBernPoly(order, order, v);
         res[2] = bernPoly(indices[ionface][0], order, u) * dBernPoly(indices[ionface][1], order, w) * bernPoly(order, order, v);
         break;
      case 4:
         res[0] = bernPoly(indices[ionface][0], order, v) * bernPoly(indices[ionface][1], order, w) * dBernPoly(0, order, u);
         res[1] = dBernPoly(indices[ionface][0], order, v) * bernPoly(indices[ionface][1], order, w) * bernPoly(0, order, u);
         res[2] = bernPoly(indices[ionface][0], order, v) * dBernPoly(indices[ionface][1], order, w) * bernPoly(0, order, u);
         break;
      case 5:
         res[0] = dBernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, v) * bernPoly(order, order, w);
         res[1] = bernPoly(indices[ionface][0], order, u) * dBernPoly(indices[ionface][1], order, v) * bernPoly(order, order, w);
         res[2] = bernPoly(indices[ionface][0], order, u) * bernPoly(indices[ionface][1], order, v) * dBernPoly(order, order, w);
         break;
      default:
         res = 0.;
         break;
   }

   return;
}

void xApproxFunctionScalarQuadBernstein::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                 xtensor::xVector<>& res) const
{
   getGradLocal(geo_appro, geo_integ, res);
   res = (geo_appro)->PushBack(res);
}

//  std::vector<std::vector<int> >  xApproxFunctionScalarHexBernstein::indices;

void xApproxFunctionScalarHexBernstein::setNodeIndices(int order, int ionhex, std::vector<std::vector<int>>& indices)
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

xApproxFunctionScalarHexBernstein::xApproxFunctionScalarHexBernstein(int _order, int _ionhex) : order(_order), ionhex(_ionhex)
{
   setNodeIndices(order, ionhex, indices);
}

void xApproxFunctionScalarHexBernstein::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   xtensor::xPoint uvwLoc = geo_appro->getUVW();
   res = bernPoly(indices[ionhex][0], order, uvwLoc(0)) * bernPoly(indices[ionhex][1], order, uvwLoc(1)) *
         bernPoly(indices[ionhex][2], order, uvwLoc(2));
   return;
}

void xApproxFunctionScalarHexBernstein::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                     xtensor::xVector<>& res) const
{
   xtensor::xPoint uvwLoc = geo_appro->getUVW();

   res[0] = dBernPoly(indices[ionhex][0], order, uvwLoc(0)) * bernPoly(indices[ionhex][1], order, uvwLoc(1)) *
            bernPoly(indices[ionhex][2], order, uvwLoc(2));
   res[1] = dBernPoly(indices[ionhex][1], order, uvwLoc(1)) * bernPoly(indices[ionhex][0], order, uvwLoc(0)) *
            bernPoly(indices[ionhex][2], order, uvwLoc(2));
   res[2] = dBernPoly(indices[ionhex][2], order, uvwLoc(2)) * bernPoly(indices[ionhex][0], order, uvwLoc(0)) *
            bernPoly(indices[ionhex][1], order, uvwLoc(1));

   return;
}

void xApproxFunctionScalarHexBernstein::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                xtensor::xVector<>& res) const
{
   getGradLocal(geo_appro, geo_integ, res);
   res = (geo_appro)->PushBack(res);
}

std::function<double(double)> bernPolys[136] = {
    [](double u) -> double { return 1.; },
    [](double u) -> double { return -0.5 * u + 0.5; },
    [](double u) -> double { return 0.5 * u + 0.5; },
    [](double u) -> double { return 0.25 * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -0.5 * u * u + 0.5; },
    [](double u) -> double { return 0.25 * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.125 * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.375 * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return 0.375 * (-u + 1.0) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.125 * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -0.25 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return 0.375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.25 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.03125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.15625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return -0.3125 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.3125 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.15625 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.03125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.015625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -0.09375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return 0.234375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.3125 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.234375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.09375 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.015625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.0078125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.0546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return -0.1640625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.2734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.2734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.1640625 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0546875 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0078125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double {
       return 0.00390625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -0.03125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return 0.109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.21875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.2734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.21875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.109375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.03125 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00390625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.001953125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.017578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -0.0703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.1640625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.24609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.24609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.1640625 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0703125 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.017578125 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.001953125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0009765625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.);
    },
    [](double u) -> double {
       return -0.009765625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.0);
    },
    [](double u) -> double {
       return 0.0439453125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -0.1171875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.205078125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -0.24609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.205078125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -0.1171875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.0439453125 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.009765625 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.0009765625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -0.00048828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.00537109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -0.02685546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.08056640625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.1611328125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.2255859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.2255859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.1611328125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.08056640625 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.02685546875 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00537109375 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00048828125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000244140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -0.0029296875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return 0.01611328125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0537109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.120849609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.193359375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.2255859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.193359375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.120849609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0537109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.01611328125 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0029296875 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000244140625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0001220703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.0015869140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -0.009521484375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.034912109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0872802734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.1571044921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.20947265625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.20947265625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.1571044921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0872802734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.034912109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.009521484375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0015869140625 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0001220703125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00006103515625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -0.0008544921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return 0.00555419921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.022216796875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.06109619140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.1221923828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.18328857421875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.20947265625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.18328857421875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.1221923828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.06109619140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.022216796875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00555419921875 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0008544921875 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00006103515625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000030517578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.000457763671875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -0.003204345703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.013885498046875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.041656494140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.091644287109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.152740478515625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.196380615234375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.196380615234375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.152740478515625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.091644287109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.041656494140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.013885498046875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.003204345703125 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000457763671875 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000030517578125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    }};

std::function<double(double)> dBernPolys[136] = {
    [](double u) -> double { return 0.; },
    [](double u) -> double { return -0.50000000000000000000; },
    [](double u) -> double { return 0.50000000000000000000; },
    [](double u) -> double { return 0.5 * u - 0.5; },
    [](double u) -> double { return -u; },
    [](double u) -> double { return 0.5 * u + 0.5; },
    [](double u) -> double { return -0.375 * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 1.125 * u * u - 0.75 * u - 0.375; },
    [](double u) -> double { return -1.125 * u * u - 0.75 * u + 0.375; },
    [](double u) -> double { return 0.375 * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.25 * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -(u - 1.) * (u - 1.) * (u + 0.5); },
    [](double u) -> double { return 1.5 * u * (u * u - 1.0); },
    [](double u) -> double { return (-u + 0.5) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.25 * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.15625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.15625 * (u - 1.) * (u - 1.) * (u - 1.) * (5.0 * u + 3.0); },
    [](double u) -> double { return -1.5625 * u * u * u * u + 1.25 * u * u * u + 1.875 * u * u - 1.25 * u - 0.3125; },
    [](double u) -> double { return 0.3125 * (u - 1.0) * (u + 1.) * (u + 1.) * (5.0 * u - 1.0); },
    [](double u) -> double { return 0.15625 * (-5.0 * u + 3.0) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.15625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.09375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -0.1875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (3.0 * u + 2.0); },
    [](double u) -> double { return 0.46875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0) * (3.0 * u + 1.0); },
    [](double u) -> double { return 1.875 * u * (-u * u * u * u + 2.0 * u * u - 1.0); },
    [](double u) -> double { return 0.46875 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (3.0 * u - 1.0); },
    [](double u) -> double { return 0.1875 * (-3.0 * u + 2.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.09375 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.0546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.0546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (7.0 * u + 5.0); },
    [](double u) -> double { return -0.1640625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0) * (7.0 * u + 3.0); },
    [](double u) -> double { return 0.2734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (7.0 * u + 1.0); },
    [](double u) -> double { return 0.2734375 * (-7.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.1640625 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (7.0 * u - 3.0); },
    [](double u) -> double { return 0.0546875 * (-7.0 * u + 5.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0546875 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.03125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double {
       return -0.0625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (4.0 * u + 3.0);
    },
    [](double u) -> double {
       return 0.4375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0) * (2.0 * u + 1.0);
    },
    [](double u) -> double {
       return -0.4375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (4.0 * u + 1.0);
    },
    [](double u) -> double { return 2.1875 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double {
       return 0.4375 * (-4.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.4375 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (2.0 * u - 1.0);
    },
    [](double u) -> double {
       return 0.0625 * (-4.0 * u + 3.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double { return 0.03125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double {
       return -0.017578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.017578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (9.0 * u + 7.0);
    },
    [](double u) -> double {
       return -0.0703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0) * (9.0 * u + 5.0);
    },
    [](double u) -> double {
       return 0.4921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (3.0 * u + 1.0);
    },
    [](double u) -> double {
       return -0.24609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (9.0 * u + 1.0);
    },
    [](double u) -> double {
       return 0.24609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (9.0 * u - 1.0);
    },
    [](double u) -> double {
       return 0.4921875 * (-3.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0703125 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (9.0 * u - 5.0);
    },
    [](double u) -> double {
       return 0.017578125 * (-9.0 * u + 7.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.017578125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.009765625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -0.01953125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (5.0 * u + 4.0);
    },
    [](double u) -> double {
       return 0.087890625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0) *
              (5.0 * u + 3.0);
    },
    [](double u) -> double {
       return -0.234375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (5.0 * u + 2.0);
    },
    [](double u) -> double {
       return 0.41015625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (5.0 * u + 1.0);
    },
    [](double u) -> double {
       return -2.4609375 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.41015625 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (5.0 * u - 1.0);
    },
    [](double u) -> double {
       return 0.234375 * (-5.0 * u + 2.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.087890625 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (5.0 * u - 3.0);
    },
    [](double u) -> double {
       return 0.01953125 * (-5.0 * u + 4.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.009765625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.00537109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.);
    },
    [](double u) -> double {
       return 0.00537109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (11.0 * u + 9.0);
    },
    [](double u) -> double {
       return -0.02685546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0) *
              (11.0 * u + 7.0);
    },
    [](double u) -> double {
       return 0.08056640625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (11.0 * u + 5.0);
    },
    [](double u) -> double {
       return -0.1611328125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (11.0 * u + 3.0);
    },
    [](double u) -> double {
       return 0.2255859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (11.0 * u + 1.0);
    },
    [](double u) -> double {
       return 0.2255859375 * (-11.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.1611328125 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (11.0 * u - 3.0);
    },
    [](double u) -> double {
       return 0.08056640625 * (-11.0 * u + 5.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.02685546875 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (11.0 * u - 7.0);
    },
    [](double u) -> double {
       return 0.00537109375 * (-11.0 * u + 9.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00537109375 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.0029296875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -0.005859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (6.0 * u + 5.0);
    },
    [](double u) -> double {
       return 0.064453125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.0) * (3.0 * u + 2.0);
    },
    [](double u) -> double {
       return -0.322265625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (2.0 * u + 1.0);
    },
    [](double u) -> double {
       return 0.4833984375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (3.0 * u + 1.0);
    },
    [](double u) -> double {
       return -0.38671875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (6.0 * u + 1.0);
    },
    [](double u) -> double {
       return 2.70703125 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.38671875 * (-6.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.4833984375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (3.0 * u - 1.0);
    },
    [](double u) -> double {
       return 0.322265625 * (-2.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.064453125 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (3.0 * u - 2.0);
    },
    [](double u) -> double {
       return 0.005859375 * (-6.0 * u + 5.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0029296875 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0015869140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.0015869140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (13.0 * u + 11.0);
    },
    [](double u) -> double {
       return -0.009521484375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.0) * (13.0 * u + 9.0);
    },
    [](double u) -> double {
       return 0.034912109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (13.0 * u + 7.0);
    },
    [](double u) -> double {
       return -0.0872802734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (13.0 * u + 5.0);
    },
    [](double u) -> double {
       return 0.1571044921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (13.0 * u + 3.0);
    },
    [](double u) -> double {
       return -0.20947265625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (13.0 * u + 1.0);
    },
    [](double u) -> double {
       return 0.20947265625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (13.0 * u - 1.0);
    },
    [](double u) -> double {
       return 0.1571044921875 * (-13.0 * u + 3.0) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0872802734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (13.0 * u - 5.0);
    },
    [](double u) -> double {
       return 0.034912109375 * (-13.0 * u + 7.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.009521484375 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (13.0 * u - 9.0);
    },
    [](double u) -> double {
       return 0.0015869140625 * (-13.0 * u + 11.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0015869140625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0008544921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -0.001708984375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (7.0 * u + 6.0);
    },
    [](double u) -> double {
       return 0.0111083984375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.0) * (7.0 * u + 5.0);
    },
    [](double u) -> double {
       return -0.04443359375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (7.0 * u + 4.0);
    },
    [](double u) -> double {
       return 0.1221923828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (7.0 * u + 3.0);
    },
    [](double u) -> double {
       return -0.244384765625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (7.0 * u + 2.0);
    },
    [](double u) -> double {
       return 0.3665771484375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (7.0 * u + 1.0);
    },
    [](double u) -> double {
       return -2.9326171875 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.3665771484375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (7.0 * u - 1.0);
    },
    [](double u) -> double {
       return 0.244384765625 * (-7.0 * u + 2.0) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.1221923828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (7.0 * u - 3.0);
    },
    [](double u) -> double {
       return 0.04443359375 * (-7.0 * u + 4.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0111083984375 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (7.0 * u - 5.0);
    },
    [](double u) -> double {
       return 0.001708984375 * (-7.0 * u + 6.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0008544921875 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000457763671875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.000457763671875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (15.0 * u + 13.0);
    },
    [](double u) -> double {
       return -0.003204345703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0) * (15.0 * u + 11.0);
    },
    [](double u) -> double {
       return 0.041656494140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (5.0 * u + 3.0);
    },
    [](double u) -> double {
       return -0.041656494140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (15.0 * u + 7.0);
    },
    [](double u) -> double {
       return 0.458221435546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (3.0 * u + 1.0);
    },
    [](double u) -> double {
       return -0.458221435546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (5.0 * u + 1.0);
    },
    [](double u) -> double {
       return 0.196380615234375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (15.0 * u + 1.0);
    },
    [](double u) -> double {
       return 0.196380615234375 * (-15.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.458221435546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (5.0 * u - 1.0);
    },
    [](double u) -> double {
       return 0.458221435546875 * (-3.0 * u + 1.0) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.041656494140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (15.0 * u - 7.0);
    },
    [](double u) -> double {
       return 0.041656494140625 * (-5.0 * u + 3.0) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.003204345703125 * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (15.0 * u - 11.0);
    },
    [](double u) -> double {
       return 0.000457763671875 * (-15.0 * u + 13.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000457763671875 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    }};

std::function<double(double)> ScaledbernPolys[136] = {
    [](double u) -> double { return 1.; },
    [](double u) -> double { return -0.25 * u + 0.25; },
    [](double u) -> double { return 0.25 * u + 0.25; },
    [](double u) -> double { return 0.0625 * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -0.125 * u * u + 0.125; },
    [](double u) -> double { return 0.0625 * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.015625 * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.046875 * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return 0.046875 * (-u + 1.0) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.015625 * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.00390625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -0.015625 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return 0.0234375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.015625 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.00390625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.0009765625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.0048828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return -0.009765625 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.009765625 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0048828125 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0009765625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.000244140625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -0.00146484375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return 0.003662109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.0048828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.003662109375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.00146484375 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.000244140625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double {
       return -0.00006103515625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 0.00042724609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -0.00128173828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00213623046875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.00213623046875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00128173828125 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00042724609375 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00006103515625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0000152587890625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -0.0001220703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return 0.00042724609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0008544921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.001068115234375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0008544921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00042724609375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0001220703125 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0000152587890625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -3.814697265625e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.);
    },
    [](double u) -> double {
       return 0.000034332275390625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.0);
    },
    [](double u) -> double {
       return -0.0001373291015625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.0003204345703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -0.00048065185546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.00048065185546875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -0.0003204345703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.0001373291015625 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 0.000034332275390625 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 3.814697265625e-6 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return 9.5367431640625e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -9.5367431640625e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return 0.00004291534423828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.00011444091796875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0002002716064453125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000240325927734375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0002002716064453125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.00011444091796875 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00004291534423828125 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 9.5367431640625e-6 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 9.5367431640625e-7 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -2.384185791015625e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 2.6226043701171875e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -0.0000131130218505859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0000393390655517578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000078678131103515625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000110149383544921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000110149383544921875 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000078678131103515625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0000393390655517578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0000131130218505859375 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 2.6226043701171875e-6 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 2.384185791015625e-7 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 5.9604644775390625e-8 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -7.152557373046875e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return 3.93390655517578125e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0000131130218505859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000029504299163818359375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000047206878662109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0000550746917724609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000047206878662109375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000029504299163818359375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.0000131130218505859375 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 3.93390655517578125e-6 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 7.152557373046875e-7 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 5.9604644775390625e-8 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -1.490116119384765625e-8 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 1.9371509552001953125e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -1.1622905731201171875e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 4.2617321014404296875e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000010654330253601074219 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000019177794456481933594 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000025570392608642578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000025570392608642578125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000019177794456481933594 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000010654330253601074219 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -4.2617321014404296875e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 1.1622905731201171875e-6 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 1.9371509552001953125e-7 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 1.490116119384765625e-8 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 3.7252902984619140625e-9 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -5.2154064178466796875e-8 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return 3.3900141716003417969e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -1.3560056686401367188e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 3.7290155887603759766e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -7.4580311775207519531e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00001118704676628112793 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000012785196304321289063 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00001118704676628112793 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -7.4580311775207519531e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 3.7290155887603759766e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -1.3560056686401367188e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 3.3900141716003417969e-7 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 5.2154064178466796875e-8 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 3.7252902984619140625e-9 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -9.3132257461547851563e-10 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return 1.3969838619232177734e-8 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -9.7788870334625244141e-8 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 4.2375177145004272461e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -1.2712553143501281738e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 2.7967616915702819824e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -4.6612694859504699707e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 5.993060767650604248e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -5.993060767650604248e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 4.6612694859504699707e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -2.7967616915702819824e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 1.2712553143501281738e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -4.2375177145004272461e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 9.7788870334625244141e-8 * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 1.3969838619232177734e-8 * (-u + 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 9.3132257461547851563e-10 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    }};

std::function<double(double)> dScaledBernPolys[136] = {
    [](double u) -> double { return 0.; },
    [](double u) -> double { return -0.25000000000000000000; },
    [](double u) -> double { return 0.25000000000000000000; },
    [](double u) -> double { return 0.125 * u - 0.125; },
    [](double u) -> double { return -0.25 * u; },
    [](double u) -> double { return 0.125 * u + 0.125; },
    [](double u) -> double { return -0.046875 * (u - 1.) * (u - 1.); },
    [](double u) -> double { return (0.140625 * u + 0.046875) * (u - 1.0); },
    [](double u) -> double { return (-0.140625 * u + 0.046875) * (u + 1.0); },
    [](double u) -> double { return 0.046875 * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.015625 * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -(0.0625 * u + 0.03125) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return 0.09375 * u * (u * u - 1.0); },
    [](double u) -> double { return (-0.0625 * u + 0.03125) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.015625 * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.0048828125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return (0.0244140625 * u + 0.0146484375) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double {
       return -0.048828125 * u * u * u * u + 0.0390625 * u * u * u + 0.05859375 * u * u - 0.0390625 * u - 0.009765625;
    },
    [](double u) -> double { return (0.048828125 * u - 0.009765625) * (u - 1.0) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return (-0.0244140625 * u + 0.0146484375) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.0048828125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.00146484375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return -(0.0087890625 * u + 0.005859375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double { return (0.02197265625 * u + 0.00732421875) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0); },
    [](double u) -> double { return u * (-0.029296875 * u * u * u * u + 0.05859375 * u * u - 0.029296875); },
    [](double u) -> double { return (0.02197265625 * u - 0.00732421875) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return (-0.0087890625 * u + 0.005859375) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return 0.00146484375 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double { return -0.00042724609375 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.); },
    [](double u) -> double {
       return (0.00299072265625 * u + 0.00213623046875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(0.00897216796875 * u + 0.00384521484375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return (0.01495361328125 * u + 0.00213623046875) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.01495361328125 * u + 0.00213623046875) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.00897216796875 * u - 0.00384521484375) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.00299072265625 * u + 0.00213623046875) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double { return 0.00042724609375 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double {
       return 0.0001220703125 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(0.0009765625 * u + 0.000732421875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (0.00341796875 * u + 0.001708984375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -(0.0068359375 * u + 0.001708984375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double { return 0.008544921875 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.); },
    [](double u) -> double {
       return (-0.0068359375 * u + 0.001708984375) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.00341796875 * u - 0.001708984375) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.0009765625 * u + 0.000732421875) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.0001220703125 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.000034332275390625 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (0.000308990478515625 * u + 0.000240325927734375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(0.0012359619140625 * u + 0.0006866455078125) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.0);
    },
    [](double u) -> double {
       return (0.0028839111328125 * u + 0.0009613037109375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -(0.00432586669921875 * u + 0.00048065185546875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return (0.00432586669921875 * u - 0.00048065185546875) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return (-0.0028839111328125 * u + 0.0009613037109375) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return (0.0012359619140625 * u - 0.0006866455078125) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return (-0.000308990478515625 * u + 0.000240325927734375) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.000034332275390625 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 9.5367431640625e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.);
    },
    [](double u) -> double {
       return -(9.5367431640625e-5 * u + 0.0000762939453125) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (0.000429153442382813 * u + 0.0002574920654296875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -(0.0011444091796875 * u + 0.000457763671875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.00200271606445313 * u + 0.000400543212890625) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.00240325927734375 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.00200271606445313 * u - 0.000400543212890625) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.0011444091796875 * u + 0.000457763671875) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000429153442382813 * u - 0.0002574920654296875) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-9.5367431640625e-5 * u + 0.0000762939453125) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 9.5367431640625e-6 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.);
    },
    [](double u) -> double {
       return -2.6226043701171875e-6 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (2.88486480712891e-5 * u + 0.0000236034393310546875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(0.000144243240356445 * u + 0.0000917911529541015625) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return (0.000432729721069336 * u + 0.0001966953277587890625) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -(0.000865459442138672 * u + 0.000236034393310546875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.00121164321899414 * u + 0.000110149383544921875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.00121164321899414 * u + 0.000110149383544921875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000865459442138672 * u - 0.000236034393310546875) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.000432729721069336 * u + 0.0001966953277587890625) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000144243240356445 * u - 0.0000917911529541015625) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-2.88486480712891e-5 * u + 0.0000236034393310546875) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 2.6226043701171875e-6 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 7.152557373046875e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(8.58306884765625e-6 * u + 7.152557373046875e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (4.72068786621094e-5 * u + 0.00003147125244140625) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -(0.000157356262207031 * u + 0.000078678131103515625) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.00035405158996582 * u + 0.0001180171966552734375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -(0.000566482543945313 * u + 0.00009441375732421875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 0.00066089630126953125 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.000566482543945313 * u + 0.00009441375732421875) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.00035405158996582 * u - 0.0001180171966552734375) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.000157356262207031 * u + 0.000078678131103515625) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (4.72068786621094e-5 * u - 0.00003147125244140625) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-8.58306884765625e-6 * u + 7.152557373046875e-6) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 7.152557373046875e-7 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -1.9371509552001953125e-7 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (2.51829624176025e-6 * u + 2.1308660507202148438e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(1.51097774505615e-5 * u + 0.000010460615158081054688) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return (5.54025173187256e-5 * u + 0.000029832124710083007813) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -(0.000138506293296814 * u + 0.000053271651268005371094) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000249311327934265 * u + 0.000057533383369445800781) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -(0.000332415103912354 * u + 0.000025570392608642578125) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000332415103912354 * u - 0.000025570392608642578125) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.000249311327934265 * u + 0.000057533383369445800781) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000138506293296814 * u - 0.000053271651268005371094) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-5.54025173187256e-5 * u + 0.000029832124710083007813) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (1.51097774505615e-5 * u - 0.000010460615158081054688) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-2.51829624176025e-6 * u + 2.1308660507202148438e-6) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 1.9371509552001953125e-7 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 5.2154064178466796875e-8 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(7.30156898498535e-7 * u + 6.258487701416015625e-7) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (4.74601984024048e-6 * u + 3.3900141716003417969e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return -(1.89840793609619e-5 * u + 0.00001084804534912109375) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (5.22062182426453e-5 * u + 0.000022374093532562255859) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -(0.000104412436485291 * u + 0.000029832124710083007813) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000156618654727936 * u + 0.000022374093532562255859) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -0.00017899274826049804688 * u * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (0.000156618654727936 * u - 0.000022374093532562255859) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-0.000104412436485291 * u + 0.000029832124710083007813) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (5.22062182426453e-5 * u - 0.000022374093532562255859) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-1.89840793609619e-5 * u + 0.00001084804534912109375) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (4.74601984024048e-6 * u - 3.3900141716003417969e-6) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-7.30156898498535e-7 * u + 6.258487701416015625e-7) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 5.2154064178466796875e-8 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -1.3969838619232177734e-8 * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return (2.09547579288483e-7 * u + 1.8160790205001831055e-7) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.);
    },
    [](double u) -> double {
       return -(1.46683305501938e-6 * u + 1.0756775736808776855e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.0);
    },
    [](double u) -> double {
       return (6.35627657175064e-6 * u + 3.8137659430503845215e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -(1.90688297152519e-5 * u + 8.8987872004508972168e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (4.19514253735542e-5 * u + 0.000013983808457851409912) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return -(6.9919042289257e-5 * u + 0.000013983808457851409912) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (8.98959115147591e-5 * u + 5.993060767650604248e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-8.98959115147591e-5 * u + 5.993060767650604248e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (6.9919042289257e-5 * u - 0.000013983808457851409912) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-4.19514253735542e-5 * u + 0.000013983808457851409912) * (u - 1.) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (1.90688297152519e-5 * u - 8.8987872004508972168e-6) * (u - 1.) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-6.35627657175064e-6 * u + 3.8137659430503845215e-6) * (u - 1.) * (u - 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (1.46683305501938e-6 * u - 1.0756775736808776855e-6) * (u - 1.0) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return (-2.09547579288483e-7 * u + 1.8160790205001831055e-7) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    },
    [](double u) -> double {
       return 1.3969838619232177734e-8 * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) *
              (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.) * (u + 1.);
    }};

inline double binomial(int n, int k) { return boost::math::binomial_coefficient<double>(n, k); }

inline double bernPoly(int v, int n, double u)
{
   //     return pow(2,-n)*binomial(n,v)*pow(1.+u,v)*pow(1.-u,n-v);
   //     return bernPolys[n*(n+1)/2 + v](u);
   return ScaledbernPolys[n * (n + 1) / 2 + v](u);
}

inline double dBernPoly(int v, int n, double u)
{
   //     if (v < 1) return -0.5*n*(bernPoly(v,n-1,u));
   //     else if (v > n-1)  return 0.5*n*(bernPoly(v-1,n-1,u));
   //     else return 0.5*n*(bernPoly(v-1,n-1,u) - bernPoly(v,n-1,u));

   //     return dBernPolys[n*(n+1)/2 + v](u);
   return dScaledBernPolys[n * (n + 1) / 2 + v](u);
}

}  // end namespace xfem
