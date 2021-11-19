/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <algorithm>
#include <iostream>
// xmapping
#include "xLagrangeMapping.h"

namespace xmapping
{
std::vector<double> GeomShapeFunctionPoint() { return std::vector<double>{1.}; }
std::vector<vector3d> GradGeomShapeFunctionPoint() { return std::vector<vector3d>{{0., 0., 0.}}; }
std::vector<double> GeomShapeFunctionLine(double u) { return std::vector<double>{(1. - u) * 0.5, (1. + u) * 0.5}; }
std::vector<vector3d> GradGeomShapeFunctionLine() { return std::vector<vector3d>{{-0.5, 0., 0.}, {0.5, 0., 0.}}; }
std::vector<double> GeomShapeFunctionTri(double u, double v) { return std::vector<double>{1. - u - v, u, v}; }
std::vector<vector3d> GradGeomShapeFunctionTri() { return std::vector<vector3d>{{-1., -1., 0.}, {1., 0., 0.}, {0., 1., 0.}}; }
std::vector<double> GeomShapeFunctionQuad(double u, double v)
{
   return std::vector<double>{0.25 * (1. - u) * (1. - v), 0.25 * (1. + u) * (1. - v),  //
                              0.25 * (1. + u) * (1. + v), 0.25 * (1. - u) * (1. + v)};
}

std::vector<vector3d> GradGeomShapeFunctionQuad(double u, double v)
{
   const double t1 = 0.25 * (v - 1.0);
   const double t2 = 0.25 * (v + 1.0);
   const double t3 = 0.25 * (u - 1.0);
   const double t4 = -0.25 * (u + 1.0);
   return std::vector<vector3d>{{t1, t3, 0.}, {-t1, t4, 0.}, {t2, -t4, 0.}, {-t2, -t3, 0.}};
}

std::vector<double> GeomShapeFunctionTet(double u, double v, double w) { return std::vector<double>{1. - u - v - w, u, v, w}; }

std::vector<vector3d> GradGeomShapeFunctionTet()
{
   return std::vector<vector3d>{{-1., -1., -1.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
}

std::vector<double> GeomShapeFunctionHex(double u, double v, double w)
{
   return std::vector<double>{0.125 * (1. - u) * (1. - v) * (1. - w), 0.125 * (1. + u) * (1. - v) * (1. - w),
                              0.125 * (1. + u) * (1. + v) * (1. - w), 0.125 * (1. - u) * (1. + v) * (1. - w),
                              0.125 * (1. - u) * (1. - v) * (1. + w), 0.125 * (1. + u) * (1. - v) * (1. + w),
                              0.125 * (1. + u) * (1. + v) * (1. + w), 0.125 * (1. - u) * (1. + v) * (1. + w)};
}

std::vector<vector3d> GradGeomShapeFunctionHex(double u, double v, double w)
{
   return std::vector<vector3d>{{-0.125 * (1. - v) * (1. - w), -0.125 * (1. - u) * (1. - w), -0.125 * (1. - u) * (1. - v)},
                                {0.125 * (1. - v) * (1. - w), -0.125 * (1. + u) * (1. - w), -0.125 * (1. + u) * (1. - v)},
                                {0.125 * (1. + v) * (1. - w), 0.125 * (1. + u) * (1. - w), -0.125 * (1. + u) * (1. + v)},
                                {-0.125 * (1. + v) * (1. - w), 0.125 * (1. - u) * (1. - w), -0.125 * (1. - u) * (1. + v)},
                                {-0.125 * (1. - v) * (1. + w), -0.125 * (1. - u) * (1. + w), 0.125 * (1. - u) * (1. - v)},
                                {0.125 * (1. - v) * (1. + w), -0.125 * (1. + u) * (1. + w), 0.125 * (1. + u) * (1. - v)},
                                {0.125 * (1. + v) * (1. + w), 0.125 * (1. + u) * (1. + w), 0.125 * (1. + u) * (1. + v)},
                                {-0.125 * (1. + v) * (1. + w), 0.125 * (1. - u) * (1. + w), 0.125 * (1. - u) * (1. + v)}};
}
std::vector<double> GeomShapeFunctionPrism(double u, double v, double w)
{
   const double k = 1. - u - v;
   return std::vector<double>{0.5 * k * (1. + w), 0.5 * u * (1. + w), 0.5 * v * (1. + w),
                              0.5 * k * (1. - w), 0.5 * u * (1. - w), 0.5 * v * (1. - w)};
}

std::vector<vector3d> GradGeomShapeFunctionPrism(double u, double v, double w)
{
   return std::vector<vector3d>{
       {-0.5 * (1. + w), -0.5 * (1. + w), 0.5 * (1. - u - v)},  {0.5 * (1. + w), 0.0, 0.5 * u}, {0.0, -0.5 * (1. + w), 0.5 * v},
       {-0.5 * (1. - w), -0.5 * (1. - w), -0.5 * (1. - u - v)}, {0.5 * (1. - w), 0., -0.5 * u}, {0.0, -0.5 * (1. - w), -0.5 * v}};
}

xLagrangeMapping::xLagrangeMapping(xReferenceElementType type, const std::vector<xtensor::xPoint> &input_knots)
    : xMapping(type), knots(input_knots)
{
}

int xLagrangeMapping::order() const { return 1; }

int xLagrangeMapping::geomOrder() const { return 1; }

void xLagrangeMapping::eval(double u, double v, double w, double &x, double &y, double &z) const
{
   x = 0.0;
   y = 0.0;
   z = 0.0;
   const auto f = GeomShapeFunction(u, v, w);
   for (size_t i = 0; i < knots.size(); i++)
   {
      xtensor::xPoint p = knots[i];
      double fct = f[i];
      x += p(0) * fct;
      y += p(1) * fct;
      z += p(2) * fct;
   }
}

void xLagrangeMapping::deval(double u, double v, double w, double &dxdu, double &dydu, double &dzdu, double &dxdv, double &dydv,
                             double &dzdv, double &dxdw, double &dydw, double &dzdw) const
{
   dxdu = 0.0;
   dydu = 0.0;
   dzdu = 0.0;
   dxdv = 0.0;
   dydv = 0.0;
   dzdv = 0.0;
   dxdw = 0.0;
   dydw = 0.0;
   dzdw = 0.0;
   std::vector<vector3d> dus = GradGeomShapeFunction(u, v, w);
   for (size_t i = 0; i < knots.size(); i++)
   {
      vector3d du(dus[i]);
      xtensor::xPoint p = knots[i];

      const double xx = p(0);
      const double yy = p(1);
      const double zz = p(2);

      dxdu += xx * du(0);
      dydu += yy * du(0);
      dzdu += zz * du(0);

      dxdv += xx * du(1);
      dydv += yy * du(1);
      dzdv += zz * du(1);

      dxdw += xx * du(2);
      dydw += yy * du(2);
      dzdw += zz * du(2);
   }
}

void xLagrangeMapping::boundingBox(xtensor::xPoint &pMin, xtensor::xPoint &pMax) const
{
   xgeom::xBoundingBox bb = boundingBox();
   pMin = bb.min;
   pMax = bb.max;
}

xgeom::xBoundingBox xLagrangeMapping::boundingBox() const
{
   xgeom::xBoundingBox bb;
   for (const auto &p : knots) bb.inclose(p);
   return bb;
}

std::vector<double> xLagrangeMapping::GeomShapeFunction(double u, double v, double w) const
{
   switch (type)
   {
      case xReferenceElementType::VERTEX:
         return GeomShapeFunctionPoint();
      case xReferenceElementType::EDGE:
         return GeomShapeFunctionLine(u);
      case xReferenceElementType::TRI:
         return GeomShapeFunctionTri(u, v);
      case xReferenceElementType::QUAD:
         return GeomShapeFunctionQuad(u, v);
      case xReferenceElementType::TET:
         return GeomShapeFunctionTet(u, v, w);
      case xReferenceElementType::HEX:
         return GeomShapeFunctionHex(u, v, w);
      case xReferenceElementType::PRISM:
         return GeomShapeFunctionPrism(u, v, w);
      default:
         std::cout << "GeomshapeFunction: Entity type not supported in" << __FILE__ << ":" << __LINE__ << std::endl;
         throw;
   }
}

std::vector<vector3d> xLagrangeMapping::GradGeomShapeFunction(double u, double v, double w) const
{
   switch (type)
   {
      case xReferenceElementType::VERTEX:
         return GradGeomShapeFunctionPoint();
      case xReferenceElementType::EDGE:
         return GradGeomShapeFunctionLine();
      case xReferenceElementType::TRI:
         return GradGeomShapeFunctionTri();
      case xReferenceElementType::QUAD:
         return GradGeomShapeFunctionQuad(u, v);
      case xReferenceElementType::TET:
         return GradGeomShapeFunctionTet();
      case xReferenceElementType::HEX:
         return GradGeomShapeFunctionHex(u, v, w);
      case xReferenceElementType::PRISM:
         return GradGeomShapeFunctionPrism(u, v, w);
      default:
         std::cout << "GradGeomshapeFunction: Entity type not supported in" << __FILE__ << ":" << __LINE__ << std::endl;
         throw;
   }
}

/*void xLagrangeMapping::normalVector(pEntity border, double u, double v, double w, vector3d &n) const
{
   // In case of a mesh mapping, the normal to an entity which is on
   // a border = sum of grad geom shape functions which node are not
   // on the border entity

   n[0] = n[1] = n[2] = 0.0;

   tensor jInv;
   double detJac = jacInverse(u, v, w, jInv);

   std::vector<pVertex> verts;
   M_GetVertices(const_cast<AOMD::mEntity *>(&ent), verts);
   std::vector<pVertex> vb;
   M_GetVertices(border, vb);

   vector3d dus[256];
   GradGeomShapeFunction(u, v, w, dus);

   for (int i = 0; i < verts.size(); i++)
   {
      std::vector<pVertex>::iterator it = std::find(vb.begin(), vb.end(), verts[i]);
      if (it == vb.end())
      {
         vector3d gr(dus[i]);
         gr *= jInv;
         n -= gr;
      }
   }
   n.norm();
}
*/

xConstantLagrangeMapping::xConstantLagrangeMapping(xReferenceElementType type, const std::vector<xtensor::xPoint> &input_knots)
    : xLagrangeMapping(type, input_knots)
{
   det = xLagrangeMapping::jacInverse(0, 0, 0, jac);
   if (det < 0) det = -det;
}

xCylindricalCoordinatesLagrangeMapping::xCylindricalCoordinatesLagrangeMapping(xReferenceElementType type,
                                                                               const std::vector<xtensor::xPoint> &input_knots)
    : xLagrangeMapping(type, input_knots)
{
}

double xConstantLagrangeMapping::detJac(double, double, double) const { return det; }

double xConstantLagrangeMapping::jacInverse(double, double, double, tensor &j) const
{
   j = jac;
   return det;
}

double xCylindricalCoordinatesLagrangeMapping::detJac(double u, double v, double w) const
{
   tensor j;
   double r, z, teta;
   eval(u, v, w, r, z, teta);
   //  printf("r = %12.5E det = %12.5E\n",r,LagrangeMapping::jacInverse(u,v,w,j) );
   return (r * xLagrangeMapping::jacInverse(u, v, w, j));
}

double xCylindricalCoordinatesLagrangeMapping::jacInverse(double u, double v, double w, tensor &j) const
{
   double r, z, teta;
   eval(u, v, w, r, z, teta);
   double det = xLagrangeMapping::jacInverse(u, v, w, j);
   return (r * det);
}

int xCylindricalCoordinatesLagrangeMapping::order() const { return 2; }

xRegularCubeLagrangeMapping::xRegularCubeLagrangeMapping(xReferenceElementType type,
                                                         const std::vector<xtensor::xPoint> &input_knots)
    : xLagrangeMapping(type, input_knots)
{
   assert(knots.size());
   double x1, y1, z1;
   xtensor::xPoint &p0 = knots[0];
   x0 = x1 = p0(0);
   y0 = y1 = p0(1);
   z0 = z1 = p0(2);
   for (size_t i = 1; i < knots.size(); i++)
   {
      xtensor::xPoint &p = knots[i];
      if (x0 > p(0)) x0 = p(0);
      if (y0 > p(1)) y0 = p(1);
      if (z0 > p(2)) z0 = p(2);
      if (x1 < p(0)) x1 = p(0);
      if (y1 < p(1)) y1 = p(1);
      if (z1 < p(2)) z1 = p(2);
   }
   dx = x1 - x0;
   dy = y1 - y0;

   //  printf("%f %f %f %f\n",x0,y0,dx,dy);

   if (type == xReferenceElementType::QUAD)
      dz = 1.0;
   else
      dz = z1 - z0;
}

void xRegularCubeLagrangeMapping::eval(double u, double v, double w, double &x, double &y, double &z) const
{
   x = x0 + .5 * (1. + u) * dx;
   y = y0 + .5 * (1. + v) * dy;
   if (type == xReferenceElementType::QUAD)
      z = 0.0;
   else
      z = z0 + .5 * (1. + w) * dz;
}

void xRegularCubeLagrangeMapping::deval([[gnu::unused]] double u, [[gnu::unused]] double v, [[gnu::unused]] double w,
                                        double &dxdu, double &dydu, double &dzdu, double &dxdv, double &dydv, double &dzdv,
                                        double &dxdw, double &dydw, double &dzdw) const
{
   dxdu = 0.5 * dx;
   dydv = 0.5 * dy;
   if (type == xReferenceElementType::QUAD)
      dzdw = 0.0;
   else
      dzdw = .5 * dz;
   dxdv = dxdw = dydu = dydw = dzdu = dzdv = 0.0;
}

bool xRegularCubeLagrangeMapping::invert(double x, double y, double z, double &u, double &v, double &w) const
{
   u = 2. * (x - x0) / dx - 1.0;
   v = 2. * (y - y0) / dy - 1.0;
   if (type == xReferenceElementType::QUAD)
      w = 0.0;
   else
      w = 2. * (z - z0) / dz - 1.0;
   return true;
}

double xRegularCubeLagrangeMapping::detJac([[gnu::unused]] double u, [[gnu::unused]] double v, [[gnu::unused]] double w) const
{
   if (type == xReferenceElementType::QUAD)
      return 0.25 * (dx * dy);
   else
      return 0.125 * (dx * dy * dz);
}

double xRegularCubeLagrangeMapping::jacInverse([[gnu::unused]] double u, [[gnu::unused]] double v, [[gnu::unused]] double w,
                                               tensor &jInv) const
{
   jInv(0, 0) = 2. / dx;
   jInv(1, 1) = 2. / dy;
   if (type == xReferenceElementType::QUAD)
      jInv(2, 2) = 1.0;
   else
      jInv(2, 2) = 2. / dz;

   jInv(1, 2) = jInv(2, 1) = jInv(2, 0) = jInv(0, 2) = jInv(0, 1) = jInv(1, 0) = 0.0;
   return detJac(0, 0, 0);
}

/*
void xRegularCubeLagrangeMapping::normalVector(pEntity border, double u, double v, double w, vector3d &n) const
{
   double u1, v1, w1;
   xtensor::xPoint p1;
   xLagrangeMapping lm(border);
   lm.COG(u1, v1, w1);
   lm.eval(u1, v1, w1, p1(0), p1(1), p1(2));

   if (fabs(p1(0) - x0) < 1.e-6 * dx)
      n = vector3d(-1, 0, 0);
   else if (fabs(p1(0) - x0 - dx) < 1.e-6 * dx)
      n = vector3d(1, 0, 0);
   else if (fabs(p1(1) - y0) < 1.e-6 * dy)
      n = vector3d(0, -1, 0);
   else if (fabs(p1(1) - y0 - dy) < 1.e-6 * dy)
      n = vector3d(0, 1, 0);
   else if (fabs(p1(2) - z0) < 1.e-6 * dz)
      n = vector3d(0, 0, -1);
   else if (fabs(p1(2) - z0 - dz) < 1.e-6 * dz)
      n = vector3d(0, 0, 1);
   else
      assert(1 == 0);
}
*/

/* double xRegularCubeLagrangeMapping::PushBack (double u, double v, double w, int vsize, vector<vector3d> &gr) const
 {
   double detJ = detJac(u,v,w);

   double xx = 2./dx;
   double yy = 2./dy;
   double zz = 2./dz;

   for(int i=0;i<vsize;i++)
     {
       gr[i](0) *= (xx);
       gr[i](1) *= (yy);
   if(M_GetElementType(ent) == HEX)
         gr[i](2) *= (zz);
     }
   return detJ;
 }*/

}  // namespace xmapping
