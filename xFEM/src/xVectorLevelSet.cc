/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include "xVectorLevelSet.h"

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

namespace xfem
{
using AOMD::AOMD_Util;
using AOMD::mEdge;
using AOMD::mVertex;
xAnalyticalVectorLevelSetCylinder::xAnalyticalVectorLevelSetCylinder(const xtensor::xPoint& p, xtensor::xVector<> d,
                                                                     const double& r)
    : point_on_axis(p), axis(d.norm()), radius(r)
{
   // assert (r > 0.);
   assert(axis.mag() > 0.);
   xtensor::xVector<> ortho1 = axis % (xtensor::xVector<>(1., 0., 0.));
   xtensor::xVector<> ortho2 = axis % (xtensor::xVector<>(0., 1., 0.));
   xtensor::xVector<> ortho3 = axis % (xtensor::xVector<>(0., 0., 1.));
   ortho = (ortho1.mag() > ortho2.mag()) ? ((ortho1.mag() > ortho3.mag()) ? ortho1 : ortho3)
                                         : ((ortho2.mag() > ortho3.mag()) ? ortho2 : ortho3);
   ortho.norm();
}

xVectorLevelSetData xAnalyticalVectorLevelSetCylinder::operator()(const xtensor::xPoint& p) const
{
   const double eps = 1.e-6;
   const xtensor::xVector<> op = xtensor::xVector<>(point_on_axis, p);
   // double d = op * axis;
   const xtensor::xVector<> ob = axis * (axis * op);
   const xtensor::xVector<> pb = -op + ob;
   const double npb = pb.mag();
   if ((npb / radius) > eps)
   {
      xtensor::xVector<> directioniso = pb * ((npb - radius)) / npb;
      int inout = ((npb - radius) > 0.) ? 1 : (((npb - radius) == 0.) ? 0 : -1);
      return xVectorLevelSetData(inout, directioniso);
   }
   else
   {
      xtensor::xVector<> directioniso = ortho * radius;
      int inout = -1;
      return xVectorLevelSetData(inout, directioniso);
   }
}

xAnalyticalVectorLevelSetSphere::xAnalyticalVectorLevelSetSphere(const xtensor::xPoint& _center, const double& _radius)
    : center(_center), radius(_radius)
{
}

xVectorLevelSetData xAnalyticalVectorLevelSetSphere::operator()(const xtensor::xPoint& p) const
{
   xtensor::xVector<> dir(center, p);
   double dist = dir.mag();
   xtensor::xPoint closestpoint(center(0) + radius, center(1), center(2));
   const double eps = 1.e-6;
   if ((dist / radius) <= eps) dist = 0.;

   int inout = ((dist / radius) > (1. + eps)) ? 1 : (((dist / radius) < (1. - eps)) ? -1 : 0);
   if (dist != 0.)
   {
      closestpoint(0) = center(0) + dir(0) / dist * radius;
      closestpoint(1) = center(1) + dir(1) / dist * radius;
      closestpoint(2) = center(2) + dir(2) / dist * radius;
   }
   xtensor::xVector<> directioniso(p, closestpoint);
   return xVectorLevelSetData(inout, directioniso);
}

xAnalyticalVectorLevelSetGrownSegment::xAnalyticalVectorLevelSetGrownSegment(const xtensor::xPoint& p1_,
                                                                             const xtensor::xPoint& p2_, const double& growth_)
    : p1(p1_),
      p2(p2_),
      sphere1(p1_, growth_),
      sphere2(p2_, growth_),
      cyl(p1_, xtensor::xVector<>(p1_, p2_), growth_),
      dir(p1_, p2_)
{
   length = dir.normValue();
   dir.norm();
}

xVectorLevelSetData xAnalyticalVectorLevelSetGrownSegment::operator()(const xtensor::xPoint& p) const
{
   xtensor::xVector<> op(p1, p);
   double proj = op * dir;
   if (proj <= 0) return sphere1(p);
   if (proj >= length) return sphere2(p);
   return cyl(p);
}
const xVectorLevelSetData* xVectorLevelSet::operator()(const AOMD::mVertex& v) const
{
   const xVectorLevelSetData* pdata = data.getData(v);
   if (pdata)
      return pdata;
   else
   {
      std::cout << "Warning : no xVectorLevelSetData attached to vertex " << &v << " file " << __FILE__ << ":" << __LINE__
                << std::endl;
      return nullptr;
   }
}

xtensor::xVector<> xVectorLevelSet::getValVector(AOMD::mEntity* e, const xtensor::xPoint& uvw) const
{
   xElement elem(e);
   elem.setUvw(uvw);
   const int nb = e->size(0);
   std::vector<xtensor::xVector<>> vects(nb);
   for (int i = 0; i < nb; ++i) vects[i] = data.at(*e->get(0, i)).getDirectionToIso();
   return elem.getInterpoVec(vects);
}

double xVectorLevelSet::getValInOut(AOMD::mEntity* e, const xtensor::xPoint& uvw) const
{
   xElement elem(e);
   elem.setUvw(uvw);
   const int nb = e->size(0);
   std::vector<double> vals(nb);
   for (int i = 0; i < nb; ++i) vals[i] = data.at(*e->get(0, i)).getInOut();
   return elem.getInterpoSca(vals);
}

xVectorLevelSet::~xVectorLevelSet() { data.clear(); }

}  // namespace xfem
