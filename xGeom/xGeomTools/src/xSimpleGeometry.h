/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XSIMPLEGEOMETRY_H
#define _XSIMPLEGEOMETRY_H

#include <cmath>
#include <functional>
#include <vector>

// xtensor
#include "xTensor2.h"
#include "xVector.h"

using xPointToDouble = std::function<double(const xtensor::xPoint&)>;

namespace xgeom
{
inline double xDistance(const xtensor::xPoint& p1, const xtensor::xPoint& p2) { return xtensor::xVector<>(p1, p2).mag(); }

inline xtensor::xPoint xPointPlusVector(const xtensor::xPoint& p, const xtensor::xVector<>& vec)
{
   return xtensor::xPoint(p(0) + vec(0), p(1) + vec(1), p(2) + vec(2));
}

class xCompl
{
  public:
   xCompl(const xPointToDouble& d) : dom(d) {}
   double operator()(const xtensor::xPoint& p) const { return (-dom(p)); }

  private:
   xPointToDouble dom;
};

class xUnion
{
  public:
   xUnion(const xPointToDouble& d1, const xPointToDouble& d2) : dom1(d1), dom2(d2) {}
   double operator()(const xtensor::xPoint& p) const { return std::min(dom1(p), dom2(p)); }

  private:
   xPointToDouble dom1;
   xPointToDouble dom2;
};

class xUnion3
{
  public:
   xUnion3(const xPointToDouble& d1, const xPointToDouble& d2, const xPointToDouble& d3)
   {
      xUnion u1(d1, d2);
      xUnion u2(u1, d3);
      result = u2;
   }
   double operator()(const xtensor::xPoint& p) const { return result(p); }

  private:
   xPointToDouble result;
};

class xUnionVector
{
  public:
   xUnionVector(const std::vector<xPointToDouble>& ds_) : ds(ds_) {}
   double operator()(const xtensor::xPoint& p) const
   {
      std::vector<double> ds_val;
      for (const auto& f : ds) ds_val.push_back(f(p));
      return *std::min_element(ds_val.begin(), ds_val.end());
   }

  private:
   std::vector<xPointToDouble> ds;
};

class xInter
{
  public:
   xInter(const xPointToDouble& d1, const xPointToDouble& d2) : dom1(d1), dom2(d2) {}
   double operator()(const xtensor::xPoint& p) const { return std::max(dom1(p), dom2(p)); }

  private:
   xPointToDouble dom1;
   xPointToDouble dom2;
};

class xInter3
{
  public:
   xInter3(const xPointToDouble& d1, const xPointToDouble& d2, const xPointToDouble& d3)
   {
      xInter u1(d1, d2);
      xInter u2(u1, d3);
      result = u2;
   }
   double operator()(const xtensor::xPoint& p) const { return result(p); }

  private:
   xPointToDouble result;
};

class xInterVector
{
  public:
   xInterVector(const std::vector<xPointToDouble>& ds_) : ds(ds_) {}
   double operator()(const xtensor::xPoint& p) const
   {
      std::vector<double> ds_val;
      for (const auto& f : ds) ds_val.push_back(f(p));
      return *std::max_element(ds_val.begin(), ds_val.end());
   }

  private:
   std::vector<xPointToDouble> ds;
};

class xMinus
{
  public:
   xMinus(const xPointToDouble& d1, const xPointToDouble& d2) : dom1(d1), dom2(d2) {}
   double operator()(const xtensor::xPoint& p) const { return std::max(dom1(p), -dom2(p)); }

  private:
   xPointToDouble dom1;
   xPointToDouble dom2;
};

class xTranslate
{
  public:
   xTranslate(const xPointToDouble& d, const xtensor::xVector<>& vec_) : dom(d), vec(vec_) {}
   double operator()(const xtensor::xPoint& p) const { return dom(xtensor::xPoint(p(0) - vec(0), p(1) - vec(1), p(2) - vec(2))); }

  private:
   xPointToDouble dom;
   xtensor::xVector<> vec;
};

class xRadialZoom
{
  public:
   xRadialZoom(const xPointToDouble& d, const xtensor::xPoint& c, const double& f) : dom(d), center(c), fac(f)
   {
      if (!(fac > 0.))
      {
         std::cerr << " zoom factor must be positive! " << std::endl;
         assert(0);
      }
   }
   double operator()(const xtensor::xPoint& p) const
   {
      return dom(xtensor::xPoint((1. / fac) * ((fac - 1.) * center(0) + p(0)), (1. / fac) * ((fac - 1.) * center(1) + p(1)),
                                 (1. / fac) * ((fac - 1.) * center(2) + p(2))));
   }

  private:
   xPointToDouble dom;
   const xtensor::xPoint center;
   const double fac;
};

class xRotate
{
  public:
   xRotate(const xPointToDouble& d, const xtensor::xPoint& center, const xtensor::xVector<>& vector, const double& angle)
       : dom(d), centre(center), axe(vector), teta(angle)
   {
      xtensor::xTensor2<> rotloc, pass, passt;
      xtensor::xVector<> vecpr;

      rotloc(0, 0) = 1.;
      rotloc(0, 1) = 0.;
      rotloc(0, 2) = 0.;
      rotloc(1, 0) = 0.;
      rotloc(1, 1) = cos(teta);
      rotloc(1, 2) = sin(teta);
      rotloc(2, 0) = 0.;
      rotloc(2, 1) = -sin(teta);
      rotloc(2, 2) = cos(teta);
      // Construction de la matrice de passage
      if (axe(0) != 0.)
      {
         vecpr(0) = axe(1) / axe(0);
         vecpr(1) = -1.;
         vecpr(2) = 0.;
      }
      else if (axe(1) != 0.)
      {
         vecpr(0) = -1.;
         vecpr(1) = axe(0) / axe(1);
         vecpr(2) = 0.;
      }
      else if (axe(2) != 0.)
      {
         vecpr(0) = 0.;
         vecpr(1) = -1.;
         vecpr(2) = axe(1) / axe(2);
      }
      else
         std::cerr << "Bad vector!\n";

      double di = sqrt((axe(0) * axe(0)) + (axe(1) * axe(1)) + (axe(2) * axe(2)));
      double dd = sqrt((vecpr(0) * vecpr(0)) + (vecpr(1) * vecpr(1)) + (vecpr(2) * vecpr(2)));
      double u = (1. / (di * dd)) * ((axe(1) * vecpr(2)) - (axe(2) * vecpr(1)));
      double v = (1. / (di * dd)) * ((axe(2) * vecpr(0)) - (axe(0) * vecpr(2)));
      double w = (1. / (di * dd)) * ((axe(0) * vecpr(1)) - (axe(1) * vecpr(0)));

      pass(0, 0) = (1 / di) * axe(0);
      pass(1, 0) = (1 / di) * axe(1);
      pass(2, 0) = (1 / di) * axe(2);
      pass(0, 1) = (1 / dd) * vecpr(0);
      pass(1, 1) = (1 / dd) * vecpr(1);
      pass(2, 1) = (1 / dd) * vecpr(2);
      pass(0, 2) = u;
      pass(1, 2) = v;
      pass(2, 2) = w;
      std::cout << "Construction de la matrice de rotation\n";
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++) passt(i, j) = pass(j, i);
      }  //	Inverse de la matrice de passage
      rot = (pass * rotloc) * passt;
   }
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> ompr, om;
      for (int i = 0; i < 3; i++)
      {
         ompr(i) = p(i) - centre(i);
      }
      om = rot * ompr;
      return dom(xtensor::xPoint(centre(0) + om(0), centre(1) + om(1), centre(2) + om(2)));
   }

  private:
   xPointToDouble dom;
   xtensor::xPoint centre;
   xtensor::xVector<> axe;
   double teta;
   xtensor::xTensor2<> rot;
};

// to compute the signed distance to a sphere
class xSphere
{
  public:
   xSphere(const xtensor::xPoint& o, const double& r) : origin(o), radius(r) {}
   double operator()(const xtensor::xPoint& p) const { return xDistance(p, origin) - radius; }

  private:
   xSphere();
   xtensor::xPoint origin;
   double radius;
};

// to compute the signed distance  to a plane
// the normal to the plane gives the positive side
class xPlane
{
  public:
   xPlane(const xtensor::xPoint& p, xtensor::xVector<> n) : point_on_plane(p), normal(n.norm()) {}
   double operator()(const xtensor::xPoint& p) const { return xtensor::xVector<>(point_on_plane, p) * normal; }

  private:
   xPlane();
   xtensor::xPoint point_on_plane;
   xtensor::xVector<> normal;
};

// Pour representer les fissures de type CCT :
class xBiPlane
{
  public:
   xBiPlane(const xtensor::xPoint& p, const double width_, xtensor::xVector<> n)
       : point_on_plane(p), normal(n.norm()), width(width_)
   {
   }
   double operator()(const xtensor::xPoint& p) const
   {
      return (std::abs(xtensor::xVector<>(point_on_plane, p) * normal) - width / 2.);
   }

  private:
   xBiPlane();
   xtensor::xPoint point_on_plane;
   xtensor::xVector<> normal;
   double width;
};

// to compute the signed distance to a cube
class xCube
{
  public:
   xCube(const xtensor::xPoint& center_, xtensor::xVector<> axis1_, xtensor::xVector<> axis2_, const double& side_)
       : center(center_), axis1(axis1_.norm()), axis2(axis2_.norm()), side(side_)
   {
      xtensor::xVector<> axis3 = axis1 % axis2;
      axis3.norm();
      std::vector<xPointToDouble> ds;
      ds.emplace_back(xPlane(xPointPlusVector(center, axis1 * side / 2.), axis1));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis1 * (-side / 2.)), axis1 * -1.));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis2 * side / 2.), axis2));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis2 * -side / 2.), axis2 * -1.));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis3 * side / 2.), axis3));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis3 * -side / 2.), axis3 * -1.));
      xInterVector inter(ds);
      total = inter;
   }
   double operator()(const xtensor::xPoint& p) const { return total(p); }

  private:
   xCube();
   xtensor::xPoint center;
   xtensor::xVector<> axis1, axis2;
   double side;
   double radius;
   mutable xPointToDouble total;
};

// to compute the signed distance to a parallelepiped rectangle
class xParallelepipedRectangle
{
  public:
   xParallelepipedRectangle(const xtensor::xPoint& center_, xtensor::xVector<> axis1_, xtensor::xVector<> axis2_,
                            const double& side1_, const double& side2_, const double& side3_)
       : center(center_), axis1(axis1_.norm()), axis2(axis2_.norm()), side1(side1_), side2(side2_), side3(side3_)
   {
      xtensor::xVector<> axis3 = axis1 % axis2;
      axis3.norm();
      std::vector<xPointToDouble> ds;
      ds.emplace_back(xPlane(xPointPlusVector(center, axis1 * side1 / 2.), axis1));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis1 * -side1 / 2.), axis1 * -1.));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis2 * side2 / 2.), axis2));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis2 * -side2 / 2.), axis2 * -1.));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis3 * side3 / 2.), axis3));
      ds.emplace_back(xPlane(xPointPlusVector(center, axis3 * -side3 / 2.), axis3 * -1.));
      xInterVector inter(ds);
      total = inter;
   }
   double operator()(const xtensor::xPoint& p) const { return total(p); }

  private:
   xParallelepipedRectangle();
   xtensor::xPoint center;
   xtensor::xVector<> axis1, axis2;
   double side1, side2, side3;
   double radius;
   mutable xPointToDouble total;
};

// to compute the signed distance to a cone
// a cone is like a ice cream cone
// the level set in negative inside
// the axis points inside the cone
// see page 43 book 1
class xCone
{
  public:
   xCone(const xtensor::xPoint& o, xtensor::xVector<> d, const double& a) : origin(o), axis(d.norm()), angle(a)
   {
      d.norm();
      // const double PI = 4. * atan(1.);
      assert(0. < angle && angle < M_PI);
   }
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> op = xtensor::xVector<>(origin, p);
      double op_axis = op * axis;
      xtensor::xVector<> op_ortho = op - (axis * (op_axis));
      double norm_op_ortho = std::sqrt(op_ortho * op_ortho);
      //    double cone_dual = norm_op_ortho + op_axis * (1.0/std::tan(angle));
      double cone_dual = std::tan(angle) * norm_op_ortho + op_axis;
      if (cone_dual <= 0.0) return std::sqrt(op * op);
      if (norm_op_ortho >= 1.e-15) op_ortho /= norm_op_ortho;
      xtensor::xVector<> e = axis * (-std::sin(angle)) + op_ortho * std::cos(angle);
      return op * e;
   }

  private:
   xCone();
   xtensor::xPoint origin;
   xtensor::xVector<> axis;
   double angle;
};

// to compute the sign distance  to a cylinder
// positive outside
class xCylinder
{
  public:
   xCylinder(const xtensor::xPoint& p, xtensor::xVector<> d, const double& r) : point_on_axis(p), axis(d.norm()), radius(r) {}
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> op = xtensor::xVector<>(point_on_axis, p);
      double d = op * axis;
      xtensor::xVector<> ortho = op - axis * d;
      return std::sqrt(ortho * ortho) - radius;
   }

  private:
   xCylinder();
   xtensor::xPoint point_on_axis;
   xtensor::xVector<> axis;
   double radius;
};

// to compute the sign distance  to a cylinder
// with hexagonal base
// negative inside the cylinder
// positive outside
// radius is the radius of the circum cicle to the hexa
// axis is the axis of the cylinder
// second_axis indicate a direction in which the axis
// is the further apart from the cylinder
class xHexaCylinder
{
  public:
   xHexaCylinder(const xtensor::xPoint& p, xtensor::xVector<> d, xtensor::xVector<> b, const double& r)
       : point_on_axis(p), axis(d.norm()), second_axis(b), radius(r)
   {
      double proj = second_axis * axis;
      second_axis = second_axis - axis * proj;
      second_axis.norm();
      double inner_radius = std::sqrt(3.) / 2. * radius;
      double rr = inner_radius;
      const double PI = 4. * atan(1.);
      assert(0);  // 6 line below are wrong
                  // to correct see the square
      xtensor::xPoint pt1 = point_on_axis + xtensor::xPoint(rr * cos(PI / 6.), rr * sin(PI / 6.), 0.);
      xtensor::xPoint pt2 = point_on_axis + xtensor::xPoint(rr * cos(3. * PI / 6.), rr * sin(3. * PI / 6.), 0.);
      xtensor::xPoint pt3 = point_on_axis + xtensor::xPoint(rr * cos(5. * PI / 6.), rr * sin(5. * PI / 6.), 0.);
      xtensor::xPoint pt4 = point_on_axis + xtensor::xPoint(rr * cos(7. * PI / 6.), rr * sin(7. * PI / 6.), 0.);
      xtensor::xPoint pt5 = point_on_axis + xtensor::xPoint(rr * cos(9. * PI / 6.), rr * sin(9. * PI / 6.), 0.);
      xtensor::xPoint pt6 = point_on_axis + xtensor::xPoint(rr * cos(11. * PI / 6.), rr * sin(11. * PI / 6.), 0.);
      xPlane pl1(pt1, xtensor::xVector<>(point_on_axis, pt1));
      xPlane pl2(pt2, xtensor::xVector<>(point_on_axis, pt2));
      xPlane pl3(pt3, xtensor::xVector<>(point_on_axis, pt3));
      xPlane pl4(pt4, xtensor::xVector<>(point_on_axis, pt4));
      xPlane pl5(pt5, xtensor::xVector<>(point_on_axis, pt5));
      xPlane pl6(pt6, xtensor::xVector<>(point_on_axis, pt6));
      xInter m1(pl1, pl2);
      xInter m2(pl3, pl4);
      xInter m3(pl5, pl6);
      xInter m4(m1, m2);
      xInter m5(m4, m3);
      total = m5;
   }
   double operator()(const xtensor::xPoint& p) const { return total(p); }

  private:
   xHexaCylinder();
   xtensor::xPoint point_on_axis;
   xtensor::xVector<> axis;
   xtensor::xVector<> second_axis;
   double radius;
   //
   double inner_radius;
   mutable xPointToDouble total;
};

// to compute the sign distance  to a cylinder
// with square base
// negative inside the cylinder
// positive outside
// radius is the radius of the inner circle to the square
// axis is the axis of the cylinder
// second_axis indicate a direction in which the axis
// is the closer apart from the cylinder
class xSquareCylinder
{
  public:
   xSquareCylinder(const xtensor::xPoint& p, xtensor::xVector<> d, xtensor::xVector<> b, const double& r)
       : point_on_axis(p), axis(d.norm()), second_axis(b), radius(r)
   {
      double proj = second_axis * axis;
      second_axis = second_axis - axis * proj;
      second_axis.norm();
      xtensor::xVector<> third_axis = axis % second_axis;
      xtensor::xPoint ax2(second_axis(0), second_axis(1), second_axis(2));
      xtensor::xPoint ax3(third_axis(0), third_axis(1), third_axis(2));
      xtensor::xPoint pt1 = point_on_axis + ax2 * radius;
      xtensor::xPoint pt2 = point_on_axis + ax3 * radius;
      xtensor::xPoint pt3 = point_on_axis - ax2 * radius;
      xtensor::xPoint pt4 = point_on_axis - ax3 * radius;
      xPlane pl1(pt1, xtensor::xVector<>(point_on_axis, pt1));
      xPlane pl2(pt2, xtensor::xVector<>(point_on_axis, pt2));
      xPlane pl3(pt3, xtensor::xVector<>(point_on_axis, pt3));
      xPlane pl4(pt4, xtensor::xVector<>(point_on_axis, pt4));
      xInter m1(pl1, pl2);
      xInter m2(pl3, pl4);
      xInter m3(m1, m2);
      total = m3;
   }
   double operator()(const xtensor::xPoint& p) const { return total(p); }

  private:
   xSquareCylinder();
   xtensor::xPoint point_on_axis;
   xtensor::xVector<> axis;
   xtensor::xVector<> second_axis;
   double radius;
   //
   mutable xPointToDouble total;
};

class xGrownHalfLine
{
  public:
   xGrownHalfLine(const xtensor::xPoint& origin_, xtensor::xVector<> dir_, const double& growth_)
       : origin(origin_), dir(dir_.norm()), growth(growth_)
   {
   }
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> op(origin, p);
      double proj = op * dir;
      double res;
      if (proj <= 0)
         res = op.mag();
      else
      {
         xtensor::xVector<> ort = op - dir * proj;
         res = ort.mag();
      }
      return res - growth;
   }

  private:
   const xtensor::xPoint origin;
   const xtensor::xVector<> dir;
   const double growth;
};

class xGrownSegment
{
  public:
   xGrownSegment(const xtensor::xPoint& p1_, const xtensor::xPoint& p2_, const double& growth_)
       : p1(p1_), p2(p2_), dir(p1_, p2_), growth(growth_)
   {
      length = dir.normValue();
      dir.norm();
   }
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> op(p1, p);
      double proj = op * dir;
      double res;
      if (proj <= 0)
         res = op.mag();
      else if (proj >= length)
      {
         xtensor::xVector<> op2(p2, p);
         res = op2.mag();
      }
      else
      {
         xtensor::xVector<> ort = op - dir * proj;
         res = ort.mag();
      }
      return res - growth;
   }

  private:
   const xtensor::xPoint p1, p2;
   xtensor::xVector<> dir;
   double length;
   const double growth;
};

class xGrownPlane
{
  public:
   xGrownPlane(const xtensor::xPoint& p1_, const xtensor::xPoint& p2_, const xtensor::xPoint& p3_, const double& growth_)
       : p1(p1_), p2(p2_), p3(p3_), segment(p1_, p2_), dir_inf(p1_, p3_), growth(growth_)
   {
      length = segment.normValue();
      dir_inf.norm();
      segment.norm();
   }

   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> x(p1, p);
      double coord1 = std::max(0., std::min(x * segment, length));
      double coord2 = x * dir_inf;
      xtensor::xVector<> x_in_plane = coord1 * segment + coord2 * dir_inf;

      xtensor::xVector<> temp = x - x_in_plane;
      double distance_x2plane = temp.normValue();

      return (distance_x2plane - growth);
   }

  private:
   const xtensor::xPoint p1, p2, p3;
   xtensor::xVector<> segment, dir_inf;
   double length;
   const double growth;
};

class xGrownEllipse
{
  public:
   xGrownEllipse(const xtensor::xPoint& p1_, const xtensor::xPoint& p2_, const xtensor::xPoint& p3_, const double& growth_)
       : p1(p1_), va(p1_, p2_), vb(p1_, p3_), a2(va.normValue()), b2(vb.normValue()), growth(growth_)
   {
      a2 *= a2;
      b2 *= b2;
   }
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> x(p1, p);
      double coord1 = x * va;
      double coord2 = x * vb;
      double eq = coord1 * coord1 / a2 + coord2 * coord2 / b2;

      if (eq > 1.)
      {
         const double lambda = sqrt(1. / eq);
         coord1 = lambda * coord1;
         coord2 = lambda * coord2;
      }

      xtensor::xVector<> x_in_disk = coord1 * va + coord2 * vb;

      xtensor::xVector<> directioniso = x_in_disk - x;

      return directioniso.normValue() - growth;
   }

  private:
   const xtensor::xPoint p1;
   xtensor::xVector<> va, vb;
   double a2, b2;
   const double growth;
};

class xGrownEllipticalRing
{
  public:
   xGrownEllipticalRing(const xtensor::xPoint& p1_, const xtensor::xPoint& p2_, const xtensor::xPoint& p3_,
                        const double& thickness, const double& growth_)
       : p1(p1_),
         va(p1_, p2_),
         vb(p1_, p3_),
         al2(va.normValue()),
         bl2(vb.normValue()),
         au2(al2 + thickness),
         bu2(bl2 + thickness),
         growth(growth_)
   {
      al2 *= al2;
      bl2 *= bl2;
      au2 *= au2;
      bu2 *= bu2;
   }
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> x(p1, p);
      double coord1 = x * va;
      double coord2 = x * vb;
      const double c12 = coord1 * coord1;
      const double c22 = coord2 * coord2;
      double eql = c12 / al2 + c22 / bl2;
      double equ = c12 / au2 + c22 / bu2;

      if (eql < 1.)
      {
         const double lambda = sqrt(1. / eql);
         coord1 *= lambda;
         coord2 *= lambda;
      }
      else if (equ > 1.)
      {
         const double lambda = sqrt(1. / equ);
         coord1 *= lambda;
         coord2 *= lambda;
      }

      xtensor::xVector<> x_in_disk = coord1 * va + coord2 * vb;

      xtensor::xVector<> directioniso = x_in_disk - x;

      return directioniso.normValue() - growth;
   }

  private:
   const xtensor::xPoint p1;
   xtensor::xVector<> va, vb;
   double al2, bl2;
   double au2, bu2;
   const double growth;
};

class xGrownPartialSphere
{
  public:
   /// Constructor :
   //!  orig_ = origin of the global sphere
   //!  radius_ = radius of the global sphere
   //!  n= vector describing cone central axe passing by orig
   //!  angle = angle of the cone (angle/2 is the angle between n and any vector made of orig and a point on the cone)
   //!  growth_= offset to impose from the geometrical surface corresponding to the intersection of the cone and the global
   //!  sphere.
   //
   xGrownPartialSphere(const xtensor::xPoint& orig_, const xtensor::xVector<>& n, const double& radius_, const double& angle,
                       const double& growth_)
       : orig(orig_),
         vn(n),
         radius(radius_),
         growth(growth_),
         cos_gamas2(cos(angle / 2.)),
         z0(cos_gamas2 * radius),
         torus_R(fabs(radius * sin(angle / 2.))),
         radiusp(radius + growth_),
         radiusm(radius - growth_)
   {
      if (growth > torus_R)
      {
         std::cout << "growth parmeter is to large. It is " << growth << " and should be less then " << torus_R << std::endl;
         std::cout << "Or change implementation to deal with torus intersection" << std::endl;
         throw -487956;
      }

      // Normalize
      vn.norm();
   }
   double operator()(const xtensor::xPoint& p) const
   {
      xtensor::xVector<> vp(orig, p);
      xtensor::xVector<> vpn(vp);
      const double dist = vpn.normValue();
      const double cos_point = vn * vpn;
      // in cone
      if (cos_point > cos_gamas2)
      {
         if (dist > radius)
            return dist - radiusp;
         else
            return radiusm - dist;
      }
      // out of cone
      else
      {
         // Use axis symmetrical information :
         const double z = z0 - dist * cos_point;
         xtensor::xVector<> vortho(vn % vp);
         const double r = vortho.mag() - torus_R;
         return sqrt(z * z + r * r) - growth;
      }
   }

  private:
   const xtensor::xPoint orig;
   xtensor::xVector<> vn;
   const double radius, growth, cos_gamas2, z0, torus_R;
   const double radiusp, radiusm;
};

class xXYZInterface
{
  public:
   xXYZInterface(const xPointToDouble& t_) : t(t_) {}
   double operator()(const double& x, const double& y, const double& z) const { return t(xtensor::xPoint(x, y, z)); }

  private:
   xPointToDouble t;
};

class xSuperEllipse
{
  public:
   // n1 and n2 are doubles because exponents can be non-integer
   xSuperEllipse(const xtensor::xPoint& p, const double rad_ = 1., const double n1_ = 2, const double n2_ = 2)
       : center(p), rad(rad_), n1(n1_), n2(n2_)
   {
   }

   double operator()(const xtensor::xPoint& p) const
   {
      double dx = p(0) - center(0);
      double dy = p(1) - center(1);

      return pow((dx + dy) * (dx + dy) * 0.5, n1) + pow((dx - dy) * (dx - dy) * 2. / 9., n2) - rad;
   }

  private:
   xtensor::xPoint center;
   const double rad, n1, n2;
};

}  // namespace xgeom

#endif
