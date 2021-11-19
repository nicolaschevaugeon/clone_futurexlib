/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef XLAGRANGEMAPPING_
#define XLAGRANGEMAPPING_

// mapping
#include "xMapping.h"

namespace xmapping
{
class xLagrangeMapping : public xMapping
{
  public:
   xLagrangeMapping(xReferenceElementType type, const std::vector<xtensor::xPoint> &input_knots);
   ~xLagrangeMapping() override = default;
   inline std::vector<xtensor::xPoint> *getKnots() { return &knots; }
   // void normalVector(AOMD::mEntity *border, double u, double v, double w, vector3d &n) const override;
   void eval(double u, double v, double w, double &x, double &y, double &z) const override;
   void deval(double u, double v, double w, double &dxdu, double &dydu, double &dzdu, double &dxdv, double &dydv, double &dzdv,
              double &dxdw, double &dydw, double &dzdw) const override;
   void boundingBox(xtensor::xPoint &min, xtensor::xPoint &max) const override;
   xgeom::xBoundingBox boundingBox() const override;
   int geomOrder() const override;
   int order() const override;

  protected:
   std::vector<double> GeomShapeFunction(double u, double v, double w) const;
   std::vector<vector3d> GradGeomShapeFunction(double u, double v, double w) const;
   std::vector<xtensor::xPoint> knots;
};

class xConstantLagrangeMapping : public xLagrangeMapping
{
  public:
   ~xConstantLagrangeMapping() override = default;
   xConstantLagrangeMapping(xReferenceElementType type, const std::vector<xtensor::xPoint> &input_knots);
   double jacInverse(double x, double y, double z, tensor &) const override;
   double detJac(double u, double v, double w) const override;

  private:
   double det;
   tensor jac;
};

class xCylindricalCoordinatesLagrangeMapping : public xLagrangeMapping
{
  public:
   ~xCylindricalCoordinatesLagrangeMapping() override = default;
   xCylindricalCoordinatesLagrangeMapping(xReferenceElementType type, const std::vector<xtensor::xPoint> &input_knots);
   double jacInverse(double x, double y, double z, tensor &) const override;
   double detJac(double u, double v, double w) const override;
   int order() const override;
};

class xRegularCubeLagrangeMapping : public xLagrangeMapping
{
  public:
   ~xRegularCubeLagrangeMapping() override = default;
   xRegularCubeLagrangeMapping(xReferenceElementType type, const std::vector<xtensor::xPoint> &input_knots);
   // void normalVector(AOMD::mEntity *border, double u, double v, double w, vector3d &n) const override;
   void eval(double u, double v, double w, double &x, double &y, double &z) const override;
   void deval(double u, double v, double w, double &dxdu, double &dydu, double &dzdu, double &dxdv, double &dydv, double &dzdv,
              double &dxdw, double &dydw, double &dzdw) const override;
   bool invert(double x, double y, double z, double &u, double &v, double &w) const override;
   double jacInverse(double x, double y, double z, tensor &) const override;
   double detJac(double u, double v, double w) const override;
   // virtual double PushBack (double u, double v, double w, int vsize, std::vector< vector3d > &vec) const override;
  private:
   double x0, y0, z0;
   double dx, dy, dz;
};
}  // end of namespace xmapping

#endif
