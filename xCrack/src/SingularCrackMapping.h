/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef H_SingularCrackMapping
#define H_SingularCrackMapping

// xmapping
#include "xLagrangeMapping.h"
// xtensor
#include "xVector.h"

namespace xcrack
{
//! A Mapping for a triangular element which has one node at the crack Tipe.
/*! The reference element is a quad between with  -1. <= u,v <= 1., with edge u =-1. map to the crack tip,
 *  with a singular mapping for u = -1.
 */
class SingularCrackMapping : public xmapping::xMapping
{
  public:
   // constructor
   SingularCrackMapping(AOMD::mEntity *e_, int singular_node_);
   // operators
   void eval(double u, double v, double w, double &x, double &y, double &z) const override;
   void deval(double u, double v, double w, double &dxdu, double &dydu, double &dzdu, double &dxdv, double &dydv, double &dzdv,
              double &dxdw, double &dydw, double &dzdw) const override;
   void boundingBox(xtensor::xPoint &pMin, xtensor::xPoint &pMax) const override;
   // void normalVector(AOMD::mEntity *, double u, double v, double w, xmapping::vector3d &n) const override;
   bool inReferenceElement(double u, double v, double w) const override;
   void COG(double &u, double &v, double &w) const override;
   int order() const override;
   int geomOrder() const override;
   double a;
   double b;
   double r0;
   double v1;
   double v2;
   double v3;
   double w1;
   double w2;
   double w3;
   double aaa;
   double f_a;
   double g_a;
   double singX;
   double singY;
   double singZ;
   xtensor::xVector<> VectNorm;
   xtensor::xVector<> vect12norm;
   double theta1, theta2;

  private:
   xmapping::xConstantLagrangeMapping mapping_regular;
   int singular_node;
   void evalprime(double u, double v, double &xp, double &yp) const;
   void devalprime(double u, double v, double &dxpdu, double &dypdu, double &dxpdv, double &dypdv) const;
   void invertprime(double up, double vp, double &xp, double &yp) const;
};

}  // namespace xcrack
#endif
