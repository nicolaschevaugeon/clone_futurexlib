/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#include <cassert>
#include <fstream>
#include <iostream>

// aomd
#include "AOMD_Internals.h"
#include "mEntity.h"
#include "mPoint.h"
#include "mVertex.h"
// xtensor
#include "xVector.h"
// xmapping
#include "xMapping.h"
// xquadradure
#include "xReferenceElement.h"
// xcrack
#include "SingularCrackMapping.h"

using std::cout;
using std::endl;

#include "xAOMDEntityUtil.h"

namespace xcrack
{
SingularCrackMapping::SingularCrackMapping(AOMD::mEntity *e_, int singular_node_)
    : xmapping::xMapping(xmapping::xReferenceElementType::QUAD),
      mapping_regular(xinterface::aomd::getRefElementType(*e_), xinterface::aomd::getPoints(*e_)),
      singular_node(singular_node_)
{
   const bool debug = false;
   if (debug) cout << "entering SingularCrackMapping.cc" << endl;
   if (mapping_regular.getType() != xmapping::xReferenceElementType::TRI)
   {
      std::cout << "SingularCrack only defined for TRI " << __LINE__ << __FILE__ << std::endl;
      throw;
   }
   // type of the reference element
   // entityType = QUAD;

   AOMD::mVertex *Ver[3];

   for (int iii = 0; iii < 3; iii++)
   {
      Ver[iii] = static_cast<AOMD::mVertex *>(e_->get(0, (singular_node_ + iii) % 3));
   }
   singular_node = 0;
   singX = Ver[0]->point()(0);
   singY = Ver[0]->point()(1);
   singZ = Ver[0]->point()(2);

   xtensor::xVector<> vect01(Ver[0]->point(), Ver[1]->point());
   xtensor::xVector<> vect02(Ver[0]->point(), Ver[2]->point());
   xtensor::xVector<> vect12(Ver[1]->point(), Ver[2]->point());
   xtensor::xVector<> orientation_vector = vect01 % vect02;
   orientation_vector.norm();

   double t1 = 0.;
   double t2 = 0.;

   double norm1 = vect01.mag();
   double norm2 = vect02.mag();
   xtensor::xVector<> vect12norm = vect12;
   vect12norm.norm();

   xtensor::xVector<> VectNorm = vect12 % orientation_vector;
   VectNorm.norm();

   xtensor::xVector<> Vtemp;

   Vtemp = (VectNorm % vect01);
   double sin1 = Vtemp * orientation_vector;
   sin1 = sin1 / norm1;
   theta1 = asin(sin1);
   Vtemp = (VectNorm % vect02);
   double sin2 = Vtemp * orientation_vector;
   sin2 = sin2 / norm2;
   theta2 = asin(sin2);

   r0 = VectNorm * vect01;
   if (sqrt(r0 * r0) > 1.e-13)
   {
      t1 = 0.5 * log((1 + sin1) / (1 - sin1));
      t2 = 0.5 * log((1 + sin2) / (1 - sin2));
      a = (t1 + t2) / 2.;
      b = (t2 - t1) / 2.;
      v1 = VectNorm(0);
      v2 = VectNorm(1);
      v3 = VectNorm(2);
      w1 = vect12norm(0);
      w2 = vect12norm(1);
      w3 = vect12norm(2);
   }
}

bool SingularCrackMapping::inReferenceElement(double u, double v, double w) const
{
   const double eps = 1.e-6;
   if (u < (-1. - eps) || u > (1 + eps) || v < (-1.0 - eps) || v > (1.0 + eps)) return false;
   return true;
}
void SingularCrackMapping::COG(double &u, double &v, double &w) const
{
   u = v = w = 0.0;
   return;
}

void SingularCrackMapping::eval(double u, double v, double w, double &x, double &y, double &z) const
{
   if (sqrt(r0 * r0) > 1.e-13)
   {
      double xl = r0 * pow((1. + u) / 2., 2);
      double yl = r0 * pow((1. + u) / 2., 2) * sinh(a + b * v);
      x = singX + xl * v1 + yl * w1;
      y = singY + xl * v2 + yl * w2;
      z = singZ + xl * v3 + yl * w3;
   }
}

void SingularCrackMapping::deval(double u, double v, double w, double &dxdu, double &dydu, double &dzdu, double &dxdv,
                                 double &dydv, double &dzdv, double &dxdw, double &dydw, double &dzdw) const
{
   if (sqrt(r0 * r0) > 1.e-13)
   {
      double dxldu = r0 * (2.) / 2. * pow((1. + u) / 2., 1.);
      double dyldu = r0 * (2.) / 2. * pow((1. + u) / 2., 1.) * sinh(a + b * v);
      double dxldv = 0.;
      double dyldv = r0 * pow((1. + u) / 2., 2) * (b * cosh(a + b * v));

      dxdu = dxldu * v1 + dyldu * w1;
      dydu = dxldu * v2 + dyldu * w2;
      dzdu = dxldu * v3 + dyldu * w3;

      dxdv = dxldv * v1 + dyldv * w1;
      dydv = dxldv * v2 + dyldv * w2;
      dzdv = dxldv * v3 + dyldv * w3;

      dxdw = dydw = dzdw = 0.;
   }
}

void SingularCrackMapping::boundingBox(xtensor::xPoint &pMin, xtensor::xPoint &pMax) const
{
   mapping_regular.boundingBox(pMin, pMax);
}

/*
void SingularCrackMapping::normalVector(AOMD::mEntity *border, double u, double v, double w, Trellis_Util::mVector &n) const
{
   mapping_regular.normalVector(border, u, v, w, n);
}*/

int SingularCrackMapping::order() const { return 1; }
int SingularCrackMapping::geomOrder() const { return 1; }

}  // namespace xcrack
