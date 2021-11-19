/*
     This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef _XCRACKGEOMETRY_H
#define _XCRACKGEOMETRY_H

#include <cmath>
#include <iostream>

#include "xPoint.h"
#include "xPointToDouble.h"
#include "xTensor2.h"
#include "xVector.h"

namespace xcrack
{
class xCenterCrack
{
  public:
   xCenterCrack(const xtensor::xPoint& p, xtensor::xVector<> n, double length_)
       : point_on_plane(p), normal(n.norm()), length(length_)
   {
   }
   double operator()(const xtensor::xPoint& p) const { return (abs(xtensor::xVector<>(point_on_plane, p) * normal) - length); }

  private:
   xCenterCrack();
   xtensor::xPoint point_on_plane;
   xtensor::xVector<> normal;
   double length;
};

}  // namespace xcrack

#endif
