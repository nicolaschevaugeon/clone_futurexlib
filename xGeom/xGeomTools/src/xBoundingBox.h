#ifndef XBOUNDINGBOX_H
#define XBOUNDINGBOX_H

#include "xPoint.h"
#include "xVector.h"

namespace xgeom
{
/// xBounding box define a cartesian bounding box by its 2 corners : min and max
class xBoundingBox
{
  public:
   using d = std::numeric_limits<double>;
   xtensor::xPoint min{d::max(), d::max(), d::max()};
   xtensor::xPoint max{d::lowest(), d::lowest(), d::lowest()};
   xBoundingBox() = default;
   xBoundingBox(const xtensor::xPoint &_min, const xtensor::xPoint &_max)
       : min{_min}, max{_max} {}  /// return true if p is inside the boundingbox
   bool contains(const xtensor::xPoint &p) const
   {
      for (int i = 0; i < 3; ++i)
         if (min[i] > p[i] || max[i] < p[i]) return false;
      return true;
   }
   /// return true if the boundingbox intersect the othe bounding box, false other wise.
   bool cover(const xBoundingBox &otherbb) const
   {
      for (int i = 0; i < 3; i++)
         if (min[i] > otherbb.max[i] || max[i] < otherbb.min[i]) return false;
      return true;
   }
   /// enlarge the bounding  box so that it include point p if needed.
   xBoundingBox &inclose(const xtensor::xPoint &p)
   {
      for (size_t i = 0; i < 3; ++i)
      {
         if (min(i) > p(i)) min(i) = p(i);
         if (max(i) < p(i)) max(i) = p(i);
      }
      return *this;
   }
   /// return the center of the bounding box
   xtensor::xPoint center() const { return (min + max) * 0.5; }
   /// return a vector from min to max
   xtensor::xVector<double> diag() const { return xtensor::xVector<double>{min, max}; }
   /// Comparaison between bounding box. return true if both bounding boxes have equal corners.
   friend bool operator==(const xBoundingBox &l, const xBoundingBox &r) { return (l.min == r.min && l.max == r.max); }
};
}  // namespace xgeom

#endif  // XBOUNDINGBOX_H
