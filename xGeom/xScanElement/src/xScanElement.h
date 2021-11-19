/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _SCAN_ELEMENT_H
#define _SCAN_ELEMENT_H

#ifdef DEBUG
#include <iostream>
#endif

#include <vector>

// AOMD
#include "mEntity.h"
#include "mVertex.h"

// xtensor
#include "xPoint.h"

namespace xgeom
{
/// \struct xAOMDToScanElement
/// \brief function object that return ScanElement type from
/// AOMD type.
struct xAOMDToScanElement
{
   typedef xtensor::xPoint PointArg_t;
   typedef AOMD::mEntity* ElementArg_t;

   typedef std::vector<double> PointRes_t;
   typedef std::vector<PointRes_t> ElementRes_t;

   PointRes_t convertPoint(PointArg_t point_)
   {
      PointRes_t point;
      point.push_back(point_(0));
      point.push_back(point_(1));
      point.push_back(point_(2));
      return point;
   }
   ElementRes_t convertElement(ElementArg_t element_)
   {
      ElementRes_t element;
      for (int i = 0; i < element_->size(0); ++i)
      {
         PointArg_t p = static_cast<AOMD::mVertex*>(element_->get(0, i))->point();
         element.push_back(convertPoint(p));
      }
      return element;
   }
};

/// \struct xScanIdentity
/// \brief function object that to nothing
struct xScanIdentity
{
   typedef std::vector<double> PointArg_t;
   typedef std::vector<PointArg_t> ElementArg_t;

   typedef PointArg_t PointRes_t;
   typedef ElementArg_t ElementRes_t;

   PointRes_t convertPoint(PointArg_t point_) { return point_; }
   ElementRes_t convertElement(ElementArg_t element_) { return element_; }
};

/// \struct xGreaterPoint2D
struct xGreaterPoint2D : std::binary_function<std::vector<double>, std::vector<double>, bool>
{
   result_type operator()(first_argument_type arg_1, first_argument_type arg_2) { return arg_1[1] > arg_2[1]; }
};

/// \struct xGreaterPoint3D
struct xGreaterPoint3D : std::binary_function<std::vector<double>, std::vector<double>, bool>
{
   result_type operator()(first_argument_type arg_1, first_argument_type arg_2) { return arg_1[2] > arg_2[2]; }
};

/// \struct xScanElement
/// \brief base abstract struct for all xScanElement classes
template <typename UNARYOP>
struct xScanElement
{
   typedef std::vector<int> IJK;
   typedef std::vector<double> XYZ;

   virtual ~xScanElement() = default;
   virtual void scan(typename UNARYOP::ElementArg_t, std::vector<int>&) = 0;
   virtual inline IJK getIJK(const int id) = 0;
   virtual inline XYZ getXYZ(const int id) = 0;
};

/// \struct xScanElement2D
/// \brief given a bounding box, a level and a tolerance to
/// the constructor, given an element to the scan(...)
/// function, return all id of points contained in the
/// element. level and id numbering are based on Octree
/// concepts. For 2D purpose.
template <typename UNARYOP>
class xScanElement2D : public xScanElement<UNARYOP>
{
  public:
   typedef std::vector<int> IJK;
   typedef std::vector<double> XYZ;
   typedef typename UNARYOP::PointRes_t Point_t;
   typedef typename UNARYOP::ElementRes_t Element_t;

   xScanElement2D(typename UNARYOP::PointArg_t p_min_, typename UNARYOP::PointArg_t p_max_, const int level_,
                  const double tolerance_);

   void scan(typename UNARYOP::ElementArg_t triangle_, std::vector<int>& active_id) override;
   inline IJK getIJK(const int id) override;
   inline XYZ getXYZ(const int id) override;

   inline double intersectionLineLine(const double x1, const double y1, const double x2, const double y2, const double y_line);
   void scanPoint(const double x, const int j, std::vector<int>& active_id);

  private:
   const Point_t p_min;
   const Point_t p_max;

   const int level;
   const double tolerance;
   const int nb;
   const double pas_x;
   const double pas_y;
#ifdef DEBUG
   const bool debug;
#endif
   UNARYOP entity_to_scan_entity;
};

/// \struct xScanElement3D
/// \brief given a bounding box, a level and a tolerance to
/// the constructor, given an element to the scan(...)
/// function, return all id of points contained in the
/// element. level and id numbering are based on Octree
/// concepts. For 3D purpose.
template <typename UNARYOP>
class xScanElement3D : public xScanElement<UNARYOP>
{
  public:
   typedef std::vector<int> IJK;
   typedef std::vector<double> XYZ;
   typedef typename UNARYOP::PointRes_t Point_t;
   typedef typename UNARYOP::ElementRes_t Element_t;

   xScanElement3D(typename UNARYOP::PointArg_t p_min_, typename UNARYOP::PointArg_t p_max_, const int level_,
                  const double tolerance_);

   void scan(typename UNARYOP::ElementArg_t tetraedre_, std::vector<int>& active_id) override;
   inline IJK getIJK(const int id) override;
   inline XYZ getXYZ(const int id) override;

  private:
   inline typename UNARYOP::PointRes_t intersectionLinePlane(const double x1, const double y1, const double z1, const double x2,
                                                             const double y2, const double z2, const double z_plane);
   void scanPoint(const double x, const double y, const int k, std::vector<int>& active_id);
   void scanSegment(const double x1, const double y1, const double x2, const double y2, const int k, std::vector<int>& active_id);

   const Point_t p_min;
   const Point_t p_max;

   const int level;
   const double tolerance;
   const int nb;
   const int nb2;
   const double pas_x;
   const double pas_y;
   const double pas_z;
#ifdef DEBUG
   const bool debug;
#endif
   xScanElement2D<xScanIdentity> scan2D;

   UNARYOP entity_to_scan_entity;
};
}  // namespace xgeom

#include "xScanElement_imp.h"

#endif
