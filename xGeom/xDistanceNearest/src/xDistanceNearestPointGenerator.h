/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _DISTANCE_NEAREST_POINT_GENERATOR_H
#define _DISTANCE_NEAREST_POINT_GENERATOR_H
#ifndef HAVE_AOMD
#define HAVE_AOMD
#endif

#include <vector>

#ifdef HAVE_CGAL
// CGAL
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Simple_cartesian.h>
// interface
#include "xAttachedDataManagerAOMD.h"
#endif

#ifdef HAVE_ANN
#include <ANN/ANN.h>
#endif

#ifdef HAVE_AOMD
// trellis
#include "mEntity.h"
#include "mVertex.h"
// xinterface
#include "xAOMDEntityUtil.h"
// xmapping
#include "xMappingBuilderHolder.h"
#endif

// geom
#include "xDistanceNearestPoint_traits.h"

/// \file xDistanceNearestPointGenerator.h
/// \brief To implement a new GENERATOR, typedef to respect are:
/// Trait_t, PointArg_t, VertexArg_t, EntityArg_t, PointRes_t and
/// EntityRes_t. Up-to-now, GENERATOR are classes that derives from
/// nothing. Nevertheless, users must respect prototypes of GENERATOR
/// member functions. They are DIFFERENT for every Trait.

#ifdef HAVE_AOMD
namespace xfem
{
class xLevelSet;
}
namespace xcut
{
class xRefCutToIsoZeroVector;
}
#endif

//#define HAVE_CGAL

namespace xgeom
{
#ifdef HAVE_AOMD
/// this template function use specific generator to generate a geometric entity if a cutter methode of the generator
//! respond some sort of true for entity e
//! It also tag node according to return code of the tagger
//! it aims to be general for generator based on cutter
//! GEN type must have the following methodes :
//!        int cutter( std::vector < double >& vals);
//!        void modifyCoord(xmapping::xMapping *mapping);
//!        void modifyCoord2(xmapping::xMapping *mapping);
//!        void addList(Container_t& entity_res_list);
//!        void addList2(Container_t& entity_res_list);
//!        bool tagTouchingIso(mEntity*e, int id);
//! GEN type must have the following object :
//!         ls (level set)
//!         tagger (function object of type void operator () (AOMD::mEntity* e) const)
//
template <typename GEN>
void genAndTagFromCut(GEN& generator, AOMD::mEntity* e, typename GEN::Container_t& container);
#ifdef HAVE_CGAL
// struct xAOMDVertexToCGALPoint
// {
//    typedef xCGALTrait Trait_t;
// 
//    typedef xtensor::xPoint PointArg_t;
//    typedef AOMD::mEntity* VertexArg_t;
//    typedef AOMD::mEntity* EntityArg_t;
// 
//    typedef CGAL::Simple_cartesian<double> Kernel_t;
//    typedef Kernel_t::Point_3 PointRes_t;
//    typedef Kernel_t::Point_3 EntityRes_t;
// 
//    typedef std::vector<EntityRes_t> Container_t;
//    typedef Container_t::iterator Iterator_t;
//    typedef CGAL::AABB_point_primitive<Kernel_t, Iterator_t> Primitive_t;
// 
//    void generateEntity(EntityArg_t entity_, Container_t& entity_res_list_) const
//    {
//       PointArg_t pt1 = static_cast<AOMD::mVertex*>(entity_)->point();
//       entity_res_list_.emplace_back(pt1(0), pt1(1), pt1(2));
//       return;
//    }
//    PointRes_t convertPoint(const PointArg_t& point_) const { return PointRes_t(point_(0), point_(1), point_(2)); }
//    PointRes_t convertVertex(VertexArg_t vertex_) const
//    {
//       PointArg_t point = static_cast<AOMD::mVertex*>(vertex_)->point();
//       return PointRes_t(point(0), point(1), point(2));
//    }
//    PointArg_t reverseConvertPoint(const PointRes_t& point_) const { return PointArg_t(point_[0], point_[1], point_[2]); }
// };

/// \struct xAOMDEdgeToCGALEdge
/// \brief A simple generator that convert AOMD mesh entity
/// types (xtensor::xPoint, mVertex, mEdge) to CGAL types.  It provides
/// as well xCGALTrait.
struct xAOMDEdgeToCGALEdge
{
   typedef xCGALTrait Trait_t;

   typedef xtensor::xPoint PointArg_t;
   typedef AOMD::mEntity* VertexArg_t;
   typedef AOMD::mEntity* EntityArg_t;

   typedef CGAL::Simple_cartesian<double> Kernel_t;
   typedef Kernel_t::Point_3 PointRes_t;
   typedef Kernel_t::Segment_3 EntityRes_t;

   typedef std::vector<EntityRes_t> Container_t;
   typedef Container_t::iterator Iterator_t;
   typedef CGAL::AABB_segment_primitive<Kernel_t, Iterator_t> Primitive_t;

   void generateEntity(EntityArg_t entity_, Container_t& entity_res_list_) const
   {
      PointArg_t pt1 = static_cast<AOMD::mVertex*>(entity_->get(0, 0))->point();
      PointArg_t pt2 = static_cast<AOMD::mVertex*>(entity_->get(0, 1))->point();
      PointRes_t ptc1(pt1(0), pt1(1), pt1(2));
      PointRes_t ptc2(pt2(0), pt2(1), pt2(2));
      entity_res_list_.emplace_back(ptc1, ptc2);
      return;
   }
   PointRes_t convertPoint(const PointArg_t& point_) const { return PointRes_t(point_(0), point_(1), point_(2)); }
   PointRes_t convertVertex(VertexArg_t vertex_) const
   {
      PointArg_t point = static_cast<AOMD::mVertex*>(vertex_)->point();
      return PointRes_t(point(0), point(1), point(2));
   }
   PointArg_t reverseConvertPoint(const PointRes_t& point_) const { return PointArg_t(point_[0], point_[1], point_[2]); }
};

/// \struct xAOMDTriToCGALTri
/// \brief A simple generator that convert AOMD mesh entity
/// types (xtensor::xPoint, mVertex, mFace) to CGAL types.  It provides
/// as well xCGALTrait.
struct xAOMDTriToCGALTri
{
   typedef xCGALTrait Trait_t;

   typedef xtensor::xPoint PointArg_t;
   typedef AOMD::mEntity* VertexArg_t;
   typedef AOMD::mEntity* EntityArg_t;

   typedef CGAL::Simple_cartesian<double> Kernel_t;
   typedef Kernel_t::Point_3 PointRes_t;
   typedef Kernel_t::Triangle_3 EntityRes_t;

   typedef std::vector<EntityRes_t> Container_t;
   typedef Container_t::iterator Iterator_t;
   typedef CGAL::AABB_triangle_primitive<Kernel_t, Iterator_t> Primitive_t;

   void generateEntity(EntityArg_t entity_, Container_t& entity_res_list_) const
   {
      PointArg_t pt1 = static_cast<AOMD::mVertex*>(entity_->get(0, 0))->point();
      PointArg_t pt2 = static_cast<AOMD::mVertex*>(entity_->get(0, 1))->point();
      PointArg_t pt3 = static_cast<AOMD::mVertex*>(entity_->get(0, 2))->point();
      PointRes_t ptc1(pt1(0), pt1(1), pt1(2));
      PointRes_t ptc2(pt2(0), pt2(1), pt2(2));
      PointRes_t ptc3(pt3(0), pt3(1), pt3(2));
      entity_res_list_.emplace_back(ptc1, ptc2, ptc3);
      return;
   }
   PointRes_t convertPoint(const PointArg_t& point_) const { return PointRes_t(point_(0), point_(1), point_(2)); }
   PointRes_t convertVertex(VertexArg_t vertex_) const
   {
      PointArg_t point = static_cast<AOMD::mVertex*>(vertex_)->point();
      return PointRes_t(point(0), point(1), point(2));
   }
   PointArg_t reverseConvertPoint(const PointRes_t& point_) const { return PointArg_t(point_[0], point_[1], point_[2]); }
};

class xLevelSetToCGALAndTag
{
  protected:
   xLevelSetToCGALAndTag(xfem::xLevelSet& ls_, const xinterface::aomd::xAttachedDataManagerAOMD<char>& tag_nodes_, bool fit,
                         double fittol);
   class tagNode
   {
     public:
      tagNode(const xinterface::aomd::xAttachedDataManagerAOMD<char>& tag_nodes_);
      void operator()(AOMD::mEntity* e) const;

     private:
      const xinterface::aomd::xAttachedDataManagerAOMD<char>& tag_nodes;
   };
   const xcut::xRefCutToIsoZeroVector& refcut;
   xfem::xLevelSet& ls;
   const tagNode tagger;
};
/// \class xLevelSetTriToCGALEdgeAndTag
/// \brief A generator that cut AOMD triangle face mesh entity,
//! using Level set values at it's node, to generate iso zero
//! CGAL edge
//! It provides as well xCGALTrait.
//! It tag node defining iso zero as well.
class xLevelSetTriToCGALEdgeAndTag : private xLevelSetToCGALAndTag
{
  public:
   typedef xCGALTrait Trait_t;

   typedef xtensor::xPoint PointArg_t;
   typedef AOMD::mEntity* VertexArg_t;
   typedef AOMD::mEntity* EntityArg_t;

   typedef CGAL::Simple_cartesian<double> Kernel_t;
   typedef Kernel_t::Point_3 PointRes_t;
   typedef Kernel_t::Segment_3 EntityRes_t;

   typedef std::vector<EntityRes_t> Container_t;
   typedef Container_t::iterator Iterator_t;
   typedef CGAL::AABB_segment_primitive<Kernel_t, Iterator_t> Primitive_t;

   xLevelSetTriToCGALEdgeAndTag(xfem::xLevelSet& ls_, const xinterface::aomd::xAttachedDataManagerAOMD<char>& tag_nodes_,
                                bool fit = true, double fittol = 1.e-2);

   // mandatory for xDistanceNearestPoint class
   void generateEntity(EntityArg_t entity_, Container_t& entity_res_list_) const;

   PointRes_t convertPoint(const PointArg_t& point_) const { return PointRes_t(point_(0), point_(1), point_(2)); }
   PointRes_t convertVertex(VertexArg_t vertex_) const
   {
      PointArg_t point = static_cast<AOMD::mVertex*>(vertex_)->point();
      return PointRes_t(point(0), point(1), point(2));
   }
   PointArg_t reverseConvertPoint(const PointRes_t& point_) const { return PointArg_t(point_[0], point_[1], point_[2]); }

  private:
   std::vector<double> cut_points;
   template <typename GEN>
   friend void genAndTagFromCut(GEN& generator, AOMD::mEntity* e, typename GEN::Container_t& container);

   // mandatory for genAndTagFromCut function
   int cutter(std::vector<double>& vals);
   void modifyCoord(xmapping::xMapping* mapping);
   void modifyCoord2(xmapping::xMapping* mapping);
   void addList(Container_t& entity_res_list);
   void addList2(Container_t& entity_res_list);
   bool tagTouchingIso(AOMD::mEntity* e, int id);
};
/// \class xLevelSetTetToCGALTriAndTag
/// \brief A generator that cut AOMD Tetraedron  mesh entity,
//! using Level set values at it's node, to generate iso zero
//! CGAL Triangle
//! It provides as well xCGALTrait.
//! It tag node defining iso zero as well.
class xLevelSetTetToCGALTriAndTag : private xLevelSetToCGALAndTag
{
  public:
   typedef xCGALTrait Trait_t;

   typedef xtensor::xPoint PointArg_t;
   typedef AOMD::mEntity* VertexArg_t;
   typedef AOMD::mEntity* EntityArg_t;

   typedef CGAL::Simple_cartesian<double> Kernel_t;
   typedef Kernel_t::Point_3 PointRes_t;
   typedef Kernel_t::Triangle_3 EntityRes_t;

   typedef std::vector<EntityRes_t> Container_t;
   typedef Container_t::iterator Iterator_t;
   typedef CGAL::AABB_triangle_primitive<Kernel_t, Iterator_t> Primitive_t;

   xLevelSetTetToCGALTriAndTag(xfem::xLevelSet& ls_, const xinterface::aomd::xAttachedDataManagerAOMD<char>& tag_nodes_,
                               bool fit = true, double fittol = 1.e-2);

   void generateEntity(EntityArg_t entity_, Container_t& entity_res_list_) const;

   PointRes_t convertPoint(const PointArg_t& point_) const { return PointRes_t(point_(0), point_(1), point_(2)); }
   PointRes_t convertVertex(VertexArg_t vertex_) const
   {
      PointArg_t point = static_cast<AOMD::mVertex*>(vertex_)->point();
      return PointRes_t(point(0), point(1), point(2));
   }
   PointArg_t reverseConvertPoint(const PointRes_t& point_) const { return PointArg_t(point_[0], point_[1], point_[2]); }

  private:
   std::vector<double> cut_points;
   template <typename GEN>
   friend void genAndTagFromCut(GEN& generator, AOMD::mEntity* e, typename GEN::Container_t& container);

   // mandatory for genAndTagFromCut function
   int cutter(std::vector<double>& vals);
   void modifyCoord(xmapping::xMapping* mapping);
   void modifyCoord2(xmapping::xMapping* mapping);
   void addList(Container_t& entity_res_list);
   void addList2(Container_t& entity_res_list);
   bool tagTouchingIso(AOMD::mEntity* e, int id);
};

#endif

#ifdef HAVE_ANN
/// \struct xAOMDToANN
/// \brief A simple generator that convert AOMD mesh entity
/// types (xtensor::xPoint, mVertex) to ANN types.  It provides
/// as well xANNTrait.
struct xAOMDToANN
{
   typedef xANNTrait Trait_t;

   typedef xtensor::xPoint PointArg_t;
   typedef AOMD::mVertex* VertexArg_t;
   typedef AOMD::mEntity* EntityArg_t;

   typedef ANNpoint PointRes_t;
   typedef AOMD::mVertex* EntityRes_t;

   void generateEntity(EntityArg_t entity_, std::vector<EntityRes_t>& entity_res_vector_) const
   {
      entity_res_vector_.push_back(static_cast<AOMD::mVertex*>(entity_));
      return;
   }
   PointRes_t convertPoint(const PointArg_t& point_) const
   {
      PointRes_t point = annAllocPt(3);
      point[0] = point_(0);
      point[1] = point_(1);
      point[2] = point_(2);
      return point;
   }
   PointRes_t convertVertex(VertexArg_t vertex_) const
   {
      PointArg_t point_ = static_cast<AOMD::mVertex*>(vertex_)->point();
      return convertPoint(point_);
   }
   PointArg_t reverseConvertVertex(VertexArg_t vertex_) const
   {
      PointArg_t point_ = static_cast<AOMD::mVertex*>(vertex_)->point();
      return point_;
   }
};
#endif

/// \struct xAOMDToBruteForce
/// \brief A simple generator that convert AOMD mesh entity
/// types (xtensor::xPoint, mVertex) to "BruteForce" types.  It
/// provides as well xBruteForceTrait. For 2D purpose only
struct xAOMDToBruteForce
{
   typedef xBruteForceTrait Trait_t;

   typedef xtensor::xPoint PointArg_t;
   typedef AOMD::mEntity* VertexArg_t;
   typedef AOMD::mEntity* EntityArg_t;

   typedef double* PointRes_t;
   typedef double EntityRes_t;

   void generateEntity(EntityArg_t entity_, std::vector<EntityRes_t>& entity_res_vector_) const
   {
      PointArg_t pt1 = static_cast<AOMD::mVertex*>(entity_->get(0, 0))->point();
      PointArg_t pt2 = static_cast<AOMD::mVertex*>(entity_->get(0, 1))->point();
      entity_res_vector_.push_back(pt1(0));
      entity_res_vector_.push_back(pt1(1));
      entity_res_vector_.push_back(pt2(0));
      entity_res_vector_.push_back(pt2(1));
      return;
   }
   void convertPoint(const PointArg_t& point_, PointRes_t point_res_) const
   {
      point_res_[0] = point_(0);
      point_res_[1] = point_(1);
      return;
   }
   void convertVertex(VertexArg_t vertex_, PointRes_t point_res_) const
   {
      PointArg_t point = static_cast<AOMD::mVertex*>(vertex_)->point();
      convertPoint(point, point_res_);
      return;
   }
   void reverseConvertPoint(PointRes_t point_, PointArg_t& point_res_) const
   {
      point_res_(0) = point_[0];
      point_res_(1) = point_[1];
      return;
   }
};
#endif
}  // namespace xgeom

#include "xDistanceNearestPointGenerator_imp.h"

#endif
