/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _DISTANCE_NEAREST_POINT_INTERNAL_H
#define _DISTANCE_NEAREST_POINT_INTERNAL_H

#include <iostream>
#include <limits>
#include <vector>

#include "mpi.h"

#ifdef HAVE_CGAL
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#endif

#ifdef HAVE_ANN
#include <ANN/ANN.h>
#endif

#include "workInProgress.h"
#include "xDistanceNearestPoint_traits.h"

namespace xgeom
{
// Internal xDistanceNearestPoint class, NOT TO BE USED.
template <typename ITER, typename GENERATOR, typename TRAIT>
class xDistanceNearestPoint_internal
{
   // empty class //
};

#ifdef HAVE_CGAL
// Partial specialisation on xCGALTrait
template <typename ITER, typename GENERATOR>
class xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>
{
  public:
   typedef CGAL::AABB_traits<typename GENERATOR::Kernel_t, typename GENERATOR::Primitive_t> CGALTraits_t;
   typedef CGAL::AABB_tree<CGALTraits_t> CGALTree_t;

   xDistanceNearestPoint_internal(ITER begin_, ITER end_, const GENERATOR& generator_, MPI_Comm world_);
   ~xDistanceNearestPoint_internal();

   // Compute from point
   double distance(const typename GENERATOR::PointArg_t& point_);
   void nearestPoint(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_);
   void nearestPointDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_);
   void nearestEntityDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_);
   void nearestVerticesInsideRadiusDistances(const typename GENERATOR::PointArg_t& point_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_)
   {
      std::cerr << "This function is not coded in CGAL context. ANN is well adapted. Use a generator working with ANN."
                << std::endl;
      throw;
      return;
   }

   // Compute from vertex
   double distance(typename GENERATOR::VertexArg_t vertex_);
   void nearestPoint(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_);
   void nearestPointDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_);
   void nearestEntityDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_);
   void nearestVerticesInsideRadiusDistances(typename GENERATOR::VertexArg_t vertex_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_)
   {
      std::cerr << "This function is not coded in CGAL context. ANN is well adapted. Use a generator working with ANN."
                << std::endl;
      throw;
      return;
   }

  private:
   MPI_Comm world;
   CGALTree_t tree;
   std::vector<typename GENERATOR::EntityArg_t> entities_vect;
   typename GENERATOR::Container_t cgal_entities;
   const GENERATOR& generator;
};
#endif

#ifdef HAVE_ANN
// Partial specialisation on xANNTrait
template <typename ITER, typename GENERATOR>
class xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>
{
  public:
   xDistanceNearestPoint_internal(ITER begin_, ITER end_, const GENERATOR& generator_, MPI_Comm world);
   ~xDistanceNearestPoint_internal();

   // Compute from point
   double distance(const typename GENERATOR::PointArg_t& point_);
   void nearestPoint(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_);
   void nearestPointDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_);
   void nearestEntityDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_);
   void nearestVerticesInsideRadiusDistances(const typename GENERATOR::PointArg_t& point_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_);

   // Compute from vertex
   double distance(typename GENERATOR::VertexArg_t vertex_);
   void nearestPoint(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_);
   void nearestPointDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_);
   void nearestEntityDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_);
   void nearestVerticesInsideRadiusDistances(typename GENERATOR::VertexArg_t vertex_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_);

  private:
   void nearestVertexDistance_internal(typename GENERATOR::PointRes_t ann_point_,
                                       typename GENERATOR::VertexArg_t& nearest_vertex_, double& distance_);
   void nearestVerticesInsideRadiusDistances_internal(
       typename GENERATOR::PointRes_t ann_point_, double radius_factor_,
       std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_);

   ANNkd_tree* tree;
   ANNpointArray data_pts;
   std::vector<typename GENERATOR::EntityRes_t> ann_vertices;
   const GENERATOR& generator;
};
#endif

// Partial specialisation on xBruteForceTrait
template <typename ITER, typename GENERATOR>
class xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>
{
  public:
   xDistanceNearestPoint_internal(ITER begin_, ITER end_, const GENERATOR& generator_, MPI_Comm world);
   ~xDistanceNearestPoint_internal();

   // Compute from point
   double distance(const typename GENERATOR::PointArg_t& point_);
   void nearestPoint(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_);
   void nearestPointDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_);
   void nearestEntityDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_)
   {
      std::cerr << "This function is not coded in BruteForce context. ANN is well adapted. Use a generator working with ANN."
                << std::endl;
      throw;
      return;
   }
   void nearestVerticesInsideRadiusDistances(const typename GENERATOR::PointArg_t& point_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_)
   {
      std::cerr << "This function is not coded in BruteForce context. ANN is well adapted. Use a generator working with ANN."
                << std::endl;
      throw;
      return;
   }

   // Compute from vertex
   double distance(typename GENERATOR::VertexArg_t vertex_);
   void nearestPoint(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_);
   void nearestPointDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_);
   void nearestEntityDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_)
   {
      std::cerr << "This function is not coded in BruteForce context. ANN is well adapted. Use a generator working with ANN."
                << std::endl;
      throw;
      return;
   }
   void nearestVerticesInsideRadiusDistances(typename GENERATOR::VertexArg_t vertex_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_)
   {
      std::cerr << "This function is not coded in BruteForce context. ANN is well adapted. Use a generator working with ANN."
                << std::endl;
      throw;
      return;
   }

  private:
   void nearestPointDistance_internal(typename GENERATOR::PointRes_t P, typename GENERATOR::PointRes_t Cp, double& distance_);

   std::vector<typename GENERATOR::EntityRes_t> brute_force_edges;
   const GENERATOR& generator;
};
}  // namespace xgeom

#include "xDistanceNearestPoint_internal_imp.h"

#endif
