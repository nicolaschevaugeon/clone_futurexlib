/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _DISTANCE_NEAREST_POINT_H
#define _DISTANCE_NEAREST_POINT_H
#ifndef HAVE_AOMD
#define HAVE_AOMD
#endif

#include "xDistanceNearestPointGenerator.h"
// <- this comment is necessery to force xDistanceNearestPointGenerator.h
// to stay first even if the file is clang-formated .. The file above must always
// come first !  indeed it is adding some overload to CGAL AABBPrimitive that need
// to be present before instanciation of some template.
#include "xDistanceNearestPoint_internal.h"

namespace xgeom
{
/// \class xDistanceNearestPoint
/// \brief This class provides a distance AND/OR nearest
/// point(nearest vertex in certain cases) computation from
/// various implementations (internal as well as external; for
/// example: ANN and CGAL).  The class is template on ITER and
/// GENERATOR.  GENERATOR classes are here to provide type
/// conversions from mesh entity types (for example AOMD), to
/// the desired types (for example ANN or CGAL). One can
/// imagine as well a more sophisticated generator (for
/// example in a CGAL context, from one AOMD tretrahedron cut
/// by the iso-zero of a level-set the GENERATOR could
/// generate several CGAL triangles corresponding to the
/// boundary mesh, so ITER would be iterator on AOMD
/// tetrahedron but distance AND/OR nearest point computation
/// would be done on CGAL boundary triangles). So the
/// GENERATOR CAN generate several entities from one.
template <typename ITER, typename GENERATOR>
class xDistanceNearestPoint
{
  public:
   xDistanceNearestPoint(ITER begin_, ITER end_, const GENERATOR& generator_ = GENERATOR(), MPI_Comm world = MPI_COMM_WORLD)
       : internal(begin_, end_, generator_, world)
   {
   }
   ~xDistanceNearestPoint() = default;

   /// Compute from point
   double distance(const typename GENERATOR::PointArg_t& point_) { return internal.distance(point_); }
   void nearestPoint(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_)
   {
      internal.nearestPoint(point_, nearest_point_);
      return;
   }
   void nearestPointDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_)
   {
      internal.nearestPointDistance(point_, nearest_point_, distance_);
      return;
   }
   void nearestEntityDistance(const typename GENERATOR::PointArg_t& point_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_)
   {
      internal.nearestEntityDistance(point_, nearest_entity_, nearest_point_, distance_);
      return;
   }
   void nearestVerticesInsideRadiusDistances(const typename GENERATOR::PointArg_t& point_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_)
   {
      internal.nearestVerticesInsideRadiusDistances(point_, radius_factor_, nearest_vertices_distances_);
      return;
   }
   /// Compute from vertex
   double distance(typename GENERATOR::VertexArg_t vertex_) { return internal.distance(vertex_); }
   void nearestPoint(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_)
   {
      internal.nearestPoint(vertex_, nearest_point_);
      return;
   }
   void nearestPointDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t& nearest_point_,
                             double& distance_)
   {
      internal.nearestPointDistance(vertex_, nearest_point_, distance_);
      return;
   }
   void nearestEntityDistance(typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::EntityArg_t& nearest_entity_,
                              typename GENERATOR::PointArg_t& nearest_point_, double& distance_)
   {
      internal.nearestEntityDistance(vertex_, nearest_entity_, nearest_point_, distance_);
      return;
   }
   void nearestVerticesInsideRadiusDistances(typename GENERATOR::VertexArg_t vertex_, double radius_factor_,
                                             std::map<typename GENERATOR::VertexArg_t, double>& nearest_vertices_distances_)
   {
      internal.nearestVerticesInsideRadiusDistances(vertex_, radius_factor_, nearest_vertices_distances_);
      return;
   }

  private:
   xDistanceNearestPoint_internal<ITER, GENERATOR, typename GENERATOR::Trait_t> internal;
};
}  // namespace xgeom

#endif
