/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _DISTANCE_NEAREST_POINT_INTERNAL_IMP_H
#define _DISTANCE_NEAREST_POINT_INTERNAL_IMP_H

#include <cmath>

namespace xgeom
{
#ifdef HAVE_CGAL
// Partial specialisation on xCGALTrait

template <typename ITER, typename GENERATOR>
xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::xDistanceNearestPoint_internal(ITER begin_, ITER end_,
                                                                                            const GENERATOR &generator_,
                                                                                            MPI_Comm world_)
    : world(world_), generator(generator_)
{
   std::copy(begin_, end_, std::inserter(entities_vect, entities_vect.begin()));
   cgal_entities.reserve(entities_vect.size());
   for (; begin_ != end_; ++begin_)
   {
      generator.generateEntity(*begin_, cgal_entities);
   }
   // create temporary type for Gatherv operation
   /* A portable implementation should have created a new type but Kernel_t::Point_3,segment_3, ... are
    * rather opaque. Even with some investigation, work would have to be set in generator. To much
    * work for now. Just use generic binary type which is not fully portable (says if communications
    * are done between computer with different OS and/or hardware encoding this will be wrong).
    *
    * For now restricted to homogeneous material/OS
    *
    */
   int size_data = sizeof(typename GENERATOR::Container_t::value_type);

   // collect and compute nb of entities and displacement used for Gatherv operation
   int nb_entities = cgal_entities.size() * size_data;
   int nb_proc;
   MPI_Comm_size(world, &nb_proc);
   std::vector<int> all_nb_entities(nb_proc);
   MPI_Allgather(&nb_entities, 1, MPI_INT, &all_nb_entities[0], 1, MPI_INT, world);
   std::vector<int> nb_entities_disp(nb_proc);
   int k = 0;
   for (int i = 0; i < nb_proc; ++i)
   {
      nb_entities_disp[i] = k;
      k += all_nb_entities[i];
   }
   unsigned char *data = new unsigned char[k];
   int nb_final_entities = k / size_data;

   // Gatherv operation
   MPI_Allgatherv(&cgal_entities[0], nb_entities, MPI_BYTE, data, &all_nb_entities[0], &nb_entities_disp[0], MPI_BYTE, world);

   cgal_entities.clear();
   cgal_entities.reserve(nb_final_entities);
   for (int i = 0; i < k; i += size_data)
      cgal_entities.push_back(*reinterpret_cast<typename GENERATOR::Container_t::value_type *>(&data[i]));

   delete[] data;
   tree.rebuild(cgal_entities.begin(), cgal_entities.end());
   tree.accelerate_distance_queries();
   return;
}

template <typename ITER, typename GENERATOR>
xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::~xDistanceNearestPoint_internal()
{
   // nothing to do here //
   return;
}

template <typename ITER, typename GENERATOR>
double xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::distance(const typename GENERATOR::PointArg_t &point_)
{
   double sqd = (double)tree.squared_distance(generator.convertPoint(point_));
   return std::sqrt(sqd);
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::nearestPoint(const typename GENERATOR::PointArg_t &point_,
                                                                               typename GENERATOR::PointArg_t &nearest_point_)
{
   typename GENERATOR::PointRes_t cgal_point = generator.convertPoint(point_);
   nearest_point_ = generator.reverseConvertPoint(tree.closest_point(cgal_point));
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::nearestPointDistance(
    const typename GENERATOR::PointArg_t &point_, typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   typename GENERATOR::PointRes_t cgal_point = generator.convertPoint(point_);
   typename GENERATOR::PointRes_t cgal_closest_point = tree.closest_point(cgal_point);
   distance_ = std::sqrt((cgal_closest_point.x() - cgal_point.x()) * (cgal_closest_point.x() - cgal_point.x()) +
                         (cgal_closest_point.y() - cgal_point.y()) * (cgal_closest_point.y() - cgal_point.y()) +
                         (cgal_closest_point.z() - cgal_point.z()) * (cgal_closest_point.z() - cgal_point.z()));
   nearest_point_ = generator.reverseConvertPoint(cgal_closest_point);
   return;
}
template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::nearestEntityDistance(
    const typename GENERATOR::PointArg_t &point_, typename GENERATOR::EntityArg_t &nearest_entity_,
    typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   assert(xtool::workInProgress());

   // point from where to find
   typename GENERATOR::PointRes_t cgal_point = generator.convertPoint(point_);

   // query the tree
   typename CGALTree_t::Point_and_primitive_id res_pair = tree.closest_point_and_primitive(cgal_point);

   // nearest point
   typename GENERATOR::PointRes_t cgal_closest_point = res_pair.first;
   nearest_point_ = generator.reverseConvertPoint(cgal_closest_point);

   // distance
   distance_ = std::sqrt((cgal_closest_point.x() - cgal_point.x()) * (cgal_closest_point.x() - cgal_point.x()) +
                         (cgal_closest_point.y() - cgal_point.y()) * (cgal_closest_point.y() - cgal_point.y()) +
                         (cgal_closest_point.z() - cgal_point.z()) * (cgal_closest_point.z() - cgal_point.z()));

   // nearest entity
   // nota : both cgal_entities and entities_vect should be vector container for this to work
   nearest_entity_ = entities_vect[res_pair.second - cgal_entities.begin()];

   return;
}

template <typename ITER, typename GENERATOR>
double xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::distance(typename GENERATOR::VertexArg_t vertex_)
{
   double sqd = (double)tree.squared_distance(generator.convertVertex(vertex_));
   return std::sqrt(sqd);
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::nearestPoint(typename GENERATOR::VertexArg_t vertex_,
                                                                               typename GENERATOR::PointArg_t &nearest_point_)
{
   typename GENERATOR::PointRes_t cgal_point = generator.convertVertex(vertex_);
   nearest_point_ = generator.reverseConvertPoint(tree.closest_point(cgal_point));
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::nearestPointDistance(
    typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   typename GENERATOR::PointRes_t cgal_point = generator.convertVertex(vertex_);
   typename GENERATOR::PointRes_t cgal_closest_point = tree.closest_point(cgal_point);
   distance_ = std::sqrt((cgal_closest_point.x() - cgal_point.x()) * (cgal_closest_point.x() - cgal_point.x()) +
                         (cgal_closest_point.y() - cgal_point.y()) * (cgal_closest_point.y() - cgal_point.y()) +
                         (cgal_closest_point.z() - cgal_point.z()) * (cgal_closest_point.z() - cgal_point.z()));
   nearest_point_ = generator.reverseConvertPoint(cgal_closest_point);
   return;
}
template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xCGALTrait>::nearestEntityDistance(
    typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::EntityArg_t &nearest_entity_,
    typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   assert(xtool::workInProgress());

   // point from where to find
   typename GENERATOR::PointRes_t cgal_point = generator.convertVertex(vertex_);

   // query the tree
   typename CGALTree_t::Point_and_primitive_id res_pair = tree.closest_point_and_primitive(cgal_point);

   // nearest point
   typename GENERATOR::PointRes_t cgal_closest_point = res_pair.first;
   nearest_point_ = generator.reverseConvertPoint(cgal_closest_point);

   // distance
   distance_ = std::sqrt((cgal_closest_point.x() - cgal_point.x()) * (cgal_closest_point.x() - cgal_point.x()) +
                         (cgal_closest_point.y() - cgal_point.y()) * (cgal_closest_point.y() - cgal_point.y()) +
                         (cgal_closest_point.z() - cgal_point.z()) * (cgal_closest_point.z() - cgal_point.z()));

   // nearest entity
   // nota : both cgal_entities and entities_vect should be vector container for this to work
   nearest_entity_ = entities_vect[res_pair.second - cgal_entities.begin()];

   return;
}
#endif

#ifdef HAVE_ANN
// Partial specialisation on xANNTrait

template <typename ITER, typename GENERATOR>
xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::xDistanceNearestPoint_internal(ITER begin_, ITER end_,
                                                                                           const GENERATOR &generator_,
                                                                                           MPI_Comm world_)
    : tree(nullptr), generator(generator_)
{
   for (; begin_ != end_; ++begin_)
   {
      generator.generateEntity(*begin_, ann_vertices);
   }
   int max_pts = ann_vertices.size();
   data_pts = annAllocPts(max_pts, 3);
   int i = 0;
   typename std::vector<typename GENERATOR::VertexArg_t>::const_iterator it = ann_vertices.begin();
   typename std::vector<typename GENERATOR::VertexArg_t>::const_iterator end = ann_vertices.end();
   for (; it != end; ++it)
   {
      typename GENERATOR::PointRes_t p = generator.convertVertex(*it);
      data_pts[i] = p;
      ++i;
   }
   tree = new ANNkd_tree(data_pts, max_pts, 3);
   return;
}

template <typename ITER, typename GENERATOR>
xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::~xDistanceNearestPoint_internal()
{
   annDeallocPts(data_pts);
   delete tree;
   tree = nullptr;
   annClose();
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestVertexDistance_internal(
    typename GENERATOR::PointRes_t ann_point_, typename GENERATOR::VertexArg_t &nearest_vertex_, double &distance_)
{
   ANNidxArray ann_idx = new ANNidx[1];
   ANNdistArray ann_dist = new ANNdist[1];

   tree->annkSearch(ann_point_, 1, ann_idx, ann_dist);
   int i = ann_idx[0];
   distance_ = std::sqrt(ann_dist[0]);

   delete[] ann_idx;
   delete[] ann_dist;
   annDeallocPt(ann_point_);

   nearest_vertex_ = ann_vertices[i];
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestVerticesInsideRadiusDistances_internal(
    typename GENERATOR::PointRes_t ann_point_, double radius_factor_,
    std::map<typename GENERATOR::VertexArg_t, double> &nearest_vertices_distances_)
{
   // first, need to know the distance to the nearest point
   ANNidxArray ann_idx = new ANNidx[1];
   ANNdistArray ann_dist = new ANNdist[1];

   tree->annkSearch(ann_point_, 1, ann_idx, ann_dist);
   double squared_radius = ann_dist[0];

   delete[] ann_idx;
   delete[] ann_dist;

   // now, computes nearest neighbors in radius
   radius_factor_ *= radius_factor_;
   double squared_distance = squared_radius * radius_factor_;
   int number_of_points = tree->annkFRSearch(ann_point_, squared_distance, 0);

   ann_idx = new ANNidx[number_of_points];
   ann_dist = new ANNdist[number_of_points];
   tree->annkFRSearch(ann_point_, squared_distance, number_of_points, ann_idx, ann_dist);
   for (int k = 0; k < number_of_points; ++k)
   {
      nearest_vertices_distances_.insert(
          std::pair<typename GENERATOR::VertexArg_t, double>(ann_vertices[ann_idx[k]], std::sqrt(ann_dist[k])));
   }

   delete[] ann_idx;
   delete[] ann_dist;
   annDeallocPt(ann_point_);
   return;
}

template <typename ITER, typename GENERATOR>
double xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::distance(const typename GENERATOR::PointArg_t &point_)
{
   double distance;
   typename GENERATOR::VertexArg_t nearest_vertex = nullptr;
   typename GENERATOR::PointRes_t ann_point = generator.convertPoint(point_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance);
   return distance;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestPoint(const typename GENERATOR::PointArg_t &point_,
                                                                              typename GENERATOR::PointArg_t &nearest_point_)
{
   double distance;
   typename GENERATOR::VertexArg_t nearest_vertex = nullptr;
   typename GENERATOR::PointRes_t ann_point = generator.convertPoint(point_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance);
   nearest_point_ = generator.reverseConvertVertex(nearest_vertex);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestPointDistance(
    const typename GENERATOR::PointArg_t &point_, typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   typename GENERATOR::VertexArg_t nearest_vertex = nullptr;
   typename GENERATOR::PointRes_t ann_point = generator.convertPoint(point_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance_);
   nearest_point_ = generator.reverseConvertVertex(nearest_vertex);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestEntityDistance(
    const typename GENERATOR::PointArg_t &point_, typename GENERATOR::EntityArg_t &nearest_entity_,
    typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   typename GENERATOR::VertexArg_t nearest_vertex = nullptr;
   typename GENERATOR::PointRes_t ann_point = generator.convertPoint(point_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance_);
   nearest_point_ = generator.reverseConvertVertex(nearest_vertex);
   nearest_entity_ = static_cast<typename GENERATOR::EntityArg_t>(nearest_vertex);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestVerticesInsideRadiusDistances(
    const typename GENERATOR::PointArg_t &point_, double radius_factor_,
    std::map<typename GENERATOR::VertexArg_t, double> &nearest_vertices_distances_)
{
   typename GENERATOR::PointRes_t ann_point = generator.convertPoint(point_);
   nearestVerticesInsideRadiusDistances_internal(ann_point, radius_factor_, nearest_vertices_distances_);
   return;
}

template <typename ITER, typename GENERATOR>
double xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::distance(typename GENERATOR::VertexArg_t vertex_)
{
   double distance;
   typename GENERATOR::VertexArg_t nearest_vertex = nullptr;
   typename GENERATOR::PointRes_t ann_point = generator.convertVertex(vertex_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance);
   return distance;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestPoint(typename GENERATOR::VertexArg_t vertex_,
                                                                              typename GENERATOR::PointArg_t &nearest_point_)
{
   double distance;
   typename GENERATOR::VertexArg_t nearest_vertex = nullptr;
   typename GENERATOR::PointRes_t ann_point = generator.convertVertex(vertex_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance);
   nearest_point_ = generator.reverseConvertVertex(nearest_vertex);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestPointDistance(
    typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   typename GENERATOR::VertexArg_t nearest_vertex = nullptr;
   typename GENERATOR::PointRes_t ann_point = generator.convertVertex(vertex_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance_);
   nearest_point_ = generator.reverseConvertVertex(nearest_vertex);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestEntityDistance(
    typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::EntityArg_t &nearest_entity_,
    typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   typename GENERATOR::VertexArg_t nearest_vertex = NULL;
   typename GENERATOR::PointRes_t ann_point = generator.convertVertex(vertex_);
   nearestVertexDistance_internal(ann_point, nearest_vertex, distance_);
   nearest_point_ = generator.reverseConvertVertex(nearest_vertex);
   nearest_entity_ = static_cast<typename GENERATOR::EntityArg_t>(nearest_vertex);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xANNTrait>::nearestVerticesInsideRadiusDistances(
    typename GENERATOR::VertexArg_t vertex_, double radius_factor_,
    std::map<typename GENERATOR::VertexArg_t, double> &nearest_vertices_distances_)
{
   typename GENERATOR::PointRes_t ann_point = generator.convertVertex(vertex_);
   nearestVerticesInsideRadiusDistances_internal(ann_point, radius_factor_, nearest_vertices_distances_);
   return;
}
#endif

// Partial specialisation on xBruteForceTrait

template <typename ITER, typename GENERATOR>
xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::xDistanceNearestPoint_internal(ITER begin_, ITER end_,
                                                                                                  const GENERATOR &generator_,
                                                                                                  MPI_Comm world_)
    : generator(generator_)
{
   for (; begin_ != end_; ++begin_)
   {
      generator.generateEntity(*begin_, brute_force_edges);
   }
   return;
}

template <typename ITER, typename GENERATOR>
xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::~xDistanceNearestPoint_internal()
{
   // nothing to do here //
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::nearestPointDistance_internal(
    typename GENERATOR::PointRes_t P, typename GENERATOR::PointRes_t Cp, double &distance_)
{
   distance_ = std::numeric_limits<double>::max();
   const double eps = 1.e-10;
   int i = 0;
   int max = brute_force_edges.size();
   for (; i != max; i += 4)
   {
      const double A[2] = {brute_force_edges[i], brute_force_edges[i + 1]};
      const double B[2] = {brute_force_edges[i + 2], brute_force_edges[i + 3]};
      const double AB[2] = {B[0] - A[0], B[1] - A[1]};
      const double ABN = AB[0] * AB[0] + AB[1] * AB[1];
      const double AP[2] = {P[0] - A[0], P[1] - A[1]};
      const double AP_AB = AP[0] * AB[0] + AP[1] * AB[1];
      double u = (ABN < eps) ? 0. : AP_AB / ABN;
      u = u < 0. ? 0 : u;
      u = u > 1. ? 1. : u;
      double Cpe[2];
      Cpe[0] = A[0] + u * AB[0];
      Cpe[1] = A[1] + u * AB[1];
      const double diste = (Cpe[0] - P[0]) * (Cpe[0] - P[0]) + (Cpe[1] - P[1]) * (Cpe[1] - P[1]);
      if (diste < distance_)
      {
         Cp[0] = Cpe[0];
         Cp[1] = Cpe[1];
         distance_ = diste;
      }
   }
   distance_ = std::sqrt(distance_);
   return;
}

template <typename ITER, typename GENERATOR>
double xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::distance(const typename GENERATOR::PointArg_t &point_)
{
   double P[2] = {0., 0.};
   double Cp[2] = {0., 0.};
   double distance;
   generator.convertPoint(point_, P);
   nearestPointDistance_internal(P, Cp, distance);
   return distance;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::nearestPoint(
    const typename GENERATOR::PointArg_t &point_, typename GENERATOR::PointArg_t &nearest_point_)
{
   double P[2] = {0., 0.};
   double Cp[2] = {0., 0.};
   double distance;
   generator.convertPoint(point_, P);
   nearestPointDistance_internal(P, Cp, distance);
   generator.reverseConvertPoint(Cp, nearest_point_);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::nearestPointDistance(
    const typename GENERATOR::PointArg_t &point_, typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   double P[2] = {0., 0.};
   double Cp[2] = {0., 0.};
   generator.convertPoint(point_, P);
   nearestPointDistance_internal(P, Cp, distance_);
   generator.reverseConvertPoint(Cp, nearest_point_);
   return;
}

template <typename ITER, typename GENERATOR>
double xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::distance(typename GENERATOR::VertexArg_t vertex_)
{
   double P[2] = {0., 0.};
   double Cp[2] = {0., 0.};
   double distance;
   generator.convertVertex(vertex_, P);
   nearestPointDistance_internal(P, Cp, distance);
   return distance;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::nearestPoint(
    typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t &nearest_point_)
{
   double P[2] = {0., 0.};
   double Cp[2] = {0., 0.};
   double distance;
   generator.convertVertex(vertex_, P);
   nearestPointDistance_internal(P, Cp, distance);
   generator.reverseConvertPoint(Cp, nearest_point_);
   return;
}

template <typename ITER, typename GENERATOR>
void xDistanceNearestPoint_internal<ITER, GENERATOR, xBruteForceTrait>::nearestPointDistance(
    typename GENERATOR::VertexArg_t vertex_, typename GENERATOR::PointArg_t &nearest_point_, double &distance_)
{
   double P[2] = {0., 0.};
   double Cp[2] = {0., 0.};
   generator.convertVertex(vertex_, P);
   nearestPointDistance_internal(P, Cp, distance_);
   generator.reverseConvertPoint(Cp, nearest_point_);
   return;
}
}  // namespace xgeom

#endif
