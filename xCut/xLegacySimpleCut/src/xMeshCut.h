#ifndef XMESHCUT_H
#define XMESHCUT_H

// xFEM
#include "xEntityToEntity.h"
#include "xLevelSet.h"

#include "xMeshToolAOMD.h"
// AOMD interface
#include "xAttachedDataManagerAOMD.h"

namespace xcut
{
template <typename T>
using datamanager_t = xinterface::aomd::xAttachedDataManagerAOMD<T>;
/// Given a xLevelSet, this function construct inside imesh the mesh of the "intersection" of the ls support with the iso-zero.
/*! Given a xLevelSet ls, this function construct inside imesh the mesh of the "intersection" of the ls support with the iso-zero
 * Each entity created in imesh point to it's parent entity in the support of of ls. This data is stored in the data mananger
 * was_created_by. Each vertex of the imesh located on an edge of the ls support point to it's parametric coordinate on the edge
 * it's located on (between 0 and 1). This data is stored in r_on_edge. It is typically used in client code to reinterpolate
 * functions defined on the parent edge to the vertex of the imesh. Optional parameter simplex_only control if we split the
 * natural quadrangle obtain during the cut of a 3d element into 2 triangles or not. Optional parameter debug control if one wants
 * to output a verbose trace of the process to help debugging process. \brief createIsoZeroMeshFromLevelSet \param ls \param imesh
 * \param was_created_by
 * \param r_on_edge
 * \param simplex_only
 * \param debug
 */
void createIsoZeroMeshFromLevelSet(const xfem::xLevelSet &ls, xfem::xMesh &imesh, datamanager_t<AOMD::mEntity *> &was_created_by,
                                   datamanager_t<double> &r_on_edge, const bool simplex_only = true, const bool debug = false);

inline void createIsoZeroMeshFromLevelSet(const xfem::xLevelSet &ls, xfem::xMesh &imesh, const bool simplex_only = true,
                                          const bool debug = false)
{
   datamanager_t<AOMD::mEntity *> was_created_by;
   datamanager_t<double> r_on_edge;
   createIsoZeroMeshFromLevelSet(ls, imesh, was_created_by, r_on_edge, simplex_only, debug);
}

/// given a mesh "interface" typically obtained with  createIsoZeroMeshFromLevelSet, that represent a surface (resp a curve)
/// embedded in a
/*! background 3D (resp 2D) mesh, and the datamanager was_created_by that link each entity of the interface to it's parent entity
 * in the background mesh, this function create a partition for each parententity in the form of a xMesh, stored in partition.
 * Each entity in a given partition know of the parent entity it belong to. This data is stored in is_in_partition_of. Each vertex
 * in a partiton is either a copy of a vertex of the interface or a copy of a vertex of a parent entity. This relation is stored
 * in was_duplicated_from. The reverse relation is stored in is_duplicated_in
 */
void cutAlongInterface(const xfem::xMesh &interface, const datamanager_t<AOMD::mEntity *> &was_created_by,
                       datamanager_t<xfem::xMesh> &partition, datamanager_t<AOMD::mEntity *> &is_duplicated_in,
                       datamanager_t<AOMD::mEntity *> &was_duplicated_from, datamanager_t<AOMD::mEntity *> &is_in_partition_of);

void cutAlongInterfaceRecursive(const xfem::xMesh &interface, const xfem::xLevelSet &ls,
                                const datamanager_t<AOMD::mEntity *> &was_created_by, datamanager_t<xfem::xMesh> &partition,
                                datamanager_t<AOMD::mEntity *> &is_duplicated_in,
                                datamanager_t<AOMD::mEntity *> &was_duplicated_from,
                                datamanager_t<AOMD::mEntity *> &is_in_partition_of);

void cutMesh(const xfem::xMesh &base_mesh, const xfem::xLevelSet &ls, xfem::xMesh &iso_zero_mesh,
             datamanager_t<AOMD::mEntity *> &was_created_by, datamanager_t<double> &r_on_edge,
             datamanager_t<AOMD::mEntity *> &is_duplicated_in, datamanager_t<AOMD::mEntity *> &was_duplicated_from,
             datamanager_t<xfem::xMesh> &partition, datamanager_t<AOMD::mEntity *> &is_in_partition_of,
             xfem::xEntityToEntity classify_in, xfem::xEntityToEntity classify_out, bool create_partition,
             bool keep_old_partition, bool recursive);

void cutMesh(const xfem::xMesh &base_mesh, const xfem::xLevelSet &ls, xfem::xMesh &iso_zero_mesh,
             xfem::xEntityToEntity classify_in, xfem::xEntityToEntity classify_out, bool create_partition,
             bool keep_old_partition, bool recursive);

/// get creator recursivelly.
AOMD::mEntity *getOriginalCreator(AOMD::mEntity *egros, const datamanager_t<AOMD::mEntity *> &was_created_by);
///                          the elements of the mesh are transformed into TRIiangles and TETrahedra elements
///                          with inheritance of tags. Division is made by cutting faces on the diagonal which
///                          contains the node with the lowest Id. In this way the mesh remains conform.
}  // namespace xcut
#endif  // XMESHCUT_H
