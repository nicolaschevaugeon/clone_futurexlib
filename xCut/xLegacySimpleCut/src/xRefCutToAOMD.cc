/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include <iostream>
#include <sstream>
// aomd
#include "mEdge.h"
#include "mEntity.h"
#include "mFace.h"
#include "mTet.h"
#include "mVertex.h"
// xinterface
#include "xAOMDEntityUtil.h"
// xfem
#include "xElement.h"
#include "xLevelSet.h"
#include "xMesh.h"
#include "xRegion.h"
// xmapping
#include "xMappingBuilderHolder.h"
// xcut
#include "xPhysSurfByTagging.h"
#include "xRefCut.h"
#include "xRefCutToAOMD.h"

namespace xcut
{
using AOMD::mEdge;
using AOMD::mEntity;
using AOMD::mFace;
using AOMD::mTet;
using AOMD::mVertex;
// static numbering
// nota : if for some reason this become opslotte all Ref class are bugged
// these are from analysing AOMD
static char tria_edge_nodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
// these are simple copie of mTet.cc modified for generality purpose
// face for prisme or ... have more then 3 node : 4
// here we put 4 nodes for tetra if in future prisme or ... are implemented
// same thing for edges
static char tetra_edge_nodes[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
static char tetra_face_nodes[4][4] = {{0, 1, 2, -1}, {0, 1, 3, -1}, {1, 2, 3, -1}, {0, 2, 3, -1}};
static char tetra_face_edges[4][4] = {{0, 1, 2, -1}, {0, 4, 3, -1}, {1, 5, 4, -1}, {2, 3, 5, -1}};
// these are edge corespondance betewen  REF and AOMD
// static char edgea_edge_edge[1] ={0};
static char tria_edge_edge[3] = {0, 1, 2};
static char tetra_edge_edge[6] = {0, 1, 2, 4, 5, 3};
// these are face corespondance betewen  REF and AOMD
static char tetra_face_face[4] = {0, 1, 2, 3};

// local function
int nodeLabel(mEntity *e) { return e->getId(); }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xRefCutToAOMD class implementation
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////// Constructor
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xRefCutToAOMD::xRefCutToAOMD(const xfem::xLevelSet &ls_,
                             xinterface::aomd::xAttachedDataManagerAOMDGeneric<AOMD::mEntity, char> &entities_tag_,
                             xinterface::aomd::xAttachedDataManagerAOMDGeneric<AOMD::mEntity, char> &support_tag_,
                             xfem::xMesh *isomesh_, xfem::xEntityToEntity in_, xfem::xEntityToEntity out_,
                             bool keep_old_partition_, bool recursive_, std::vector<xPhysSurfByTagging *> &promotors_,
                             bool only_higher_dim_partition_, datamanager_t<AOMD::mEntity *> &_was_created_by,
                             datamanager_t<double> &_r_on_edge, datamanager_t<xfem::xMesh> &_partition,
                             datamanager_t<AOMD::mEntity *> &_was_duplicated_from,
                             datamanager_t<AOMD::mEntity *> &_is_in_partition_of)
    : ls(ls_),
      region(ls_.getSupport()),
      dim_region(region.dim()),
      entities_tag(entities_tag_),
      support_tag(support_tag_),
      nb_promotors(promotors_.size()),
      do_promote(false),
      promotors(promotors_),
      isomesh(isomesh_),
      classify_in(in_),
      classify_out(out_),
      keep_old_partition(keep_old_partition_),
      recursive(recursive_),
      only_higher_dim_partition(only_higher_dim_partition_),
      was_created_by(_was_created_by),
      r_on_edge(_r_on_edge),
      partition(_partition),
      was_duplicated_from(_was_duplicated_from),
      is_in_partition_of(_is_in_partition_of)

{
}
/////////////////////////////////////// End constructor
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Destructor
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xRefCutToAOMD::~xRefCutToAOMD()
{
   // cleaning iso_zero node tagging (xfem::xMesh::get_is_the_creator_of_tag()) as recursive or other object may have to do with
   // this tagging this permit createIsoZero function to work from a clean mesh
   // not really needed since is_the_creator of is cleaned up at destruction ...
   /*
   // loop on degre to clean all subentity
   for (int k = 0; k < dim_region; ++k)
   {
       for (mEntity *e : region.range(k))
       {
           // clean if tag attached
           is_the_creator_of.deleteData(*e);
           // for now sub subentities are not cleaned has only 2 level of recursion is possible => tag is_the_creator_of_tag may
   be present
           // in sub subentities only if a entity have been cut 2 times. Looking for tag in sub subentities means that it's the
   3trd times we try to cut !
           // to be checked

       }
   }
   */
}
/////////////////////////////////////// End Destructor
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Private methode
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void xRefCutToAOMD::setFunctionForEntity(const mEntity *e)
{
   ptfunc_cutElemRefByLevelSet = nullptr;
   ptfunc_createSubMeshFromRefCut = nullptr;
   ptfunc_createSubEntitiesSubMeshFromRefCut = nullptr;
   ptfunc_createSubEntitiesSubMeshFromRefCutByLoop = nullptr;
   ptfunc_createIsoZeroFromRefCut = nullptr;
   ptfunc_createIsoZeroFromFrontier = nullptr;
   ptfunc_cleanSubEntitiesClassification = nullptr;
   ptfunc_tagEntitiesOfElement = nullptr;

   switch (e->getType())
   {
         // 1D elements
      case mEntity::EDGE:
      {
         ptfunc_cutElemRefByLevelSet = &cutEdgeRefByLevelSet;
         ptfunc_createSubMeshFromRefCut = &xRefCutToAOMD::createSubMeshEdgeFromRef;
         ptfunc_createIsoZeroFromRefCut = &xRefCutToAOMD::createIsoZeroFromEdgeRefCut;
         ptfunc_cleanSubEntitiesClassification = &xRefCutToAOMD::cleanSubEntitiesClassificationEdge;
         ptfunc_createSubEntitiesSubMeshFromRefCut = &xRefCutToAOMD::createSubEntitiesSubMeshForEdgeFromRef;
         ptfunc_tagEntitiesOfElement = &xRefCutToAOMD::tagEntitiesOfElementEdge;
         ptfunc_createIsoZeroFromFrontier = &xRefCutToAOMD::createIsoZeroFromEdgeFrontier;
         ptfunc_createSubEntitiesSubMeshFromRefCutByLoop = &xRefCutToAOMD::createSubEntitiesSubMeshForEdgeFromRefCutByLoop;
         break;
      }
         // 2D elements
      case mEntity::TRI:
      {
         ptfunc_cutElemRefByLevelSet = &cutTriRefByLevelSet;
         ptfunc_createSubMeshFromRefCut = &xRefCutToAOMD::createSubMeshTriFromRef;
         ptfunc_createSubEntitiesSubMeshFromRefCut = &xRefCutToAOMD::createSubEntitiesSubMeshForTriFromRef;
         ptfunc_createSubEntitiesSubMeshFromRefCutByLoop = &xRefCutToAOMD::createSubEntitiesSubMeshForTriFromRefCutByLoop;
         ptfunc_createIsoZeroFromRefCut = &xRefCutToAOMD::createIsoZeroFromTriRefCut;
         ptfunc_createIsoZeroFromFrontier = &xRefCutToAOMD::createIsoZeroFromTriFrontier;
         ptfunc_cleanSubEntitiesClassification = &xRefCutToAOMD::cleanSubEntitiesClassificationTri;
         ptfunc_tagEntitiesOfElement = &xRefCutToAOMD::tagEntitiesOfElementTri;
         break;
      }
         // case mEntity::QUAD :
         //   break;
         // 3D elements
      case mEntity::TET:
      {
         ptfunc_cutElemRefByLevelSet = &cutTetRefByLevelSet;
         ptfunc_createSubMeshFromRefCut = &xRefCutToAOMD::createSubMeshTetFromRef;
         ptfunc_createSubEntitiesSubMeshFromRefCut = &xRefCutToAOMD::createSubEntitiesSubMeshForTetFromRef;
         ptfunc_createSubEntitiesSubMeshFromRefCutByLoop = &xRefCutToAOMD::createSubEntitiesSubMeshForTetFromRefCutByLoop;
         ptfunc_createIsoZeroFromRefCut = &xRefCutToAOMD::createIsoZeroFromTetRefCut;
         ptfunc_createIsoZeroFromFrontier = &xRefCutToAOMD::createIsoZeroFromTetFrontier;
         ptfunc_cleanSubEntitiesClassification = &xRefCutToAOMD::cleanSubEntitiesClassificationTet;
         ptfunc_tagEntitiesOfElement = &xRefCutToAOMD::tagEntitiesOfElementTet;
         break;
      }
         // case mEntity::PRISM :
         //    break;
         // case mEntity::HEX :
         //    break;
         // not implemented
      default:
         std::ostringstream oss;
         oss << " Wrong entity type : " << e->getType() << ". Not covered by this class.\n";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
         break;
   }
}

bool xRefCutToAOMD::getSubmesh(mEntity *e, xfem::xMesh **submesh)
{
   // by default it is considered as the first cut
   bool first_cut = true;
   // looking for a existing submesh
   // if a submesh exist 2 cases : recursive is true => cutting existing submesh
   //                              recursive is false => erasing and do as if it was the first cut

   if (!recursive)
   {
      xfem::xMesh *pMesh = partition.getData(*e);
      if (pMesh)
      {
         for (int d = 0; d <= 3; ++d)
         {
            for (auto e : pMesh->range(d)) cleanEltPartition(e);
         }
      }
      partition.deleteData(*e);
   }
   else if (partition.getData(*e))
      first_cut = false;
   *submesh = &partition.setData(*e);
   return first_cut;
}

void xRefCutToAOMD::createSubMeshEdgeFromRef(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *submesh,
                                             std::vector<mEntity *> &subentities_vect)
{
   // local
   int i;
   xRefMesh::elem_it it, itend;
   xRefMesh::point_it itp, itpend;
   mVertex *nodes[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};

   // Creation of nodes (vector of pointer to node on which sub elements are created )
   // they are duplicated from iso-zero node and nodes of e
   // iso-zero nodes are already created in iso_nodes
   itp = cutmesh.pointBegin();
   itpend = cutmesh.pointEnd();
   for (; itp != itpend; ++itp)
   {
      // node index
      i = itp->first;

      // if node is not from iso-zero nodes
      if (!iso_nodes[i])
         duplicateNode((mVertex *)subentities_vect[i], &nodes[i], submesh);
      else
         duplicateNode(iso_nodes[i], &nodes[i], submesh);
   }

   // if promotion of other level set tagging is to be donne
   if (nb_promotors) setPromotorsStatus(e);

   createSubEdge(e, cutmesh, submesh, &nodes[0]);
}

void xRefCutToAOMD::createSubMeshTriFromRef(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *submesh,
                                            std::vector<mEntity *> &subentities_vect)
{
   // nota ///////////////////////////////////////////////////////////
   // axiom : edge numbering betewn RefCut and AOMD are equivalente
   //         edge connectivity betewn RefCut and AOMD differs :
   //               edge 2 start from verice 2 to vertice 0 in AOMD
   //               edge 2 start from verice 0 to vertice 2 in RefCut
   // nota ///////////////////////////////////////////////////////////

   // local
   int i;
   xRefMesh::elem_it it, itend;
   xRefMesh::point_it itp, itpend;
   mVertex *nodes[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};

   // Creation of nodes (vector of pointer to node on which sub elements are created )
   // they are duplicated from iso-zero node and nodes of e
   // iso-zero nodes are already created in iso_nodes
   itp = cutmesh.pointBegin();
   itpend = cutmesh.pointEnd();
   for (; itp != itpend; ++itp)
   {
      // node index
      i = itp->first;

      // if node is not from iso-zero nodes
      if (!iso_nodes[i])
         duplicateNode((mVertex *)subentities_vect[i], &nodes[i], submesh);
      else
         duplicateNode(iso_nodes[i], &nodes[i], submesh);
   }

   // if promotion of other level set tagging is to be donne
   if (nb_promotors) setPromotorsStatus(e);

   createSubFace(e, cutmesh, submesh, &nodes[0]);
}

void xRefCutToAOMD::createSubMeshTetFromRef(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *submesh,
                                            std::vector<mEntity *> &subentities_vect)
{
   // local
   const bool debug = false;
   int i;
   double subvol = 0.0;
   mTet *tet_new;
   xRefMesh::elem_it it, itend;
   xRefMesh::point_it itp, itpend;
   mVertex *nodes[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
   AOMD::mMesh &submmesh = submesh->getMesh();
   // Creation of nodes (vector of pointer to node on which sub elements are created )
   // they are duplicated from iso-zero node and nodes of e
   // iso-zero nodes are already created in iso_nodes
   itp = cutmesh.pointBegin();
   itpend = cutmesh.pointEnd();
   for (; itp != itpend; ++itp)
   {
      // node index
      i = itp->first;

      // if node is not from iso-zero nodes
      if (!iso_nodes[i])
         duplicateNode((mVertex *)subentities_vect[i], &nodes[i], submesh);
      else
         duplicateNode(iso_nodes[i], &nodes[i], submesh);
   }

   // if promotion of other level set tagging is to be donne
   if (nb_promotors) setPromotorsStatus(e);

   // loop on element "in" to create sub elements
   it = cutmesh.elemBegin();
   itend = cutmesh.elemLim();
   for (; it != itend; it++)
   {
      // get connectivity
      xRefMesh::elemdef_t &conectivity = it->second;
      // create tetra
      tet_new = submmesh.createTetWithVertices(nodes[conectivity[0]], nodes[conectivity[1]], nodes[conectivity[2]],
                                               nodes[conectivity[3]], e->getClassification());
      is_in_partition_of.setData(*tet_new) = e;
      // classify "in"
      classify_in(tet_new);
      // tag the new entity
      entities_tag.setData(*tet_new) = LOOS_IN_STAT;
      //  promotion of other level set tagging
      if (do_promote) promoteStatusOfPromoters(tet_new);
      // sub volume calculation
      if (debug) subvol += xfem::xElement(tet_new).getVolume();
   }

   // loop on element "out" to create sub elements
   itend = cutmesh.elemIsozero();
   for (; it != itend; it++)
   {
      // get connectivity
      xRefMesh::elemdef_t &conectivity = it->second;
      // create tetra
      tet_new = submmesh.createTetWithVertices(nodes[conectivity[0]], nodes[conectivity[1]], nodes[conectivity[2]],
                                               nodes[conectivity[3]], e->getClassification());
      is_in_partition_of.setData(*tet_new) = e;
      // classify "out"
      classify_out(tet_new);
      // tag the new entity
      entities_tag.setData(*tet_new) = LOOS_OUT_STAT;
      //  promotion of other level set tagging
      if (do_promote) promoteStatusOfPromoters(tet_new);
      // sub volume calculation
      if (debug) subvol += xfem::xElement(tet_new).getVolume();
   }

   do_promote = false;
   // volume checking
   if (debug)
   {
      const double vol = xfem::xElement(e).getVolume();
      if (fabs(vol - subvol) / vol > 1.e-10)
      {
         throw;
      }
   }
}

void xRefCutToAOMD::createSubEntitiesSubMeshForEdgeFromRef(xRefMesh &cutmesh, std::vector<mEntity *> &subentities_vect)
{
   //
   // for Edge subentities are Nodes => creat nothing
   //
}

void xRefCutToAOMD::createSubEntitiesSubMeshForTriFromRef(xRefMesh &cutmesh, std::vector<mEntity *> &subentities_vect)
{
   //
   // for TRI subentities are edges => create edges
   //
   createSubEdgeFromRef(cutmesh, subentities_vect);
}

void xRefCutToAOMD::createSubEntitiesSubMeshForTetFromRef(xRefMesh &cutmesh, std::vector<mEntity *> &subentities_vect)
{
   //
   // for Tetrahedron subentities are faces and edges
   //
   // local
   int i, k;
   xRefMesh::elem_it it, itend;
   mFace *face_new;
   mVertex *nodes[4][8] = {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
                           {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
                           {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
                           {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};

   //  => create edges
   createSubEdgeFromRef(cutmesh, subentities_vect);

   //  => create faces
   //

   // if promotion of other level set tagging is to be donne
   if (nb_promotors) do_promote = true;

   // loop on face "in"
   it = cutmesh.faceBegin();
   itend = cutmesh.faceLim();
   for (; it != itend; ++it)
   {
      // get creator face number
      i = it->first % 10;

      // if submesh for this face have to be constructed
      if (face_submesh[i])
      {
         // get connectivity
         xRefMesh::elemdef_t &conectivity = it->second;
         //
         // sub face are always dicribed from cut node to original node (original=node of tetrahedron).
         // The last node is always from tetrahedron, the second may or may not commes from tetrahedron and first
         // is always from iso-zero.
         //
         // last
         k = conectivity[2];

         // if not already created in this submesh
         if (!nodes[i][k]) duplicateNode((mVertex *)subentities_vect[k], &nodes[i][k], face_submesh[i]);

         // midlle
         k = conectivity[1];

         // if not already created in this submesh
         if (!nodes[i][k])
         {
            // if it is a iso-zero
            if (k > 3) duplicateNode(iso_nodes[k], &nodes[i][k], face_submesh[i]);
            // if it's not a iso-zero
            else
               duplicateNode((mVertex *)subentities_vect[k], &nodes[i][k], face_submesh[i]);
         }

         // first
         k = conectivity[0];

         // if not already created in this submesh
         if (!nodes[i][k]) duplicateNode(iso_nodes[k], &nodes[i][k], face_submesh[i]);

         // create sub face
         face_new = face_submesh[i]->getMesh().createFaceWithVertices(nodes[i][conectivity[0]], nodes[i][conectivity[1]],
                                                                      nodes[i][conectivity[2]], faces[i]->getClassification());
         is_in_partition_of.setData(*face_new) = faces[i];

         // classify "in"
         classify_in(face_new);

         // tagging
         entities_tag.setData(*face_new) = LOOS_IN_STAT;

         //  promotion of other level set tagging
         if (do_promote)
         {
            setPromotorsStatus(faces[i]);
            promoteStatusOfPromoters(face_new);
         }
      }
   }

   // loop on face "out"
   itend = cutmesh.faceIsozero();
   for (; it != itend; ++it)
   {
      // get creator face number
      i = it->first % 10;

      // if submesh for this face have to be constructed
      if (face_submesh[i])
      {
         // get connectivity
         xRefMesh::elemdef_t &conectivity = it->second;
         //
         // sub face are always dicribed from cut node to original node (original=node of tetrahedron).
         // The last node is always from tetrahedron, the second may or may not commes from tetrahedron and first
         // is always from iso-zero.
         //
         // The "in" pass have construct duplicated iso-zero. This pass don't have to take care of them
         // only first and eventualy midlle node have to be duplicated
         //
         // last
         k = conectivity[2];

         // if not already created in this submesh
         if (!nodes[i][k]) duplicateNode((mVertex *)subentities_vect[k], &nodes[i][k], face_submesh[i]);

         // midlle
         k = conectivity[1];

         // if not iso and not already created in this submesh
         if (k < 4 && !nodes[i][k]) duplicateNode((mVertex *)subentities_vect[k], &nodes[i][k], face_submesh[i]);

         // create sub face
         face_new = face_submesh[i]->getMesh().createFaceWithVertices(nodes[i][conectivity[0]], nodes[i][conectivity[1]],
                                                                      nodes[i][conectivity[2]], faces[i]->getClassification());

         is_in_partition_of.setData(*face_new) = faces[i];
         // classify "out"
         classify_out(face_new);

         // tagging
         entities_tag.setData(*face_new) = LOOS_OUT_STAT;

         //  promotion of other level set tagging
         if (do_promote)
         {
            setPromotorsStatus(faces[i]);
            promoteStatusOfPromoters(face_new);
         }
      }
   }
}
void xRefCutToAOMD::createSubEdgeFromRef(xRefMesh &cutmesh, std::vector<mEntity *> &subentities_vect)
{
   //
   // local
   int i, k;
   xRefMesh::elem_it it, itend;
   mVertex *pnv;
   mEdge *edge_new;
   mVertex *nodes[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

   // if promotion of other level set tagging is to be donne
   if (nb_promotors) do_promote = true;

   // loop on edge "in"
   it = cutmesh.edgeBegin();
   itend = cutmesh.edgeLim();
   for (; it != itend; ++it)
   {
      // get creator edge number
      i = it->first % 10;

      // if submesh for this edge have to be constructed
      if (edge_submesh[i])
      {
         // get connectivity
         xRefMesh::elemdef_t &conectivity = it->second;
         //
         // sub edge are always dicribed from cut node to original node (original=node of entity e)
         // original apearse only one time per submesh so they are added directely
         duplicateNode((mVertex *)subentities_vect[conectivity[1]], &pnv, edge_submesh[i]);

         // cut node apearse 2 times per submesh in the general description of sub edges.
         // but it have to be created only one time
         // it is done during this "in" pass
         k = conectivity[0];

         // duplication from iso-zero
         duplicateNode(iso_nodes[k], &nodes[k], edge_submesh[i]);

         // create sub edge
         edge_new = edge_submesh[i]->getMesh().createEdge(nodes[k], pnv, edges[i]->getClassification());
         is_in_partition_of.setData(*edge_new) = edges[i];
         // classify "in"
         classify_in(edge_new);

         // tagging of submesh
         entities_tag.setData(*edge_new) = LOOS_IN_STAT;

         //  promotion of other level set tagging
         if (do_promote)
         {
            setPromotorsStatus(edges[i]);
            promoteStatusOfPromoters(edge_new);
         }
      }
   }

   // loop on edge "out"
   itend = cutmesh.edgeEnd();
   for (; it != itend; ++it)
   {
      // get creator edge number
      i = it->first % 10;

      // if submesh for this edge have to be constructed
      if (edge_submesh[i])
      {
         // get connectivity
         xRefMesh::elemdef_t &conectivity = it->second;
         //
         // sub edge are always dicribed from cut node to original node (original=node of entity e)
         // original apearse only one time per submesh so they are added directely
         duplicateNode((mVertex *)subentities_vect[conectivity[1]], &pnv, edge_submesh[i]);

         // cut node was created by the "in" pass so it is used directly

         // create sub edge
         edge_new = edge_submesh[i]->getMesh().createEdge(nodes[conectivity[0]], pnv, edges[i]->getClassification());
         is_in_partition_of.setData(*edge_new) = edges[i];

         // classify "out"
         classify_out(edge_new);

         // tagging of submesh
         entities_tag.setData(*edge_new) = LOOS_OUT_STAT;

         //  promotion of other level set tagging
         if (do_promote)
         {
            setPromotorsStatus(edges[i]);
            promoteStatusOfPromoters(edge_new);
         }
      }
   }

   do_promote = false;
}
void xRefCutToAOMD::createSubEntitiesSubMeshForEdgeFromRefCutByLoop(mEntity *e, xRefMesh &cutmesh, const xfem::xLevelSet &ls)
{
   //
   // for EDGE subentities are nodes => nothing to do
}
void xRefCutToAOMD::createSubEntitiesSubMeshForTriFromRefCutByLoop(mEntity *e, xRefMesh &cutmesh, const xfem::xLevelSet &ls)
{
   //
   // for TRI subentities are edges => create edges
   createSubEntitiesFromRefCutByLoop(e, cutmesh.edgeBegin(), cutmesh.edgeLim(), &edge_submesh_recursive[0], ls, 1,
                                     &tria_edge_edge[0], &cutEdgeRefByLevelSet, &xRefCutToAOMD::createSubEdge);

   // now set iso-zero reference for sub iso-mesh creation for entity as subenties have been treated
   iso_entity_ref = is_the_creator_of.at(*e);  // <- check ... may be a nullptr is acceptable.
}
void xRefCutToAOMD::createSubEntitiesSubMeshForTetFromRefCutByLoop(mEntity *e, xRefMesh &cutmesh, const xfem::xLevelSet &ls)
{
   //
   // for Tetrahedron subentities are edges and faces
   //
   // => create edges
   createSubEntitiesFromRefCutByLoop(e, cutmesh.edgeBegin(), cutmesh.edgeLim(), &edge_submesh_recursive[0], ls, 1,
                                     &tetra_edge_edge[0], &cutEdgeRefByLevelSet, &xRefCutToAOMD::createSubEdge);

   // => create faces
   createSubEntitiesFromRefCutByLoop(e, cutmesh.faceBegin(), cutmesh.faceLim(), &face_submesh_recursive[0], ls, 2,
                                     &tetra_face_face[0], &cutTriRefByLevelSet, &xRefCutToAOMD::createSubFace);

   // now set iso-zero reference for sub iso-mesh creation for entity as subenties have been treated
   iso_entity_ref = is_the_creator_of.at(*e);  // <- check ... may be a nullptr is acceptable.
}
void xRefCutToAOMD::createSubEntitiesFromRefCutByLoop(
    mEntity *e, xRefMesh::elem_it it, xRefMesh::elem_it itend, xfem::xMesh **entity_submesh_recursive, const xfem::xLevelSet &ls,
    int dim_sub, char *idx, int (*ptfunc_cutter)(const std::vector<double> &, const std::vector<int> &, xRefMesh &),
    void (xRefCutToAOMD::*ptfunc_createSubEntites)(mEntity *, xRefMesh &, xfem::xMesh *, mVertex **))
{
   //
   //  loop on cut entity to
   //            get parent ( entity of e)
   //            see if this parent entity have been treated or not by the normal procedure (no recursion)
   //            treat sub entities of parent entity if sub entity of sub entity is needed
   //
   //
   // here cut entities comes from cutmesh wich is a cut a the level of e. It don't give the right cut of sub entities of entities
   // of e. It just give a quike way to find what to do.
   //
   //
   // local
   int i;

   // loop on sub entities "in" of e
   // here loop on sub entities "out" of e is of no intrest : it will give the same sub entities of e to be treated
   // passing only on "in" sub entities leads to unique pass thru sub entities of e for edge
   // for face see coment below
   for (; it != itend; ++it)
   {
      // get creator entity number
      i = it->first % 10;

      // if submesh for this entity have not been treated by normal procedure
      if (entity_submesh_recursive[i])
      {
         // set iso-zero reference for sub iso-mesh creation for sub entity  i
         iso_entity_ref = is_the_creator_of.at(*e->get(dim_sub, idx[i]));  // <- check ... may be a nullptr is acceptable.

         // create subentities of submesh entities
         createSubSubEntitiesFromMultiRefCut(entity_submesh_recursive[i], ls, dim_sub, ptfunc_cutter, ptfunc_createSubEntites);

         // reset pointer as now sub entity "i" of e have been treated
         // this is for the case of face treatement where more then 1 sub faces may appears in "in" part
         entity_submesh_recursive[i] = nullptr;

      }  // endif entity have to be treated

   }  // end loop on sub entities "in" of e
}
void xRefCutToAOMD::createSubSubEntitiesFromMultiRefCut(xfem::xMesh *submesh, const xfem::xLevelSet &ls, int dim_sub,
                                                        int (*ptfunc_cutter)(const std::vector<double> &,
                                                                             const std::vector<int> &, xRefMesh &),
                                                        void (xRefCutToAOMD::*ptfunc_createSubEntites)(mEntity *, xRefMesh &,
                                                                                                       xfem::xMesh *, mVertex **))
{
   //
   // local
   char res;
   int i, k, d;
   double x, y, z, min, max;
   mEntity *se;
   mVertex *pviso = nullptr;
   // highest numbers of nodes is for subelement of type triangle cut on 2 edges : 3+2=5
   // this is for dim_sub>0
   mVertex *nodes[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
   xRefMesh::point_it itp, itpend;
   xRefMesh::elem_it it, itend;
   xfem::xMesh *sub_submesh;

   const double zero = 0.0;

   // list of values of level set at node
   std::vector<double> vals;

   // here we do a durty thing :
   //   we use a RefCut to cut the sub entities and create them (it look compliquate to use cutmesh to cut sub entities)
   //   it was done somehow like this in the old version
   //
   //
   // when called from main loop on element
   if (dim_sub < 1)
   {
      d = dim_region;
   }
   else
   {
      d = dim_sub;
   }
   submesh->getMesh().modifyState(3, 2, true);
   submesh->getMesh().modifyState(2, 3, true);
   submesh->getMesh().modifyState(2, 1, true);
   submesh->getMesh().modifyState(1, 2, true);
   submesh->getMesh().modifyState(1, 0, true);
   submesh->getMesh().modifyState(0, 1, true);

   //   creation of a mesh container to store a iso-zero mesh of the submesh : now this submesh is keepped (old version remove it)
   //   It is now a attachable mesh  attached to the corresponding entity of the global iso-zero mesh of this xPhysSurfByTagging
   xfem::xMesh *sub_isomesh;
   if (!(getSubmesh(iso_entity_ref, &sub_isomesh)))
   {
      std::ostringstream oss;
      oss << " Bugg : more then 1 level of recursion is impossible\nA sub iso mesh have already been attached to iso-zero entity "
             "choosen";
      throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set level set on submesh
   xfem::xLevelSet ls_submesh;
   ls.interpolateTo(*submesh, was_duplicated_from, was_created_by, ls_submesh);

   // to have unicity of node id generation by duplicateNode sub_isomesh node id generator have to be offseted
   // by simply using max nod id of submesh. Here no need to have unicity of node id in betewn sub-isomesh of 2
   // entity (tetra mainely) as only 2 level of recursion is possible
   // Unicity is here to avoid durty message from AOMD
   sub_isomesh->getMesh().setNewVertexId(submesh->getMesh().getNewVertexId());

   // loop on sub entities of curent entity via its submesh
   for (auto se : submesh->range(d))
   {
      // create a data container to retrive cut mesh definition in reference element
      // by declaring inside the loop this "empty" the container
      xRefMesh sub_cutmesh;
      // vector of subentities entity
      // by declaring inside the loop this "empty" the vector
      std::vector<mEntity *> sub_subentities_vect;

      // getting constitutve node of entity se
      xinterface::aomd::getSubentities(se, 0, sub_subentities_vect);

      // get value of level set on nodes of the elements
      vals = ls_submesh.getVals(sub_subentities_vect);

      // set node label vector of se
      std::vector<int> node_label(sub_subentities_vect.size());
      std::transform(sub_subentities_vect.begin(), sub_subentities_vect.end(), node_label.begin(), nodeLabel);

      // cut reference element of se and check if something was done
      if (ptfunc_cutter(vals, node_label, sub_cutmesh))
      {
         // force new tagging
         entities_tag.setData(*se) = CUT_STAT;

         // when called from a loop on subentities of reference element
         if (dim_sub > 0)
         {
            // get the submesh of the sub entity
            if (!(getSubmesh(se, &sub_submesh)))
            {
               std::ostringstream oss;
               oss << " Bugg : more then 1 level of recursion is impossible\nA sub entity have already a sub mesh";
               throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
            }

            // define a mapping to declare new nodes in physical coordinate
            xmapping::xMapping *mapping = xfem::xMappingBuilderHolderSingleton::instance().buildMapping(*se);

            // Minimal creation of "iso-zero" : mainly nodes are created as they may be used to retrive tag information of
            // their creators
            // Now was_duplicated_from_tag will always give a attached data correct : a node of sub-iso mesh wich is keeped as
            // attached mesh of entity of global iso-zero mesh
            //
            // for edge submesh
            if (dim_sub < 2)
            {
               itp = sub_cutmesh.pointBegin();
               itpend = sub_cutmesh.pointEnd();
               for (; itp != itpend; ++itp)
               {
                  // node index
                  k = itp->first;

                  // "iso-node" on edges
                  if (k > dim_sub)
                  {
                     // get new point coordinate in physical space
                     xPoint &point = itp->second;
                     mapping->eval(point.u, point.v, point.w, x, y, z);
                     // create point in iso-zero
                     pviso = sub_isomesh->getMesh().createVertex(x, y, z, se->getClassification());
                     was_created_by.setData(*pviso) = se;
                     entities_tag.setData(*pviso) = ISO_ZERO_STAT;
                     // duplicate node in sub submesh
                     duplicateNode(pviso, &nodes[k], sub_submesh);
                  }
                  // nodes from subentities
                  else
                  {
                     // force status, no check. for at most one node of the submesh this will be done twice
                     if (vals[k] < zero)
                        entities_tag.setData(*sub_subentities_vect[k]) = STRICT_IN_STAT;
                     else
                        entities_tag.setData(*sub_subentities_vect[k]) = STRICT_OUT_STAT;

                     // duplicate node in sub submesh
                     duplicateNode((mVertex *)sub_subentities_vect[k], &nodes[k], sub_submesh);
                  }
               }
            }
            // for face submesh
            else
            {
               itp = sub_cutmesh.pointBegin();
               itpend = sub_cutmesh.pointEnd();
               for (; itp != itpend; ++itp)
               {
                  // node index
                  i = itp->first;

                  // "iso-node" on edges
                  if (i > dim_sub)
                  {
                     // here 3 things are done with this loop on cut edge "in" :
                     //        * identify parent edge having this node
                     //        * create if necesary the node
                     //        * tag appartenance
                     // there is no need to loop on subedge "out" as they will give the same set of parent edge and iso-node
                     it = sub_cutmesh.edgeBegin();
                     itend = sub_cutmesh.edgeLim();
                     for (; it != itend; ++it)
                     {
                        // get edge conectivity
                        xRefMesh::elemdef_t &edge_conectivity = it->second;

                        // does the node considered is part of this edge
                        // sub edge are always dicribed from cut to original so only one test is needed
                        if (edge_conectivity[0] == i)
                        {
                           // get parent edge
                           k = it->first % 10;
                           // here tria_edge_edge is used directly has subface are considered to be triangle only. If quad are
                           // introduced the table have to become at least a parameter of the methode. Not realy importante here
                           // as REF and AOMD have the same definition for edge order of a triangle face. It is used to don't
                           // forget the eventuel mis correspondance in future use.
                           mEdge *pedge = (mEdge *)se->get(1, tria_edge_edge[k]);

                           // check if this edge is not already treated
                           // the presence of the iso-zero node attached to it indicate that this egde have already
                           // been cut by an other cut pass on adjancy sub element
                           // Here use of is_the_creator_of_tag is safe has edges were constructed just before this loop
                           // on sub faces by submesh->modifyState(2,1,true)
                           AOMD::mEntity *const *ppviso = is_the_creator_of.getData(*pedge);
                           if (ppviso)
                              pviso = static_cast<mVertex *>(*ppviso);
                           else  // if not already created the node is create in the sub iso-zero mesh
                           {
                              // node creation
                              xPoint &point = itp->second;
                              mapping->eval(point.u, point.v, point.w, x, y, z);
                              pviso = sub_isomesh->getMesh().createVertex(x, y, z, pedge->getClassification());

                              // tagging
                              is_the_creator_of.setData(*pedge) = pviso;
                              was_created_by.setData(*pviso) = pedge;

                              // as node has never been created befor it's parent edge was not tagged
                              // => tag it cut by iso-zero
                              entities_tag.setData(*pedge) = CUT_STAT;
                              // itself have to be taged to promote it's status to sub submesh naturaly with duplicateNode
                              entities_tag.setData(*pviso) = ISO_ZERO_STAT;
                           }

                           // stop loop on cut edges as the node have been treated
                           break;

                        }  // endif the node considered is part of this edge

                     }  // end loop on cut edges

                     // duplicate node in sub submesh
                     duplicateNode(pviso, &nodes[i], sub_submesh);
                  }
                  // nodes from subentities
                  else
                  {
                     // force status, no check. somme status will be set more then once but we don't test all node status
                     if (vals[i] < zero)
                        entities_tag.setData(*sub_subentities_vect[i]) = STRICT_IN_STAT;
                     else if (vals[i] > zero)
                        entities_tag.setData(*sub_subentities_vect[i]) = STRICT_OUT_STAT;
                     else
                        entities_tag.setData(*sub_subentities_vect[i]) = ISO_ZERO_STAT;

                     // duplicate node in sub submesh
                     duplicateNode((mVertex *)sub_subentities_vect[i], &nodes[i], sub_submesh);
                  }
               }
            }

            // clean
            delete mapping;

            // for faces, tagging of edges have to be finished
            // loop above only treat edges cut by iso zero. Not cutted edge have to be tagged
            if (dim_sub > 1)
            {
               // loop on edges of subfaces se
               for (i = 0; i < 3; ++i)
               {
                  mEdge *pedge = (mEdge *)se->get(1, i);
                  res = getTag(*pedge);
                  if (res == NO_STAT) setTagForTaggingEdgesOfElement(sub_subentities_vect, pedge, &tria_edge_nodes[i][0]);
               }
            }

            // if promotion of other level set tagging is to be donne
            if (nb_promotors) setPromotorsStatus(se);

            // creation of sub entities
            (this->*ptfunc_createSubEntites)(se, sub_cutmesh, sub_submesh, &nodes[0]);
         }
         // when called from main loop on element
         else
         {
            // create part of tempory iso-zero mesh from sub_cutmesh of the element se
            // only used to have consitant nodes creation and tagging for now
            // it give way to tag sub sub mesh when needed
            (*this.*ptfunc_createIsoZeroFromRefCut)(se, sub_cutmesh, sub_isomesh, sub_subentities_vect, vals);

            // get the submesh of the sub entity
            if (!(getSubmesh(se, &sub_submesh)))
            {
               std::ostringstream oss;
               oss << " Bugg : more then 1 level of recursion is impossible\nA sub entity have already a sub mesh";
               throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
            }

            // create sub element of se and classify from sub_cutmesh
            (*this.*ptfunc_createSubMeshFromRefCut)(se, sub_cutmesh, sub_submesh, sub_subentities_vect);
         }
      }
      // element se is not cut
      else
      {
         //    se is completely in or out
         //    or the new level set is passing on the frontier of se
         //
         //
         // Only classification and tagging have to be checked
         // for Classification test min on vals is simple and gives the right classification in any case :
         //    se is completely in or out : min gives a value strictely positive or negative
         //
         min = 1.0;
         max = -1.0;
         double val;

         // loop on nodes to force tag, it erase previous tag of nodes
         std::vector<mEntity *>::iterator itnbeg = sub_subentities_vect.begin();
         std::vector<mEntity *>::iterator itnend = sub_subentities_vect.end();
         std::vector<mEntity *>::iterator itn = itnbeg;
         for (; itn != itnend; ++itn)
         {
            val = vals[itn - itnbeg];
            if (val < min) min = val;
            if (val > max) max = val;
            if (val < zero)
               entities_tag.setData(*(*itn)) = STRICT_IN_STAT;
            else if (val > zero)
               entities_tag.setData(*(*itn)) = STRICT_OUT_STAT;
            else
               entities_tag.setData(*(*itn)) = ISO_ZERO_STAT;
         }

         // se is in
         if (min < zero)
         {
            classify_in(se);
            // se is strictly in
            if (max < zero)
            {
               entities_tag.setData(*se) = STRICT_IN_STAT;
            }
            // se is touching iso zero from in (max=zero no test as max>=0 and element not cut so max >0 imposible with min<0)
            else
               entities_tag.setData(*se) = LOOS_IN_STAT;
         }
         // se is strictly out
         else if (min > zero)
         {
            entities_tag.setData(*se) = STRICT_OUT_STAT;
            classify_out(se);
         }
         // se is touching iso zero from out (min=zero)
         else
         {
            entities_tag.setData(*se) = LOOS_OUT_STAT;
            classify_out(se);
         }

         // for faces, tagging of edges have to be done
         if (dim_sub > 1)
         {
            // loop on edges of subfaces se
            for (i = 0; i < 3; ++i)
            {
               mEdge *pedge = (mEdge *)se->get(1, i);
               res = getTag(*pedge);
               if (res == NO_STAT) setTagForTaggingEdgesOfElement(sub_subentities_vect, pedge, &tria_edge_nodes[i][0]);
            }
         }
      }

   }  // end loop on sub entities
}

void xRefCutToAOMD::createIsoZeroFromEdgeRefCut(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *iso_zero_mesh,
                                                std::vector<mEntity *> &subentities_vect, const std::vector<double> &vals)
{
   // nota ////////////////////////////////////////////////////////////
   // axiom : edge numbering betewn RefCut and AOMD are equivalente
   //         edge connectivity betewn RefCut and AOMD differs :
   //               edge 2 start from verice 2 to vertice 0 in AOMD
   //               edge 2 start from verice 0 to vertice 2 in RefCut
   // nota ///////////////////////////////////////////////////////////

   // local
   int i;
   xRefMesh::elem_it it;

   // define a mapping to declare new nodes in physical coordinate
   xmapping::xMapping *mapping = xfem::xMappingBuilderHolderSingleton::instance().buildMapping(*e);

   // reset comon members
   // nota : here 3 is ok as max edg... is 2
   for (i = 0; i < 6; ++i)
   {
      iso_nodes[i] = nullptr;
      edge_submesh[i] = nullptr;
      edge_submesh_recursive[i] = nullptr;
      edges[i] = nullptr;
   }

   // for edge only one iso-zero element : a point
   it = cutmesh.elemIsozero();

   // get connectivity
   xRefMesh::elemdef_t &conectivity = it->second;

   // create/retrive iso-zero node in iso_nodes member and init edge submes
   // cout << " EDGE case getIso" << endl;
   // A mesh node is cut
   if (conectivity[0] < 2)
   {
   }
   // The Edge is trully cut
   else
   {
      int i = conectivity[0];
      double x, y, z;
      xPoint &point = cutmesh.getPoint(i);
      mapping->eval(point.u, point.v, point.w, x, y, z);
      mVertex *pviso = iso_nodes[i] = iso_zero_mesh->getMesh().createVertex(x, y, z, e->getClassification());
      is_the_creator_of.setData(*e) = pviso;
      was_created_by.setData(*pviso) = e;
      r_on_edge.setData(*pviso) = point.u;
      entities_tag.setData(*e) = CUT_STAT;
      entities_tag.setData(*pviso) = ISO_ZERO_STAT;
   }
   //    getIsoNodesFromContext(e,conectivity,cutmesh,1,2,iso_zero_mesh, subentities_vect, mapping, &edgea_edge_edge[0]);

   // clean
   delete mapping;

   // finish tagging of nodes : getIsoNodesFromContext has taged iso-zero nodes
   tagNodes(vals, subentities_vect, cutmesh);
}

void xRefCutToAOMD::createIsoZeroFromTriRefCut(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *iso_zero_mesh,
                                               std::vector<mEntity *> &subentities_vect, const std::vector<double> &vals)
{
   // nota ////////////////////////////////////////////////////////////
   // axiom : edge numbering betewn RefCut and AOMD are equivalente
   //         edge connectivity betewn RefCut and AOMD differs :
   //               edge 2 start from verice 2 to vertice 0 in AOMD
   //               edge 2 start from verice 0 to vertice 2 in RefCut
   // nota ///////////////////////////////////////////////////////////

   // local
   int i;
   xRefMesh::elem_it it;

   // define a mapping to declare new nodes in physical coordinate
   xmapping::xMapping *mapping = xfem::xMappingBuilderHolderSingleton::instance().buildMapping(*e);

   // reset comon members
   // nota : here 5 is ok as max edg... is 6
   for (i = 0; i < 5; ++i)
   {
      iso_nodes[i] = nullptr;
      edge_submesh[i] = nullptr;
      edge_submesh_recursive[i] = nullptr;
      edges[i] = nullptr;
   }

   // for triangle only one iso-zero element : a edge
   it = cutmesh.elemIsozero();

   // get connectivity
   xRefMesh::elemdef_t &conectivity = it->second;

   // create/retrive iso-zero node in iso_nodes member and init edge submesh
   getIsoNodesFromContext(e, conectivity, cutmesh, 2, 3, iso_zero_mesh, subentities_vect, mapping, &tria_edge_edge[0]);
   // create iso-zero element in iso-zero mesh
   mEdge *edge_new =
       iso_zero_mesh->getMesh().createEdge(iso_nodes[conectivity[0]], iso_nodes[conectivity[1]], e->getClassification());
   was_created_by.setData(*edge_new) = e;
   // the edge is attached to the triangle for is_the_creator_of_tag
   // This is only used for sub iso-zero mesh attachment
   is_the_creator_of.setData(*e) = edge_new;
   // clean
   delete mapping;
   // finish tagging of nodes : getIsoNodesFromContext has taged iso-zero nodes
   tagNodes(vals, subentities_vect, cutmesh);
}
void xRefCutToAOMD::createIsoZeroFromTetRefCut(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *iso_zero_mesh,
                                               std::vector<mEntity *> &subentities_vect, const std::vector<double> &vals)
{
   // nota ////////////////////////////////////////////////////////////
   // axiom : face numbering betewn RefCut and AOMD are equivalent
   // nota ///////////////////////////////////////////////////////////

   // local
   int i, k;
   xRefMesh::elem_it it, itend;
   mFace *face_new = nullptr;
   mEdge *peiso = nullptr;

   // define a mapping to declare new nodes in physical coordinate
   xmapping::xMapping *mapping = xfem::xMappingBuilderHolderSingleton::instance().buildMapping(*e);

   // reset comon members
   iso_nodes.fill(nullptr);
   edge_submesh.fill(nullptr);
   edge_submesh_recursive.fill(nullptr);
   edges.fill(nullptr);
   face_submesh.fill(nullptr);
   face_submesh_recursive.fill(nullptr);
   faces.fill(nullptr);
   // for tetraedron one or two iso-zero element :  triangle(s)
   // loop on them
   it = cutmesh.elemIsozero();
   itend = cutmesh.elemEnd();
   for (i = 0; it != itend; ++it, ++i)
   {
      // get connectivity
      xRefMesh::elemdef_t &conectivity = it->second;

      // create/retrive iso-zero node in iso_nodes member and init edge submesh
      getIsoNodesFromContext(e, conectivity, cutmesh, 3, 4, iso_zero_mesh, subentities_vect, mapping, &tetra_edge_edge[0]);

      // create iso-zero element in iso-zero mesh
      face_new = iso_zero_mesh->getMesh().createFaceWithVertices(iso_nodes[conectivity[0]], iso_nodes[conectivity[1]],
                                                                 iso_nodes[conectivity[2]], e->getClassification());
      was_created_by.setData(*face_new) = e;

      // 2 triangles iso-zero mesh have to be completed by edge intersection
      if (i > 0)
      {
         // ref cut gives this edge directely : it's given by the first to nodes of the second triangle
         peiso =
             iso_zero_mesh->getMesh().createEdge(iso_nodes[conectivity[0]], iso_nodes[conectivity[1]], e->getClassification());
         was_created_by.setData(*peiso) = e;
      }
      // Arbitrary the first face is attached to the tetraedron for is_the_creator_of_tag
      // This is only used for sub iso-zero mesh attachment
      // As iso-zero have at least one face this procedure is safe
      // this is arbitrary but it respect the fact that attachable mesh of iso-zero face of the global iso-zero mesh are
      // sub iso-zero mesh of cut of sub tetrahedron
      // a alternative yould have been to use the edge in betwen them but in this case  attachable mesh of edge yould be of 2
      // differents nature : iso-zero mesh of cut of sub faces or iso-zero mesh of cut of sub tetrahedron ... rather confusing,
      // no?!
      else
         is_the_creator_of.setData(*e) = face_new;
   }

   // clean
   delete mapping;

   // finish tagging of nodes : getIsoNodesFromContext has taged iso-zero nodes
   tagNodes(vals, subentities_vect, cutmesh);

   // for tetraedron iso-zero edge have to be constructed
   // loop on them
   it = cutmesh.faceIsozero();
   itend = cutmesh.faceEnd();
   for (i = 0; it != itend; ++it, ++i)
   {
      // get parent number
      k = it->first % 10;

      // parent is a face
      if (k < 4)
      {
         // get parent face
         faces[k] = static_cast<mFace *>(e->get(2, k));

         // check if this face is not already treated
         // the presence of the iso-zero edge attached to it indicate that this egde have already
         // been cut by an other cut pass on adjancy element
         AOMD::mEntity *const *ppeiso = is_the_creator_of.getData(*faces[k]);
         if (ppeiso) peiso = static_cast<AOMD::mEdge *>(*ppeiso);
         // if not already created, the edge is create in the iso-zero mesh
         // it means also that submesh creation for this face have to be considered
         else
         {
            // get edge conectivity
            xRefMesh::elemdef_t &edge_conectivity = it->second;

            // edge creation
            AOMD::mEdge *peiso = iso_zero_mesh->getMesh().createEdge(iso_nodes[edge_conectivity[0]],
                                                                     iso_nodes[edge_conectivity[1]], e->getClassification());

            // tagging
            is_the_creator_of.setData(*faces[k]) = peiso;
            was_created_by.setData(*peiso) = faces[k];

            // as edge has never been created befor it's parent face was not tagged
            // => tag it cut by iso-zero
            entities_tag.setData(*faces[k]) = CUT_STAT;

            if (!only_higher_dim_partition)
            {
               // create submesh
               bool first_cut_face = getSubmesh(faces[k], &face_submesh[k]);

               // if a old submesh exist for this face
               if (!first_cut_face)
               {
                  // in recursive context it simply means that finaly this face don't have to be treated by
                  // the standard procedure as it have already been cut in a previous call to xPhysSurfByRef
                  if (recursive)
                  {
                     face_submesh_recursive[k] = face_submesh[k];
                     face_submesh[k] = nullptr;
                  }
                  // almost imposible as check on peiso above should avoid such case when not recursive
                  // leaved here for security
                  else
                  {
                     std::ostringstream oss;
                     oss << " Bugg : face entity is already cut says getSubmesh but it's not the case says iso-zero mesh";
                     throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
                  }
               }  // endif a old submesh exist for this face
            }     // endif only_higher_dim_partition
         }        // endif edge have to be created
      }
      // parent is a edge
      else
      {
         k -= 4;
         mEdge *parent_edge = static_cast<mEdge *>(e->get(1, tetra_edge_edge[k]));

         // check if this edge is not already treated
         // the presence of the iso-zero edge attached to it indicate that this egde have already
         // been created
         mEntity *const *ppeiso = is_the_creator_of.getData(*parent_edge);
         if (ppeiso) peiso = static_cast<mEdge *>(*ppeiso);
         // if not already created, the edge is create in the iso-zero mesh
         else
         {
            // get edge conectivity
            xRefMesh::elemdef_t &edge_conectivity = it->second;
            // edge creation
            mEdge *peiso = iso_zero_mesh->getMesh().createEdge(iso_nodes[edge_conectivity[0]], iso_nodes[edge_conectivity[1]],
                                                               e->getClassification());
            // tagging
            is_the_creator_of.setData(*parent_edge) = peiso;
            was_created_by.setData(*peiso) = parent_edge;
            // as it's the first time this parent edge is used it have to be tagged also
            entities_tag.setData(*parent_edge) = ISO_ZERO_STAT;
         }  // endif edge have to be created

      }  // end parent type

   }  // loop on iso-zero edge
}
void xRefCutToAOMD::createIsoZeroFromEdgeFrontier(mEntity *e, const bool is_out, const std::vector<double> &vals,
                                                  xfem::xMesh *iso_zero_mesh, std::vector<mEntity *> &subentities_vect)
{
   // not used except as fictuous argument
   char zero[1][4] = {{0, 0, 0}};
   char zero2[1][2] = {{0}};

   // call general methode with correct argument for triangle
   createIsoZeroFromEltFrontier(e, is_out, vals, iso_zero_mesh, subentities_vect, 2, 0, 0, 0, zero2, zero, zero);
}
void xRefCutToAOMD::createIsoZeroFromTriFrontier(mEntity *e, const bool is_out, const std::vector<double> &vals,
                                                 xfem::xMesh *iso_zero_mesh, std::vector<mEntity *> &subentities_vect)
{
   // not used except as fictuous argument
   char zero[1][4] = {{0, 0, 0}};

   // call general methode with correct argument for triangle
   createIsoZeroFromEltFrontier(e, is_out, vals, iso_zero_mesh, subentities_vect, 3, 3, 0, 0, tria_edge_nodes, zero, zero);
}
void xRefCutToAOMD::createIsoZeroFromTetFrontier(mEntity *e, const bool is_out, const std::vector<double> &vals,
                                                 xfem::xMesh *iso_zero_mesh, std::vector<mEntity *> &subentities_vect)
{
   // call general methode with correct argument for tetrahedron
   createIsoZeroFromEltFrontier(e, is_out, vals, iso_zero_mesh, subentities_vect, 4, 6, 4, 3, tetra_edge_nodes, tetra_face_nodes,
                                tetra_face_edges);
}
void xRefCutToAOMD::createIsoZeroFromEltFrontier(mEntity *e, const bool is_out, const std::vector<double> &vals,
                                                 xfem::xMesh *iso_zero_mesh, std::vector<mEntity *> &subentities_vect,
                                                 const char nb_nodes, const char nb_edges, const char nb_faces,
                                                 const char nb_edges_per_face, char(Tev)[][2], char(Tfv)[][4], char(Tfe)[][4])
{
   // local
   int i;
   char r = 0, res, new_tag;
   // char ress;
   mFace *face;
   mEdge *edge;
   mVertex *pv, *pviso;
   const double zero = 0.0;
   const char strict_tag = ((is_out) ? STRICT_OUT_STAT : STRICT_IN_STAT);
   const char loos_tag = ((is_out) ? LOOS_OUT_STAT : LOOS_IN_STAT);
   const char srel_tag = ((is_out) ? SREL_OUT_STAT : SREL_IN_STAT);
   AOMD::mMesh &iso_zero_mmesh = iso_zero_mesh->getMesh();
   // reinit
   for (i = 0; i < nb_nodes; ++i) iso_nodes[i] = nullptr;

   // loop on nodes ///////////////////////////////
   for (i = 0; i < nb_nodes; ++i)
   {
      pv = (mVertex *)subentities_vect[i];

      // the node is in iso-zero
      if (vals[i] == zero)
      {
         // track already created node
         mEntity *const *ppviso = is_the_creator_of.getData(*pv);
         if (ppviso)
         {
            pviso = static_cast<mVertex *>(*ppviso);
            iso_nodes[i] = pviso;
         }
         // if not created do so
         else
         {
            pviso = iso_nodes[i] = iso_zero_mmesh.createVertex(pv->getId(), pv->point(), pv->getClassification());
            was_created_by.setData(*pviso) = pv;
            is_the_creator_of.setData(*pv) = pviso;
            // as it has never been created befor it was not tagged
            // => tag it in iso-zero
            entities_tag.setData(*pv) = ISO_ZERO_STAT;
            // itself have to be taged to promote it's status to sub submesh naturaly with duplicateNode
            entities_tag.setData(*pviso) = ISO_ZERO_STAT;
         }
         // track number of nodes in iso-zero
         ++r;
      }
      // the node is not in iso-zero
      else
      {
         res = getTag(*pv);

         // if it is not already tagged
         if (res == NO_STAT) entities_tag.setData(*pv) = strict_tag;
         // check and if ok nothing to do
         else if (res != strict_tag)
         {
            std::ostringstream oss;
            oss << " Bugg : Inconsitant tagging , a node already taged as " << (int)res << " is now to be " << (int)strict_tag
                << " !";
            throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
         }
      }

      // In all cases (support in or out) its support is "related to" iso zero
      // if it is not already tagged as "cut by iso-zero Elt wise"    (SCUTEW_STAT win)
      // or if it is not already tagged as "cut by iso-zero"          (SCUT_STAT win)
      // or if it is not already tagged as "in or out touching iso-zero"
      // => it is tagged as "in or out touching iso-zero"
      if ((res = getSupportTag(*pv)) == NO_STAT)
      {
         support_tag.setData(*pv) = srel_tag;
         ;
      }
      // if it is allready tagged as "in or out touching iso-zero" check if
      // it is the same tag. If not it means that this entity have a support
      // with at least 2 elements having SREL_IN_STAT and SREL_OUT_STAT status
      // => it have to be tagged "cut by iso-zero"
      else if (res & TS_RELATED)
      {
         if (res != srel_tag) support_tag.setData(*pv) = SCUT_STAT;
         ;
      }

   }  // end loop on nodes //////////////////////////

   // check if element is not completely included in iso-zero
   if (r == nb_nodes)
   {
      std::ostringstream oss;
      oss << " Bugg : A element completely inside iso-zero is considered as a bug for now ! split element so that at least one "
             "of it's node is not in iso-zero";
      throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   if (nb_edges < 1) return;

   // loop on edges ///////////////////////////////
   mEntity **pe = &subentities_vect[nb_nodes];
   for (i = 0; i < nb_edges; ++i)
   {
      edge = (mEdge *)pe[i];

      res = getTag(*edge);
      // ress = xfem::getAttachedChar(edge,tag_support);

      const int i1 = Tev[i][0];
      const int i2 = Tev[i][1];

      // the edge is included in the iso-zero
      if (iso_nodes[i1] && iso_nodes[i2])
      {
         // if it was not already created
         // here use of tagging to take the decision : if it is not already
         // tagged it imply that it was not created befor
         if (res == NO_STAT)
         {
            // create iso-zero element in iso-zero mesh
            mEdge *edge_new = iso_zero_mmesh.createEdge(iso_nodes[i1], iso_nodes[i2], edge->getClassification());
            // tagging
            was_created_by.setData(*edge_new) = edge;
            new_tag = ISO_ZERO_STAT;
         }
         else
            new_tag = res;
      }
      // the edge is touching the iso-zero
      else if (iso_nodes[i1] || iso_nodes[i2])
      {
         new_tag = loos_tag;
      }
      // the edge is not touching nor in iso-zero
      // possible when only 1 (for tri and 2 for tetraedron) node(s) is in iso-zero
      else
      {
         new_tag = strict_tag;
      }

      // if it is not already tagged
      if (res == NO_STAT) entities_tag.setData(*edge) = new_tag;
      // check and if ok nothing to do
      else if (res != new_tag)
      {
         std::ostringstream oss;
         oss << " Bugg : Inconsitant tagging , a edge already taged as " << (int)res << " is now to be " << (int)new_tag << " !";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // In all cases (support in or out) its support is "related to" iso zero
      // if it is not already tagged as "cut by iso-zero Elt wise"    (SCUTEW_STAT win)
      // or if it is not already tagged as "cut by iso-zero"          (SCUT_STAT win)
      // or if it is not already tagged as "in or out touching iso-zero"
      // => it is tagged as "in or out touching iso-zero"
      if ((res = getSupportTag(*edge)) == NO_STAT)
      {
         support_tag.setData(*edge) = srel_tag;
         ;
      }
      // if it is allready tagged as "in or out touching iso-zero" check if
      // it is the same tag. If not it means that this entity have a support
      // with at least 2 elements having SREL_IN_STAT and SREL_OUT_STAT status
      // => it have to be tagged "cut by iso-zero"
      else if (res & TS_RELATED)
      {
         if (res != srel_tag) support_tag.setData(*edge) = SCUT_STAT;
         ;
      }

   }  // end loop on edges //////////////////////////

   if (nb_faces < 1) return;

   // loop on faces ///////////////////////////////
   pe = &subentities_vect[nb_nodes + nb_edges];
   for (i = 0; i < nb_faces; ++i)
   {
      face = (mFace *)pe[i];

      res = getTag(*face);
      // ress = xfem::getAttachedChar(edge,tag_support);

      const int i1 = Tfv[i][0];
      const int i2 = Tfv[i][1];
      const int i3 = Tfv[i][2];

      // the face is included in the iso-zero
      if (iso_nodes[i1] && iso_nodes[i2] && iso_nodes[i3])
      {
         // if it was not already created
         // here use of tagging to take the decision : if it is not already
         // tagged it imply that it was not created befor
         if (res == NO_STAT)
         {
            // create iso-zero element in iso-zero mesh
            mFace *face_new =
                iso_zero_mmesh.createFaceWithVertices(iso_nodes[i1], iso_nodes[i2], iso_nodes[i3], face->getClassification());
            // tagging
            was_created_by.setData(*face_new) = face;
            new_tag = ISO_ZERO_STAT;
         }
         else
            new_tag = res;
      }
      // the face is touching the iso-zero
      else if (iso_nodes[i1] || iso_nodes[i2] || iso_nodes[i3])
      {
         new_tag = loos_tag;
      }
      // the face is not touching nor in iso-zero
      // possible when only 1 (for tetraedron) node is in iso-zero
      else
      {
         new_tag = strict_tag;
      }

      // if entity is not already tagged
      if (res == NO_STAT) entities_tag.setData(*face) = new_tag;
      // check and if ok nothing to do
      else if (res != new_tag)
      {
         std::ostringstream oss;
         oss << " Bugg : Inconsitant tagging , a face already taged as " << (int)res << " is now to be " << (int)new_tag << " !";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // In all cases (support in or out) its support is "related to" iso zero
      // if it is not already tagged as "cut by iso-zero Elt wise"    (SCUTEW_STAT win)
      // or if it is not already tagged as "cut by iso-zero"          (SCUT_STAT win)
      // or if it is not already tagged as "in or out touching iso-zero"
      // => it is tagged as "in or out touching iso-zero"
      if ((res = getSupportTag(*face)) == NO_STAT)
      {
         support_tag.setData(*face) = srel_tag;
         ;
      }
      // if it is allready tagged as "in or out touching iso-zero" check if
      // it is the same tag. If not it means that this entity have a support
      // with at least 2 elements having SREL_IN_STAT and SREL_OUT_STAT status
      // => it have to be tagged "cut by iso-zero"
      else if (res & TS_RELATED)
      {
         if (res != srel_tag) support_tag.setData(*face) = SCUT_STAT;
         ;
      }

   }  // end loop on faces //////////////////////////

   return;
}
inline void xRefCutToAOMD::tagNodes(const std::vector<double> &vals, const std::vector<mEntity *> &subentities_vect,
                                    xRefMesh &cutmesh)
{
   int i;
   char res;
   const double zero = 0.0;

   // set iterator
   xRefMesh::point_it itp = cutmesh.pointBegin();
   xRefMesh::point_it itpend = cutmesh.pointEnd();

   // loop on nodes
   for (; itp != itpend; ++itp)
   {
      // node index
      i = itp->first;

      // if node is not from iso-zero nodes
      if (!iso_nodes[i])
      {
         char new_tag;
         res = getTag(*subentities_vect[i]);
         if (vals[i] < zero)
            new_tag = STRICT_IN_STAT;
         else
            new_tag = STRICT_OUT_STAT;
         if (res != new_tag)
         {
            if (res == NO_STAT)
               entities_tag.setData(*subentities_vect[i]) = new_tag;
            else
            {
               std::ostringstream oss;
               oss << " Bugg : Inconsitant tagging , a node already taged as " << (int)res << " is now to be " << (int)new_tag
                   << " !";
               throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
            }
         }
         // => tag it's support as cut
         support_tag.setData(*subentities_vect[i]) = SCUTEW_STAT;
         ;
      }
   }
}
inline void xRefCutToAOMD::tagSubSubEntities(mEntity *se, int d, char res, xfem::xEntityToEntity classify)
{
   const xfem::xMesh *sam = partition.getData(*se);
   if (sam)
   {
      xfem::xIter sseit = sam->begin(d);
      xfem::xIter sseitend = sam->end(d);
      for (; sseit != sseitend; ++sseit)
      {
         classify(*sseit);
         entities_tag.setData(**sseit) = res;
      }
   }
}
inline void xRefCutToAOMD::tagSubSubEntitiesTouching(mEntity *se, bool is_out, int n, int d, char res)
{
   std::vector<mEntity *> se_nodes_vect;
   xinterface::aomd::getSubentities(se, 0, se_nodes_vect);
   assert((int)(se_nodes_vect.size()) == n);
   bool is_touching = false;
   char res_new = res;
   char rese;
   int i;

   // loop on nodes to check touching or not
   // if one at least of the nodes of the sub entity is in iso-zero
   // it is touching. Otherwise this sub entity is not touching iso-zero
   for (i = 0; i < n; ++i)
   {
      // first get node wich was duplicate to create this one
      mEntity **pp = was_duplicated_from.getData(*se_nodes_vect[i]);
      if (pp)
      {
         mEntity *p = *pp;
         // second get the entity attached to this node by the tag was_created_by_tag
         mEntity *const *pp2 = was_created_by.getData(*p);

         // if something is attached to p it is a node from "iso-zero"
         if (pp2)
         {
            mEntity *p2 = *pp2;
            // if p2 is a node  it is from entity observed and its tagging is correct
            //   => his status is in iso-zero, the sub entitie is touching
            //   => his status is not iso-zero, the sub entitie se may be of the same status (strictly in or out), keep status in
            //   res_new
            if (p2->getLevel() == 0)
            {
               if ((res_new = getTag(*p2)) == ISO_ZERO_STAT)
               {
                  is_touching = true;
                  // tag the node according to the actuel level set
                  entities_tag.setData(*se_nodes_vect[i]) = res_new;
               }
            }
            // otherwise p2 is a edge of entity observed
            else
            {
               // the edge is included in iso-zero
               // => its node is
               if ((rese = getTag(*p2)) == ISO_ZERO_STAT)
               {
                  is_touching = true;
               }
               // the edge touching iso-zero
               // => its node won't touch it
               // rese have to be modified from loss to strict
               else if (rese & TE_LOOS)
               {
                  rese >>= 1;
               }
               // other wise
               // the edge is not touching iso-zero
               // => its node won't touch it
               // rese is tagged strict

               // rese is tagged strict and can give its value to res_new
               res_new = rese;

               // tag the node according to the actuel level set
               entities_tag.setData(*se_nodes_vect[i]) = res_new;
            }
         }
         // if nothing is attached to p it is a node from entity observed
         //   => his status is in iso-zero, the sub entitie is touching
         //   => his status is not iso-zero, the sub entitie se may be of the same status (strictly in or out), keep status in
         //   res_new
         else
         {
            if ((res_new = getTag(*p)) == ISO_ZERO_STAT)
            {
               is_touching = true;
               // tag the node according to the actuel level set
               entities_tag.setData(*se_nodes_vect[i]) = ISO_ZERO_STAT;
            }
         }
      }
      else
      {
         std::ostringstream oss;
         oss << " Bugg : nothing attached to node of submesh !";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   // after loop on node is_touching is correct and tagging can be done
   // if it is touching iso zero
   if (is_touching)
   {
      entities_tag.setData(*se) = res;

      // here we have to check sub sub entities for the extremely rare case of 2 level set cutting on element
      // and the 3th one (the one treated here) is touching this element
      // or 1 level set cut the element and the level-set treated here recursivly cut and touch it
      // interpolation have to be done
      const xfem::xMesh *sam = partition.getData(*se);
      if (sam)
      {
         xfem::xIter sseit = sam->begin(d);
         xfem::xIter sseitend = sam->end(d);
         for (; sseit != sseitend; ++sseit)
         {
            mEntity *sse = *sseit;
            if (is_out)
               classify_out(sse);
            else
               classify_in(sse);
            std::vector<mEntity *> sse_nodes_vect;
            xinterface::aomd::getSubentities(sse, 0, sse_nodes_vect);
            is_touching = false;
            res_new = res;
            for (i = 0; i < n; ++i)
            {
               // first get node wich was duplicate to create this one
               mEntity **pp = was_duplicated_from.getData(*sse_nodes_vect[i]);
               if (pp)
               {
                  mEntity *p = *pp;
                  // second get the entity attached to this node by the tag was_created_by_tag
                  mEntity *const *pp2 = was_created_by.getData(*p);

                  // if something is attached to p it is a node from "iso-zero"
                  if (pp2)
                  {
                     mEntity *p2 = *pp2;
                     // if p2 is a node  it is from se : it is tagged correctly from loop above on nodes of se
                     //   => his status is in iso-zero, the sub entitie is touching
                     //   => his status is not iso-zero, the sub entitie se may be of the same status (strictly in or out), keep
                     //   status in res_new
                     if (p2->getLevel() == 0)
                     {
                        if ((res_new = getTag(*p2)) == ISO_ZERO_STAT)
                        {
                           is_touching = true;
                           // tag the node according to the actuel level set
                           entities_tag.setData(*sse_nodes_vect[i]) = ISO_ZERO_STAT;
                        }
                     }
                     // otherwise p2 is a edge of se
                     else
                     {
                        rese = res;
                        // get edge node status
                        char res1 = getTag(*p2->get(0, 0));
                        char res2 = getTag(*p2->get(0, 1));
                        // the edge is included in iso-zero
                        // => its node is
                        if ((res1 & res2) & ISO_ZERO_STAT)
                        {
                           is_touching = true;
                        }
                        // the edge touching iso-zero
                        // => its node won't touch it
                        // or the edge is not touching iso-zero
                        // => its node won't touch it
                        // rese have to be modified from loss to strict in all case
                        else
                        {
                           rese >>= 1;
                        }
                        // rese is tagged strict and can give its value to res_new
                        res_new = rese;
                        // tag the node according to the actuel level set
                        entities_tag.setData(*sse_nodes_vect[i]) = rese;
                     }
                  }
                  // if nothing is attached to p it is a node from se : it is tagged correctly from loop above on nodes of se
                  //   => his status is in iso-zero, the sub entitie is touching
                  //   => his status is not iso-zero, the sub entitie sse may be of the same status (strictly in or out), keep
                  //   status in res_new
                  else
                  {
                     if ((res_new = getTag(*p)) == ISO_ZERO_STAT)
                     {
                        is_touching = true;
                        // tag the node according to the actuel level set
                        entities_tag.setData(*sse_nodes_vect[i]) = ISO_ZERO_STAT;
                     }
                  }
               }
               else
               {
                  std::ostringstream oss;
                  oss << " Bugg : nothing attached to node of submesh !";
                  throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
               }
            }
            // after loop on node is_touching is corect and tagging can be done
            // if it is touching iso zero
            if (is_touching) entities_tag.setData(*sse) = res;
            // if it is not touching iso zero
            else
               entities_tag.setData(*sse) = res_new;
         }  // end loop on sub sub entity
      }     // end if sub submesh exist
   }
   // if it is not touching iso zero
   else
   {
      entities_tag.setData(*se) = res_new;
      // here we have to check sub sub entities for the extremely rare case of 2 level set cutting on element
      // and the 3th one (the one treated here) is touching this element
      // or 1 level set cut the element and the level-set treated here recursivly cut and touch it
      // it is a more simple case then the is_touching case has here all sub sub entity are with the same
      // status res_new strictlly in or out
      if (is_out)
         tagSubSubEntities(se, d, res_new, classify_out);
      else
         tagSubSubEntities(se, d, res_new, classify_in);
   }
}

void xRefCutToAOMD::tagEntitiesOfElementEdge(mEntity *e, std::vector<mEntity *> &subentities_vect)
{
   // for Edge subenties are nodes
   // => tag nothing
}

void xRefCutToAOMD::tagEntitiesOfElementTri(mEntity *e, std::vector<mEntity *> &subentities_vect)
{
   // for triangle subenties are edges
   // => tag edges
   tagEdgesOfElement(e, subentities_vect, 3, tria_edge_edge, tria_edge_nodes);
}

void xRefCutToAOMD::tagEntitiesOfElementTet(mEntity *e, std::vector<mEntity *> &subentities_vect)
{
   // for tetrahedron subenties are egdges and faces
   // => tag edges
   tagEdgesOfElement(e, subentities_vect, 6, tetra_edge_edge, tetra_edge_nodes);
   // => tag faces
   tagFacesOfElement(e, subentities_vect, 4, tetra_face_face, tetra_face_nodes);
}
void xRefCutToAOMD::tagEdgesOfElement(mEntity *e, std::vector<mEntity *> &subentities_vect, int nb_edge, char(edge_edge)[],
                                      char(edge_nodes)[][2])
{
   // local
   int i, k;
   char rese;
   mEdge *edge;

   // loop on edge (folowing REF order to be consistante with edge_submesh and edge_submesh_recursive)
   for (i = 0; i < nb_edge; ++i)
   {
      // get edge index folowing AOMD order
      k = edge_edge[i];
      // get edge
      edge = (mEdge *)e->get(1, k);
      // get it's status
      rese = getTag(*edge);

      // if not tagged or to be tagged by a other procedure
      // and don't have a status yet
      if ((rese == NO_STAT) && (!edge_submesh[i]) && (!edge_submesh_recursive[i]))
      {
         setTagForTaggingEdgesOfElement(subentities_vect, edge, &edge_nodes[k][0]);
      }

      // in all case support tag have to be set
      support_tag.setData(*edge) = SCUTEW_STAT;
      ;
   }
}
inline void xRefCutToAOMD::setTagForTaggingEdgesOfElement(std::vector<mEntity *> &nodes_vect, mEntity *edge, char *edge_nodes)
{
   char res, new_tag;
   const int n1 = (int)(*edge_nodes);
   const int n2 = (int)(*(edge_nodes + 1));

   // use static numbering to get node tag
   res = getTag(*nodes_vect[n1]);
   if (res == ISO_ZERO_STAT)
   {
      res = getTag(*nodes_vect[n2]);
      if (res == STRICT_OUT_STAT)
         new_tag = LOOS_OUT_STAT;
      else if (res == STRICT_IN_STAT)
         new_tag = LOOS_IN_STAT;
      else if (res == ISO_ZERO_STAT)
         new_tag = ISO_ZERO_STAT;
      else
      {
         std::ostringstream oss;
         oss << " Bugg : A node is not tagged correctly has it sould be, it's tag value is " << (int)res << " ?";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }
   else
   {
      if (res == STRICT_OUT_STAT)
         new_tag = STRICT_OUT_STAT;
      else if (res == STRICT_IN_STAT)
         new_tag = STRICT_IN_STAT;
      else
      {
         std::ostringstream oss;
         oss << " Bugg : A node is not tagged correctly has it sould be, it's tag value is " << (int)res << " ?";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      res = getTag(*nodes_vect[n2]);
      if (res == ISO_ZERO_STAT) new_tag <<= 1;  // use the fact that loos=2*strict
   }

   // tag entity
   entities_tag.setData(*edge) = new_tag;
}

void xRefCutToAOMD::tagFacesOfElement(mEntity *e, std::vector<mEntity *> &subentities_vect, int nb_face, char(face_face)[],
                                      char(face_nodes)[][4])
{
   // local
   int i, k;
   char res, rese, new_tag;
   mFace *face;

   // loop on faces (folowing REF order to be consistante with face_submesh and face_submesh_recursive)
   for (i = 0; i < nb_face; ++i)
   {
      // get face index folowing AOMD order
      k = face_face[i];
      // get face
      face = (mFace *)e->get(2, k);

      // get it's status
      rese = getTag(*face);

      // if not tagged or to be tagged by a other procedure
      // and don't have a status yet
      if ((rese == NO_STAT) && (!face_submesh[i]) && (!face_submesh_recursive[i]))
      {
         // use static numbering to get node tag
         res = getTag(*subentities_vect[face_nodes[k][0]]);
         if (res == ISO_ZERO_STAT)
         {
            res = getTag(*subentities_vect[face_nodes[k][1]]);
            if (res == STRICT_OUT_STAT)
               new_tag = LOOS_OUT_STAT;
            else if (res == STRICT_IN_STAT)
               new_tag = LOOS_IN_STAT;
            else
            {
               res = getTag(*subentities_vect[face_nodes[k][2]]);
               if (res == STRICT_OUT_STAT)
                  new_tag = LOOS_OUT_STAT;
               else if (res == STRICT_IN_STAT)
                  new_tag = LOOS_IN_STAT;
               else
               {
                  std::ostringstream oss;
                  oss << " Bugg : Inconsitant tagging , a face not cut and not in iso-zero have 2 nodes tagged in iso-zero and "
                         "the other have "
                      << (int)res << " which is not strictly in or out !";
                  throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
               }
            }
         }
         else
         {
            if (res == STRICT_OUT_STAT)
               new_tag = STRICT_OUT_STAT;
            else if (res == STRICT_IN_STAT)
               new_tag = STRICT_IN_STAT;
            else
            {
               std::ostringstream oss;
               oss << " Bugg : A node is not tagged correctly has it sould be, it's tag value is " << (int)res << " ?";
               throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
            }
            res = getTag(*subentities_vect[face_nodes[k][1]]);
            if (res == ISO_ZERO_STAT)
               new_tag <<= 1;  // use the fact that loos=2*strict
            else
            {
               // no check on res : res sould be the same as for node 0

               res = getTag(*subentities_vect[face_nodes[k][2]]);
               if (res == ISO_ZERO_STAT) new_tag <<= 1;  // use the fact that loos=2*strict
            }
         }

         // tag entity
         entities_tag.setData(*face) = new_tag;

      }  // end if tagging to be done

      // in all case support tag have to be set
      support_tag.setData(*face) = SCUTEW_STAT;
      ;

   }  // end loop on face
}

void xRefCutToAOMD::cleanSubEntitiesClassificationEdge(mEntity *e, std::vector<mEntity *> &e_nodes_vect, const bool is_out,
                                                       const bool subcheck)
{
   // local
   char res;

   // if subentity exist for e clean classification
   const xfem::xMesh *am = partition.getData(*e);
   if (am)
   {
      // loop on sub entities of e via the submesh of e
      xfem::xIter seit = am->begin(1);
      xfem::xIter seitend = am->end(1);
      // For sub element when existing submesh exist it have to be re-tagged. But in the case
      // of e touching iso-zero sub-element have to be separated in not touching and touching iso.
      // For that use of tag of nodes
      if (subcheck)
      {
         if (is_out)
         {
            // set status of e from is_out
            res = LOOS_OUT_STAT;

            // loop on sub element
            for (; seit != seitend; ++seit)
            {
               // classify
               classify_out(*seit);
               // tag against curent level-set
               tagSubSubEntitiesTouching(*seit, is_out, 2, 1, res);
            }
         }
         else
         {
            // set status of e from is_out
            res = LOOS_IN_STAT;

            // loop on sub element
            for (; seit != seitend; ++seit)
            {
               // classify
               classify_in(*seit);
               // tag against curent level-set
               tagSubSubEntitiesTouching(*seit, is_out, 2, 1, res);
            }
         }

         // treat nodes of e as they may have a submesh with entity to be classify
         // Not done !!
         //            cleanSubEntitiesClassificationCheck(e,is_out,3,1,2);
      }
      // no checking as e is not touching iso-zero. Tag may be promoted as is
      else
      {
         // get status of e
         res = getTag(*e);

         if (is_out)
         {
            for (; seit != seitend; ++seit)
            {
               classify_out(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, 1, res, classify_out);
            }
         }
         else
         {
            for (; seit != seitend; ++seit)
            {
               classify_in(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, 1, res, classify_in);
            }
         }

         // treat nodes of e as they may have a submesh with entity to be classify
         // Not done !!
         //            cleanSubEntitiesClassificationCheck(e,is_out,3,1,2);
         // cleanSubEntitiesClassificationNoCheck(e,is_out,3,1);

      }  // end else no checking

   }  // end if subentity exist for e clean classification
}

void xRefCutToAOMD::cleanSubEntitiesClassificationTri(mEntity *e, std::vector<mEntity *> &e_nodes_vect, const bool is_out,
                                                      const bool subcheck)
{
   // local
   char res;

   // if subentity exist for e clean classification
   const xfem::xMesh *am = partition.getData(*e);
   if (am)
   {
      // loop on sub entities of e via the submesh of e
      xfem::xIter seit = am->begin(2);
      xfem::xIter seitend = am->end(2);
      // For sub element when existing submesh exist it have to be re-tagged. But in the case
      // of e touching iso-zero sub-element have to be separated in not touching and touching iso.
      // For that use of tag of nodes
      if (subcheck)
      {
         if (is_out)
         {
            // set status of e from is_out
            res = LOOS_OUT_STAT;

            // loop on sub element
            for (; seit != seitend; ++seit)
            {
               // classify
               classify_out(*seit);
               // tag against curent level-set
               tagSubSubEntitiesTouching(*seit, is_out, 3, 2, res);
            }
         }
         else
         {
            // set status of e from is_out
            res = LOOS_IN_STAT;

            // loop on sub element
            for (; seit != seitend; ++seit)
            {
               // classify
               classify_in(*seit);
               // tag against curent level-set
               tagSubSubEntitiesTouching(*seit, is_out, 3, 2, res);
            }
         }

         // treat edges of e as they may have a submesh with entity to be classify
         cleanSubEntitiesClassificationCheck(e, is_out, 3, 1, 2);
      }
      // no checking as e is not touching iso-zero. Tag may be promoted as is
      else
      {
         // get status of e
         res = getTag(*e);

         if (is_out)
         {
            for (; seit != seitend; ++seit)
            {
               classify_out(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, 2, res, classify_out);
            }
         }
         else
         {
            for (; seit != seitend; ++seit)
            {
               classify_in(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, 2, res, classify_in);
            }
         }

         // treat edges of e as they may have a submesh with entity to be classify
         cleanSubEntitiesClassificationNoCheck(e, is_out, 3, 1);

      }  // end else no checking

   }  // end if subentity exist for e clean classification
}

void xRefCutToAOMD::cleanSubEntitiesClassificationTet(mEntity *e, std::vector<mEntity *> &e_nodes_vect, const bool is_out,
                                                      const bool subcheck)
{
   // local
   char res;

   // if subentity exist for e clean classification

   const xfem::xMesh *am = partition.getData(*e);
   if (am)
   {
      // loop on sub entities of e via the submesh of e
      xfem::xIter seit = am->begin(3);
      xfem::xIter seitend = am->end(3);

      // For sub element when existing submesh exist it have to be re-tagged. But in the case
      // of e touching iso-zero sub-element have to be separated in stricly not touching and touching iso.
      // For that use of tag of nodes
      if (subcheck)
      {
         if (is_out)
         {
            // set status of e from is_out
            res = LOOS_OUT_STAT;

            // loop on sub element
            for (; seit != seitend; ++seit)
            {
               // classify
               classify_out(*seit);
               // tag against curent level-set
               tagSubSubEntitiesTouching(*seit, is_out, 4, 3, res);
            }
         }
         else
         {
            // set status of e from is_out
            res = LOOS_IN_STAT;

            // loop on sub element
            for (; seit != seitend; ++seit)
            {
               // classify
               classify_in(*seit);
               // tag against curent level-set
               tagSubSubEntitiesTouching(*seit, is_out, 4, 3, res);
            }
         }

         // treat faces of e as they may have a submesh with entity to be classify
         cleanSubEntitiesClassificationCheck(e, is_out, 4, 2, 3);

         // treat edges of e as they may have a submesh with entity to be classify
         cleanSubEntitiesClassificationCheck(e, is_out, 6, 1, 2);
      }
      // no checking as e is not touching iso-zero. Tag may be promoted as is
      else
      {
         // get status of e
         res = getTag(*e);

         if (is_out)
         {
            for (; seit != seitend; ++seit)
            {
               classify_out(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, 3, res, classify_out);
            }
         }
         else
         {
            for (; seit != seitend; ++seit)
            {
               classify_in(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, 3, res, classify_in);
            }
         }

         // treat faces of e as they may have a submesh with entity to be classify
         cleanSubEntitiesClassificationNoCheck(e, is_out, 4, 2);

         // treat edges of e as they may have a submesh with entity to be classify
         cleanSubEntitiesClassificationNoCheck(e, is_out, 6, 1);
      }

   }  // end if subentity exist for e clean classification
}
void xRefCutToAOMD::cleanSubEntitiesClassificationNoCheck(mEntity *e, const bool is_out, const int nb_subentities,
                                                          const int dim_sub)
{
   // local
   int i;
   char res;
   mEntity *sub_entity;

   // loop on sub entities (edge/face) of e as they may have a submesh with entity to be classify
   for (i = 0; i < nb_subentities; ++i)
   {
      sub_entity = (mEntity *)e->get(dim_sub, i);

      // if submesh exist for sub entity clean classification
      const xfem::xMesh *am = partition.getData(*sub_entity);
      if (am)
      {
         // get status of sub entity i
         res = getTag(*sub_entity);

         // loop on sub entities of sub entity via the submesh of sub entity
         xfem::xIter seit = am->begin(dim_sub);
         xfem::xIter seitend = am->end(dim_sub);
         if (is_out)
         {
            for (; seit != seitend; ++seit)
            {
               classify_out(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, dim_sub, res, classify_out);
            }
         }
         else
         {
            for (; seit != seitend; ++seit)
            {
               classify_in(*seit);
               entities_tag.setData(**seit) = res;
               // here we have to check sub sub entities !
               tagSubSubEntities(*seit, dim_sub, res, classify_in);
            }
         }
      }
   }
}
void xRefCutToAOMD::cleanSubEntitiesClassificationCheck(mEntity *e, const bool is_out, const int nb_subentities,
                                                        const int dim_sub, const int nb_nodes_per_subentities)
{
   // local
   int i;
   char res;
   mEntity *sub_entity;

   // loop on sub entities (edge/face) of e as they may have a submesh with entity to be classify
   for (i = 0; i < nb_subentities; ++i)
   {
      sub_entity = (mEntity *)e->get(dim_sub, i);

      // if submesh exist for sub entity clean classification
      const xfem::xMesh *am = partition.getData(*sub_entity);
      if (am)
      {
         // get status of sub entity i
         res = getTag(*sub_entity);

         // loop on sub entities of sub entity via the submesh of sub entity
         xfem::xIter seit = am->begin(dim_sub);
         xfem::xIter seitend = am->end(dim_sub);
         if (is_out)
         {
            for (; seit != seitend; ++seit)
            {
               classify_out(*seit);
               tagSubSubEntitiesTouching(*seit, is_out, nb_nodes_per_subentities, dim_sub, res);
            }
         }
         else
         {
            for (; seit != seitend; ++seit)
            {
               classify_in(*seit);
               tagSubSubEntitiesTouching(*seit, is_out, nb_nodes_per_subentities, dim_sub, res);
            }
         }
      }
   }
}
void xRefCutToAOMD::setPromotorsStatus(mEntity *e)
{
   for (int i = 0; i < nb_promotors; ++i)
   {
      (promotors[i])->setPromotorStatus(e);
   }

   do_promote = true;
}
void xRefCutToAOMD::promoteStatusOfPromoters(mEntity *e)
{
   for (int i = 0; i < nb_promotors; ++i)
   {
      (promotors[i])->promoteStatus(e);
   }
}
void xRefCutToAOMD::getIsoNodesFromContext(const mEntity *e, const xRefMesh::elemdef_t &iso_elem_conectivity, xRefMesh &cutmesh,
                                           const char &nb_node_per_iso_elem, const char &nb_node_per_elem,
                                           xfem::xMesh *iso_zero_mesh, const std::vector<mEntity *> &subentities_vect,
                                           xmapping::xMapping *mapping, const char *edge_ref_to_AOMD)
{
   // local
   int i, j, k;
   double x, y, z, one = 1.;
   xRefMesh::elem_it it, itend;
   mVertex *pv, *pviso;
   AOMD::mMesh &iso_zero_mmesh = iso_zero_mesh->getMesh();
   // loop on iso-zero node of curent iso-zero element
   for (j = 0; j < nb_node_per_iso_elem; ++j)
   {
      // if it's a node from  element
      if ((i = iso_elem_conectivity[j]) < nb_node_per_elem)
      {
         pv = static_cast<mVertex *>(subentities_vect[i]);
         // track already created node
         mEntity *const *ppviso = is_the_creator_of.getData(*pv);
         if (ppviso)
         {
            pviso = static_cast<mVertex *>(*ppviso);
            iso_nodes[i] = static_cast<mVertex *>(*ppviso);
         }
         // if not already created the node is create in the iso-zero mesh
         if (!ppviso)
         {
            mVertex *pviso = iso_nodes[i] = iso_zero_mmesh.createVertex(pv->getId(), pv->point(), pv->getClassification());
            was_created_by.setData(*pviso) = pv;
            is_the_creator_of.setData(*pv) = pviso;
            // as it has never been created befor element node was not tagged
            // => tag it in iso-zero
            entities_tag.setData(*pv) = ISO_ZERO_STAT;
            // itself have to be taged to promote it's status to submesh naturaly with duplicateNode
            entities_tag.setData(*pviso) = ISO_ZERO_STAT;
         }
         // For node element they may be SREL_IN/OUT_STAT  and already created in iso_zero mesh ....
         // In all case it's support have to be forced to SCUTEW_STAT
         // no preliminary test, just force. See if pre test give better performance
         support_tag.setData(*pv) = SCUTEW_STAT;
         ;
      }
      // if it's a node "on edge"
      else
      {
         // here 4 things are done with this loop on cut edge :
         //        * identify parent edge having this node
         //        * create if necesary the node
         //        * prepare submesh generation for edges
         //        * tag appartenance
         it = cutmesh.edgeBegin();
         itend = cutmesh.edgeEnd();
         for (; it != itend; ++it)
         {
            // get edge conectivity
            xRefMesh::elemdef_t &edge_conectivity = it->second;
            // does the node considered is part of this edge
            // sub edge are always dicribed from cut to original so only one test is needed
            if (edge_conectivity[0] == i)
            {
               // get parent edge
               k = it->first % 10;
               edges[k] = (mEdge *)e->get(1, edge_ref_to_AOMD[k]);
               // check if this edge is not already treated
               // the presence of the iso-zero node attached to it indicate that this edge have already
               // been cut by an other cut pass on adjancy element
               mEntity *const *ppviso = is_the_creator_of.getData(*edges[k]);
               if (ppviso)
               {
                  pviso = static_cast<mVertex *>(*ppviso);
                  iso_nodes[i] = static_cast<mVertex *>(pviso);
               }
               // if not already created the node is create in the iso-zero mesh
               // it means also that submesh creation for this edge have to be considered
               else
               {
                  // node creation
                  xPoint &point = cutmesh.getPoint(i);
                  mapping->eval(point.u, point.v, point.w, x, y, z);
                  pviso = iso_nodes[i] = iso_zero_mmesh.createVertex(x, y, z, edges[k]->getClassification());

                  // tagging
                  is_the_creator_of.setData(*edges[k]) = pviso;
                  was_created_by.setData(*pviso) = edges[k];
                  // as node has never been created befor it's parent edge was not tagged
                  // => tag it cut by iso-zero
                  entities_tag.setData(*edges[k]) = CUT_STAT;
                  // itself have to be taged to promote it's status to submesh naturaly with duplicateNode
                  entities_tag.setData(*pviso) = ISO_ZERO_STAT;

                  // store r (notation from xPhysSurf) which is here coordinate in reference element
                  // nota :
                  //   * for edge 2 reverse r values as conetivity differs betwen AOMD and RefCut
                  //            see tables Tev from mTet.cc for AOMD
                  //            and coment in xRefCut.cc for RefCut
                  // TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO ---
                  // TODO --- TODO --- TODO --- TODO --- TODO
                  //   * here 2D/3D are mixed together as onely triangle and tetraedron are treated so far
                  //      if other kind of element are implemented this switch have to be encapsulated in a
                  //      function given as a pointeur parameter or a policy class template parameter
                  //      If other element are compatible it may stays like that
                  // TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO ---
                  // TODO --- TODO --- TODO --- TODO --- TODO
                  //
                  switch (k)
                  {
                     case 0:
                     {
                        r_on_edge.setData(*pviso) = point.u;
                        break;
                     }
                     case 1:
                     {
                        r_on_edge.setData(*pviso) = point.v;
                        break;
                     }
                     case 2:
                     {
                        r_on_edge.setData(*pviso) = one - point.v;
                        break;
                     }
                     case 3:
                     case 4:
                     case 5:
                     {
                        r_on_edge.setData(*pviso) = point.w;
                        break;
                     }
                  }

                  // if the edge_submesh is not already created do so
                  if (!edge_submesh[k] && !only_higher_dim_partition)
                  {
                     bool first_cut_edge = getSubmesh(edges[k], &edge_submesh[k]);

                     // if a old submesh exist for this edge
                     if (!first_cut_edge)
                     {
                        // in recursive context it simply means that finaly this edge don't have to be treated by
                        // the standard procedure as it have already been cut in a previous call to xPhysSurfBytagging
                        if (recursive)
                        {
                           edge_submesh_recursive[k] = edge_submesh[k];
                           edge_submesh[k] = nullptr;
                        }
                        // almost imposible as check on pviso above should avoid such case when not recursive
                        // leaved here for security
                        else
                        {
                           std::ostringstream oss;
                           oss << " Bugg : edge entity is already cut says getSubmesh but it's not the case says iso-zero mesh";
                           throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
                        }
                     }
                  }
               }

               // stop loop on cut edges as the node have been treated
               break;

            }  // endif the node considered is part of this edge

         }  // end loop on cut edges

      }  // endif it's a node on edge

   }  // end loop on iso-zero nodes of the curent iso-zero element
}
inline void xRefCutToAOMD::duplicateNode(mVertex *pv, mVertex **pnv, xfem::xMesh *mesh)
{
   char res = getTag(*pv);
   // id is propagated to have a conforming mesh
   // using pv->getId() is ok as node may be of 3 kinds :
   //   - a node of a element : it's id is wat we are looking for
   //   - a iso-zero node duplicate from a element node : it's id is normaly comming from the duplicated node (see
   //   createIsoZeroFromEltFrontier,getIsoNodesFromContext)
   //   - a iso-zero node on edge : it's id is comming from iso-zero mesh or submesh. No choice here
   //
   mVertex *pnv_loc = *pnv = mesh->getMesh().createVertex(pv->getId(), pv->point(), pv->getClassification());
   was_duplicated_from.setData(*pnv_loc) = pv;
   if (res) entities_tag.setData(*pnv_loc) = res;
   /*if (nb_promotors){
     promoteStatusOfPromoters(pv);
     promoteStatusOfPromoters(pnv_loc);
     }*/
}
void xRefCutToAOMD::createSubFace(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *submesh, mVertex **nodes)
{
   // local
   const bool debug = false;
   mFace *face_new;
   xRefMesh::elem_it it, itend;
   double subvol = 0.0;

   // loop on element "in" to create sub elements
   it = cutmesh.elemBegin();
   itend = cutmesh.elemLim();
   for (; it != itend; it++)
   {
      // get connectivity
      xRefMesh::elemdef_t &conectivity = it->second;
      // create face
      face_new = submesh->getMesh().createFaceWithVertices(nodes[conectivity[0]], nodes[conectivity[1]], nodes[conectivity[2]],
                                                           e->getClassification());
      is_in_partition_of.setData(*face_new) = e;
      // classify "in"
      classify_in(face_new);
      // tag the new entity
      entities_tag.setData(*face_new) = LOOS_IN_STAT;
      //  promotion of other level set tagging
      if (do_promote) promoteStatusOfPromoters(face_new);
      // sub volume calculation
      if (debug) subvol += xfem::xElement(face_new).getVolume();
   }

   // loop on element "out" to create sub elements
   itend = cutmesh.elemIsozero();
   for (; it != itend; it++)
   {
      // get connectivity
      xRefMesh::elemdef_t &conectivity = it->second;
      // create face
      face_new = submesh->getMesh().createFaceWithVertices(nodes[conectivity[0]], nodes[conectivity[1]], nodes[conectivity[2]],
                                                           e->getClassification());

      is_in_partition_of.setData(*face_new) = e;
      // classify "out"
      classify_out(face_new);
      // tag the new entity
      entities_tag.setData(*face_new) = LOOS_OUT_STAT;
      //  promotion of other level set tagging
      if (do_promote) promoteStatusOfPromoters(face_new);
      // sub volume calculation
      if (debug) subvol += xfem::xElement(face_new).getVolume();
   }

   do_promote = false;
   // surface checking
   if (debug)
   {
      const double vol = xfem::xElement(e).getVolume();
      if (fabs(vol - subvol) / vol > 1.e-10) throw;
   }
}
void xRefCutToAOMD::createSubEdge(mEntity *e, xRefMesh &cutmesh, xfem::xMesh *submesh, mVertex **nodes)
{
   // local
   mEdge *edge_new_in, *edge_new_out;

   // creation of edge "in"
   {
      xRefMesh::elemdef_t &conectivity = cutmesh.getElem(0);
      edge_new_in = submesh->getMesh().createEdge(nodes[conectivity[0]], nodes[conectivity[1]], e->getClassification());
   }
   is_in_partition_of.setData(*edge_new_in) = e;
   classify_in(edge_new_in);
   entities_tag.setData(*edge_new_in) = LOOS_IN_STAT;

   // creation of edge "out"
   {
      xRefMesh::elemdef_t &conectivity = cutmesh.getElem(1);
      edge_new_out = submesh->getMesh().createEdge(nodes[conectivity[0]], nodes[conectivity[1]], e->getClassification());
   }
   is_in_partition_of.setData(*edge_new_out) = e;
   classify_out(edge_new_out);
   entities_tag.setData(*edge_new_out) = LOOS_OUT_STAT;
   //  promotion of other level set tagging
   if (do_promote)
   {
      promoteStatusOfPromoters(edge_new_in);
      promoteStatusOfPromoters(edge_new_out);
      do_promote = false;
   }
}
void xRefCutToAOMD::cleanAllPartition(xfem::xMesh *mesh)
{
   for (int i = mesh->dim(); i >= 0; --i)
      for (auto e : mesh->range(i)) cleanEltPartition(e);
}
void xRefCutToAOMD::cleanEltPartition(mEntity *e)
{
   if (xfem::xMesh *am = partition.getData(*e))
   {
      // loop on sub entities of e via the submesh of e
      for (int i = 0; i <= 3; ++i)
         for (auto e : am->range(i)) cleanEltPartition(e);
      // clean the sub mesh as now it have been cleaned
      partition.deleteData(*e);
   }
   // remove the creator
   was_created_by.deleteData(*e);
   // remove entity status associated to e
   entities_tag.deleteData(*e);
   // remove  support status associated to e
   support_tag.deleteData(*e);
}
/////////////////////////////////////// End Private methode
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Public methode
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void xRefCutToAOMD::cutMeshByRef()
{
   // local variable
   mEntity *e;
   xfem::xMesh *submesh;
   int i;
   char res;
   double min, max;
   const double zero = 0.0;
   bool is_out;

   // mesh from region
   xfem::xMesh *mesh = region.getMesh();

   // initialize node id generator of iso-zero mesh to a value greater then the maximum
   // node id of the mesh => this is mandatory in duplicateNode to generate unique node
   // number for submeshs
   // see remarks for recursive generation in createSubSubEntitiesFromMultiRefCut
   //
   // nota : here mesh->getNewVertexId increment the id generator of mesh which is not problematic excepte
   // if a huge number of call to this function is donne (huge here is 2147483647(INT_MAX 26/10/2011) - nb_node of the mesh )
   isomesh->getMesh().setNewVertexId(mesh->getMesh().getNewVertexId());

   // iterator
   xfem::xIter eit = region.begin(dim_region);
   xfem::xIter eitend = region.end(dim_region);

   // reinite tags if they existe
   // nota : sub entities have their tag modyfied, normaly, by folowing methode called in main loop
   // element are treated in the main loop
   // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---
   // This task wil be no more necessary as tag at this level will be temporary in multilevelset context.
   // At the end of there use they are removed so they don't exist in a recursive/update context.
   // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---
   for (i = 0; i < dim_region; ++i)
   {
      eitend = region.end(i);
      for (eit = region.begin(i); eit != eitend; ++eit)
      {
         mEntity *e = *eit;
         entities_tag.setData(*e) = NO_STAT;
         support_tag.deleteData(*e);
      }
   }

   // see if it canot be included in the loop above
   // clean classification
   mesh->cleanClassification();

   // One can ask to keep/destroy the old partition
   //
   // This is done for full compatibility with old xPhysSurf
   // it should be linked with recursive as from a logical point of view it is irelevant to ask for recursion and clean partition
   // at the same time !
   //
   if (!keep_old_partition) cleanAllPartition(mesh);

   // list of values of level set at node
   std::vector<double> vals;

   // ===============================================================
   // main loop on elements
   eitend = region.end(dim_region);
   for (eit = region.begin(dim_region); eit != eitend; ++eit)
   {
      // set element pointeur
      e = *eit;

      // vector of subentities entity
      // by declaring inside the loop this "empty" the vector
      std::vector<mEntity *> subentities_vect;

      // getting constitutve node of entity e
      xinterface::aomd::getSubentities(e, 0, subentities_vect);

      // getting value of the level set ls for all nodes of entity  e
      // vals is declared outside the loop and is normaly crached by object returned by getVals methode
      vals = ls.getVals(subentities_vect);

      // min max
      min = *min_element(vals.begin(), vals.end());
      max = *max_element(vals.begin(), vals.end());

      // if the element is related to iso-zero it have to be inspected
      // -------------------------------------------------------------
      if (min * max <= zero)
      {
         // set pointeur function for this type of entity
         setFunctionForEntity(e);

         // create a data container to retrive cut mesh definition in reference element
         // by declaring inside the loop this "empty" the container
         xRefMesh cutmesh;

         // set node label vector
         std::vector<int> node_label(subentities_vect.size());
         std::transform(subentities_vect.begin(), subentities_vect.end(), node_label.begin(), nodeLabel);

         // cut reference element of e and check if something was done
         if (ptfunc_cutElemRefByLevelSet(vals, node_label, cutmesh))
         {
            // the element is strictly cut
            // tag the entity
            entities_tag.setData(*e) = CUT_STAT;
            support_tag.setData(*e) = SCUTEW_STAT;
            ;

            // get the submesh and see if it is the first time this element is cut
            bool first_cut = getSubmesh(e, &submesh);

            // create part of iso-zero mesh from cutmesh of the element e
            (*this.*ptfunc_createIsoZeroFromRefCut)(e, cutmesh, isomesh, subentities_vect, vals);

            //                if(!only_higher_dim_partition){
            // create sub element of sub entities of e and classify, from cutmesh
            // nota : ptfunc_createIsoZeroFromRefCut have done the job of initializing information needed to
            // decides if a sub entities have to be treated or not. In recursive context only sub entities
            // that where not cut, are treated here. For allready cut entity a special procedure is done in
            // the recursive part below
            (*this.*ptfunc_createSubEntitiesSubMeshFromRefCut)(cutmesh, subentities_vect);

            // finalise tagging of subenties of e
            if (!only_higher_dim_partition) (*this.*ptfunc_tagEntitiesOfElement)(e, subentities_vect);
            //                }

            // if it is the first time this element is cut
            if (first_cut)
            {
               // create sub element of e and classify, from cutmesh
               (*this.*ptfunc_createSubMeshFromRefCut)(e, cutmesh, submesh, subentities_vect);
            }
            // recursive : it's the second time this element is cut
            else
            {
               // special treatement of sub entities of e :
               // if sub entity of e is already cut it have to be cut again.
               // Here we talk about sub entities of e and not sub entities of sub entities of e.
               // This confusion is a big temptation as recursion would have been then more simple/clear => same shemas for entity
               // and sub entity why not ? well from a boundary condition point of view we start from a edge and not from the
               // element having this edge. It is then difficulte to look down for partition if we don't start from the edge
               // excepte if we attache edge from sub element to edge from  parent element .... to think about
               (*this.*ptfunc_createSubEntitiesSubMeshFromRefCutByLoop)(e, cutmesh, ls);

               // create sub element from sub element of e
               createSubSubEntitiesFromMultiRefCut(submesh, ls, 0, ptfunc_cutElemRefByLevelSet, nullptr);

            }  // endif recursive or not
         }
         // element e is not cut but touch iso-zero
         else
         {
            // nothing done :  it's a element with one at least of it's sub entity in the iso-zero
            //
            // Classification have to be checked anyway

            // if the element is outside and touch the iso-zero
            if (min == zero)
            {
               is_out = true;

               // Classification
               classify_out(e);

               // tag the entity
               entities_tag.setData(*e) = LOOS_OUT_STAT;

               // The support of this entity which is itself is touching iso-zero
               // => it have to be tagged as "out touching iso zero"
               support_tag.setData(*e) = SREL_OUT_STAT;
               ;
            }
            // if the element is inside and touch the iso-zero
            else if (max == zero)
            {
               is_out = false;

               // Classification
               classify_in(e);

               // tag the entity
               entities_tag.setData(*e) = LOOS_IN_STAT;

               // The support of this entity which is itself is touching iso-zero
               // => it have to be tagged as "in touching iso zero"
               support_tag.setData(*e) = SREL_IN_STAT;
               ;
            }
            else
            {
               std::ostringstream oss;
               oss << " Bugg :: Entity  number " << e->getId()
                   << " is not cut and don't have any node in iso-zero which is imposible as min/max test leads to inspect this "
                      "element.\n";
               throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
            }

            // store nodes for cleanSubEntitiesClassification
            std::vector<mEntity *> nodes_vect;
            if (keep_old_partition) nodes_vect = subentities_vect;

            // get all subentities entity of e
            // nodes already donne
            // and itself not needed
            for (i = 1; i < dim_region; ++i) xinterface::aomd::getSubentities(e, i, subentities_vect);

            // Create iso-zero mesh  part from sub entity in the iso-zero
            (*this.*ptfunc_createIsoZeroFromFrontier)(e, is_out, vals, isomesh, subentities_vect);

            // depending on value of keep_old_partition cleaning classification of sub entities
            // is needed :  is_out give the way to classify sub entity
            if (keep_old_partition)
            {
               // clean
               (*this.*ptfunc_cleanSubEntitiesClassification)(e, nodes_vect, is_out, true);
            }
         }
      }
      // e is not related to iso-zero
      // ----------------------------
      else
      {
         // store nodes for cleanSubEntitiesClassification
         std::vector<mEntity *> nodes_vect;
         if (keep_old_partition) nodes_vect = subentities_vect;

         // get all subentities entity of e
         // nodes already donne
         // and itself not needed
         for (i = 1; i < dim_region; ++i) xinterface::aomd::getSubentities(e, i, subentities_vect);

         // set iterator
         std::vector<mEntity *>::iterator itc;
         std::vector<mEntity *>::iterator itcb = subentities_vect.begin();
         std::vector<mEntity *>::const_iterator itce = subentities_vect.end();

         // if the element is completely outside
         if (min > zero)
         {
            // tag the entity
            entities_tag.setData(*e) = STRICT_OUT_STAT;

            // Classification
            classify_out(e);

            // set bool for classification
            is_out = true;

            // support tag have to be reset
            // element is the support of itself and min>0 so NO_STAT
            support_tag.deleteData(*e);

            // loop on subentities entity
            for (itc = itcb; itc != itce; ++itc)
            {
               // get mEntity pointeur
               mEntity *c = *itc;

               res = getTag(*c);

               // entity allready marked as out :
               //   correctely tagged => nothing to do
               if (res & STRICT_OUT_STAT)
                  ;
               // entity allready marked but not strictly out
               // => bug !
               else if (res)
               {
                  std::ostringstream oss;
                  oss << " Bugg : entity is strictly out but it is already tagged as " << (int)res
                      << "  without having it's support touched or cut by iso-zero";
                  throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
               }
               // not tagged yet
               else
                  entities_tag.setData(*c) = STRICT_OUT_STAT;
            }
         }
         // if the element is completely inside
         // no other test as !(min.max<=0) & !(min>0)  ==> max < 0
         else
         {
            // tag the entity
            entities_tag.setData(*e) = STRICT_IN_STAT;

            // Classification
            classify_in(e);

            // set bool for classification
            is_out = false;

            // boundary tag have to be reset
            support_tag.deleteData(*e);

            // loop on subentities entity
            for (itc = itcb; itc != itce; ++itc)
            {
               // get mEntity pointeur
               mEntity *c = *itc;

               res = getTag(*c);

               // entity allready marked as in :
               //   correctely tagged => nothing to do
               if (res & STRICT_IN_STAT)
                  ;
               // entity allready marked but not strictly in
               // => bug !
               else if (res)
               {
                  std::ostringstream oss;
                  oss << " Bugg : entity is strictly in but it is already tagged as " << (int)res
                      << "  without having it's support touched or cut by iso-zero";
                  throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
               }
               // not tagged yet
               else
                  entities_tag.setData(*c) = STRICT_IN_STAT;
            }
         }

         // depending on value of keep_old_partition cleaning classification of sub entities
         // todo : see if it canot be mixed with loop above
         if (keep_old_partition)
         {
            // set pointeur function for this type of entity
            setFunctionForEntity(e);

            // clean
            (*this.*ptfunc_cleanSubEntitiesClassification)(e, nodes_vect, is_out, false);
         }
      }

   }  // end loop on element
   // ===============================================================

   return;
}
/////////////////////////////////////// End Public methode
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xRefCutToAOMDException class implementation
///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// general exception used for all methodes of xRefCutToAOMDException throw
xRefCutToAOMDException::xRefCutToAOMDException(std::string info, std::string file, int Line, std::string date, std::string time)
{
   std::ostringstream oss;
   oss << "In file " << file << " line " << Line << " compiled " << date << " at " << time << std::endl;
   oss << "xRefCutToAOMDException : " << info << std::endl;
   msg = oss.str();
}
/// general exception object : destructor
xRefCutToAOMDException ::~xRefCutToAOMDException() throw() = default;

/// mandatory what method
const char *xRefCutToAOMDException::what() const throw() { return this->msg.c_str(); }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xRefCutToAOMDException class implementation
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace xcut
