/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef ___XREFCUTTOAOMD__H
#define ___XREFCUTTOAOMD__H
// xinterface
#include "xAttachedDataManagerAOMD.h"
// xfem
#include "xEntityToEntity.h"
// xcut
#include "xRefMesh.h"
// xmapping
#include "xMapping.h"

namespace AOMD
{
class mEntity;
class mVertex;
class mEdge;
class mFace;
}  // namespace AOMD
namespace xfem
{
class xLevelSet;
}

namespace xcut
{
class xPhysSurfByTagging;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xRefCutToAOMD class
///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class xRefCutToAOMD
{
   template <typename T>
   using datamanager_t = xinterface::aomd::xAttachedDataManagerAOMD<T>;

  public:
   xRefCutToAOMD(const xfem::xLevelSet &ls_,
                 xinterface::aomd::xAttachedDataManagerAOMDGeneric<AOMD::mEntity, char> &entities_tag_,
                 xinterface::aomd::xAttachedDataManagerAOMDGeneric<AOMD::mEntity, char> &support_tag_, xfem::xMesh *isomesh_,
                 xfem::xEntityToEntity in_, xfem::xEntityToEntity out_, bool keep_old_partition_, bool recursive_,
                 std::vector<xPhysSurfByTagging *> &promotors_, bool only_higher_dim_partition,
                 datamanager_t<AOMD::mEntity *> &was_created_by, datamanager_t<double> &r_on_edge,
                 datamanager_t<xfem::xMesh> &partition, datamanager_t<AOMD::mEntity *> &was_duplicated_from,
                 datamanager_t<AOMD::mEntity *> &is_in_partition_of);
   virtual ~xRefCutToAOMD();

   // public methodes ///////////////
   /// this is the equivalent of xfem::xMesh::cutMesh
   // it does the same things but differentelly :
   //    creat the iso-zero mesh
   //    cut all entities of the mesh cutted by the iso-zero and attache a submesh to it
   //    accepte one level of recursion for the cut
   //    classify all cutted entities
   void cutMeshByRef();

   // public members ////////////////
  protected:
   // private members ////////////////

   // levelSet mesh and dimension
   const xfem::xLevelSet &ls;
   const xfem::xRegion &region;
   int dim_region;

   // xPhysSurfByTagging tag
   xinterface::aomd::xAttachedDataManagerAOMDGeneric<AOMD::mEntity, char> &entities_tag;
   // short cut to get the entities tag. If not in the Data Manager return NO_STAT (0)
   char getTag(AOMD::mEntity &e) const
   {
      const char *ptag = entities_tag.getData(e);
      return ptag ? *ptag : 0;
   }

   xinterface::aomd::xAttachedDataManagerAOMDGeneric<AOMD::mEntity, char> &support_tag;
   // short cut to get the supprot tag. If not in the Data Manager return NO_STAT (0)
   char getSupportTag(AOMD::mEntity &e) const
   {
      const char *ptag = support_tag.getData(e);
      return ptag ? *ptag : 0;
   }

   // xPhysSurfByTagging promotos
   int nb_promotors;
   bool do_promote;
   std::vector<xPhysSurfByTagging *> &promotors;

   // iso-zero submesh
   xfem::xMesh *isomesh;

   // classifyers
   xfem::xEntityToEntity classify_in;
   xfem::xEntityToEntity classify_out;

   // boleen argument
   bool keep_old_partition;
   bool recursive;
   bool only_higher_dim_partition;
   datamanager_t<AOMD::mEntity *> &was_created_by;
   datamanager_t<double> &r_on_edge;
   datamanager_t<xfem::xMesh> &partition;
   datamanager_t<AOMD::mEntity *> &was_duplicated_from;
   datamanager_t<AOMD::mEntity *> &is_in_partition_of;

   // data structure to keep track of the entities created by those from the cut region
   datamanager_t<AOMD::mEntity *> is_the_creator_of;
   // pointer to function
   int (*ptfunc_cutElemRefByLevelSet)(const std::vector<double> &, const std::vector<int> &, xRefMesh &);
   void (xRefCutToAOMD::*ptfunc_createSubMeshFromRefCut)(AOMD::mEntity *, xRefMesh &, xfem::xMesh *,
                                                         std::vector<AOMD::mEntity *> &);
   void (xRefCutToAOMD::*ptfunc_createSubEntitiesSubMeshFromRefCut)(xRefMesh &, std::vector<AOMD::mEntity *> &);
   void (xRefCutToAOMD::*ptfunc_createSubEntitiesSubMeshFromRefCutByLoop)(AOMD::mEntity *, xRefMesh &, const xfem::xLevelSet &);
   void (xRefCutToAOMD::*ptfunc_createIsoZeroFromRefCut)(AOMD::mEntity *, xRefMesh &, xfem::xMesh *,
                                                         std::vector<AOMD::mEntity *> &, const std::vector<double> &);
   void (xRefCutToAOMD::*ptfunc_createIsoZeroFromFrontier)(AOMD::mEntity *, const bool, const std::vector<double> &,
                                                           xfem::xMesh *, std::vector<AOMD::mEntity *> &);
   void (xRefCutToAOMD::*ptfunc_cleanSubEntitiesClassification)(AOMD::mEntity *, std::vector<AOMD::mEntity *> &, const bool,
                                                                const bool);
   void (xRefCutToAOMD::*ptfunc_tagEntitiesOfElement)(AOMD::mEntity *, std::vector<AOMD::mEntity *> &);

   // iso-zero node pointer for one element
   // size = maximum number of nodes for all kind of element treated so far + maximum of cut edge possible
   // here :  4 nodes for tetraedon + 4 cut edge => 8
   std::array<AOMD::mVertex *, 8> iso_nodes;
   // submesh pointer for one edge element
   // size = maximum number of edge for all kind of element treated so far
   // here :  6 edges for tetraedon
   std::array<xfem::xMesh *, 6> edge_submesh;
   // submesh pointer for one edge element in case of recursivity
   // size = maximum number of edge for all kind of element treated so far
   // here :  6 edges for tetraedon
   std::array<xfem::xMesh *, 6> edge_submesh_recursive;
   // constitutive edge pointer for one element
   // size = maximum number of edge for all kind of element treated so far
   // here :  6 edges for tetraedon
   std::array<AOMD::mEdge *, 6> edges;
   // constitutive face pointer for one element
   // size = maximum number of face for all kind of element treated so far
   // here :  4 faces for tetraedon
   std::array<AOMD::mFace *, 4> faces;
   // submesh pointer for one face element
   // size = maximum number of face for all kind of element treated so far
   // here :  4 faces for tetraedon
   std::array<xfem::xMesh *, 4> face_submesh;
   // submesh pointer for one face element in case of recursivity
   // size = maximum number of face for all kind of element treated so far
   // here :  4 faces for tetraedon
   std::array<xfem::xMesh *, 4> face_submesh_recursive;
   // iso-zero reference entity  pointer for attachement of sub iso-zero mesh
   AOMD::mEntity *iso_entity_ref;

   // private methodes ///////////////
   void setFunctionForEntity(const AOMD::mEntity *);
   bool getSubmesh(AOMD::mEntity *, xfem::xMesh **);
   void getIsoNodesFromContext(const AOMD::mEntity *, const xRefMesh::elemdef_t &, xRefMesh &, const char &, const char &,
                               xfem::xMesh *, const std::vector<AOMD::mEntity *> &, xmapping::xMapping *, const char *);
   void createSubEdgeFromRef(xRefMesh &, std::vector<AOMD::mEntity *> &);
   inline void duplicateNode(AOMD::mVertex *, AOMD::mVertex **, xfem::xMesh *);
   void createSubEntitiesFromRefCutByLoop(AOMD::mEntity *, xRefMesh::elem_it it, xRefMesh::elem_it itend,
                                          xfem::xMesh **entity_submesh_recursive, const xfem::xLevelSet &ls, int d, char *,
                                          int (*)(const std::vector<double> &, const std::vector<int> &, xRefMesh &),
                                          void (xRefCutToAOMD::*)(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, AOMD::mVertex **));
   void createSubSubEntitiesFromMultiRefCut(xfem::xMesh *, const xfem::xLevelSet &, int,
                                            int (*)(const std::vector<double> &, const std::vector<int> &, xRefMesh &),
                                            void (xRefCutToAOMD::*)(AOMD::mEntity *, xRefMesh &, xfem::xMesh *,
                                                                    AOMD::mVertex **));
   void createSubFace(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, AOMD::mVertex **);
   void createSubEdge(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, AOMD::mVertex **);
   void createIsoZeroFromEltFrontier(AOMD::mEntity *, const bool, const std::vector<double> &, xfem::xMesh *,
                                     std::vector<AOMD::mEntity *> &, const char, const char, const char, const char, char[][2],
                                     char[][4], char[][4]);
   void tagEdgesOfElement(AOMD::mEntity *, std::vector<AOMD::mEntity *> &, int, char[], char[][2]);
   void tagFacesOfElement(AOMD::mEntity *, std::vector<AOMD::mEntity *> &, int, char[], char[][4]);
   inline void tagNodes(const std::vector<double> &, const std::vector<AOMD::mEntity *> &, xRefMesh &);
   inline void setTagForTaggingEdgesOfElement(std::vector<AOMD::mEntity *> &, AOMD::mEntity *, char *);

   void setPromotorsStatus(AOMD::mEntity *);
   void promoteStatusOfPromoters(AOMD::mEntity *);
   void cleanAllPartition(xfem::xMesh *);
   void cleanEltPartition(AOMD::mEntity *);
   inline void tagSubSubEntities(AOMD::mEntity *, int, char, xfem::xEntityToEntity);
   inline void tagSubSubEntitiesTouching(AOMD::mEntity *, bool, int, int, char);
   inline void cleanSubEntitiesClassificationNoCheck(AOMD::mEntity *, const bool, const int, const int);
   inline void cleanSubEntitiesClassificationCheck(AOMD::mEntity *, const bool, const int, const int, const int);
   // for Edge
   void createSubMeshEdgeFromRef(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, std::vector<AOMD::mEntity *> &);
   void createIsoZeroFromEdgeRefCut(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, std::vector<AOMD::mEntity *> &,
                                    const std::vector<double> &);
   void createSubEntitiesSubMeshForEdgeFromRef(xRefMesh &, std::vector<AOMD::mEntity *> &);
   void cleanSubEntitiesClassificationEdge(AOMD::mEntity *, std::vector<AOMD::mEntity *> &, const bool, const bool);
   void tagEntitiesOfElementEdge(AOMD::mEntity *, std::vector<AOMD::mEntity *> &);
   void createIsoZeroFromEdgeFrontier(AOMD::mEntity *, const bool, const std::vector<double> &, xfem::xMesh *,
                                      std::vector<AOMD::mEntity *> &);
   void createSubEntitiesSubMeshForEdgeFromRefCutByLoop(AOMD::mEntity *, xRefMesh &, const xfem::xLevelSet &);
   // for triangle
   void createSubMeshTriFromRef(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, std::vector<AOMD::mEntity *> &);
   void createSubEntitiesSubMeshForTriFromRef(xRefMesh &, std::vector<AOMD::mEntity *> &);
   void createSubEntitiesSubMeshForTriFromRefCutByLoop(AOMD::mEntity *, xRefMesh &, const xfem::xLevelSet &);
   void createIsoZeroFromTriRefCut(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, std::vector<AOMD::mEntity *> &,
                                   const std::vector<double> &);
   void createIsoZeroFromTriFrontier(AOMD::mEntity *, const bool is_out, const std::vector<double> &, xfem::xMesh *,
                                     std::vector<AOMD::mEntity *> &);
   void cleanSubEntitiesClassificationTri(AOMD::mEntity *, std::vector<AOMD::mEntity *> &, const bool, const bool);
   void tagEntitiesOfElementTri(AOMD::mEntity *, std::vector<AOMD::mEntity *> &);

   // for tetraedron
   void createSubMeshTetFromRef(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, std::vector<AOMD::mEntity *> &);
   void createSubEntitiesSubMeshForTetFromRef(xRefMesh &, std::vector<AOMD::mEntity *> &);
   void createSubEntitiesSubMeshForTetFromRefCutByLoop(AOMD::mEntity *, xRefMesh &, const xfem::xLevelSet &);
   void createIsoZeroFromTetRefCut(AOMD::mEntity *, xRefMesh &, xfem::xMesh *, std::vector<AOMD::mEntity *> &,
                                   const std::vector<double> &);
   void createIsoZeroFromTetFrontier(AOMD::mEntity *, const bool, const std::vector<double> &, xfem::xMesh *,
                                     std::vector<AOMD::mEntity *> &);
   void cleanSubEntitiesClassificationTet(AOMD::mEntity *, std::vector<AOMD::mEntity *> &, const bool, const bool);
   void tagEntitiesOfElementTet(AOMD::mEntity *, std::vector<AOMD::mEntity *> &);
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xRefCutToAOMD class
///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xRefCutToAOMDException class
//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// interface derived class of standart exception for xsurf
class xRefCutToAOMDException : public std::exception
{
  public:
   xRefCutToAOMDException(std::string, std::string, int, std::string, std::string);
   ~xRefCutToAOMDException() throw() override;
   const char *what() const throw() override;

  private:
   std::string msg;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xRefCutToAOMDException  class
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace xcut

#endif
