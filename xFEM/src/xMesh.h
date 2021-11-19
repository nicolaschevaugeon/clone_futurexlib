/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _XMESH_H_
#define _XMESH_H_

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

// mpi
#include "mpi.h"
// AOMD includes
#include "GEntity.h"
#include "mAOMD.h"
#include "mVertex.h"
// xfem
#include "xDataExchangerTools.h"
#include "xEntityFilter.h"
#include "xEntityToEntity.h"
#include "xIter.h"
#include "xOctreeGrid.h"
#include "xPartition.h"
#include "xRegularGrid.h"
#include "xUnorderedMapDataManager.h"
// xinterface aomd
#include "xAttachedDataManagerAOMD.h"
// xtensor
#include "xPoint.h"
// xgeom
#include "xBoundingBox.h"

namespace xfem
{
class xSubMeshCreator;
class xRegion;
class xLevelSet;
class xSubMesh;
typedef std::function<void(AOMD::mEntity*, xPartition&, xEntityFilter)> xGetPartition;

/// This function return the "source " of e :
/*!
  It check for the tag was_duplicated_from.
  If it exist it call it self on the    entity pointed by the tag.
  else it return the input entity.
  In short, the function return e, if e is not the copy of an other entity,
  other wise, it goes up the hierachies of copies to return the "source" entity.
*/
AOMD::mEntity* getSource(AOMD::mEntity* e);

xgeom::xBoundingBox compute_bounding_box(const AOMD::mEntity& e);

/// the xMesh adds a functionality to the mMesh class of the Aomd package.
/*!
  It adds the functionnality to create subset of a mesh.
  a subset is a group of AOMD::mEntity *, for which a name is associated.
  xMesh offer the possibility to iterate on a subset.
  !*/
class xMesh  //: public AOMD::mMesh
{
  public:
   AOMD::mMesh& getMesh();
   const AOMD::mMesh& getMesh() const;
   xIter beginVertex() const;
   xIter endVertex() const;
   xIter beginEdge() const;
   xIter endEdge() const;
   xIter beginFace() const;
   xIter endFace() const;
   xIter beginSolid() const;
   xIter endSolid() const;
   xIter begin(int what) const;
   xIter end(int what) const;
   template <typename T>
   using datamanager_t = xinterface::aomd::xAttachedDataManagerAOMD<T>;

   using partman_t = partmanAOMD_t;

   xMesh(MPI_Comm world = MPI_COMM_WORLD, int id = 1);
   xMesh(const string& filename, MPI_Comm world = MPI_COMM_WORLD);
   xMesh(const string& filename, xinterface::aomd::xAttachedDataManagerAOMD<int>& entities_id, MPI_Comm world = MPI_COMM_WORLD);
   xMesh(xMesh&&) = default;
   xMesh& operator=(xMesh&&) = default;
   xMesh& operator=(const xMesh&) = delete;
   /// return an iterator range on entity of level what.
   xtool::xRange<xIter> range(int what) const { return xtool::make_range(begin(what), end(what)); }
   //[[deprecated]] void modifyAllState();
   //[[deprecated]] void modifyAllStateFalse();
   /// remove the specific data of an xMesh relative to e. To fully Delete an entity from the xMesh, see the comment below :
   void clear(AOMD::mEntity* e);
   /// change the del function of the base class (call clear(e) the mMesh::del(e))
   /// Note : the del function was initially thought as an overload to the AOMD::mMesh del function.
   /// del in AOMD do not destroy the elment but cremove it from the mMesh. In AOMD, to do the previous action and delete the
   /// entity, one was expected to call DEL(e), which itsef called del(e). But mMesh is not any more a base of xMesh. The intent I
   /// think, was to forbid to modify as much as possible the mesh contained by xMesh. Therefore, the del function as it was need
   /// to be removed. To delete an entity from xMesh, one should first clear the specifique data associated to the entity
   /// controlled by the xMesh class, then get the mMesh via getMesh and calle DEL on the mMesh.
   /// void del( AOMD::mEntity* e );

   /// clean up the xMesh by calling clear on each entity, then clean on the mMesh after that the mesh is empty.
   void clear();
   // xMesh destructor
   ~xMesh();

   void copyMesh(const xMesh& other);
   void copyMesh(const xMesh& other, xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& associated_new_entity);

   /// return 2 points forming the bounding box of the mesh
   xgeom::xBoundingBox compute_bounding_box() const;
   //[[deprecated]]
   void compute_bounding_box(xtensor::xPoint& min, xtensor::xPoint& max) const;

   /// return the dimension (0, 1, 2, 3 ) of the highest level entity in the mesh (local to the proc)
   int dim() const;

   /// return the dimension (0, 1, 2, 3 ) of the highest level entity in the mesh (global  : max for the mesh on all proc)
   int dim_global() const;

   /// Mirror vertex are used for periodic boundary condition.
   // bool isMirror(AOMD::mEntity*);
   // AOMD::mMirrorVertex* getMirrorVertex(AOMD::mEntity*);

   // new functions to deal with subsets
   xSubMesh& createSubMesh(const string& subMeshName) const;
   xSubMesh& createSubMesh(const string& subMeshName, xSubMeshCreator& creator) const;
   void deleteSubMesh(const string& subset) const;
   void deleteAllSubMesh() const;
   xSubMesh& getSubMesh(const string& sub) const;

   int size(int) const;
   AOMD::mEntity* find(AOMD::mEntity* e) const;
   // return the support of the entity e
   // the level of all the elements in l
   // is the same as dim(), all the elements in l are leaves
   // lookupSupport add things to l
   // So make sure it is empty before calling it
   // lookupSupport takes into account the fact that the entity can be a mirror
   void lookupSupport(AOMD::mEntity* e, std::set<AOMD::mEntity*>& l);

   /// check mesh : calculate bound and average "volume" value of elements.
   //! if there is a problem mim value are fare less then average.
   //! fact is a thresold to print entity with a "volume" less then 1/fact times the average "volume"
   std::pair<double, double> checkMesh(const int fact = 100000000) const;

   // Finding the hanging nodes
   // map vertex to edge or face
   typedef std::map<AOMD::mEntity*, AOMD::mEntity*> hanging_t;
   void findHangingNodes();
   // void periodicAssociation();
   /// fill up partition with the entities of same level as e found in the mesh associated to e in get_const_partition(),
   /// recursivelly. If no mesh or mesh empty fill with just e.
   static void getPartition(AOMD::mEntity* e, xPartition& partition, xEntityFilter filter = xAcceptAll());
   /// get the partition attached to entity e
   static void getPartitionEntity(AOMD::mEntity* e, std::set<AOMD::mEntity*>& partition);
   // Where is the implementation ??? not in Xfem thow
   static void createSubSimplices(xRegion& support, xMesh* x_interface,
                                  std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& ls);

   /// Partition manager
   partman_t& getPartitionManager();
   const partman_t& getPartitionManager() const;

   void cleanPartition();
   void cleanClassification();
   /// create octree data structure to locate elements (xOctreeGrid).
   void createOctree() const;
   /// clear the octree data structure to locate elements (xOctreeGrid).
   void clearOctree() const;
   /// return the octree data structure to locate elements (xOctreeGrid).
   const xgeom::xOctreeGrid& getOctreeGrid() const;
   /// create regular grid data structure to locate elements (xRegularGrid).
   void createGrid() const;
   /// clear the regular grid data structure to locate elements (xRegularGrid).
   void clearGrid() const;
   /// return the regular grid data structure to locate elements (xRegularGrid).
   const xgeom::xRegularGrid& getRegularGrid() const;
   std::vector<std::pair<const AOMD::mEntity*, const xtensor::xPoint>> locateElementOctree(const xtensor::xPoint& p) const;
   std::vector<std::pair<const AOMD::mEntity*, const xtensor::xPoint>> locateElement(const xtensor::xPoint& p) const;

   static datamanager_t<AOMD::mEntity*>& get_was_created_by();
   static const datamanager_t<AOMD::mEntity*>& get_const_was_created_by();
   static datamanager_t<double>& get_r_on_edge();
   static const datamanager_t<double>& get_const_r_on_edge();
   static datamanager_t<AOMD::mEntity*>& get_is_duplicated_in();
   static const datamanager_t<AOMD::mEntity*>& get_const_is_duplicated_in();
   static datamanager_t<AOMD::mEntity*>& get_was_duplicated_from();
   static const datamanager_t<AOMD::mEntity*>& get_const_was_duplicated_from();
   /// return  entity to partition mesh datamanager
   static datamanager_t<xMesh>& get_partition();
   /// return a const  entity to partition mesh datamanager
   static const datamanager_t<xMesh>& get_const_partition();
   static datamanager_t<AOMD::mEntity*>& get_is_in_partition_of();
   static const datamanager_t<AOMD::mEntity*>& get_const_is_in_partition_of();
   // for octree meshes created xinterface/octree.
   /// return return the entity to hanging entity datamanager
   static datamanager_t<AOMD::mEntity*>& get_is_hanging_on();
   static const datamanager_t<AOMD::mEntity*>& get_const_is_hanging_on();

   static datamanager_t<AOMD::mEntity*>& get_is_hanging_by();
   static const datamanager_t<AOMD::mEntity*>& get_const_is_hanging_by();

   static datamanager_t<std::vector<AOMD::mEntity*>>& get_down_group();
   static const datamanager_t<std::vector<AOMD::mEntity*>>& get_const_down_group();

   static datamanager_t<std::vector<AOMD::mEntity*>>& get_bnd_group();
   static const datamanager_t<std::vector<AOMD::mEntity*>>& get_const_bnd_group();

   static datamanager_t<int>& get_octree_level();
   static const datamanager_t<int>& get_const_octree_level();

   int getOctreeLevelMax() const { return level_max; }
   void setOctreeLevelMax(int l) { level_max = l; }
   int level_max;

   /// Refine the mesh n times (only for 2D triangles)
   xMesh* refineInNewMesh(int n, bool debut = true);
   /// Return if an element is on the boundary
   bool isOnBoundary(AOMD::mEntity* e);
   /// Return a subset with the neighbours vertices
   set<AOMD::mVertex*> getNeighbors(AOMD::mVertex* v);
   /// Return the boundary
   xSubMesh* setBoundary();
   xSubMesh* getBoundary();

  private:
   partman_t part_man;
   mutable std::unique_ptr<xgeom::xRegularGrid> regular_grid = nullptr;
   mutable std::unique_ptr<xgeom::xOctreeGrid> octree_grid = nullptr;
   AOMD::mMesh mesh;
   // access functions to the subsets
   // Subset are implemented as a map between the string that named the subset and an
   // Entity Container. Note :  Could be more efficient too implement it as a tag
   // associated to the string attached to each Entity that is in the sub
   mutable map<string, xSubMesh*> subsetEntities;
   typedef map<string, xSubMesh*>::const_iterator const_iter_subsets;
   typedef map<string, xSubMesh*>::iterator iter_subsets;

   void lookupSupportBasic(AOMD::mEntity* e, std::set<AOMD::mEntity*>& l);

   void copyMeshInternal(const xMesh& other, xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& associated_new_entity);

   /// cutmesh utilities
   void createSubSimplices(const xLevelSet& lst, xMesh* m_interface, xEntityToEntity& class_in, xEntityToEntity& class_out);
   //  AOMD::mVertex* vertexInBetween(AOMD::mEntity* v0, AOMD::mEntity* v1, AOMD::mEntity* q);
   /*! note internal function fo xMeshCut.
      prerequist : vertices in q have entity associated to them via was_created_by and is_duplicated_in
  */
   AOMD::mVertex* vertexInBetween(AOMD::mEntity* v0, AOMD::mEntity* v1);
   void addNewTo(AOMD::mEntity* e, const std::string& side, xEntityToEntity& class_in, xEntityToEntity& class_out);
   static void setDataManagers();

   xSubMesh* boundary;
};
////////////////////////////////// end class xMesh //////////////////////////////////////////////////

}  // namespace xfem

#endif
