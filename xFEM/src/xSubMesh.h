/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xSubMesh_
#define _xSubMesh_

#include "mAOMD.h"
#include "xEntityFilter.h"
#include "xIteratorTools.h"
#include "xMesh.h"


namespace xfem
{
/// xSubMesh Concept:
/*! A xSubMesh is a collection of entity belonging to a xMesh.
    The dimension of a xSubMesh is the dimension of the AOMD::mEntity of highest dimension in the xSubMesh.
!*/

class xSubMesh
{
  public:
   /// xSubMesh Constructor. Construct a submesh, with name _submeshname off mesh mesh.
   xSubMesh(const std::string submeshname, const xMesh& mesh);
   /// xSubMesh destructor
   ~xSubMesh();
   /// return the dimension (0, 1, 2, 3) of the part of the submesh on the current processor.
   int dim() const;
   /// return the dimension of the submesh, on the current processor.
   /*! Warning : this function involve communication accross all processors. !*/
   int dim_all_process() const;

   int size(int what) const;

   /// return the number of entities of dim what off the subMesh
   /*! Warning : this function involve communication accross all processors. !*/
   //  int  size_all_process(int what) const;

   /// return the begin iterator on entities of dimension what on the submesh.
   xIter begin() const { return begin(dim()); }
   /// return the end iterator on entities of dimension what on the submesh.
   xIter end() const { return end(dim()); }
   xtool::xRange<xIter> range() const { return range(dim()); }

   /// return the begin iterator on entities of dimension what on the submesh.
   xIter begin(int what) const;
   /// return the end iterator on entities of dimension what on the submesh.
   xIter end(int what) const;
   xtool::xRange<xIter> range(int what) const;

   xIterall beginall(int what) const;  //
   xIterall endall(int what) const;    //
   /// return the begin iterator on entities of the bnd
   xIterall beginBnd();
   /// return the end iterator on entities of the bnd
   xIterall endBnd();
   /// return a pointer to the entity given in entry if the entity is in the subset, else return 0
   AOMD::mEntity* find(AOMD::mEntity*) const;
   /// return true if the entity is on the boundary of the subset.
   bool isOnBoundary(AOMD::mEntity*) const;
   void modifyState(int from, int to, bool state, int with = 0);  //
   void modifyAllState();
   /// export the subset in gmsh format
   void exportGmsh(const std::string& filename, int formatoption = 0, int merging = 0) const;

   void add(AOMD::mEntity*);  //
   void del(AOMD::mEntity*);  //
   /// make the submesh cooherent accross partion boundary. This is collective across all process.
   void updatePartitionBoundary();
   /// construct the boundary of the domain. For the result to be correct, it must call updatePartitionBoundary, if parallel is
   /// not uptodate. Therefore it can be collective across all processors.
   void updateBoundary();
   bool stateUpToDate() const;
   bool boundaryUpToDate() const;
   bool parallelUpToDate() const;
//    const xPartitionBoundaryInfo& getPartitionBoundaryInfo() const
//    {
//       (const_cast<xSubMesh&>(*this)).updatePartitionBoundary();
//       return ownermanager;
//    };
   std::string getName() const;
   const xMesh& getMesh() const;
   /// Copy the elements of the submesh in a new mesh
   xMesh* copyInNewMesh();
   /// Compute bounding box:
   void compute_bounding_box(xtensor::xPoint& min, xtensor::xPoint& max) const;
   /// get partition manager
   //! For now it is not clear what should be given
   //! A priori it may a subset of the mesh partition manager on which it rely ... or not
   //! In sequential it make no difference as the partition manager is empty ...
   //! A assert is checking that it is not used on more then one proc
   //! TODO: port xSubMesh properly and fix this choice
   const xMesh::partman_t& getPartitionManager() const;

  private:
   std::string submeshname;
   const xMesh& mesh;
   xinterface::aomd::xAttachedDataManagerAOMD<int> tagged;
   mutable bool _stateuptodate;
   mutable bool _boundaryuptodate;
   mutable bool _paralleluptodate;

   // int dim;
   AOMD::mMeshEntityContainer container;
   AOMD::mMeshEntityContainer::CONTAINER boundary;
//    mutable xPartitionBoundaryInfo ownermanager;

   struct xCreateDownwardFunctor
   {
     public:
      xCreateDownwardFunctor(int i, int j, xSubMesh& submesh, bool f = true, bool up = false);  //
      void operator()(AOMD::mEntity* e);                                                        //
     private:
      int i_dim, j_dim;
      xSubMesh& submesh;
      bool force_create;
      bool create_upward_too;
   };
};

class xSubMeshCreator
{
  public:
   xSubMeshCreator() = default;
   virtual ~xSubMeshCreator() = default;
   virtual void create(const xMesh&, const std::string& sub) = 0;
};

/// create a new subset of the mesh containing all entities that pass the filter
class xSubMeshCreatorFilter : public xSubMeshCreator
{
  public:
   xSubMeshCreatorFilter(const xEntityFilter& filter);
   ~xSubMeshCreatorFilter() override = default;
   void create(const xMesh&, const std::string& sub) override;

  private:
   xEntityFilter filter;
};

/// create a new subset which is the union of all subset passed at construction.
class xUnionCreator : public xSubMeshCreator
{
  public:
   xUnionCreator(const std::vector<std::string>& a) : xSubMeshCreator(), all(a) {}
   xUnionCreator(const std::string& s1, const std::string& s2)
   {
      all.push_back(s1);
      all.push_back(s2);
   }
   void create(const xMesh&, const std::string& name) override;

  private:
   std::vector<std::string> all;
};

// Intersection
class xIntersectionCreator : public xSubMeshCreator
{
  public:
   xIntersectionCreator(const std::vector<std::string>& a) : xSubMeshCreator(), all(a) {}
   void create(const xMesh&, const std::string& name) override;

  private:
   const std::vector<std::string>& all;
};

/// Create a xSubMesh made off all entities from first minus all entities from second
class xFirstMinusSecondCreator : public xSubMeshCreator
{
  public:
   xFirstMinusSecondCreator(const std::string& f, const std::string& s) : xSubMeshCreator(), first(f), second(s) {}
   void create(const xMesh&, const std::string& name) override;

  private:
   std::string first;
   std::string second;
};

/// Create a xSubMesh made off all entity in the mesh.
class xAllCreator : public xSubMeshCreator
{
  public:
   xAllCreator();
   void create(const xMesh&, const std::string& name) override;
};

/// Create a xSubMesh made off all entity in  a xSubMesh.
class xCopyCreator : public xSubMeshCreator
{
  public:
   xCopyCreator(const std::string& source);
   void create(const xMesh&, const std::string& name) override;

  private:
   const std::string source;
};

class xAddLayerCreator : public xSubMeshCreator
{
  public:
   xAddLayerCreator(int nblayers, const std::string& i, int d = 0, xEntityFilter _filter = xAcceptAll())
       : xSubMeshCreator(), layers(nblayers), initial(i), dim_growth(d), filter(_filter)
   {
   }
   void create(const xMesh& m, const std::string& name) override;

  private:
   const int layers;
   const std::string initial;
   const int dim_growth;
   xEntityFilter filter;
   // meaning of dim_growth
   // if it is zero, we add all the elements connected to the nodes of the current subset
   // it it is one, we add all the elements connected to the edges of the current subset
   // it it is two, we add all the elements connected to the faces of the current subset
};

class xSubMeshModifier
{
  public:
   virtual ~xSubMeshModifier() = default;
   ;
   virtual void modify(const xMesh&, const std::string& name) const = 0;
};

class xAddLayerModifier : public xSubMeshModifier
{
  public:
   xAddLayerModifier(int nb_layers, int d = 0, xEntityFilter _filter = xAcceptAll())
       : layers(nb_layers), dim_growth(d), filter(_filter)
   {
   }
   typedef AOMD::mMeshEntityContainer::CONTAINER Container;
   // implementation of base class xSubSet member
   void modify(const xMesh&, const std::string& name) const override;

  private:
   const int layers;
   const int dim_growth;
   xEntityFilter filter;
   // mutable Container received;
};

/*!
  Unify the submesh accross process :
    By calling this modifiyer, one insure that  if an entity exist in a  submesh at least on one proc and is on the partition of
the background mesh, this entity will be in the submesh of all the proc that share this entity.
!*/
class xUnifySubMeshAccrossProcess : public xSubMeshModifier
{
  public:
   void modify(const xMesh&, const std::string& name) const override;
};

}  // namespace xfem
#endif
