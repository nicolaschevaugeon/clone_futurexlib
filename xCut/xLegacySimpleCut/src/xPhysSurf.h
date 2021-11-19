/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef ___XPHYSURF__H
#define ___XPHYSURF__H

#include "AOMDfwd.h"
#include "xDebug.h"
#include "xEntityToEntity.h"
#include "xLevelSet.h"
#include "xMesh.h"
#include "xPointToDouble.h"
// parallel
//#include "xParallel.h"
#include "xSubMeshManager.h"

/// 17-06-09
#include "mAOMD.h"
#include "mAttachableDataContainer.h"
#include "mEdge.h"
#include "mEntity.h"
#include "mFace.h"
#include "mRegion.h"
#include "mTet.h"
#include "mVertex.h"
#include "xElement.h"

namespace xfem
{
class xGeomElem;
}
// A domain is defined by the negative part of a level set over a mesh
// The domain owns the level set representing it
// To create a domian, a mesh is given as well as a functor producing the
// level set value at each node.
namespace xcut
{
//! domaine

class xPhysSurfBase /*: public xfem::xExchangeBool*/
{
   // definition du support
   // ensemble des ements connect au noeud (si un parent est connect il est
   // remplac par son enfant)
   // UNION le support des hanging des ements ci-dessus
   // IN   : support covering  the domain
   // OUT  : support completely outside the domain
   // dans le maillage, les sous-domaines correspondant
   // sont not domain_in, domain_out, domain_inout

  public:
   xPhysSurfBase(const xfem::xEntityToEntity &in, const xfem::xEntityToEntity &out, bool fit_ = true, double fittol_ = 1.e-2)
       : mesh(nullptr), mesh_bnd(nullptr), classify_in(in), classify_out(out), fit(fit_), fittol(fittol_), region()
   {
   }

//   ~xPhysSurfBase() override = default;
   xfem::xMesh &getMesh() { return *mesh; }
   const xfem::xMesh &getMesh() const { return *mesh; }

   xfem::xMesh *getMesh_bnd() { return mesh_bnd; }
   const xfem::xMesh *getMesh_bnd() const { return mesh_bnd; }
   xfem::xMesh *getMeshPtr() { return mesh; }
   const xfem::xMesh *getMeshPtr() const { return mesh; }
   bool covers(AOMD::mEntity *e) const;  //! indique si le support de e
   //! est au moins en un point strictement
   //! negatif

   bool boundary(AOMD::mEntity *e) const;  //! indique si le support de e
   //! passe
   //! d'une valeur strictement negative
   //! a une valeur strictement positive

   bool boundary_strict(AOMD::mEntity *e) const;  //! indique si au moins un des elements
   //! du support de e passe
   //! d'une valeur strictement negative
   //! a une valeur strictement positive

   void createSupportInfo();
   // virtual void  update(bool fit=true, bool keep_old_partition_flag = false) = 0;
   xfem::xIter begin();  // beginning of the elements covered by the
   xfem::xIter end();    // domain
   int size();
   int dim() const;

   xfem::xEntityToEntity &getClassifyerIn() { return classify_in; }
   xfem::xEntityToEntity &getClassifyerOut() { return classify_out; }

   xfem::xSubMeshManager *getSubMeshManager() { return &submesh_manager; }

  protected:
   enum inout
   {
      in,
      out,
      on
   };
   virtual inout side(const AOMD::mEntity &e) const = 0;

   xfem::xMesh *mesh;
   xfem::xMesh *mesh_bnd;
   xfem::xEntityToEntity classify_in;
   xfem::xEntityToEntity classify_out;
   bool fit = true;
   double fittol;
   xfem::xRegion region;  // utile dans createSupportInfo()

   xfem::xSubMeshManager submesh_manager;
   typedef std::map<AOMD::mEntity*,bool>   bucket_type;
   bucket_type map_in, map_out, map_bdry;
   void pre_exchange() /*override*/;   // initialisation du sac des noeuds
   void post_exchange() /*override*/;  // mise a jour des sacs
};

class xPhysSurf : public xPhysSurfBase
{
  public:
   template <class T>
   using datamanager_t = xfem::xMesh::datamanager_t<T>;
   xPhysSurf(xfem::xLevelSet &ls1, const xfem::xEntityToEntity &in, const xfem::xEntityToEntity &out, bool fit = true,
             double fittol_ = 1.e-2, bool keep_old_partition = false, bool recursive = false,
             datamanager_t<AOMD::mEntity *> &was_created_by_ = xfem::xMesh::get_was_created_by(),
             datamanager_t<double> &r_on_edge_ = xfem::xMesh::get_r_on_edge(),
             datamanager_t<AOMD::mEntity *> &is_duplicated_in_ = xfem::xMesh::get_is_duplicated_in(),
             datamanager_t<AOMD::mEntity *> &was_duplicated_from_ = xfem::xMesh::get_was_duplicated_from(),
             datamanager_t<xfem::xMesh> &partition_ = xfem::xMesh::get_partition(),
             datamanager_t<AOMD::mEntity *> &is_in_partition_of_ = xfem::xMesh::get_is_in_partition_of())
       : xPhysSurfBase(in, out, fit, fittol_),
         ls(ls1),
         was_created_by(was_created_by_),
         r_on_edge(r_on_edge_),
         is_duplicated_in(is_duplicated_in_),
         was_duplicated_from(was_duplicated_from_),
         partition(partition_),
         is_in_partition_of(is_in_partition_of_)
   {
      mesh = ls1.getSupport().getMesh();
      region = ls.getSupport();
      construct_(fit, keep_old_partition, recursive);
   }

   ~xPhysSurf() /*override*/;
   virtual void update(bool fit = true, bool keep_old_partition_flag = false);
   xfem::xLevelSet &getLevelSet();
   const xfem::xLevelSet &getLevelSet() const;

  private:
   void construct_(bool fit = true, bool keep_old_partition = false, bool recursive = false);

   int side_of(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ) const;

  protected:
   inout side(const AOMD::mEntity &e) const override;
   xfem::xLevelSet &ls;
   datamanager_t<AOMD::mEntity *> &was_created_by;       //      = xfem::xMesh::get_was_created_by();
   datamanager_t<double> &r_on_edge;                     //= xfem::xMesh::get_r_on_edge();
   datamanager_t<AOMD::mEntity *> &is_duplicated_in;     //   = xfem::xMesh::get_is_duplicated_in();
   datamanager_t<AOMD::mEntity *> &was_duplicated_from;  //= xfem::xMesh::get_was_duplicated_from();
   datamanager_t<xfem::xMesh> &partition;                //          = xfem::xMesh::get_partition();
   datamanager_t<AOMD::mEntity *> &is_in_partition_of;   // = xfem::xMesh::get_is_in_partition_of();
};

class xPhysSurfOctree : public xPhysSurfBase
{
  public:
   template <class T>
   using datamanager_t = xfem::xMesh::datamanager_t<T>;
   /// Constructor.
   /*!
     matint: true for material interface . In This case Refinment is not recursiv.
     n : level of refinment.
     !*/
   xPhysSurfOctree(xfem::xMesh &_mesh, const xfem::xPointToDouble &_lsv, const xfem::xEntityToEntity &in,
                   const xfem::xEntityToEntity &out, bool fit, double fittol_, bool matint, int n_,
                   xfem::xEntityFilter filter = xfem::xAcceptAll(), xfem::xEntityFilter filter_split = xfem::xAcceptAll(),
                   datamanager_t<AOMD::mEntity *> &was_created_by_ = xfem::xMesh::get_was_created_by(),
                   datamanager_t<xfem::xMesh> &partition_ = xfem::xMesh::get_partition());

//   ~xPhysSurfOctree() override = default;
   const datamanager_t<AOMD::mEntity *> &getWasCreatedBy() const { return was_created_by; }
   const datamanager_t<double> &getLevelSetValue() const { return levelset_value; }
   const datamanager_t<int> &getCutEdge() const { return cut_edge; }
   const datamanager_t<xfem::xMesh> &getRefinedElements() const { return refined_elements; }
   const datamanager_t<AOMD::mEntity *> &getRootEntity() const { return root_entity; }

  private:
   void construct_(bool fit = true);
   xfem::xPointToDouble lsv;
   int refineE;
   inout side(const AOMD::mEntity &e) const override;
   void setLevelSetValue(const AOMD::mVertex &v, double);
   std::vector<double> getLevelSetValues(const AOMD::mEntity &e) const;
   double getLevelSetValue(const AOMD::mVertex &v) const;
   void refineElementSelectionMatInt(const AOMD::mEntity &e, const xfem::xPointToDouble &lsv, int n, int &refineE);
   void refineElementSelection(const AOMD::mEntity &e, int n, int &refineE);
   void doIrefineEdge(const AOMD::mEdge &q, int n, int &refineE) const;
   void nRefineElement(const AOMD::mEntity &eorigine, const xfem::xPointToDouble &lsv, bool matint, int n, double fittol);
   void nFitToVertices(AOMD::mEntity *e, double fittol);
   void classifyElementOctreeNew(AOMD::mEntity *egros, AOMD::mEntity *epetit);
   // Warning only work for 2d triangle mesh for now
   void createElementPartitionAndInterfaceFromOneLevelSet(AOMD::mEntity *e);

   void createEdgePartitionFromOneLevelSet(AOMD::mEdge *e);

   datamanager_t<AOMD::mEntity *> &was_created_by;     //      = xfem::xMesh::get_was_created_by();
   datamanager_t<xfem::xMesh> &partition;              //          = xfem::xMesh::get_partition();
   datamanager_t<AOMD::mEntity *> is_in_partition_of;  // = xfem::xMesh::get_is_in_partition_of();

   datamanager_t<xfem::xMesh> refined_elements;
   datamanager_t<double> levelset_value;
   datamanager_t<int> cut_edge;
   datamanager_t<int> cut_element;

   datamanager_t<AOMD::mEntity *> is_the_creator_of2;
   datamanager_t<AOMD::mEntity *> was_created_by2;
   datamanager_t<AOMD::mEntity *> root_entity;
   datamanager_t<AOMD::mEntity *> is_the_copy_of;
};

void TagSidePositive(xfem::xMesh *mesh, const char *side_tag_name = "side_tag", const bool PositiveSide = true);

void TagSideNegative(xfem::xMesh *mesh, const char *side_tag_name = "side_tag", const bool PositiveSide = false);

void UnTagSide(xfem::xMesh *mesh, const char *side_tag_name = "side_tag");

}  // namespace xcut

#endif
