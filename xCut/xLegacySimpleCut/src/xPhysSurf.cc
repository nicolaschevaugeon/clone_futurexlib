/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include <algorithm>
#include <iostream>
#include <sstream>

// AOMD
#include "mAOMD.h"
#include "mEdge.h"
#include "mFace.h"
#include "mRegion.h"
#include "mTet.h"
#include "mVertex.h"
// xfem
#include "xDebug.h"
#include "xElement.h"
#include "xEntityFilter.h"
#include "xGeomElem.h"
#include "xLevelSet.h"
#include "xLevelSetOperators.h"
#include "xMesh.h"
#include "xPhysSurf.h"
#include "xPointToDouble.h"
#include "xSimpleGeometry.h"
#include "xVector.h"
// xcut
#include "xMeshCut.h"

#ifdef PARALLEL
#include "ParUtil.h"
using AOMD::ParUtil;
#endif

namespace xcut
{
using namespace std;
using AOMD::mEdge;
using AOMD::mEntity;
using AOMD::mFace;
using AOMD::mTet;
using AOMD::mVertex;

int xPhysSurfBase::dim() const { return mesh->dim(); }
xPhysSurf::~xPhysSurf() { delete mesh_bnd; }
void xPhysSurf::update(bool fit, bool keep_old_partition_flag) { construct_(fit, keep_old_partition_flag); }

void xPhysSurf::construct_(bool fit, bool keep_old_partition_flag, bool recursive)
{
   if (fit)
   {
      xfem::xFitToVertices fit(fittol);
      ls.accept(fit);
   }

   if (mesh_bnd) delete mesh_bnd;
   mesh_bnd = new xfem::xMesh;
   const bool createpartition = true;

   xcut::cutMesh(*mesh, ls, *mesh_bnd, was_created_by, r_on_edge, is_duplicated_in, was_duplicated_from, partition,
                 is_in_partition_of, classify_in, classify_out, createpartition, keep_old_partition_flag, recursive);

   // mesh->cutMesh(ls, mesh_bnd, classify_in, classify_out, createpartition, keep_old_partition_flag, recursive);  //
   // keep_old_partition_flag  indicates that we keep the old partition.

   createSupportInfo();

   for (AOMD::mEntity *pe : submesh_manager.range_sub(mesh->dim(), "in_domain")) classify_in(pe);
   for (AOMD::mEntity *pe : submesh_manager.range_sub(mesh->dim(), "out_domain")) classify_out(pe);
   return;
}

int xPhysSurf::side_of(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ) const
{
   return ls.side_of(geo_appro, geo_integ);
}

xfem::xLevelSet &xPhysSurf::getLevelSet() { return ls; }
const xfem::xLevelSet &xPhysSurf::getLevelSet() const { return ls; }

void xPhysSurfBase::createSupportInfo()
{
   // mat - modif for para
   const bool debug_csi = false;
   if (debug_csi) cout << "xphys ds csi - 1" << endl;
   pre_exchange();
   if (debug_csi) cout << "xphys ds csi - 2" << endl;
//   exchange(mesh);
   if (debug_csi) cout << "xphys ds csi - 3" << endl;
   post_exchange();
   if (debug_csi) cout << "xphys ds csi - 4" << endl;
   return;
}

xPhysSurfBase::inout xPhysSurf::side(const AOMD::mEntity &e) const
{
   std::vector<double> vals = ls.getVals(e);
   const double min = *min_element(vals.begin(), vals.end());
   const double max = *max_element(vals.begin(), vals.end());
   if (min >= 0.0)
      return out;
   else if (max <= 0.0)
      return in;
   else
      return on;
}

xPhysSurfBase::inout xPhysSurfOctree::side(const AOMD::mEntity &e) const
{
   if (const int *prefined = cut_element.getData(e))
      if (*prefined == 1)
      {
         return on;
      }

   std::vector<double> vals = getLevelSetValues(e);
   const double min = *min_element(vals.begin(), vals.end());
   const double max = *max_element(vals.begin(), vals.end());
   if (min >= 0.0)
      return out;
   else if (max <= 0.0)
      return in;
   else
      return on;
}

void xPhysSurfBase::pre_exchange()
{
   // cout << "xPhysSurf::pre_exchange() --------------------" << endl; ///
   const bool debug_pre = false;
   const bool debug_id_pre = false;
   // const bool debug = xfem::xdebug_flag;
   string out = "out_domain";  //! indicates if the support of e
   //! has not a single strictly negative value

   string in = "in_domain";  //! indicates if the support of e
   //! has at least one strictly negative value

   string boundary = "boundary";  //! indique si le support de e
   //! passe
   //! d'une valeur strictement negative
   //! a une valeur strictement positive

   string boundary_strict = "strict_boundary";  //! indique si au moins un des ��ents
   //! du support de e passe
   //! d'une valeur strictement negative
   //! a une valeur strictement positive
   //
   // boundary_strict est un sous-ensemble de boundary
   // pour un elts les notions de boundary et boundary_strict
   //      sont identiques

   //////////////////////////////////////////////////////////////////////////////
   // initialisation des maps
   //
   // map_in[e] est true si au moins un des elements de son support appartient a in_domain
   // map_out[e] est true si au moins un des elements de son support appartient a out_domain
   // attention map_in[e] at map_out[e] peuvent etre tous les deux true alors que
   //   in_domain et out_domain sont exclusifs

   // first clean things (important for an evolving level set)
   region.getMesh()->cleanClassification();
   submesh_manager.removeAllSubsets();
   map_in.clear();
   map_out.clear();
   map_bdry.clear();

   const int d = mesh->dim();
   if (debug_pre) cout << "XPHYS ds pre_exchange - 1.1" << endl;
   for (int i = 0; i <= d; ++i)
   {
      for (mEntity *pe : region.range(i))
      {
         map_in[pe] = false;
         map_out[pe] = false;
         map_bdry[pe] = false;
      }
   }

   ////////////////////////////////////////////////////////////////////////////
   submesh_manager.allocateSubsetEntities(out);
   submesh_manager.allocateSubsetEntities(in);
   submesh_manager.allocateSubsetEntities(boundary);
   submesh_manager.allocateSubsetEntities(boundary_strict);

   for (mEntity *e : region.range(d))
   {
      switch (side(*e))
      {
         case inout::out:
            break;
         case inout::in:
         {
            submesh_manager.add_sub(e, in);
            break;
         }
         case inout::on:
         {
            submesh_manager.add_sub(e, in);
            submesh_manager.add_sub(e, boundary);
            submesh_manager.add_sub(e, boundary_strict);
            break;
         }
      }
   }

   // loop on the elements
   // loop on the faces of each elt
   // select the elt touching the face
   //
   // bug possible
   //  pourquoi on boucle jusque d et pas d-1 bug ??
   // v�ifier que les elts sont bien mis �out
   //
   //
   for (int i = 0; i <= d; ++i)
   {
      for (mEntity *pe : region.range(i))
      {
         std::set<mEntity *> support_e;
         mesh->lookupSupport(pe, support_e);
         if (debug_id_pre)
         {
            cout << "the support for the entity\n";
            pe->print();
            cout << " is " << endl;
            for (mEntity *pes : support_e)
            {
               pes->print();
               cout << endl;
            }
         }
         //
         // e est in si au moins un elt du support est in
         //
         //
         // e est out si il n'est pas in
         //
         // e est  boundary
         //  si un des elts du support est boundary
         //  ou si il existe des in et des out
         //
         // e est boundary_strict
         //  si un des elts du support est boundary
         //

         ///////////////////////////////////////////////////////////////////////////
         // remplissage des maps
         // voir these Cloirec
         std::set<mEntity *>::iterator itc = support_e.begin();
         for (mEntity *pes : support_e)
         {
            if (submesh_manager.find_sub(pes, in))
               map_in[pe] = true;
            else
               map_out[pe] = true;
            if (submesh_manager.find_sub(pes, boundary)) map_bdry[pe] = true;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////
   // remplissage du vecteur de pointeurs de map

//   add_bucket(&map_in);
//   add_bucket(&map_out);
//   add_bucket(&map_bdry);

   ///////////////////////////////////////////////////////////////////////////

   if (debug_id_pre)
   {
      cout << "before exchange" << endl;
      cout << "====================================================" << endl;
      submesh_manager.printSubsetEntities();
      cout << "====================================================" << endl;
      cout << "before exchange the maps are" << endl;
      cout << "----------------------------------------------------" << endl;
      cout << "map in: " << endl;
      cout << "----------------------------------------------------" << endl;
      for (const auto &pe_bool : map_in)
      {
         pe_bool.first->print();
         cout << pe_bool.second << endl;
      }
      cout << "----------------------------------------------------" << endl;
      cout << "map out: " << endl;
      cout << "----------------------------------------------------" << endl;
      for (const auto &pe_bool : map_out)
      {
         pe_bool.first->print();
         cout << pe_bool.second << endl;
      }
      cout << "----------------------------------------------------" << endl;
      cout << "map bdry: " << endl;
      cout << "----------------------------------------------------" << endl;
      for (const auto &pe_bool : map_bdry)
      {
         pe_bool.first->print();
         cout << pe_bool.second << endl;
      }
   }
}

void xPhysSurfBase::post_exchange()
{
   const bool debug_post_id = false;
   string out = "out_domain";
   string in = "in_domain";
   string boundary = "boundary";
   string boundary_strict = "strict_boundary";

   // attention l'ordre des *map stocker ds buckets
   // vector<bucket_type*>::iterator it = buckets.begin();
//   bucket_type &map_post_in = *buckets[0];
//   bucket_type &map_post_out = *buckets[1];
//   bucket_type &map_post_bdry = *buckets[2];

   //GREG: Pour faire fonctionner sans xParallel.h (menage pre-future)
   bucket_type &map_post_in = map_in;
   bucket_type &map_post_out = map_out;
   bucket_type &map_post_bdry = map_bdry;



   /////////////////////////////////////////////////////////////////////////////////////////////////
   if (debug_post_id)
   {
      cout << "====================================================" << endl;
      cout << "after exchange the maps are" << endl;
      cout << "----------------------------------------------------" << endl;
      cout << "map post_in: " << endl;
      cout << "----------------------------------------------------" << endl;
      for (const auto &pe_bool : map_post_in)
      {
         pe_bool.first->print();
         cout << pe_bool.second << endl;
      }
      cout << "----------------------------------------------------" << endl;
      cout << "map post_out: " << endl;
      cout << "----------------------------------------------------" << endl;
      for (const auto &pe_bool : map_post_out)
      {
         pe_bool.first->print();
         cout << pe_bool.second << endl;
      }
      cout << "----------------------------------------------------" << endl;
      cout << "map post_bdry: " << endl;
      cout << "----------------------------------------------------" << endl;
      for (const auto &pe_bool : map_post_bdry)
      {
         pe_bool.first->print();
         cout << pe_bool.second << endl;
      }
      cout << "====================================================" << endl;
   }

   //////////////////////////////////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////////////////////////////////
   // ajout des entity ds les mesh containers en fonction des maps

   const int d = mesh->dim();
   for (int i = 0; i <= d; ++i)
   {
      for (mEntity *e : region.range(i))
      {
         if (map_post_in[e] == true)
            submesh_manager.add_sub(e, in);
         else
            submesh_manager.add_sub(e, out);
         if ((map_post_bdry[e] == true) || ((map_post_in[e] == true) && (map_post_out[e] == true)))
         {
            submesh_manager.add_sub(e, boundary);
         }
         if (map_post_bdry[e] == true)
         {
            submesh_manager.add_sub(e, boundary_strict);
         }
      }
   }

   //////////////////////////////////////////////////////////////////////////////////////////////////

   if (debug_post_id)
   {
      cout << "after exchange " << endl;
      cout << "====================================================" << endl;
      submesh_manager.printSubsetEntities();
      cout << "====================================================" << endl;
   }
}

bool xPhysSurfBase::covers(mEntity *e) const { return (submesh_manager.find_sub(e, "in_domain")); }
bool xPhysSurfBase::boundary(mEntity *e) const { return (submesh_manager.find_sub(e, "boundary")); }
bool xPhysSurfBase::boundary_strict(mEntity *e) const
{
   const bool debug = xfem::xdebug_flag;
   if (debug)
   {
      cout << " boundary_strict called for " << endl;
      e->print();
      cout << "result is " << (submesh_manager.find_sub(e, "strict_boundary")) << endl;
   }

   return (submesh_manager.find_sub(e, "strict_boundary"));
}

xfem::xIter xPhysSurfBase::begin()
{
   // cout << "size of in_domain " << mesh->size_sub(mesh->dim(), "in_domain");
   return submesh_manager.begin_sub(mesh->dim(), "in_domain");
}
xfem::xIter xPhysSurfBase::end() { return submesh_manager.end_sub(mesh->dim(), "in_domain"); }
int xPhysSurfBase::size() { return submesh_manager.size_sub(mesh->dim(), "in_domain"); }

void TagSidePositive(xfem::xMesh *mesh, const char *side_tag_name, const bool PositiveSide)
{
   unsigned int side_tag;
   side_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name);
   const bool debug = false;
   if (debug) cout << "In TagsidePositive : DIM=" << mesh->dim() << endl;

   for (auto e : mesh->range(mesh->dim()))
   {
      if (debug)
      {
         cout << "Entity in tagside:";
         e->print();
         cout << endl;
         cout << "In Tagside : Type=" << e->getType() << endl;
      }

      if (PositiveSide)
         e->attachInt(side_tag, 1);
      else
         e->attachInt(side_tag, -1);
      if (debug) cout << "in tagside : tag=" << e->getAttachedInt(side_tag) << endl;
   }
   return;
}

void TagSideNegative(xfem::xMesh *mesh, const char *side_tag_name, const bool PositiveSide)
{
   TagSidePositive(mesh, side_tag_name, PositiveSide);
}

void UnTagSide(xfem::xMesh *mesh, const char *side_tag_name)
{
   unsigned int side_tag;
   side_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name);
   const bool debug = false;
   if (debug) cout << "In UnTagside : DIM=" << mesh->dim() << endl;

   for (mEntity *pe : mesh->range(mesh->dim()))
   {
      if (debug)
      {
         cout << "Entity in UnTagside:";
         pe->print();
         cout << endl;
         cout << "In UnTagside : Type=" << pe->getType() << endl;
      }
      pe->deleteData(side_tag);
      if (debug) cout << "in UnTagside : tag=" << pe->getAttachedInt(side_tag) << endl;
   }
   return;
}

xPhysSurfOctree::xPhysSurfOctree(xfem::xMesh &_mesh, const xfem::xPointToDouble &_lsv, const xfem::xEntityToEntity &in,
                                 const xfem::xEntityToEntity &out, bool fit, double fittol_, bool matint, int n_,
                                 xfem::xEntityFilter filter, xfem::xEntityFilter filter_split,
                                 xfem::xMesh::datamanager_t<AOMD::mEntity *> &was_created_by_,
                                 xfem::xMesh::datamanager_t<xfem::xMesh> &partition_)
    : xPhysSurfBase(in, out, fit, fittol_), lsv(_lsv), was_created_by(was_created_by_), partition(partition_)
{
   // const bool debug = xfem::xdebug_flag;

   mesh = &_mesh;
   region = xfem::xRegion(mesh);

   mesh_bnd = new xfem::xMesh();
   int dim = mesh->dim();
   for (AOMD::mEntity *ev : mesh->range(0))
   {
      mVertex &v = static_cast<mVertex &>(*ev);
      setLevelSetValue(v, lsv(v.point()));
   }
   for (AOMD::mEntity *eorigine : mesh->range(dim))
   {
      if (filter(eorigine))
      {
         int n = n_;
         if (!filter_split(eorigine)) n = 1;
         int refineE = 0;
         if (matint)
            refineElementSelectionMatInt(*eorigine, lsv, n, refineE);  // homogeneous refinement
         else
            refineElementSelection(*eorigine, n, refineE);  // adaptative refinement

         if (refineE == 1)
         {
            nRefineElement(*eorigine, lsv, matint, n, fittol);
         }
         else
         {
            nFitToVertices(eorigine, fittol);
         }
         // creation of partitions tied to mesh elements, and interface mesh
         createElementPartitionAndInterfaceFromOneLevelSet(eorigine);
         for (int i = 0; i < eorigine->size(1); ++i)
         {
            createEdgePartitionFromOneLevelSet((AOMD::mEdge *)eorigine->get(1, i));
         }
      }
   }

   xinterface::aomd::modifyAllState(mesh_bnd->getMesh());

   AOMD::AOMD_Util::Instance()->ex_port("interface.msh", &mesh_bnd->getMesh(), 0);

   // Print partitions in a "mesh"
   // xfem::xMesh partitionMesh;
   // mesh->getPartitionMesh(partitionMesh);
   // AOMD::AOMD_Util::Instance()->ex_port("testpartition.msh",&partitionMesh,0);

   createSupportInfo();  /// mesh : classification of elements in subdomains
   for (mEntity *pe : submesh_manager.range_sub(mesh->dim(), "in_domain")) classify_in(pe);
   for (mEntity *pe : submesh_manager.range_sub(mesh->dim(), "out_domain")) classify_out(pe);
}

// Function used for octree cells
void xPhysSurfOctree::nFitToVertices(mEntity *e, double fittol)
{
   if (e->size(1) < 2) throw;
   double nfit = fittol;
   for (int i = 0; i < e->size(1); i++)
   {
      mEdge *q = (mEdge *)e->get(1, i);
      std::vector<double> lsq = getLevelSetValues(*q);
      if (lsq[0] * lsq[1] <= 0.0 && lsq[0] != 0.0 && lsq[1] != 0.0)
      {
         double r = lsq[0] / (lsq[1] - lsq[0]);
         // node too close to the interface
         if (fabs(r) <= nfit)
            setLevelSetValue(static_cast<mVertex &>(*q->get(0, 0)), 0.);
         else if (fabs(1 - r) <= nfit)
            setLevelSetValue(static_cast<mVertex &>(*q->get(0, 1)), 0.);
      }
   }
   return;
}

void xPhysSurfOctree::setLevelSetValue(const mVertex &v, double val) { levelset_value.setData(v) = val; }

std::vector<double> xPhysSurfOctree::getLevelSetValues(const AOMD::mEntity &e) const
{
   const size_t nv = size_t(e.size(0));
   if (nv < 2) throw;
   std::vector<double> lsvalues(nv, 0.);
   for (size_t i = 0; i < nv; ++i) lsvalues[i] = levelset_value.at(*e.get(0, int(i)));
   return lsvalues;
}
double xPhysSurfOctree::getLevelSetValue(const AOMD::mVertex &v) const { return levelset_value.at(v); }

// Indicates if the element has to be refined or not
// Case of material interface : indicates if edges are crossed by the iso-zero value of the level set
//    1 : edge not cut
//    2 : edge cut
//    3 : interface only crossed the edge by two nodes
void xPhysSurfOctree::refineElementSelectionMatInt(const mEntity &e, const xfem::xPointToDouble &lsv, int n, int &refineE)
{
   const size_t nnode = size_t(e.size(0));
   if (!cut_element.getData(e))
   {
      for (int i = 0; i < e.size(1); i++)
      {
         const mEdge &q = *static_cast<mEdge *>(e.get(1, int(i)));
         if (!cut_edge.getData(q)) cut_edge.setData(q) = 1;
         if (cut_edge.at(q) == 2)
         {
            refineE = 1;
            cut_element.setData(e) = 1;
         }
         else
         {
            std::vector<double> lsq = getLevelSetValues(q);
            if (lsq[0] * lsq[1] < 0.0 && lsq[0] != 0. && lsq[1] != 0.)
            {
               cut_edge.setData(q) = 2;
               cut_element.setData(e) = 1;
               refineE = 1;
            }
            else
            {
               const size_t ne = size_t(std::pow(2., n));  // number of refined edges
               const size_t np = ne + 1;                   // number of "refined nodes" on mesh edge
               const mVertex &v0 = *static_cast<mVertex *>(q.get(0, 0));
               const mVertex &v1 = *static_cast<mVertex *>(q.get(0, 1));
               std::vector<xtensor::xPoint> p(np);
               p[0] = v0.point();
               p[np - 1] = v1.point();
               for (size_t j = 1; j < np - 1; j++)
               {
                  p[j] = p[0] + (p[np - 1] - p[0]) * (j / ne);
               }
               for (size_t k = 0; k < np - 2; k++)
               {
                  const double ls0 = lsv(p[k]);
                  const double ls1 = lsv(p[k + 1]);
                  if (ls0 * ls1 <= 0.0 && lsq[0] != 0. &&
                      lsq[1] != 0.)  // an old bug ? should it be ls0 end ls1 instead of lsq0, lsq1 ??
                  {
                     cut_edge.setData(q) = 2;
                     cut_element.setData(e) = 1;
                     refineE = 1;
                  }
               }
               if (cut_edge.at(q) != 2)
               {
                  if (lsq[0] == 0. && lsq[1] == 0.)
                  {
                     cut_edge.setData(q) = 3;
                     cut_element.setData(e) = 1;
                     refineE = 1;
                  }
               }
            }
         }
      }
   }
   else
   {
      cut_element.setData(e) = 0;  // iso-zero value of the level set does not cross the element
   }
   return;
}

void xPhysSurfOctree::doIrefineEdge(const mEdge &q, int n, int &refineE) const
{
   // const std::vector <double> lsq = getLevelSetValues(q);
   const mVertex &v0 = *static_cast<mVertex *>(q.get(0, 0));
   const mVertex &v1 = *static_cast<mVertex *>(q.get(0, 1));
   const double lsq0 = lsv(v0.point());
   const double lsq1 = lsv(v1.point());

   if (lsq0 * lsq1 < 0.0 && lsq0 != 0. && lsq1 != 0.)
   {
      refineE = 1;
      return;
   }
   else if (lsq0 == 0. && lsq1 == 0.)
   {
      refineE = 1;
      return;
   }
   else
   {
      // look at finer scale if the edge is cut ...
      const size_t ne = size_t(std::pow(2., n));  // number of refined edges
      const size_t np = ne + 1;                   // number of "refined nodes" on mesh edge
      const mVertex &v0 = *static_cast<mVertex *>(q.get(0, 0));
      const mVertex &v1 = *static_cast<mVertex *>(q.get(0, 1));
      std::vector<xtensor::xPoint> p(np);
      p[0] = v0.point();
      p[np - 1] = v1.point();
      for (size_t j = 1; j < np - 1; j++)
      {
         p[j] = p[0] + (p[np - 1] - p[0]) * (j / ne);
      }
      for (size_t k = 0; k < np - 2; k++)
      {
         const double ls0 = lsv(p[k]);
         const double ls1 = lsv(p[k + 1]);
         //      if(ls0*ls1 <= 0.0 && lsq[0]!=0. && lsq[1]!=0.) {
         if (ls0 * ls1 <= 0.0)
         {
            refineE = 1;
            return;
         }
      }
   }
}

// Indicates if the element has to be refined or not
void xPhysSurfOctree::refineElementSelection(const mEntity &e, int n, int &refineE)
{
   if (cut_element.getData(e))
   {
      cut_element.setData(e) = 0;
   }
   else
   {
      if (e.size(0) == 2)
      {  // e is an edge ..
         doIrefineEdge(static_cast<const mEdge &>(e), n, refineE);
      }
      else
      {
         for (int i = 0; i < e.size(1); i++)
         {
            const mEdge &q = static_cast<const mEdge &>(*e.get(1, i));
            doIrefineEdge(q, n, refineE);
         }
      }
      if (refineE) cut_element.setData(e) = 1;
   }
   return;
}

//-----------------------------------------------------

// Creates a mesh of refined elements, tied to the mEntity* eorigine
void xPhysSurfOctree::nRefineElement(const mEntity &eorigine, const xfem::xPointToDouble &lsv, bool matint, int n, double fittol)
{
   xfem::xMesh *octree = &refined_elements.setData(eorigine);
   AOMD::mMesh &moctree = octree->getMesh();
   int dim = mesh->dim();
   if (dim == 2)
   {
      const mVertex &v0 = *static_cast<mVertex *>(eorigine.get(0, 0));
      const mVertex &v1 = *static_cast<mVertex *>(eorigine.get(0, 1));
      const mVertex &v2 = *static_cast<mVertex *>(eorigine.get(0, 2));

      mVertex &ov0 = *moctree.createVertex(v0.point(), v0.getClassification());
      mVertex &ov1 = *moctree.createVertex(v1.point(), v1.getClassification());
      mVertex &ov2 = *moctree.createVertex(v2.point(), v2.getClassification());

      root_entity.setData(ov0) = const_cast<AOMD::mVertex *>(&v0);
      root_entity.setData(ov1) = const_cast<AOMD::mVertex *>(&v1);
      root_entity.setData(ov2) = const_cast<AOMD::mVertex *>(&v2);

      const mEdge &e0 = *static_cast<mEdge *>(eorigine.get(1, 0));
      const mEdge &e1 = *static_cast<mEdge *>(eorigine.get(1, 1));
      const mEdge &e2 = *static_cast<mEdge *>(eorigine.get(1, 2));

      mEdge &oe0 = *moctree.createEdge_oneLevel(&ov0, &ov1, e0.getClassification());
      mEdge &oe1 = *moctree.createEdge_oneLevel(&ov1, &ov2, e1.getClassification());
      mEdge &oe2 = *moctree.createEdge_oneLevel(&ov2, &ov0, e2.getClassification());

      root_entity.setData(oe0) = const_cast<AOMD::mEdge *>(&e0);
      root_entity.setData(oe1) = const_cast<AOMD::mEdge *>(&e1);
      root_entity.setData(oe2) = const_cast<AOMD::mEdge *>(&e2);

      mFace &of = *moctree.createFaceWithEdges_oneLevel(&oe0, &oe1, &oe2, eorigine.getClassification());
      root_entity.setData(of) = const_cast<AOMD::mEntity *>(&eorigine);
   }
   if (dim == 3) throw;  // octree.createTetWithVertices(v[0], v[1], v[2], v[3], eorigine.getClassification());

   if (n > 0)
   {
      for (int i = 0; i < n; ++i)
      {
         std::list<mEntity *> elementstodel;
         std::list<mEntity *> elementstotreat;
         for (mEntity *pe : octree->range(dim)) elementstotreat.push_back(pe);
         for (mEntity *pe : elementstotreat)
         {
            mEntity &e = *pe;
            int refineE = 0;
            if (matint)
               refineE = 1;  // case material interface
            else
               refineElementSelection(e, n, refineE);
            if (refineE)
            {
               elementstodel.push_back(&e);
               for (int j = 0; j < e.size(1); j++)
               {
                  mEntity &q = *e.get(1, j);
                  // delete big edges
                  if (!is_the_creator_of2.getData(q))
                  {
                     mVertex &v0 = *static_cast<mVertex *>(q.get(0, 0));
                     mVertex &v1 = *static_cast<mVertex *>(q.get(0, 1));
                     const Trellis_Util::mPoint middle = (v0.point() + v1.point()) * 0.5;  // node : middle of the edge
                     mVertex &vmid = *moctree.createVertex(middle, e.getClassification());
                     mEdge &oe0 = *moctree.createEdge_oneLevel(&v0, &vmid, e.getClassification());
                     mEdge &oe1 = *moctree.createEdge_oneLevel(&vmid, &v1, e.getClassification());
                     root_entity.setData(oe0) = root_entity.at(q);
                     root_entity.setData(oe1) = root_entity.at(q);
                     root_entity.setData(vmid) = root_entity.at(q);
                     is_the_creator_of2.setData(q) = &vmid;
                     was_created_by2.setData(vmid) = &q;
                  }
               }

               if (dim == 2)  // creation of 4 triangles
               {
                  // for (int k = 0; k < e.size(1); k++) elementstodel.push_back(e.get(1, k));
                  mVertex *v0 = static_cast<mVertex *>(e.get(0, 0));
                  mVertex *v1 = static_cast<mVertex *>(e.get(0, 1));
                  mVertex *v2 = static_cast<mVertex *>(e.get(0, 2));
                  mEntity *q0 = e.get(1, 0);
                  mVertex *v3 = static_cast<mVertex *>(is_the_creator_of2.at(*q0));
                  mEntity *q1 = e.get(1, 1);
                  mVertex *v4 = static_cast<mVertex *>(is_the_creator_of2.at(*q1));
                  mEntity *q2 = e.get(1, 2);
                  mVertex *v5 = static_cast<mVertex *>(is_the_creator_of2.at(*q2));

                  mEdge &oe0 = *moctree.createEdge_oneLevel(v3, v5, e.getClassification());
                  root_entity.setData(oe0) = const_cast<AOMD::mEntity *>(&eorigine);
                  mFace *fac_sse0 = moctree.createFaceWithEdges_oneLevel(moctree.getEdge(v0, v3), moctree.getEdge(v3, v5),
                                                                         moctree.getEdge(v5, v0), e.getClassification());
                  root_entity.setData(*fac_sse0) = const_cast<AOMD::mEntity *>(&eorigine);

                  mEdge &oe1 = *moctree.createEdge_oneLevel(v4, v3, e.getClassification());
                  root_entity.setData(oe1) = const_cast<AOMD::mEntity *>(&eorigine);
                  mFace *fac_sse1 = moctree.createFaceWithEdges_oneLevel(moctree.getEdge(v1, v4), moctree.getEdge(v4, v3),
                                                                         moctree.getEdge(v3, v1), e.getClassification());
                  root_entity.setData(*fac_sse1) = const_cast<AOMD::mEntity *>(&eorigine);

                  mEdge &oe2 = *moctree.createEdge_oneLevel(v5, v4, e.getClassification());
                  root_entity.setData(oe2) = const_cast<AOMD::mEntity *>(&eorigine);
                  mFace *fac_sse2 = moctree.createFaceWithEdges_oneLevel(moctree.getEdge(v2, v5), moctree.getEdge(v5, v4),
                                                                         moctree.getEdge(v4, v2), e.getClassification());
                  root_entity.setData(*fac_sse2) = const_cast<AOMD::mEntity *>(&eorigine);

                  mFace *fac_sse3 = moctree.createFaceWithEdges_oneLevel(&oe1, &oe2, &oe0, e.getClassification());
                  root_entity.setData(*fac_sse3) = const_cast<AOMD::mEntity *>(&eorigine);
               }

               else if (dim == 3)  // creation of 8 tetrahedrons
               {                   /*
                                    for (int k = 0; k < e.size(1); k++) elementstodel.push_back(e.get(1, k));
                                    for (int k = 0; k < e.size(2); k++) elementstodel.push_back(e.get(2, k));
                                    mVertex *v0 = (mVertex *)e.get(0, 0);
                                    mVertex *v1 = (mVertex *)e.get(0, 1);
                                    mVertex *v2 = (mVertex *)e.get(0, 2);
                                    mVertex *v3 = (mVertex *)e.get(0, 3);
                  
                                    mEntity *q0 = e.get(1, 0);  // 6 edges
                                    mVertex *v4 = (mVertex *)is_the_creator_of2.getData(*q0);
                                    mEntity *q1 = e.get(1, 1);
                                    mVertex *v5 = (mVertex *)is_the_creator_of2.getData(*q1);
                                    mEntity *q2 = e.get(1, 2);
                                    mVertex *v6 = (mVertex *)is_the_creator_of2.getData(*q2);
                                    mEntity *q3 = e.get(1, 3);
                                    mVertex *v7 = (mVertex *)is_the_creator_of2.getData(*q3);
                                    mEntity *q4 = e.get(1, 4);
                                    mVertex *v8 = (mVertex *)is_the_creator_of2.getData(*q4);
                                    mEntity *q5 = e.get(1, 5);
                                    mVertex *v9 = (mVertex *)is_the_creator_of2.getData(*q5);
                  
                                      // look for the mesh element eorigine corresponding to e
                  
                  
                  
                                    mTet *tet_sse0 = octree.createTetWithVertices(v0, v4, v6, v7, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse0) = const_cast< AOMD::mEntity *>(&eorigine);
                  
                                    mTet *tet_sse1 = octree.createTetWithVertices(v1, v5, v8, v4, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse1) = const_cast< AOMD::mEntity *>(&eorigine);
                  
                                    mTet *tet_sse2 = octree.createTetWithVertices(v2, v9, v6, v5, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse2) = const_cast< AOMD::mEntity *>(&eorigine);
                  
                                    mTet *tet_sse3 = octree.createTetWithVertices(v3, v7, v9, v8, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse3) = const_cast< AOMD::mEntity *>(&eorigine);
                  
                                    mTet *tet_sse4 = octree.createTetWithVertices(v4, v5, v6, v8, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse4) = const_cast< AOMD::mEntity *>(&eorigine);
                  
                                    mTet *tet_sse5 = octree.createTetWithVertices(v7, v4, v6, v8, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse5) = const_cast< AOMD::mEntity *>(&eorigine);
                  
                                    mTet *tet_sse6 = octree.createTetWithVertices(v5, v9, v6, v8, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse6) = const_cast< AOMD::mEntity *>(&eorigine);
                  
                                    mTet *tet_sse7 = octree.createTetWithVertices(v6, v8, v9, v7, e.getClassification());
                                    is_the_element0_of.setData(*tet_sse7) = const_cast< AOMD::mEntity *>(&eorigine);
                                 }
                              }
                               */
               }
            }
         }
         std::set<mEntity *> edge_todel;
         for (mEntity *e : elementstodel)
         {
            root_entity.deleteData(*e);
            cut_element.deleteData(*e);
            for (int i = 0; i < e->size(1); ++i) edge_todel.insert(e->get(1, i));
            moctree.DEL_oneLevel(e);
         }
         for (mEntity *e : edge_todel)
         {
            if (!e->size(2))
            {
               AOMD::mEntity *v = is_the_creator_of2.at(*e);
               is_the_creator_of2.deleteData(*e);
               was_created_by2.deleteData(*v);
               root_entity.deleteData(*e);
               moctree.DEL_oneLevel(e);
            }
         }
      }
      for (mEntity *eoct : octree->range(0))
      {
         mVertex *v = (mVertex *)eoct;
         const double lsve = lsv(v->point());
         setLevelSetValue(*v, lsve);
      }
      for (mEntity *eoct : octree->range(2)) nFitToVertices(eoct, fittol);
   }
   xinterface::aomd::modifyAllState(moctree);
   return;
}

void xPhysSurfOctree::classifyElementOctreeNew(mEntity *egros, mEntity *epetit)
{
   const bool debug = xfem::xdebug_flag;
   if (debug)
   {
      std::cout << "classifyElement... : ENTREE \n\n";
      std::cout << "classifyElement... : entite courante e petit is :\n\n";
      epetit->print();
      std::cout << "\n";
      std::cout << "classifyElement... : entite courante e gros  is :\n\n";
      egros->print();
      std::cout << "\n";
   }
   std::vector<double> vals = getLevelSetValues(*egros);
   xfem::xElement elemgros(egros);
   elemgros.xyz2uvw(epetit->getCentroid());
   double val = elemgros.getInterpoSca(vals);
   if (debug) std::cout << " val is " << val << std::endl;
   if (val < 0.)
   {
      classify_in(epetit);
      if (debug)
         std::cout
             << "classifyElementAccordingToOneLS : ajoute lelement d'interface courant cote negatif de lautre interface\n\n";
   }
   else if (val > 0.)
   {
      classify_out(epetit);
      if (debug)
         std::cout
             << "classifyElementAccordingToOneLS : ajoute lelement interface courant cote positif de lautre interface\n\n\n\n";
   }
   return;
}

void xPhysSurfOctree::createEdgePartitionFromOneLevelSet(mEdge *edge)
{
   // xfem::xAttachableMesh *octree_att = (xfem::xAttachableMesh *)edge->getData(xfem::xMesh::get_refined_elements_tag());
   xfem::xMesh *octree = refined_elements.getData(*edge);
   if (octree)
   {
      xfem::xMesh *Mesh_partition = &partition.setData(*edge);
      AOMD::mMesh &mMesh_partition = Mesh_partition->getMesh();
      std::vector<mEdge *> tosplit;
      tosplit.resize(octree->size(1));
      std::transform(octree->begin(1), octree->end(1), tosplit.begin(), [](mEntity *e) { return static_cast<mEdge *>(e); });

      // else tosplit.push_back(edge);

      for (mEdge *q : tosplit)
      {
         std::vector<double> lsq = getLevelSetValues(*q);
         Trellis_Util::mPoint p[] = {static_cast<mVertex *>(q->get(0, 0))->point(),
                                     static_cast<mVertex *>(q->get(0, 1))->point()};
         mVertex *cv[2];
         for (int k = 0; k < 2; ++k) cv[k] = mMesh_partition.createVertex(p[k], edge->get(0, k)->getClassification());
         if ((lsq[0] * lsq[1] >= 0.0))
         {
            mEdge *e = mMesh_partition.createEdge(cv[0], cv[1], edge->getClassification());
            was_created_by.setData(*e) = edge;
            if (lsq[0] < 0.)
               classify_in(e);
            else
               classify_out(e);
         }
         else
         {
            double r = -lsq[0] / (lsq[1] - lsq[0]);
            Trellis_Util::mPoint inter = p[0] + (p[1] - p[0]) * r;
            mVertex *vinter = mMesh_partition.createVertex(inter, edge->getClassification());
            mEdge *edges[2];
            for (int k = 0; k < 2; ++k)
            {
               edges[k] = mMesh_partition.createEdge(cv[k], vinter, edge->getClassification());
               was_created_by.setData(*edges[k]) = edge;
            }
            if (lsq[0] < 0)
            {
               classify_in(edges[0]);
               classify_out(edges[1]);
            }
            else
            {
               classify_in(edges[1]);
               classify_out(edges[0]);
            }
         }
      }
   }
}

void xPhysSurfOctree::createElementPartitionAndInterfaceFromOneLevelSet(mEntity *ep)
{
   int dim = ep->getLevel();

   std::vector<mEntity *> er;  // vector of octree cells created from e
   int nelem;
   xfem::xMesh *octree = refined_elements.getData(*ep);
   if (octree)
   {
      nelem = octree->size(dim);
      er.reserve(nelem);
      for (mEntity *ec : octree->range(dim)) er.push_back(ec);
   }
   else
   {
      nelem = 1;
      er.resize(1, ep);
   }

   xfem::xMesh *Mesh_partition_m = &partition.setData(*ep);
   AOMD::mMesh &mMesh_partition_m = Mesh_partition_m->getMesh();

   for (int i = 0; i < nelem; i++)
   {
      mEntity *e = er[i];
      int nnode = e->size(0);
      int nedge = e->size(1);
      std::vector<mVertex *> v(nnode, nullptr);
      std::vector<mVertex *> v_part(nnode, nullptr);

      std::vector<mVertex *> interface_vertices_node;
      std::vector<mVertex *> interface_vertices_edge;
      std::vector<mVertex *> nodes_of_interface;

      std::vector<double> lse = getLevelSetValues(*e);
      for (int i = 0; i < nnode; ++i)
      {
         v[i] = dynamic_cast<mVertex *>(e->get(0, i));
         v_part[i] = mMesh_partition_m.createVertex(v[i]->point(), e->getClassification());
      }
      std::vector<mEdge *> q(nedge, (mEdge *)nullptr);
      for (int i = 0; i < nedge; ++i) q[i] = dynamic_cast<mEdge *>(e->get(1, i));
      int i0 = -1;  // last node / ls(v)=0
      int j0 = -1;  // isolated node (one side of ls)
      for (int i = 0; i < nnode; ++i)
      {
         if (lse[i] == 0.0)
         {
            i0 = i;
            mVertex *v_interface = mesh_bnd->getMesh().createVertex(v[i]->point(), e->getClassification());
            nodes_of_interface.push_back(v_interface);
            interface_vertices_node.push_back(v_part[i]);
         }
         if (lse[(i + 1) % nnode] * lse[i] <= 0.0 && lse[(i + 1) % nnode] * lse[(i + 2) % nnode] <= 0.0 &&
             lse[(i + 1) % nnode] != 0.0)
         {
            j0 = (i + 1) % nnode;
         }
      }
      for (int j = 0; j < nedge; ++j)
      {
         mEdge *q = (mEdge *)e->get(1, j);
         std::vector<double> lsq = getLevelSetValues(*q);
         mVertex *v0 = (mVertex *)q->get(0, 0);
         mVertex *v1 = (mVertex *)q->get(0, 1);
         if ((lsq[0] * lsq[1] <= 0.0 && lsq[0] != 0.0 && lsq[1] != 0.0))  // cut edge
         {
            double r = -lsq[0] / (lsq[1] - lsq[0]);
            Trellis_Util::mPoint p0 = v0->point();
            Trellis_Util::mPoint p1 = v1->point();
            Trellis_Util::mPoint inter = p0 + (p1 - p0) * r;
            mVertex *v_new = mMesh_partition_m.createVertex(inter, e->getClassification());
            was_created_by.setData(*v_new) = q;
            mVertex *v_interface = mesh_bnd->getMesh().createVertex(inter, e->getClassification());
            interface_vertices_edge.push_back(v_new);
            nodes_of_interface.push_back(v_interface);
         }
      }

      if (dim == 2)  // interface : 1 edge + 2 nodes
      {
         std::cout << "#### " << interface_vertices_node.size() << " " << interface_vertices_edge.size() << std::endl;
         if (interface_vertices_node.size() < 2 && interface_vertices_edge.size() == 0)
         {
            mFace *face_new = mMesh_partition_m.createFaceWithVertices(v_part[0], v_part[1], v_part[2], e->getClassification());
            was_created_by.setData(*face_new) = e;
            classifyElementOctreeNew(e, face_new);
         }
         else if (interface_vertices_node.size() == 2)
         {
            mFace *face_new = mMesh_partition_m.createFaceWithVertices(v_part[0], v_part[1], v_part[2], e->getClassification());
            mEdge *edge_int =
                mesh_bnd->getMesh().createEdge(nodes_of_interface[0], nodes_of_interface[1], e->getClassification());
            was_created_by.setData(*edge_int) = e;
            // classify_interface(edge_int);
            was_created_by.setData(*face_new) = e;
            classifyElementOctreeNew(e, face_new);
         }
         else if (interface_vertices_node.size() == 1 && interface_vertices_edge.size() == 1)
         {
            mFace *face_new1 = mMesh_partition_m.createFaceWithVertices(v_part[i0], v_part[(i0 + 1) % nnode],
                                                                        interface_vertices_edge[0], e->getClassification());
            mFace *face_new2 = mMesh_partition_m.createFaceWithVertices(v_part[i0], interface_vertices_edge[0],
                                                                        v_part[(i0 + 2) % nnode], e->getClassification());

            mEdge *edge_int =
                mesh_bnd->getMesh().createEdge(nodes_of_interface[0], nodes_of_interface[1], e->getClassification());

            was_created_by.setData(*edge_int) = e;
            was_created_by.setData(*face_new1) = e;
            was_created_by.setData(*face_new2) = e;
            classifyElementOctreeNew(e, face_new1);
            classifyElementOctreeNew(e, face_new2);
         }
         else if (interface_vertices_node.size() == 0 && interface_vertices_edge.size() == 2)
         {
            //  j0 : last cut edge &  number of isolated node if 2 cut edges
            // take care of faces orientation : if j0=0, opposite way
            mFace *face_new1;
            mFace *face_new2;
            mFace *face_new3;
            if (j0 == 0)
            {
               face_new1 = mMesh_partition_m.createFaceWithVertices(v_part[j0], interface_vertices_edge[0],
                                                                    interface_vertices_edge[1], e->getClassification());
               face_new2 = mMesh_partition_m.createFaceWithVertices(v_part[(j0 + 1) % nnode], interface_vertices_edge[1],
                                                                    interface_vertices_edge[0], e->getClassification());
               face_new3 = mMesh_partition_m.createFaceWithVertices(v_part[(j0 + 2) % nnode], interface_vertices_edge[1],
                                                                    v_part[(j0 + 1) % nnode], e->getClassification());
            }
            else
            {
               face_new1 = mMesh_partition_m.createFaceWithVertices(v_part[j0], interface_vertices_edge[1],
                                                                    interface_vertices_edge[0], e->getClassification());
               face_new2 = mMesh_partition_m.createFaceWithVertices(v_part[(j0 + 1) % nnode], interface_vertices_edge[0],
                                                                    interface_vertices_edge[1], e->getClassification());
               face_new3 = mMesh_partition_m.createFaceWithVertices(v_part[(j0 + 2) % nnode], interface_vertices_edge[0],
                                                                    v_part[(j0 + 1) % nnode], e->getClassification());
            }
            was_created_by.setData(*face_new1) = e;
            was_created_by.setData(*face_new2) = e;
            was_created_by.setData(*face_new3) = e;

            mEdge *edge_int =
                mesh_bnd->getMesh().createEdge(nodes_of_interface[0], nodes_of_interface[1], e->getClassification());
            // classify_interface(edge_int);
            was_created_by.setData(*edge_int) = e;
            classifyElementOctreeNew(e, face_new1);
            classifyElementOctreeNew(e, face_new2);
            classifyElementOctreeNew(e, face_new3);

            xinterface::aomd::modifyAllState(mMesh_partition_m);
            xinterface::aomd::modifyAllState(mesh_bnd->getMesh());
         }
         else
            throw;
      }
      // °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°  A FINIR
      if (dim == 3)
      {
         // The code below is very strange to me ... It seems that it does nothing .. And my compiler aggrees with me. I commented
         // it out Nicolas C 9 mai 2016. cout << " °° interface_vertices_node.size() is : "<< interface_vertices_node.size() <<
         // endl; cout << " °° interface_vertices_edge.size() is : "<< interface_vertices_edge.size() << endl;
         if (interface_vertices_node.size() == 0 && interface_vertices_edge.size() == 3)
         {
            //	     mFace* face_inter = interface.createFaceWithVertices( nodes_of_interface[0], nodes_of_interface[1] ,
            // nodes_of_interface[2] , e->getClassification());

            //	     mTet* tet_1 = Mesh_partition_m.createTetWithVertices( interface_vertices_edge[0], interface_vertices_edge[1]
            //, interface_vertices_edge[2], v_part[j0], e->getClassification());

            // int j1=-1;
            // int j2=-1;
            // int j3 = -1;
            // mEdge* q0 = (mEdge *) v[j0]->get(1,0); // car v_part[j0] n'est associé a aucune arete
            // cout << " q0 is : "<< q0 << endl;
            // mVertex* vq0 = (mVertex*) q0->get(0,0);
            // cout << " vq0 is : "<< vq0 << endl;
            // vq0->print();
            /*
              mVertex* vj1;
              if( vq0 != v[j0] )  vj1 = vq0;
              else
              {
              vj1 = (mVertex*) q0->get(0,1);
              //cout << " ** vj1 : " << vj1 << endl;
              //vj1->print();
              }
            */
            // cout << " vj1 is : "<< vj1 << endl;
            // vj1->print();
            //	     mVertex* vinter0 = (mVertex *) q0->getAttachedEntity(xfem::xMesh::is_the_creator_of_tag);
            // cout << " vinter0 is : "<< vinter0<< endl;

            /*for(int k=0; k<nnode; k++)
              {
              if(v[k]==vj1)  j1 = k; // v_part[] ou v[] ?
              }*/

            // cout << " j1 is : "<< j1 << endl;
            // cout << " v[j1] is : "<<  v[j1] << endl;
            // cout << " v_part[j1] is : "<<  v_part[j1] << endl;

            // mEdge* q1 = (mEdge *) v[j0]->get(1,1);
            // cout << " q1 is : "<< q1 << endl;
            // mVertex* vq1 = (mVertex*) q1->get(0,0);
            // cout << " vq1 is : "<< vq1 << endl;
            // vq1->print();
            /* mVertex* vj2;
               if( vq1 != v[j0] )  vj2 = vq1;
               else vj2 = (mVertex*) q0->get(0,1);
            */
            // cout << " vj2 is : "<< vj2 << endl;
            // vj2->print();
            //	     mVertex* vinter1 = (mVertex *) q1->getAttachedEntity(xfem::xMesh::is_the_creator_of_tag);
            // cout << " vinter1 is : "<< vinter1<< endl;
            /*for(int k=0; k<nnode; k++)
              {
              if(v[k]==vj2)  j2 = k;
              }*/
            // cout << " j2 is : "<< j2 << endl;
            // cout << " v[j2] is : "<<  v[j2] << endl;
            // cout << " v_part[j2] is : "<<  v_part[j2] << endl;

            // mEdge* q2 = (mEdge *) v[j0]->get(1,2);
            // cout << " q2 is : "<< q2 << endl;
            // mVertex* vq2 = (mVertex*) q2->get(0,0);
            // cout << " vq2 is : "<< vq2 << endl;
            // vq2->print();
            /* mVertex* vj3;
               if( vq2 != v[j0] )  vj3 = vq2;
               else vj3 = (mVertex*) q0->get(0,1);
            */
            // cout << " vj3 is : "<< vj3 << endl;
            // vj3->print();
            //	     mVertex* vinter2 = (mVertex *) q2->getAttachedEntity(xfem::xMesh::is_the_creator_of_tag);
            // cout << " vinter2 is : "<< vinter2<< endl;
            /*for(int k=0; k<nnode; k++)
              {
              if(v[k]==vj3)  j3 = k;
              }
            */
            // cout << " j3 is : "<< j3 << endl;
            // cout << " v[j3] is : "<<  v[j3] << endl;
            // v[j3]->print();
            // cout << " v_part[j3] is : "<<  v_part[j3] << endl;
            // v_part[j3]->print();
         }
      }

      // °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°  A FINIR !!  24-07-09

      //      v.resize(1,0);
      //      interface_vertices_node.resize(1,0);
      //      interface_vertices_edge.resize(1,0);
   }

   // cout << "avant mdfAS " << endl;

   xinterface::aomd::modifyAllState(mMesh_partition_m);
   // cout << "Mesh_partition_m->size(2) 2  " << Mesh_partition_m->size(2) << endl;
   // AOMD_Util::Instance()->ex_port("mesh_partition_m.msh", Mesh_partition_m);
   // cout << "apres mdfAS 1 " << endl;
   xinterface::aomd::modifyAllState(mesh_bnd->getMesh());
   // cout << "apres mdfAS 2" << endl;
   return;
}

// ---------------------------------------------------------------

}  // namespace xcut
