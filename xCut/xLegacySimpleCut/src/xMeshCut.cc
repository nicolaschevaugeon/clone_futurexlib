/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

// AOMD
#include "AOMD_OwnerManager.h"
#include "mAttachableDataContainer.h"
#include "mEdge.h"
#include "mEntity.h"
#include "mFace.h"
#include "mTet.h"
#include "mVertex.h"

// xinterface/aomd
#include "xAOMDEntityUtil.h"
#include "xMeshToolAOMD.h"

// xfem
#include "xDebug.h"
#include "xElement.h"
#include "xLevelSet.h"
#include "xMesh.h"

// xcut
#include "xMeshCut.h"

using namespace AOMD;
using AOMD::mEdge;
using std::cout;
using std::endl;
using std::swap;

namespace xcut
{
using namespace xfem;

mEntity* getOriginalCreator(mEntity* egros, const datamanager_t<AOMD::mEntity*>& was_created_by)
{
   mEntity* const* pcreator = was_created_by.getData(*egros);
   if (!pcreator)
      return egros;
   else
      return getOriginalCreator(*pcreator, was_created_by);
}

// works only for planar convex quads. The vertex p point to are not modifyed, they are only reordered.
void _orientQuad(std::vector<mVertex*>& p)
{
   if (p.size() != 4) throw;
   xtensor::xVector<> vec1 = xtensor::xVector<>(p[0]->point(), p[1]->point()) % xtensor::xVector<>(p[1]->point(), p[2]->point());
   xtensor::xVector<> vec2 = xtensor::xVector<>(p[1]->point(), p[2]->point()) % xtensor::xVector<>(p[2]->point(), p[3]->point());
   xtensor::xVector<> vec3 = xtensor::xVector<>(p[0]->point(), p[1]->point()) % xtensor::xVector<>(p[1]->point(), p[3]->point());
   if (vec1 * vec2 < 0.)
   {
      std::swap(p[2], p[3]);
   }
   else
   {
      if (vec1 * vec3 < 0.)
      {
         std::swap(p[1], p[2]);
         std::swap(p[1], p[3]);
      }
   }
}

void _classifyElement(mEntity* egros, mEntity* epetit, const xLevelSet& ls, const datamanager_t<AOMD::mEntity*>& was_created_by,
                      xEntityToEntity& classify_in, xEntityToEntity& classify_out)
{
   const bool debug = xdebug_flag;
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
   xtensor::xPoint g = xinterface::aomd::centroid(*epetit);
   if (debug)
   {
      std::cout << "centroid is : " << g(0) << " " << g(1) << " " << g(2) << std::endl;
      std::cout << " egros is : " << std::endl;
      egros->print();
   }
   double val = 0.;
   mEntity* creator = xcut::getOriginalCreator(egros, was_created_by);

   if (debug)
   {
      std::cout << " creator is " << std::endl;
      creator->print();
      for (int i = 0; i <= 2; i++)
      {
         std::cout << "egros->size(" << i << ") " << egros->size(i) << std::endl;
         for (int j = 0; j < egros->size(i); j++)
         {
            std::cout << "sous-element " << i << " " << j << " de egros " << std::endl;
            (egros->get(i, j))->print();
            mEntity* testEnti = egros->get(i, j);
            //	      std::cout<<"dont le creator est (->getAttachedEntity(xMesh::was_created_by_tag) "<<std::endl;
            mEntity* const* p_was_created_by = was_created_by.getData(*testEnti);
            if (p_was_created_by)
            {
               std::cout << "dont le creator est (->getAttachedEntity(xMesh::was_created_by_tag) " << std::endl;
               (*p_was_created_by)->print();
               // mEntity* attEnti = testEnti->getAttachedEntity(xMesh::was_created_by_tag);
               //		  std::cout<<"test de la LS (avec ls.getVals(attEnti))"<<std::endl;
               //		  ls.getVals(attEnti); std::cout<<"la LS testee est definie"<<std::endl;
            }
            else
            {
               std::cout << " n'a pas de creator" << std::endl;
            }
            std::cout << std::endl << std::endl;
         }
      }
      for (int i = 0; i <= 2; i++)
      {
         for (int j = 0; j < creator->size(i); j++)
         {
            std::cout << "sous-element " << i << " " << j << " du creator " << std::endl;
            (creator->get(i, j))->print();
         }
      }
   }
   if (debug) std::cout << " creator->getLevel() is " << creator->getLevel() << std::endl;
   if (creator->getLevel() == 0)
      val = ls(static_cast<mVertex*>(creator));
   else if (creator->getLevel() == 1)
   {
      xGeomElem edge(creator);
      edge.setUVWForXYZ(g);
      double u = edge.getUVW()(0);
      mVertex* v0 = static_cast<mVertex*>(creator->get(0, 0));
      mVertex* v1 = static_cast<mVertex*>(creator->get(0, 1));
      val = 0.5 * (1. - u) * ls(v0) + 0.5 * (1. + u) * ls(v1);
   }
   else
   {
      xElement el(creator);
      if (debug)
      {
         std::cout << std::endl << "creator before getVal fatidique is " << std::endl;
         creator->print();
         for (int i = 0; i <= 2; i++)
         {
            for (int j = 0; j < creator->size(i); j++)
            {
               std::cout << "sous-element " << i << " " << j << " du creator " << std::endl;
               (creator->get(i, j))->print();
            }
         }
         std::cout << "std::vector<double> valsl=ls.getVals(creator);" << std::endl;
         creator->get(0, 0)->print();
         creator->get(0, 1)->print();
         creator->get(0, 2)->print();
         std::vector<double> valsl = ls.getVals(creator);
         std::cout << "ls on creator " << std::endl << "vals 0 = " << valsl[0] << "\n";
         creator->get(0, 0)->print();
         std::cout << "vals 1 = " << valsl[1] << "\n";
         creator->get(0, 1)->print();
         std::cout << "vals 2 = " << valsl[2] << "\n";
         creator->get(0, 2)->print();
      }
      el.xyz2uvw(g);
      if (debug)
      {
         std::cout << std::endl
                   << "el.getUvw is : " << el.getUvw()(0) << " " << el.getUvw()(1) << " " << el.getUvw()(2) << std::endl;
         std::cout << "ls.getVals(creator) size :  " << ls.getVals(creator).size() << std::endl;

         std::cout << "ls.getVals(creator) :  " << ls.getVals(creator)[0] << " " << ls.getVals(creator)[1] << " "
                   << ls.getVals(creator)[2] << std::endl;
      }
      val = ls.getVal(creator, el.getUvw());
   }
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
}

void _processVertex(xfem::xMesh& imesh, mVertex& v, const xfem::xLevelSet& ls, datamanager_t<AOMD::mEntity*>& is_the_creator_of,
                    datamanager_t<AOMD::mEntity*>& was_created_by)
{
   const bool debug = false;
   const double val = ls(&v);
   if (val == 0.0)
   {
      mVertex& v_new = *imesh.getMesh().createVertex(v.point(), v.getClassification());
      is_the_creator_of.setData(v) = &v_new;
      was_created_by.setData(v_new) = &v;
      if (debug) std::cout << " creating an vertex with address " << &v_new << std::endl;
      if (debug) v_new.print();
   }
}

void _processEdge(xfem::xMesh& imesh, mEdge& e, const xfem::xLevelSet& ls, datamanager_t<mEntity*>& is_the_creator_of,
                  datamanager_t<mEntity*>& was_created_by, datamanager_t<double>& r_on_edge)
{
   const bool debug = false;
   if (debug) cout << "In process Edge\n";
   mVertex* v0 = static_cast<mVertex*>(e.get(0, 0));
   mVertex* v1 = static_cast<mVertex*>(e.get(0, 1));
   mEntity** pcv0 = is_the_creator_of.getData(*v0);
   mEntity** pcv1 = is_the_creator_of.getData(*v1);
   if (pcv0 && pcv1)  // e is an iso_zeo edge. clone it in the interface mesh.
   {
      mEdge& edge_new =
          *imesh.getMesh().createEdge(static_cast<mVertex*>(*pcv0), static_cast<mVertex*>(*pcv1), e.getClassification());
      was_created_by.setData(edge_new) = &e;
      if (debug) std::cout << " creating an edge with address " << &edge_new << std::endl;
      if (debug) edge_new.print();
      return;
   }

   const double lsv0 = ls(v0);
   const double lsv1 = ls(v1);

   if (lsv0 * lsv1 <= 0.0 && lsv0 != 0.0 && lsv1 != 0.0)  // edge is cut (strictly)
   {
      const double r = -lsv0 / (lsv1 - lsv0);
      const Trellis_Util::mPoint p0 = v0->point();
      const Trellis_Util::mPoint p1 = v1->point();
      const Trellis_Util::mPoint inter = p0 + (p1 - p0) * r;
      mVertex& v_new = *imesh.getMesh().createVertex(inter, e.getClassification());
      is_the_creator_of.setData(e) = &v_new;
      was_created_by.setData(v_new) = &e;
      r_on_edge.setData(v_new) = r;
      if (debug)
      {
         std::cout << " creating a vertex with address " << &v_new << std::endl;
         v_new.print();
         cout << "Attached on\n";
         e.print();
      }
   }
}

/// createIsoZeroMesh:
/// Takes a mesh (dim d) and the level-set ls that will cut it
/// Gives the "interface": mesh (of dim d-1) where ls=0
/// For each entity in the interface: info concerning its creator is attached
/// The initial mesh and relative attached data are not modified
/// // notes : maybe it should be a xMesh Constructor ...or a function that return a pointer to a new xMesh ...
// What if the level set has no support ...
void createIsoZeroMeshFromLevelSet(const xfem::xLevelSet& ls, xfem::xMesh& imesh, datamanager_t<AOMD::mEntity*>& was_created_by,
                                   datamanager_t<double>& r_on_edge, const bool simplex_only, const bool debug)
{
   const xfem::xRegion& region = ls.getSupport();

   int dim = region.dim();
   datamanager_t<AOMD::mEntity*> is_the_creator_of;
   if (debug) std::cout << "createInterface : examen des noeuds qui sont au nombre de " << region.size(0) << std::endl;
   for (mEntity* pv : region.range(0)) _processVertex(imesh, static_cast<mVertex&>(*pv), ls, is_the_creator_of, was_created_by);
   if (debug) std::cout << "createInterface : examen des segments qui sont au nombre de " << region.size(1) << std::endl;
   for (mEntity* pe : region.range(1))
      _processEdge(imesh, static_cast<mEdge&>(*pe), ls, is_the_creator_of, was_created_by, r_on_edge);

   if (debug) std::cout << "createInterface : examen des faces qui sont au nombre de " << region.size(2) << std::endl;
   for (mEntity* e : region.range(2))
   {
      if (debug) std::cout << " e->size(0) " << e->size(0) << " e->size(1) " << e->size(1) << std::endl;
      std::vector<mVertex*> zero_on_vertex, zero_on_edge;
      for (int j = 0; j < e->size(0); j++)
      {
         mEntity* q = e->get(0, j);
         if (mEntity** pv = is_the_creator_of.getData(*q)) zero_on_vertex.push_back(static_cast<mVertex*>(*pv));
      }

      for (int j = 0; j < e->size(1); j++)
      {
         mEntity* q = e->get(1, j);
         if (mEntity** pv = is_the_creator_of.getData(*q)) zero_on_edge.push_back(static_cast<mVertex*>(*pv));
      }
      if (debug)
         std::cout << " zero_on_vertex.size() " << zero_on_vertex.size() << " zero_on_edge.size()   " << zero_on_edge.size()
                   << std::endl;
      assert(zero_on_vertex.size() + zero_on_edge.size() <= 3);

      if ((zero_on_edge.size() || zero_on_vertex.size() == 2) && e->getType() == mEntity::QUAD)
      {
         // throw;
         // attachSimplicies(e);  //<- NOTE ... need to check what is the impact of this ..
         if (zero_on_vertex.size() == 2)
         {
            mEdge* edge_new = imesh.getMesh().createEdge(zero_on_vertex[0], zero_on_vertex[1], e->getClassification());
            was_created_by.setData(*edge_new) = e;
            is_the_creator_of.setData(*e) = edge_new;
            if (debug) std::cout << " creating an edge with address " << edge_new << std::endl;
            if (debug) edge_new->print();
         }
      }

      if (zero_on_edge.size() == 2)
      {
         mEdge* edge_new = imesh.getMesh().createEdge(zero_on_edge[0], zero_on_edge[1], e->getClassification());
         was_created_by.setData(*edge_new) = e;
         is_the_creator_of.setData(*e) = edge_new;
         if (debug) std::cout << " creating an edge with address " << edge_new << std::endl;
         if (debug) edge_new->print();
      }
      else if (zero_on_vertex.size() == 1 && zero_on_edge.size() == 1)
      {
         mEdge* edge_new = imesh.getMesh().createEdge(zero_on_edge[0], zero_on_vertex[0], e->getClassification());
         was_created_by.setData(*edge_new) = e;
         is_the_creator_of.setData(*e) = edge_new;
         if (debug) std::cout << " creating an edge with address " << edge_new << std::endl;
         if (debug) edge_new->print();
      }
      else if (zero_on_vertex.size() == 3 && dim == 3)
      {
         mFace* face_new = imesh.getMesh().createFaceWithVertices(zero_on_vertex[0], zero_on_vertex[1], zero_on_vertex[2],
                                                                  e->getClassification());
         was_created_by.setData(*face_new) = e;
         is_the_creator_of.setData(*e) = face_new;
         if (debug) std::cout << " creating a face on tri (all vertices at zero) with address " << face_new << std::endl;
         if (debug) face_new->print();
      }
   }

   if (debug) std::cout << "createInterface : examen des elements 3D\n";
   for (mEntity* e : region.range(3))
   {
      std::vector<mVertex*> ends;
      int node = 0, edge = 0;
      for (int lower = 0; lower < 2; lower++)
      {
         for (int j = 0; j < e->size(lower); j++)
         {
            mEntity* q = e->get(lower, j);
            if (mEntity** pv = is_the_creator_of.getData(*q)) ends.push_back(static_cast<mVertex*>(*pv));
         }
         if (lower == 0) node = int(ends.size());
         if (lower == 1) edge = int(ends.size()) - node;
      }

      if ((edge > 0 || node == 4) && e->getType() == mEntity::HEX)
      {
         // throw;
         // attachSimplicies(e);  //<- WARNING check this.
      }
      if (node + edge == 3 && edge > 0)
      {
         mFace* fac_new = imesh.getMesh().createFaceWithVertices(ends[0], ends[1], ends[2], e->getClassification());
         was_created_by.setData(*fac_new) = e;
         is_the_creator_of.setData(*e) = fac_new;
         if (debug) std::cout << " creating a face on tri with address " << fac_new << std::endl;
         if (debug) fac_new->print();
      }
      else if (node + edge == 4 && edge >= 0)  //>= added by greg (quad/hex)
      {
         if (ends.size() != 4) throw;
         _orientQuad(ends);
         if (simplex_only)
         {
            // decomposition en 2 triangles orientes
            mFace* fac1_new = imesh.getMesh().createFaceWithVertices(ends[0], ends[1], ends[2], e->getClassification());
            was_created_by.setData(*fac1_new) = e;
            if (debug) std::cout << " creating a face on quad with address " << fac1_new << std::endl;
            if (debug) fac1_new->print();
            mFace* fac2_new = imesh.getMesh().createFaceWithVertices(ends[2], ends[3], ends[0], e->getClassification());
            was_created_by.setData(*fac2_new) = e;
            if (debug) std::cout << " creating a face on quad with address " << fac2_new << std::endl;
            if (debug) fac2_new->print();
            // creation d un segment entre les deux triangles
            mEdge* edge_new = imesh.getMesh().createEdge(ends[0], ends[2], e->getClassification());
            if (debug) std::cout << " creating an edge diag to a quad with address " << edge_new << std::endl;
            if (debug) edge_new->print();
            was_created_by.setData(*edge_new) = e;
            is_the_creator_of.setData(*e) = edge_new;
         }
         else
         {
            if (debug) std::cout << " creation of a quadrilateral " << std::endl;
            // pas de decoupage ici (decomposition)
            // on cree simplement la face a partir des 4 noeuds
            mFace* face_new = imesh.getMesh().createFaceWithVertices(ends[0], ends[1], ends[2], ends[3], e->getClassification());
            was_created_by.setData(*face_new) = e;
            is_the_creator_of.setData(*e) = face_new;
         }
      }
   }  // en loop of 3d entities

   // Start Parallel Link between elem of the interface mesh // Note : this need to be rewritten in terms of dataexchange


   imesh.getMesh().modifyState(2, 1, true);
   imesh.getMesh().modifyState(1, 2, true);
   imesh.getMesh().modifyState(0, 1, true);
   if (debug) std::cout << "createInterface : SORTIE\n\n";
   // return imesh;
}

// add simplicies derived from e to the xMesh submesh
void _addSimplicies(const mEntity& e, xMesh& submesh, std::function<mVertex&(const mVertex&)> cloneVertex)
{
   AOMD::mMesh& submmesh = submesh.getMesh();
   const bool debug = false;
   if (debug) cout << "Attach simplicies on :\n";
   if (debug) e.print();
   auto eclass = e.getClassification();
   switch (e.getType())
   {
      case mEntity::QUAD:
      {
         std::array<mVertex*, 4> v_new;
         for (size_t i = 0; i < 4; ++i)
         {
            mVertex* v = static_cast<mVertex*>(e.get(0, int(i)));
            v_new[i] = &cloneVertex(*v);
         }
         // 012 - 023
         mEdge* ed01 = submmesh.createEdge(v_new[0], v_new[1], eclass);
         mEdge* ed03 = submmesh.createEdge(v_new[0], v_new[3], eclass);
         mEdge* ed02 = submmesh.createEdge(v_new[0], v_new[2], eclass);
         mEdge* ed12 = submmesh.createEdge(v_new[1], v_new[2], eclass);
         mEdge* ed23 = submmesh.createEdge(v_new[2], v_new[3], eclass);
         submmesh.createFaceWithEdges(ed01, ed12, ed02, eclass);
         submmesh.createFaceWithEdges(ed02, ed23, ed03, eclass);
         xinterface::aomd::modifyAllState(submmesh);
         return;
      }
      case mEntity::HEX:
      {
         std::array<mVertex*, 8> v_new;
         for (size_t i = 0; i < 8; ++i)
         {
            mVertex* v = static_cast<mVertex*>(e.get(0, int(i)));
            v_new[i] = &cloneVertex(*v);
         }
         submmesh.createTetWithVertices(v_new[3], v_new[7], v_new[5], v_new[6], eclass);
         submmesh.createTetWithVertices(v_new[3], v_new[6], v_new[5], v_new[2], eclass);
         submmesh.createTetWithVertices(v_new[3], v_new[4], v_new[5], v_new[7], eclass);
         submmesh.createTetWithVertices(v_new[3], v_new[2], v_new[5], v_new[1], eclass);
         submmesh.createTetWithVertices(v_new[3], v_new[1], v_new[5], v_new[0], eclass);
         submmesh.createTetWithVertices(v_new[3], v_new[0], v_new[5], v_new[4], eclass);
         xinterface::aomd::modifyAllState(submmesh);
         return;
      }
      default:
         throw;
   }
}

// This version assume vertex -> edge and edge->vertex
const mEdge* _getEdge(const mVertex& v0, const mVertex& v1)
{
   for (int i = 0; i < v0.size(1); ++i)
   {
      const mEdge& e = static_cast<mEdge&>(*v0.get(1, i));
      const mVertex& otherv = *static_cast<const mVertex*>((e.get(0, 0) == &v0) ? e.get(0, 1) : e.get(0, 0));
      if (&otherv == &v1) return &e;
   }
   return nullptr;
}

void _cutEdge(const mEdge& e, const mVertex& v_interface, std::function<mVertex&(const mVertex&)> cloneVertex, xMesh& submesh,
              xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& is_in_partition_of)
{
   AOMD::mMesh& submmesh = submesh.getMesh();
   mVertex& v0 = static_cast<mVertex&>(*e.get(0, 0));
   mVertex& v1 = static_cast<mVertex&>(*e.get(0, 1));
   mVertex& vnewi = cloneVertex(v_interface);  //*submesh.createVertex(v_interface.point(), v_interface.getClassification());
   mVertex& vnew0 = cloneVertex(v0);           //*submesh.createVertex(v0.point(), v0.getClassification());
   mVertex& vnew1 = cloneVertex(v1);           //*submesh.createVertex(v1.point(), v1.getClassification());
   mEdge& enew0 = *submmesh.createEdge(&vnew0, &vnewi, e.getClassification());
   mEdge& enew1 = *submmesh.createEdge(&vnewi, &vnew1, e.getClassification());
   is_in_partition_of.setData(enew0) = const_cast<mEdge*>(&e);
   is_in_partition_of.setData(enew1) = const_cast<mEdge*>(&e);
   xinterface::aomd::modifyAllState(submmesh);
}

/*
void _cutFace(const mFace &f, const mEdge &e_interface, xMesh &submesh,
              xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>&is_duplicated_in,
              xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>&was_duplicated_from,
              xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>&is_in_partition_of){
    if (f.getType() != TRI) throw;


}*/

// cutElement:
// From an interface resulting from an element cut away by a level-set,
// generates the consequent elementary submesh (elementary_sub)
// precondition
// e_interface  must be of the level of e -1.
// clean up if any partition attached to e.
// create a new one
void _cutElement(const mEntity& e, const mEntity& e_interface,
                 const xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& c_was_created_by, xMesh& submesh,
                 xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& is_duplicated_in,
                 xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& was_duplicated_from,
                 xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& is_in_partition_of)
{
   // !note, the baseMesh is needed because we need a "non-local" Version of ge(v0,v1) ... since we can not guaranties for now
   // that v0 and v1 know there edges.
   // !
   const bool debug = xdebug_flag;
   const xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& c_is_duplicated_in = is_duplicated_in;
   const xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& c_was_duplicated_from = was_duplicated_from;
   AOMD::mMesh& submmesh = submesh.getMesh();

   auto vertexInBetween = [&c_was_created_by, &c_is_duplicated_in, &c_was_duplicated_from](
                              const mVertex& v0, const mVertex& v1, const std::vector<mVertex*>& vertices) -> mVertex* {
      const mVertex& va = *static_cast<mVertex*>(c_was_duplicated_from.at(v0));
      const mVertex& vb = *static_cast<mVertex*>(c_was_duplicated_from.at(v1));
      const mEdge& ed = *_getEdge(va, vb);
      for (mVertex* v : vertices)
      {
         mEdge* pe_v = dynamic_cast<mEdge*>(c_was_created_by.at(*v));
         if (&ed == pe_v) return static_cast<mVertex*>(c_is_duplicated_in.at(*v));
      }
      return nullptr;
   };

   auto cloneVertex = [&submesh, &is_duplicated_in, &was_duplicated_from](const mVertex& v) -> mVertex& {
      mVertex& vnew = *submesh.getMesh().createVertex(v.point(), v.getClassification());
      is_duplicated_in.setData(v) = &vnew;
      was_duplicated_from.setData(vnew) = const_cast<mVertex*>(&v);
      return vnew;
   };

   static int nbCutElementCall = 0;

   submesh.getMesh().setNewVertexId(20 * nbCutElementCall);
   int memo = -1;
   int dim = e.getLevel();
   if (dim == 1)
   {
      if (e_interface.getLevel() != 0) throw;
      return _cutEdge(static_cast<const mEdge&>(e), static_cast<const mVertex&>(e_interface), cloneVertex, submesh,
                      is_in_partition_of);
   }
   // we gather the vertoces of the interface in  e_interface_vertices
   // if e_interface is a triangle, it can be one of two triangles
   // composing the cut inside a tetrahedron
   // (in the case where a tetrahedron is cut along 4 of its edges).
   std::vector<mVertex*> e_interface_vertices;
   auto tomvertexp = [](const mEntity* e) { return const_cast<mVertex*>((static_cast<const mVertex*>(e))); };

   for (int i = 0; i < e_interface.size(0); ++i)
   {
      e_interface_vertices.push_back(tomvertexp(e_interface.get(0, i)));
   }
   if (e_interface.getLevel() == 2)
   {
      for (int i = 0; i < e_interface.size(1); ++i)
      {
         mEntity* ed = e_interface.get(1, i);
         mEntity* const* p_was_created_by = c_was_created_by.getData(*ed);
         if ((p_was_created_by) && (*p_was_created_by == &e))
         {
            std::set<mVertex*> e_tmp_set(e_interface_vertices.begin(), e_interface_vertices.end());
            for (int j = 0; j < ed->size(2); ++j)
            {
               mEntity* f = ed->get(2, j);
               for (int k = 0; k < f->size(0); ++k) e_tmp_set.insert(tomvertexp(f->get(0, k)));
            }
            if (e_tmp_set.size() != 4) throw;
            e_interface_vertices.resize(4);
            std::copy(e_tmp_set.begin(), e_tmp_set.end(), e_interface_vertices.begin());
            _orientQuad(e_interface_vertices);
         }
      }
   }

   std::vector<mVertex*> interface_duplicated_vertices;
   std::vector<mVertex*> duplicated_side1_vertices;
   std::vector<mVertex*> duplicated_side2_vertices;
   //
   // duplication des noeuds de l'interface
   //
   // std::cout << "#############: " <<  e_interface_vertices.size() << std::endl;
   for (mVertex* v : e_interface_vertices)
   {
      mVertex& v_new = cloneVertex(*v);
      interface_duplicated_vertices.push_back(&v_new);
   }

   xtensor::xVector<> n;
   // int numdupvert = interface_duplicated_vertices.size();
   // if (numdupvert < 2) std::cout << "Warning in xCutFunction 1"<< std ::std::endl;
   mVertex* v1 = interface_duplicated_vertices[0];

   if (interface_duplicated_vertices.size() < 2) throw;
   mVertex* v2 = interface_duplicated_vertices[1];
   if (dim == 3)
   {
      //  if (numdupvert < 3) std::cout << "Warning in xCutFunction 2"<< std ::endl;
      mVertex* v3 = interface_duplicated_vertices[2];
      n = (xtensor::xVector<>(v1->point(), v2->point())) % (xtensor::xVector<>(v1->point(), v3->point()));
   }

   if (dim == 2)
   {
      mVertex* v1 = interface_duplicated_vertices[0];
      n = (xtensor::xVector<>(tomvertexp(e.get(0, 0))->point(), tomvertexp(e.get(0, 1))->point()) %
           xtensor::xVector<>(tomvertexp(e.get(0, 0))->point(), tomvertexp(e.get(0, 2))->point())) %
          xtensor::xVector<>(v1->point(), v2->point());
   }

   // duplication ET TRI des noeuds de l'element ayant cree qqch

   for (int j = 0; j < e.size(0); j++)
   {
      int compt = 0;
      mVertex* v = tomvertexp(e.get(0, j));
      if (debug) std::cout << "\n\ncutElement : noeud v de e: " << v->point() << " **********\n\n";
#if 1
      for (int k = 0; k < e_interface.size(0); k++)
      {
         mVertex* vq = tomvertexp(e_interface.get(0, k));
         if (debug) std::cout << "\n\ncutElement : noeud v de e_interface: " << vq->point() << " **********\n\n";
         mEntity* const* p_vqc = c_was_created_by.getData(*vq);
         if ((p_vqc) && (*p_vqc == v))
         {
            compt += 1;
         }
      }
#endif
      if (debug) std::cout << "\ncutElement: valeur de compt = " << compt << "\n\n";
      if (compt == 0)
      {
         if (debug) std::cout << "\ncutElement: le noeud est distinct de son fils compt = 0\n\n";
         double s = xtensor::xVector<>(v1->point(), v->point()) * n;
         if (debug) std::cout << "********** valeur de s = " << s << " **********\n\n";
         if (s > 0.)
         {
            if (debug) std::cout << "cutElement : je duplique le noeud de e dans side1 **********\n\n";
            mVertex& v_new = cloneVertex(*v);
            duplicated_side1_vertices.push_back(&v_new);
         }
         else if (s < 0.)
         {
            if (debug) std::cout << "cutElement : je duplique le noeud de e dans side 2 **********\n\n";
            mVertex& v_new = cloneVertex(*v);
            duplicated_side2_vertices.push_back(&v_new);
         }
      }
      if (debug) std::cout << "\ncutElement: size list1 = " << duplicated_side1_vertices.size() << "\n\n";
      if (debug) std::cout << "\ncutElement: size list2 = " << duplicated_side2_vertices.size() << "\n\n";
   }
   if (duplicated_side1_vertices.size() != 1)
   {
      swap(duplicated_side1_vertices, duplicated_side2_vertices);
      swap(interface_duplicated_vertices[1], interface_duplicated_vertices[interface_duplicated_vertices.size() - 1]);
      if (debug) std::cout << "\n cutElement: listes arrangees size list1 !=1\n\n";
   }
   if (duplicated_side2_vertices.size() == 3)
   {
      xtensor::xVector<> nn = xtensor::xVector<>(duplicated_side2_vertices[1]->point(), duplicated_side2_vertices[0]->point()) %
                              xtensor::xVector<>(duplicated_side2_vertices[2]->point(), duplicated_side2_vertices[1]->point());
      xtensor::xVector<> nbis(duplicated_side2_vertices[0]->point(), duplicated_side1_vertices[0]->point());
      if (nn * nbis < 0.) swap(duplicated_side2_vertices[0], duplicated_side2_vertices[1]);
      if (debug) std::cout << "\n cutElement: listes arrangees size list2 = 3\n\n";
   }
   // fin de duplication triee cdes noeuds de e, en dim 2 et 3

   if (debug) std::cout << "\n cutElement: entre dans la boucle de decoupage effectif en dim >1 \n\n";
   double subvol = 0.0;
   if (duplicated_side1_vertices.size() == 1)
   {
      if (dim == 2)  // triangle ayant engendre un triangle d'un cote construction du ss-triangle
      {
         if (debug) std::cout << "\n cutElement: (dim 2) je construis le 1er sous triangle... \n\n";
         mFace* fac1_new = submmesh.createFaceWithVertices(interface_duplicated_vertices[0], interface_duplicated_vertices[1],
                                                           duplicated_side1_vertices[0], e.getClassification());

         is_in_partition_of.setData(*fac1_new) = const_cast<mEntity*>(&e);
         subvol += xElement(fac1_new).getVolume();
      }
      else if (dim == 3)  // tetraedre ayant engendre un tetraedre d'un cote construction du ss-tet
      {
         if (debug) std::cout << "\n cutElement: (dim 3) je construis le 1er sous_tetraedre...\n\n";

         mTet* tet1_new = submmesh.createTetWithVertices(interface_duplicated_vertices[0], interface_duplicated_vertices[1],
                                                         interface_duplicated_vertices[2], duplicated_side1_vertices[0],
                                                         e.getClassification());
         is_in_partition_of.setData(*tet1_new) = const_cast<mEntity*>(&e);
         xElement elem(tet1_new);
         subvol += elem.getVolume();
      }

      if (duplicated_side2_vertices.size() == 1)
      {
         if (dim == 2)  // triangle decoupe en deux triangles, construction du 2eme
         {
            mFace* fac2_new = submmesh.createFaceWithVertices(duplicated_side2_vertices[0], interface_duplicated_vertices[1],
                                                              interface_duplicated_vertices[0], e.getClassification());

            is_in_partition_of.setData(*fac2_new) = const_cast<mEntity*>(&e);
            subvol += xElement(fac2_new).getVolume();
            if (debug)
               std::cout << "\n cutElement: (dim 2) fin de construction du 2eme triangle en decoupage d'un triangle en 2 "
                            "triangles \n\n";
         }
         if (dim == 3)  // tetraedre decoupe en deux tetraedres, construction du 2eme
         {
            mTet* tet2_new = submmesh.createTetWithVertices(duplicated_side2_vertices[0], interface_duplicated_vertices[2],
                                                            interface_duplicated_vertices[0], interface_duplicated_vertices[1],
                                                            e.getClassification());
            is_in_partition_of.setData(*tet2_new) = const_cast<mEntity*>(&e);
            xElement elem(tet2_new);
            subvol += elem.getVolume();
            if (debug) std::cout << "\n cutElement: (dim 3) fin de construction du 2eme tet en decoupage d'un tet en 2 tets \n\n";
         }
      }  // ferme if (duplicated_side2_vertices.size() == 1)

      else if (duplicated_side2_vertices.size() == 2)
      {
         if (dim == 2)  // triangle decoupe en un triangle et un quadilatere
         {
            if (debug) std::cout << "\n \n cutElement: debut de decoupage du triangle en un triangle et un quad \n \n";
            // decomposition du quadrilatere
            std::vector<mVertex*> quad_vertices(interface_duplicated_vertices.begin(), interface_duplicated_vertices.end());
            quad_vertices.insert(quad_vertices.begin(), duplicated_side2_vertices.begin(), duplicated_side2_vertices.end());
            _orientQuad(quad_vertices);
            mFace* fac1_new =
                submmesh.createFaceWithVertices(quad_vertices[0], quad_vertices[1], quad_vertices[2], e.getClassification());

            is_in_partition_of.setData(*fac1_new) = const_cast<mEntity*>(&e);
            mFace* fac2_new =
                submmesh.createFaceWithVertices(quad_vertices[0], quad_vertices[2], quad_vertices[3], e.getClassification());
            is_in_partition_of.setData(*fac2_new) = const_cast<mEntity*>(&e);
            subvol += xElement(fac1_new).getVolume() + xElement(fac2_new).getVolume();
            if (debug) std::cout << "\n \n cutElement ): fin de decoupage du triangle en un triangle et un quad \n \n";
         }

         if (dim == 3)  // tetraedre decoupe en un tetraedre et une pyramide
         {
            if (debug) std::cout << "\n \n cutElement: debut de decoupage du tet en un tet et 1 pyr \n \n";
            if (debug) std::cout << "duplicated_side2_vertices[0] " << duplicated_side2_vertices[0]->point() << "\n\n";
            if (debug) std::cout << "duplicated_side2_vertices[1] " << duplicated_side2_vertices[1]->point() << "\n\n";
            if (debug) std::cout << "interface_duplicated_vertices[0] " << interface_duplicated_vertices[0]->point() << "\n\n";
            if (debug) std::cout << "interface_duplicated_vertices[1] " << interface_duplicated_vertices[1]->point() << "\n\n";
            if (debug) std::cout << "interface_duplicated_vertices[2] " << interface_duplicated_vertices[2]->point() << "\n\n";

            mTet* tet1_new = submmesh.createTetWithVertices(duplicated_side2_vertices[0], interface_duplicated_vertices[2],
                                                            interface_duplicated_vertices[0], interface_duplicated_vertices[1],
                                                            e.getClassification());
            mVertex* vv1 = vertexInBetween(*duplicated_side2_vertices[0], *duplicated_side1_vertices[0], e_interface_vertices);

            if (debug) std::cout << "vertex in between vv1 " << vv1->point() << "\n\n";

            for (size_t i = 0; i < 3; i++)
            {
               if (vv1 == interface_duplicated_vertices[i])
               {
                  memo = int(i);
                  if (debug) std::cout << "valeur de memo int " << memo << "**********\n\n";
               }
            }
            if (debug) std::cout << "duplicated_side1_vertices[0] " << duplicated_side1_vertices[0]->point() << "\n\n";
            if (debug) std::cout << "valeur de memo " << memo << "**********\n\n";

            if (debug)
               std::cout << "interface_duplicated_vertices[(memo+1)%3] " << interface_duplicated_vertices[(memo + 1) % 3]->point()
                         << "\n\n";
            if (debug)
               std::cout << "interface_duplicated_vertices[(memo+2)%3] " << interface_duplicated_vertices[(memo + 2) % 3]->point()
                         << "\n\n";

            // if(debug) std::cout <<"interface_duplicated_vertices[1] "<<interface_duplicated_vertices[0]->point()<<"\n\n";
            // if(debug) std::cout <<"interface_duplicated_vertices[2] "<<interface_duplicated_vertices[0]->point()<<"\n\n";
            mTet* tet2_new = submmesh.createTetWithVertices(
                duplicated_side2_vertices[1], interface_duplicated_vertices[(memo + 1) % 3],
                interface_duplicated_vertices[(memo + 2) % 3], duplicated_side2_vertices[0], e.getClassification());

            is_in_partition_of.setData(*tet1_new) = const_cast<mEntity*>(&e);
            is_in_partition_of.setData(*tet2_new) = const_cast<mEntity*>(&e);
            xElement elem1(tet1_new);
            xElement elem2(tet2_new);
            subvol += elem1.getVolume() + elem2.getVolume();
            if (debug) std::cout << "\n \n cutElement: fin de decoupage du tet en un tet et 1 pyr \n \n";
         }
      }

      else
      {
         if (debug) std::cout << "\n \n cutElement: entre dans le cas de decoupage de tet en tet et prisme \n \n";

         if (duplicated_side2_vertices.size() == 3)
         {
            mVertex* vv2 = vertexInBetween(*duplicated_side2_vertices[0], *duplicated_side1_vertices[0], e_interface_vertices);
            for (size_t i = 0; i < 3; i++)
            {
               if (vv2 == interface_duplicated_vertices[i]) memo = int(i);
            }
            if (memo < 0) throw;
            if (debug) std::cout << " noeud vv2 = " << vv2->point() << "\n";
            if (debug) std::cout << " noeud ds2[0] = " << duplicated_side2_vertices[0]->point() << "\n";
            if (debug) std::cout << " noeud ds1[0] = " << duplicated_side1_vertices[0]->point() << "\n";
            if (debug) std::cout << " \n\n valeur de memo = " << memo << "\n\n\n";
            if (debug)
               std::cout << " \n\n nbre de neouds dupliqu� ds e_interface = " << interface_duplicated_vertices.size() << "\n\n\n";
            mTet* tet3_new =
                submmesh.createTetWithVertices(duplicated_side2_vertices[0], interface_duplicated_vertices[size_t(memo)],
                                               interface_duplicated_vertices[(memo + 1) % 3],
                                               interface_duplicated_vertices[(memo + 2) % 3], e.getClassification());

            mTet* tet4_new = submmesh.createTetWithVertices(
                duplicated_side2_vertices[0], interface_duplicated_vertices[(memo + 1) % 3], duplicated_side2_vertices[2],
                interface_duplicated_vertices[(memo + 2) % 3], e.getClassification());

            mTet* tet5_new =
                submmesh.createTetWithVertices(duplicated_side2_vertices[0], interface_duplicated_vertices[(memo + 1) % 3],
                                               duplicated_side2_vertices[1], duplicated_side2_vertices[2], e.getClassification());
            is_in_partition_of.setData(*tet3_new) = const_cast<mEntity*>(&e);
            is_in_partition_of.setData(*tet4_new) = const_cast<mEntity*>(&e);
            is_in_partition_of.setData(*tet5_new) = const_cast<mEntity*>(&e);
            xElement elem3(tet3_new);
            xElement elem4(tet4_new);
            xElement elem5(tet4_new);
            subvol += elem3.getVolume() + elem4.getVolume() + elem5.getVolume();
         }
      }
   }

   else /*if (duplicated_side1_vertices.size() == 2)*/  // Tetraedre decoupe en 2 prismes
   {
      if (debug) std::cout << "\n \n cutElement: entre dans le cas de decoupage de tet en 2 prismes \n \n";

      _orientQuad(interface_duplicated_vertices);
      xtensor::xVector<> vec1bis =
          xtensor::xVector<>(interface_duplicated_vertices[0]->point(), interface_duplicated_vertices[1]->point()) %
          xtensor::xVector<>(interface_duplicated_vertices[1]->point(), interface_duplicated_vertices[2]->point());

      if (xtensor::xVector<>(interface_duplicated_vertices[0]->point(), duplicated_side1_vertices[0]->point()) * vec1bis < 0.)
      {
         swap(interface_duplicated_vertices[1], interface_duplicated_vertices[3]);
      }

      mVertex* vv4 = vertexInBetween(*duplicated_side2_vertices[0], *duplicated_side1_vertices[0], e_interface_vertices);

      if (debug) std::cout << " noeud vv4 = " << vv4->point() << "\n";
      if (debug) std::cout << " noeud ds2[0] = " << duplicated_side2_vertices[0]->point() << "\n";
      if (debug) std::cout << " noeud ds1[0] = " << duplicated_side1_vertices[0]->point() << "\n";

      for (size_t i = 0; i < 4; i++)
      {
         if (vv4 == interface_duplicated_vertices[i]) memo = int(i);
      }
      if (debug) std::cout << " \n\n valeur de memo1 = " << memo << "\n\n\n";

      mVertex* vv5 = vertexInBetween(*duplicated_side2_vertices[0], *duplicated_side1_vertices[1], e_interface_vertices);

      if (debug) std::cout << " noeud vv5 = " << vv5->point() << "\n";
      if (debug) std::cout << " noeud ds2[0] = " << duplicated_side2_vertices[0]->point() << "\n";
      if (debug) std::cout << " noeud ds1[1] = " << duplicated_side1_vertices[1]->point() << "\n";

      if (vv5 != interface_duplicated_vertices[(memo + 1) % 4]) memo = (memo + 1) % 4;

      if (debug) std::cout << " \n\n valeur de memo2 = " << memo << "\n\n\n";
      if (debug)
         std::cout << " \n\n nbre de neouds dupliqu� de la interface = " << interface_duplicated_vertices.size() << "\n\n\n";

      // for (int i= 0; i < interface_duplicated_vertices.size(); i++)
      //{std::cout<<" noeud sdv["<<i<<"] = "<<interface_duplicated_vertices[i]->point()<<"\n";}

      mTet* tet6_new = submmesh.createTetWithVertices(duplicated_side2_vertices[0], interface_duplicated_vertices[size_t(memo)],
                                                      interface_duplicated_vertices[(memo + 1) % 4],
                                                      interface_duplicated_vertices[(memo + 2) % 4], e.getClassification());

      mTet* tet7_new = submmesh.createTetWithVertices(duplicated_side2_vertices[0], interface_duplicated_vertices[size_t(memo)],
                                                      interface_duplicated_vertices[(memo + 2) % 4],
                                                      interface_duplicated_vertices[(memo + 3) % 4], e.getClassification());
      is_in_partition_of.setData(*tet6_new) = const_cast<mEntity*>(&e);
      is_in_partition_of.setData(*tet7_new) = const_cast<mEntity*>(&e);

      xtensor::xVector<> vec4 =
          xtensor::xVector<>(duplicated_side2_vertices[0]->point(), interface_duplicated_vertices[(memo + 3) % 4]->point()) %
          xtensor::xVector<>(interface_duplicated_vertices[(memo + 3) % 4]->point(),
                             interface_duplicated_vertices[(memo + 2) % 4]->point());
      if (xtensor::xVector<>(duplicated_side2_vertices[0]->point(), duplicated_side2_vertices[1]->point()) * vec4 < 0.)
      {
         mTet* tet8_new =
             submmesh.createTetWithVertices(duplicated_side2_vertices[1], interface_duplicated_vertices[size_t(memo)],
                                            interface_duplicated_vertices[(memo + 1) % 4],

                                            duplicated_side2_vertices[0], e.getClassification());
         is_in_partition_of.setData(*tet8_new) = const_cast<mEntity*>(&e);
         subvol = xElement(tet8_new).getVolume();
      }
      else
      {
         mTet* tet9_new = submmesh.createTetWithVertices(
             duplicated_side2_vertices[0], interface_duplicated_vertices[(memo + 3) % 4],
             interface_duplicated_vertices[(memo + 2) % 4], duplicated_side2_vertices[1], e.getClassification());
         is_in_partition_of.setData(*tet9_new) = const_cast<mEntity*>(&e);
         subvol = xElement(tet9_new).getVolume();
      }
      subvol += xElement(tet6_new).getVolume() + xElement(tet7_new).getVolume();
      swap(interface_duplicated_vertices[1], interface_duplicated_vertices[3]);
      mVertex* vv6 = vertexInBetween(*duplicated_side1_vertices[0], *duplicated_side2_vertices[0], e_interface_vertices);
      for (size_t i = 0; i < 4; i++)
      {
         if (vv6 == interface_duplicated_vertices[i]) memo = int(i);
      }
      mVertex* vv7 = vertexInBetween(*duplicated_side1_vertices[0], *duplicated_side2_vertices[1], e_interface_vertices);
      if (vv7 != interface_duplicated_vertices[(memo + 1) % 4]) memo = (memo + 1) % 4;
      mTet* tet10_new = submmesh.createTetWithVertices(duplicated_side1_vertices[0], interface_duplicated_vertices[size_t(memo)],
                                                       interface_duplicated_vertices[(memo + 1) % 4],
                                                       interface_duplicated_vertices[(memo + 2) % 4], e.getClassification());
      is_in_partition_of.setData(*tet10_new) = const_cast<mEntity*>(&e);

      mTet* tet11_new = submmesh.createTetWithVertices(duplicated_side1_vertices[0], interface_duplicated_vertices[size_t(memo)],
                                                       interface_duplicated_vertices[(memo + 2) % 4],
                                                       interface_duplicated_vertices[(memo + 3) % 4], e.getClassification());
      is_in_partition_of.setData(*tet11_new) = const_cast<mEntity*>(&e);
      xtensor::xVector<> vec8 =
          xtensor::xVector<>(duplicated_side1_vertices[0]->point(), interface_duplicated_vertices[(memo + 3) % 4]->point()) %
          xtensor::xVector<>(interface_duplicated_vertices[(memo + 3) % 4]->point(),
                             interface_duplicated_vertices[(memo + 2) % 4]->point());
      if (xtensor::xVector<>(duplicated_side1_vertices[0]->point(), duplicated_side1_vertices[1]->point()) * vec8 < 0.)
      {
         mTet* tet12_new =
             submmesh.createTetWithVertices(duplicated_side1_vertices[1], interface_duplicated_vertices[size_t(memo)],
                                            interface_duplicated_vertices[(memo + 1) % 4],

                                            duplicated_side1_vertices[0], e.getClassification());
         is_in_partition_of.setData(*tet12_new) = const_cast<mEntity*>(&e);
         subvol += xElement(tet12_new).getVolume();
      }
      else
      {
         mTet* tet13_new = submmesh.createTetWithVertices(
             duplicated_side1_vertices[0], interface_duplicated_vertices[(memo + 3) % 4],
             interface_duplicated_vertices[(memo + 2) % 4], duplicated_side1_vertices[1], e.getClassification());
         is_in_partition_of.setData(*tet13_new) = const_cast<mEntity*>(&e);
         subvol += xElement(tet13_new).getVolume();
      }
      subvol += xElement(tet10_new).getVolume() + xElement(tet11_new).getVolume();

      double vol = xElement(&e).getVolume();
      assert(fabs(vol - subvol) / vol <= 1.e-10);

      if (debug) std::cout << "cutElement: somme volumes elementaires = " << subvol << "\n\n";
      if (debug) std::cout << "cutElement: volume de l'element = " << vol << "\n\n";

      if (debug) std::cout << "\n \n cutElement: sort du decoupage de tet en 2 prismes \n \n";
   }

   xinterface::aomd::modifyAllState(submmesh);

   if (debug) std::cout << "cutElement : SORTIE\n\n";
}

void cutAlongInterface(const xMesh& interface, const datamanager_t<mEntity*>& was_created_by, datamanager_t<xMesh>& partition,
                       datamanager_t<mEntity*>& is_duplicated_in, datamanager_t<mEntity*>& was_duplicated_from,
                       datamanager_t<mEntity*>& is_in_partition_of)
{
   const bool debug = xdebug_flag;
   std::unordered_set<mEntity*, EntityHashKey, EntityEqualKey> already_cut;

   for (int d = 0; d <= interface.dim(); ++d)
   {
      for (mEntity* e_interf : interface.range(d))
      {
         mEntity* e = was_created_by.at(*e_interf);
         // we need to cut the element e if all condition below are fullfilled
         //  e has not already been cut
         //  e->getLevel() = e_interf->getLevel() + 1
         //  e has no partition
         if ((e->getLevel() == e_interf->getLevel() + 1) && (!partition.getData(*e)) && !already_cut.count(e))
         {
            xMesh& submesh = partition.setData(*e);
            _cutElement(*e, *e_interf, was_created_by, submesh, is_duplicated_in, was_duplicated_from, is_in_partition_of);
            already_cut.insert(e);
            if (debug) AOMD_Util::Instance()->ex_port("submesh.msh", &submesh.getMesh());
         }
      }
   }
   return;
}

void cutAlongInterfaceRecursive(const xMesh& interface, const xLevelSet& ls, const datamanager_t<mEntity*>& was_created_by,
                                datamanager_t<xMesh>& partition, datamanager_t<mEntity*>& is_duplicated_in,
                                datamanager_t<mEntity*>& was_duplicated_from, datamanager_t<mEntity*>& is_in_partition_of)
{
   const bool debug = xdebug_flag;
   std::unordered_set<mEntity*, EntityHashKey, EntityEqualKey> already_cut;

   for (int d = 0; d <= interface.dim(); ++d)
   {
      if (debug) std::cout << " d in cutAlongInterfaceRecursive is " << d << std::endl;
      for (mEntity* e_interf : interface.range(d))
      {
         mEntity* const* pe = was_created_by.getData(*e_interf);
         if (!pe) continue;  // added by greg (Hex/Quad)
         mEntity* e = *pe;
         if (debug)
         {
            std::cout << "e_interf is " << std::endl;
            e_interf->print();
            std::cout << "e_interf was created by " << std::endl;
            e->print();
         }

         if ((e->getLevel() == e_interf->getLevel() + 1) && !already_cut.count(e))
         {
            xfem::xMesh* am = partition.getData(*e);
            if (!am && ((e->getType() == mEntity::QUAD) || (e->getType() == mEntity::HEX)))
            {
               xMesh& submesh = partition.setData(*e);
               auto cloneVertex = [&submesh, &is_duplicated_in, &was_duplicated_from](const mVertex& v) -> mVertex& {
                  mVertex& vnew = *submesh.getMesh().createVertex(v.point(), v.getClassification());
                  is_duplicated_in.setData(v) = &vnew;
                  was_duplicated_from.setData(vnew) = const_cast<mVertex*>(&v);
                  return vnew;
               };
               _addSimplicies(*e, submesh, cloneVertex);
               am = &submesh;
            }
            if (!am)
            {
               if (debug) std::cout << " no existing partition " << std::endl;
               xMesh& submesh = partition.setData(*e);
               _cutElement(*e, *e_interf, was_created_by, submesh, is_duplicated_in, was_duplicated_from, is_in_partition_of);
               already_cut.insert(e);
               if (debug) AOMD_Util::Instance()->ex_port("submesh.msh", &submesh.getMesh());
            }
            else  // initial partition exists
            {
               if (debug) std::cout << "    existing partition " << std::endl;
               xMesh& submesh = *am;
               if (debug) AOMD_Util::Instance()->ex_port("submesh_of_the_existing_parition.msh", &submesh.getMesh());
               xLevelSet ls_submesh;
               // const datamanagerxMesh_t<AOMD::mEntity*> &was_duplicated_from,
               //        const datamanagerxMesh_t<AOMD::mEntity*> &was_created_by,
               ls.interpolateTo(submesh, was_duplicated_from, was_created_by, ls_submesh);
               xMesh interface_on_submesh;
               datamanager_t<double> r_on_edge_sub;
               datamanager_t<mEntity*> was_created_by_sub;
               createIsoZeroMeshFromLevelSet(ls_submesh, interface_on_submesh, was_created_by_sub, r_on_edge_sub);
               r_on_edge_sub.clear();
               if (debug)
               {
                  std::cout << "interface_on_submesh details " << std::endl;
                  interface_on_submesh.getMesh().printAll();
                  AOMD_Util::Instance()->ex_port("interface_on_submesh.msh", &interface_on_submesh.getMesh());
                  for (mEntity* e_interfe : interface_on_submesh.range(1))
                  {
                     std::cout << "address of e_interf before (d=1) : " << e_interfe << std::endl;
                     mEntity* const* ee = was_created_by.getData(*e_interf);
                     if (!ee)
                     {
                        std::cout << " l'edge e_interfe n'a pas de createur" << std::endl;
                        e_interfe->print();
                     }
                     assert(ee);
                  }
                  for (mEntity* e_interfe : interface_on_submesh.range(2))
                  {
                     std::cout << "address of e_interf before (d=2) : " << e_interfe << std::endl;
                     std::cout << "nb nodes on it : " << e_interfe->size(0) << std::endl;
                     mEntity* const* ee = was_created_by.getData(*e_interfe);
                     assert(ee);
                     std::cout << " e_interfe->getAttachedEntity(was_created_by_tag) " << ee << std::endl;
                  }
               }
               cutAlongInterface(interface_on_submesh, was_created_by_sub, partition, is_duplicated_in, was_duplicated_from,
                                 is_in_partition_of);
               already_cut.insert(e);
               was_created_by_sub.clear();
            }
         }
      }
   }
   return;
}

// The mesh resulting from the level set partition is always composed of simplex
// elements : triangles in 2D and edges in 1D.
// Each mesh subentity has a creator whose level is the same, one more or two more.
// The case "two more" happens when a tetrahedra generates for the iso-zero two
// triangles. The common edge of the two triangles has the tetrahedron as creator.
void _createPartitionFromOneLevelSet(const xMesh& basemesh, const xLevelSet& ls, const xMesh& interface,
                                     const datamanager_t<mEntity*>& was_created_by, xEntityToEntity classify_in,
                                     xEntityToEntity classify_out, xEntityToEntity classify_interface,
                                     datamanager_t<mEntity*>& is_duplicated_in, datamanager_t<mEntity*>& was_duplicated_from,
                                     datamanager_t<xMesh>& partition, datamanager_t<mEntity*>& is_in_partition_of, bool recursive)
{
   const bool debug = xdebug_flag;

   if (debug) std::cout << "createPartitionFromOneLevelSet : ENTREE \n\n";
   if (debug) AOMD_Util::Instance()->ex_port("interface.msh", &interface.getMesh());
   int dim_interface = interface.dim();
   if (!recursive)
      cutAlongInterface(interface, was_created_by, partition, is_duplicated_in, was_duplicated_from, is_in_partition_of);
   else
      cutAlongInterfaceRecursive(interface, ls, was_created_by, partition, is_duplicated_in, was_duplicated_from,
                                 is_in_partition_of);
   int dim_ambient = basemesh.dim();
   if (debug)
   {
      std::cout << "basicMesh is " << std::endl;
      basemesh.getMesh().printAll();
      for (mEntity* e : basemesh.range(2))
      {
         e->print();
         std::vector<double> vals = ls.getVals(e);
         std::cout << "ls " << std::endl << "vals 0 = " << vals[0] << "\n";
         std::cout << "vals 1 = " << vals[1] << "\n";
         std::cout << "vals 2 = " << vals[2] << "\n";
      }
   }
   for (int d = 1; d <= dim_ambient; ++d)
   {
      for (mEntity* egros : basemesh.range(d))
      {
         xPartition partition;
         xMesh::getPartition(egros, partition);
         if (debug)
         {
            std::cout << "taille de la partition : " << partition.size() << std::endl;
            std::cout << "in create partitionfromone levelset, dim is : " << d << std::endl;
         }
         for (mEntity* epetit : partition)
         {
            if (debug)
            {
               epetit->print();
               egros->print();
            }
            _classifyElement(egros, epetit, ls, was_created_by, classify_in, classify_out);
         }
      }
   }
   for (mEntity* egros : interface.range(dim_interface))
   {
      xPartition partition;
      xMesh::getPartition(egros, partition);
      for (mEntity* epetit : partition)
      {
         classify_interface(epetit);
         if (debug) std::cout << "classifie lelement courant sur linterface\n";
      }
   }
   if (debug) std::cout << "createPartitionFromInterface : SORTIE\n";
   return;
}
// loop over the nodes
//   if ls == 0.0 create a node and give creator
// loop over the edges
//   if all vertices cut,      create edge and give creator
//   if edge cut,              create node and give creator
// loop over face
//   if all vertices cut,            create face and give creator
//   if two edges cut
//   or
//   if one edge and one vertex cut, create edge from the 2 points and give creator
// loop over regions
//   if #edge cut and #vertex cut >=3 create face with 3 or 4 points and give creator
// remark an vertex is "is_the_creator_of" if ls = 0.0
//       an edge   is "is_the_creator_of" if ls = 0.0 strictly in between it's vertices
// passer le filtre
// note filter is not used !!!! should be removed.
void cutMesh(const xMesh& base_mesh, const xLevelSet& ls, xMesh& iso_zero_mesh, datamanager_t<mEntity*>& was_created_by,
             datamanager_t<double>& r_on_edge, datamanager_t<mEntity*>& is_duplicated_in,
             datamanager_t<mEntity*>& was_duplicated_from, datamanager_t<xMesh>& partition,
             datamanager_t<mEntity*>& is_in_partition_of, xEntityToEntity classify_in, xEntityToEntity classify_out,
             bool create_partition, bool keep_old_partition, bool recursive)
{
   xcut::createIsoZeroMeshFromLevelSet(ls, iso_zero_mesh, was_created_by, r_on_edge);
   if (!keep_old_partition)
      for (int i = base_mesh.dim(); i >= 0; i--)
         for (mEntity* e : base_mesh.range(i)) partition.deleteData(*e);
   if (create_partition)
      xcut::_createPartitionFromOneLevelSet(base_mesh, ls, iso_zero_mesh, was_created_by,             //
                                            classify_in, classify_out, xtool::xIdentity<mEntity*>(),  //
                                            is_duplicated_in, was_duplicated_from,                    //
                                            partition, is_in_partition_of,                            //
                                            recursive);
   return;
}
void cutMesh(const xfem::xMesh& base_mesh, const xfem::xLevelSet& ls, xfem::xMesh& iso_zero_mesh,  //
             xfem::xEntityToEntity classify_in, xfem::xEntityToEntity classify_out,                //
             bool create_partition, bool keep_old_partition, bool recursive)
{
   datamanager_t<mEntity*>& was_created_by = xMesh::get_was_created_by();
   datamanager_t<double>& r_on_edge = xMesh::get_r_on_edge();
   datamanager_t<mEntity*>& is_duplicated_in = xMesh::get_is_duplicated_in();
   datamanager_t<mEntity*>& was_duplicated_from = xMesh::get_was_duplicated_from();
   datamanager_t<xMesh>& partition = xMesh::get_partition();
   datamanager_t<mEntity*>& is_in_partition_of = xMesh::get_is_in_partition_of();
   cutMesh(base_mesh, ls, iso_zero_mesh, was_created_by, r_on_edge,  //
           is_duplicated_in, was_duplicated_from,                    //
           partition, is_in_partition_of,                            //
           classify_in, classify_out, create_partition, keep_old_partition, recursive);
}

}  // end namespace xcut
