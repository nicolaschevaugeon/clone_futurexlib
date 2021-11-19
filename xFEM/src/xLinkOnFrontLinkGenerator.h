#ifndef _LINKONFRONTLINKGENERATOR__H
#define _LINKONFRONTLINKGENERATOR__H

#include "xEntityFilter.h"
#include "xEntityToEntity.h"
#include "xMesh.h"
#include "xVector.h"

namespace xfem
{
// the following class selects node whose dofs are to be linked around a front
// on the front support according to the paper entitled
// "A stable Lagrange multiplier space for stiff interface conditions within the extended finite element method", IJNME,
// Bechet, Moes and Wohlmuth DOI: 10.1002/nme.2515

class xLinkOnFrontLinkGenerator
{
   //    typedef xfem::xValueManagerDist<double>  xValueManagerDist<double>;
   // typedef map<mEdge*,double> orthoinfo;// the higher is the double, the more orthogonal to the front is the edge... double in
   // [0,1]

   // the following map indicates which nodes will be linked
   // if the second AOMD::mVertex* is not NULL, the two vertices must be linked (linear combination state for the first vertex)
   // if not, the first vertex has to be declared as new dof ("double" state)
   typedef map<AOMD::mVertex*, AOMD::mVertex*> creationinfo;

  public:
   // front is given by iterators on edges (1D front) or faces (2D front)
   // this constructor is more generic since one can filter the front
   template <typename ITERFRONT>
   xLinkOnFrontLinkGenerator(ITERFRONT it, ITERFRONT end, bool link_isolated = true, xEntityFilter filter = xAcceptAll(),
                             xEntityToEntity interf2appro = xCreator());

   // front is given by the mesh
   xLinkOnFrontLinkGenerator(const xMesh* m, bool link_isolated = true, xEntityFilter filt = xAcceptAll(),
                             xEntityToEntity interf2appro = xCreator());
   ~xLinkOnFrontLinkGenerator();
   // gives iterators on vertices around the front
   // should use iterators directly on the map, instead of creating a copy in a vector... ?
   vector<AOMD::mVertex*>::iterator beginVertexIter();
   vector<AOMD::mVertex*>::iterator endVertexIter();

   AOMD::mVertex* getLinkedVertex(AOMD::mVertex* v) const;

  private:
   void buildVertexVector();
   template <typename ITERFRONT>
   void buildNormalAtNodes(ITERFRONT it, ITERFRONT end);
   void buildTab(bool link_isolated, xEntityFilter filter, xEntityToEntity interf2appro);
   // double compute_orthogonality(mEdge *e, xtensor::xVector<> vector1, AOMD::mVertex *v);
   // orthoinfo orthotab;
   vector<AOMD::mVertex*> vertexvector;
   creationinfo tab;
   typedef std::map<AOMD::mEntity*, xtensor::xVector<>> normal_at_nodes_t;
   typedef normal_at_nodes_t::iterator normal_at_nodes_it_t;
   normal_at_nodes_t normal_at_nodes;

   // private class to sort AOMD::mEntity viewed as AOMD::mVertex with criteria on edge.normale value and coordinate acending
   // order
   class nodeSortingCriteria
   {
     public:
      nodeSortingCriteria(normal_at_nodes_t* norms_);
      bool operator()(AOMD::mEntity* ent1, AOMD::mEntity* ent2) const;

     private:
      normal_at_nodes_t* norms;
   };

   typedef set<AOMD::mEntity*, nodeSortingCriteria> node_set_t;

   // normal computation
   // nota :
   //    normForElem2D is dedicated to edge
   //    normForElem3D is dedicated to triangle
   //
   std::function<void(AOMD::mEntity*, xtensor::xVector<>&)> normForElem;
   void normForElem2D(AOMD::mEntity* elem, xtensor::xVector<>& norm);
   void normForElem3D(AOMD::mEntity* elem, xtensor::xVector<>& norm);
};

//----------------------------------------------------------------------------------------------------------------

template <typename ITERFRONT>
xLinkOnFrontLinkGenerator::xLinkOnFrontLinkGenerator(ITERFRONT it, ITERFRONT end, bool link_isolated, xEntityFilter filter,
                                                     xEntityToEntity interf2appro)
{
   if (it == end) return;

   // nota :
   //  sorting is base on 2 criteria :
   //    * node of front having related edge localy perpendicular to front mesh are treated first
   //    * node of front having same edge orthogonality criteria are sorted by theire coordinates
   //
   buildNormalAtNodes(it, end);
   buildTab(link_isolated, filter, interf2appro);
}

template <typename ITERFRONT>
void xLinkOnFrontLinkGenerator::buildNormalAtNodes(ITERFRONT it, ITERFRONT end)
{
   // choose the right normal compute function
   if ((*it)->getType() == AOMD::mEntity::TRI)
      normForElem = std::bind(&xLinkOnFrontLinkGenerator::normForElem3D, this, std::placeholders::_1, std::placeholders::_2);
   else if ((*it)->getType() == AOMD::mEntity::EDGE)
      normForElem = std::bind(&xLinkOnFrontLinkGenerator::normForElem2D, this, std::placeholders::_1, std::placeholders::_2);
   else
      throw;

   std::pair<normal_at_nodes_it_t, bool> found;
   xtensor::xVector<> norm;
   const double zero = 0.;

   // loop on elements of front to compute and accumultate normal at nodes
   for (; it != end; ++it)
   {
      AOMD::mEntity* elem = *it;
      normForElem(elem, norm);
      AOMD::mAdjacencyContainer::iter itn = elem->begin(0);
      AOMD::mAdjacencyContainer::iter endn = elem->end(0);
      for (; itn != endn; ++itn)
      {
         found = normal_at_nodes.insert(make_pair((*itn), norm));
         if (!found.second)
         {
            if (((*found.first).second * norm) < zero)
               (*found.first).second -= norm;
            else
               (*found.first).second += norm;
         }
      }
   }
}
}

#endif
