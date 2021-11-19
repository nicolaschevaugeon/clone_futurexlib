#include "xLinkOnFrontLinkGenerator.h"

using namespace AOMD;

namespace xfem
{
xLinkOnFrontLinkGenerator::~xLinkOnFrontLinkGenerator() { vertexvector.clear(); }
void xLinkOnFrontLinkGenerator::buildVertexVector()
{
   creationinfo::iterator it = tab.begin();
   // cout << "building vector vextex " << endl;
   for (; it != tab.end(); it++)
   {
      vertexvector.push_back(it->first);
      /*	cout << it->first->getId() << " " << it->first->point() << " ";
         if (it->second)
                  cout << it->second->getId() << " " << it->second->point() << endl;
         else
          cout << "NULL"<<endl;
       */
   }
   return;
}

vector<AOMD::mVertex *>::iterator xLinkOnFrontLinkGenerator::beginVertexIter()
{
   if (!vertexvector.size()) buildVertexVector();
   return vertexvector.begin();
}

vector<AOMD::mVertex *>::iterator xLinkOnFrontLinkGenerator::endVertexIter()
{
   if (!vertexvector.size()) buildVertexVector();
   return vertexvector.end();
}

AOMD::mVertex *xLinkOnFrontLinkGenerator::getLinkedVertex(AOMD::mVertex *v) const
{
   creationinfo::const_iterator it = tab.find(v);
   if (it == tab.end())
      throw;
   else
      return it->second;
}

xLinkOnFrontLinkGenerator::nodeSortingCriteria::nodeSortingCriteria(normal_at_nodes_t *norms_) : norms(norms_) {}
void xLinkOnFrontLinkGenerator::normForElem2D(AOMD::mEntity *e, xtensor::xVector<> &norm)
{
   xtensor::xVector<> edge_vect((static_cast<AOMD::mVertex *>(e->get(0, 0)))->point(),
                                (static_cast<AOMD::mVertex *>(e->get(0, 1)))->point());
   xtensor::xVector<> plan_norm(0., 0., 1.);
   norm = plan_norm % edge_vect;
   norm.norm();
   return;
}

void xLinkOnFrontLinkGenerator::normForElem3D(AOMD::mEntity *e, xtensor::xVector<> &norm)
{
   xtensor::xPoint p0 = (static_cast<AOMD::mVertex *>(e->get(0, 0)))->point();
   xtensor::xVector<> edge_1(p0, (static_cast<AOMD::mVertex *>(e->get(0, 1)))->point());
   xtensor::xVector<> edge_2(p0, (static_cast<AOMD::mVertex *>(e->get(0, 2)))->point());
   norm = edge_1 % edge_2;
   norm.norm();
   return;
}

bool xLinkOnFrontLinkGenerator::nodeSortingCriteria::operator()(AOMD::mEntity *ent1, AOMD::mEntity *ent2) const
{
   normal_at_nodes_it_t f1 = norms->find(ent1);
   assert(f1 != norms->end());
   normal_at_nodes_it_t f2 = norms->find(ent2);
   assert(f2 != norms->end());

   if ((*f1).second(0) > (*f2).second(0)) return true;
   if ((*f1).second(0) < (*f2).second(0)) return false;
   return ((AOMD::mVertex *)ent1)->point().lexicographicLessThan(((AOMD::mVertex *)ent2)->point(), 0.);
}

xLinkOnFrontLinkGenerator::xLinkOnFrontLinkGenerator(const xMesh *m, bool link_isolated, xEntityFilter filter,
                                                     xEntityToEntity interf2appro)
{
   if (m->size(0) == 0) return;

   // nota :
   //  sorting is base on 2 criteria :
   //    * node of front having related edge localy perpendicular to front mesh are treated first
   //    * node of front having same edge orthogonality criteria are sorted by theire coordinates
   //
   buildNormalAtNodes(m->begin(m->dim()), m->end(m->dim()));
   buildTab(link_isolated, filter, interf2appro);
}

void xLinkOnFrontLinkGenerator::buildTab(bool link_isolated, xEntityFilter filter, xEntityToEntity interf2appro)
{
   const double zero = 0.;
   double prod;

   // second, compute edge normal product at the nodes of front
   // and store it according to sorting criteria
   // =========================================================

   // generate container to store ordered nodes
   nodeSortingCriteria sc(&normal_at_nodes);
   node_set_t front_nodes(sc);

   // loop on container to avoid finding operation
   normal_at_nodes_it_t itvn = normal_at_nodes.begin();
   normal_at_nodes_it_t itvne = normal_at_nodes.end();
   for (; itvn != itvne; ++itvn)
   {
      AOMD::mEntity *node = (*itvn).first;
      AOMD::mEntity *e = interf2appro(node);
      if (e && filter(e))
      {
         (*itvn).second.norm();
         if (e->getLevel() == 1)
         {
            xtensor::xVector<> edge_vect((static_cast<AOMD::mVertex *>(e->get(0, 0)))->point(),
                                         (static_cast<AOMD::mVertex *>(e->get(0, 1)))->point());
            edge_vect.norm();
            prod = ((*itvn).second) * edge_vect;
            prod = fabs(prod);
         }
         else
            prod = zero;

         // to minimize storage consumption use of xtensor::xVector<> first component to store product
         (*itvn).second(0) = prod;

         // this is inserting and sorting
         front_nodes.insert(node);
      }
   }

   // map is no more used => clear it
   normal_at_nodes.clear();

   // set of nodes already visited
   AOMD::mMeshEntityContainer visited_nodes;

   // third loop to select the vital edge
   // ====================================
   typename node_set_t::iterator beg_nodes = front_nodes.begin();
   typename node_set_t::iterator end_nodes = front_nodes.end();
   typename node_set_t::iterator it_nodes;
   for (it_nodes = beg_nodes; it_nodes != end_nodes; ++it_nodes)
   {
      AOMD::mEntity *e_bnd = *it_nodes;
      AOMD::mEntity *e = interf2appro(e_bnd);

      if (e->getLevel() == 0)
      {
         // regular dof will be created
         tab.insert(pair<AOMD::mVertex *, AOMD::mVertex *>((AOMD::mVertex *)e, (AOMD::mVertex *)nullptr));
         visited_nodes.add(e);
      }
      else if (e->getLevel() == 1)
      {
         // check if 0, 1 or 2 nodes have dofs
         // if 0, will have to create two dofs and link them
         AOMD::mEntity *v1 = e->get(0, 0);
         AOMD::mEntity *v2 = e->get(0, 1);
         if (!visited_nodes.find(v1) && !visited_nodes.find(v2))
         {
            tab.insert(pair<AOMD::mVertex *, AOMD::mVertex *>((AOMD::mVertex *)v1, (AOMD::mVertex *)nullptr));
            tab.insert(pair<AOMD::mVertex *, AOMD::mVertex *>((AOMD::mVertex *)v2, (AOMD::mVertex *)v1));
            visited_nodes.add(v1);
            visited_nodes.add(v2);
         }
      }
   }
   // fourth loop for remaining vertices
   // ==================================
   for (it_nodes = beg_nodes; it_nodes != end_nodes; ++it_nodes)
   {
      AOMD::mEntity *e_bnd = *it_nodes;
      AOMD::mEntity *e = interf2appro(e_bnd);
      if (e->getLevel() == 1)
      {
         // check if 0, 1 or 2 nodes have dofs
         // if 1, create the other dof and link it
         AOMD::mEntity *v1 = e->get(0, 0);
         AOMD::mEntity *v2 = e->get(0, 1);
         bool b1 = (visited_nodes.find(v1));
         bool b2 = (visited_nodes.find(v2));
         assert(b1 || b2);
         if (!b1 || !b2)
         {
            if (b2) swap(v1, v2);
            if (link_isolated)
               tab.insert(pair<AOMD::mVertex *, AOMD::mVertex *>((AOMD::mVertex *)v2, (AOMD::mVertex *)v1));
            else
               tab.insert(pair<AOMD::mVertex *, AOMD::mVertex *>((AOMD::mVertex *)v2, (AOMD::mVertex *)nullptr));
            visited_nodes.add(v2);
         }
      }
   }
}

}  // namespace xfem
