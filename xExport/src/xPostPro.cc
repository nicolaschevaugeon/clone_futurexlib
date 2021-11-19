/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xPostPro.h"

#include "xMesh.h"

using namespace std;
using namespace AOMD;

namespace xexport
{
xPostPro::xPostPro(xfem::xMesh *m)
    : begin(m->begin(m->dim())), end(m->end(m->dim())), filter(xfem::xAcceptAll()), lsn(nullptr), lst(nullptr)
{
   AOMD::mMesh &mm = m->getMesh();
   // Fill the data structure according to the export in Trellis
   // Assumptions :
   // * contiguously numbered nodes
   // * Simplex only !

   // Check if node numbering begins at zero or one
   int nodesOffset{0};
   if (!mm.getVertex(0)) nodesOffset = 1;

   // can't use mesh iterators as they won't be ordered
   for (int i = 0 + nodesOffset; i < m->size(0) + nodesOffset; ++i)
   {
      mVertex *v = mm.getVertex(i);
      nodesList.push_back(v);
      approNodeList.push_back(v->get(mm.getDim(), 0));
      integNodeList.push_back(v->get(mm.getDim(), 0));
   }

   // Elements
   for (xfem::xIter it = begin; it != end; ++it)
   {
      mEntity *e = *it;

      std::vector<int> elementConnectivity = {e->get(0, 0)->getId() - nodesOffset, e->get(0, 1)->getId() - nodesOffset,
                                              e->get(0, 2)->getId() - nodesOffset};

      if (mm.getDim() == 3) elementConnectivity.push_back(e->get(0, 3)->getId() - nodesOffset);

      element_type.push_back(e->getType());
      connectivity.push_back(elementConnectivity);
      approElemList.push_back(e);
      integElemList.push_back(e);
   }
}

xPostPro::xPostPro(const xfem::xIter &_begin, const xfem::xIter &_end, const xfem::xEntityFilter _filter,
                   const xfem::xLevelSet *_lsn, const xfem::xLevelSet *_lst)
    : begin(_begin), end(_end), filter(_filter), lsn(_lsn), lst(_lst)
{
   if ((lsn == nullptr) && (lst == nullptr))
   {
      mergeMesh();
   }
   else if ((lsn != nullptr) && (lst == nullptr))
   {
      mergeMesh(lsn);
   }
   else if ((lsn != nullptr) && (lst != nullptr))
   {
      mergeMesh(lsn, lst);
   }
}

void xPostPro::mergeMesh()
{
   unsigned int new_id_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("new_id_tag");
   for (xfem::xIter itElem = begin; itElem != end; ++itElem)
   {
      mEntity *e = *itElem;
      if (filter(e))
      {
         xfem::xPartition partition;
         xfem::xMesh::getPartition(e, partition);
         // Loop on subelements
         for (std::set<mEntity *>::iterator itSubElem = partition.begin(); itSubElem != partition.end(); ++itSubElem)
         {
            mEntity *se = *itSubElem;
            if (filter(se))
            {
               std::vector<int> elementConnectivity;
               // Loop on nodes
               for (int n = 0; n < se->size(0); ++n)
               {
                  mVertex *v = (mVertex *)se->get(0, n);
                  v = (mVertex *)xfem::getSource(v);
                  if (!v->getAttachedInt(new_id_tag))
                  {
                     // renumber the node if not already renumbered
                     nodesList.push_back(v);
                     v->attachInt(new_id_tag, nodesList.size());
                     approNodeList.push_back(e);
                     integNodeList.push_back(se);
                  }
                  elementConnectivity.push_back(v->getAttachedInt(new_id_tag) - 1);
               }
               element_type.push_back(e->getType());
               connectivity.push_back(elementConnectivity);
               approElemList.push_back(e);
               integElemList.push_back(se);
            }
         }
      }
   }
   for (vector<mVertex *>::iterator it = nodesList.begin(); it != nodesList.end(); ++it)
   {
      (*it)->deleteData(new_id_tag);
   }
}

void xPostPro::mergeMesh(const xfem::xLevelSet *lsn)
{
   unsigned int new_id_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("new_id_tag");
   unsigned int side_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
   xfem::xEvalLevelSet<xtool::xIdentity<double>> eval_lsn(*lsn);
   double lsnv;
   // Loop on elements
   for (xfem::xIter itElem = begin; itElem != end; ++itElem)
   {
      mEntity *e = *itElem;
      if (filter(e))
      {
         xfem::xPartition partition;
         xfem::xMesh::getPartition(e, partition);
         // Loop on subelements
         for (std::set<mEntity *>::iterator itSubElem = partition.begin(); itSubElem != partition.end(); ++itSubElem)
         {
            mEntity *se = *itSubElem;
            if (filter(se))
            {
               std::vector<int> elementConnectivity;
               xfem::xGeomElem geo_appro(e);
               xfem::xGeomElem geo_integ(se);
               // Loop on nodes
               for (int n = 0; n < se->size(0); ++n)
               {
                  mVertex *v = (mVertex *)se->get(0, n);
                  v = (mVertex *)xfem::getSource(v);
                  geo_appro.setUVWForXYZ(v->point());
                  geo_integ.setUVWForXYZ(v->point());
                  eval_lsn(&geo_appro, &geo_integ, lsnv);
                  double precision = 1e-12;
                  bool lsnNull = (lsnv > -precision) && (lsnv < precision);
                  if (!lsnNull)
                  {
                     ////////////////// Node not on crack  //////////////////
                     if (!v->getAttachedInt(new_id_tag))
                     {
                        nodesList.push_back(v);
                        v->attachInt(new_id_tag, nodesList.size());
                        approNodeList.push_back(e);
                        integNodeList.push_back(se);
                     }
                     elementConnectivity.push_back(v->getAttachedInt(new_id_tag) - 1);
                  }
                  else
                  {
                     //////////////////// Node on crack  ////////////////////
                     int side = lsn->side_of(&geo_appro, &geo_integ);
                     if (!v->getAttachedInt(new_id_tag))
                     {  // nodes not in nodesList
                        nodesList.push_back(v);
                        nodesList.push_back(v);
                        approNodeList.push_back(e);
                        approNodeList.push_back(e);
                        integNodeList.push_back(se);
                        integNodeList.push_back(se);
                        int id = 0;
                        if (side == -1)
                        {
                           id = nodesList.size() - 1;
                        }
                        else if (side == 1)
                        {
                           id = nodesList.size();
                        }
                        v->attachInt(new_id_tag, id);
                        v->attachInt(side_tag, side);
                        elementConnectivity.push_back(v->getAttachedInt(new_id_tag) - 1);
                     }
                     else
                     {  // root node in list
                        int id;
                        if (side == v->getAttachedInt(side_tag))
                        {  // actual node and root node on the same side
                           id = v->getAttachedInt(new_id_tag);
                        }
                        else
                        {  // actual node and root node on opposite side
                           id = v->getAttachedInt(new_id_tag);
                           id = v->getAttachedInt(new_id_tag);
                           id += side;
                           approNodeList[id - 1] = e;
                           integNodeList[id - 1] = se;
                        }
                        elementConnectivity.push_back(id - 1);
                     }
                  }
               }
               element_type.push_back(e->getType());
               connectivity.push_back(elementConnectivity);
               approElemList.push_back(e);
               integElemList.push_back(se);
            }
         }
      }
   }
   for (vector<mVertex *>::iterator it = nodesList.begin(); it != nodesList.end(); ++it)
   {
      (*it)->deleteData(new_id_tag);
      if ((*it)->getAttachedInt(side_tag)) (*it)->deleteData(side_tag);
   }
}

void xPostPro::mergeMesh(const xfem::xLevelSet *lsn, const xfem::xLevelSet *lst)
{
   unsigned int new_id_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("new_id_tag");
   unsigned int side_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
   xfem::xEvalLevelSet<xtool::xIdentity<double>> eval_lsn(*lsn);
   xfem::xEvalLevelSet<xtool::xIdentity<double>> eval_lst(*lst);
   double lsnv;
   double lstv;

   // Loop on elements
   for (xfem::xIter itElem = begin; itElem != end; ++itElem)
   {
      mEntity *e = *itElem;
      if (filter(e))
      {
         xfem::xPartition partition;
         xfem::xMesh::getPartition(e, partition);
         // Loop on subelements
         for (std::set<mEntity *>::iterator itSubElem = partition.begin(); itSubElem != partition.end(); ++itSubElem)
         {
            mEntity *se = *itSubElem;
            if (filter(se))
            {
               std::vector<int> elementConnectivity;
               xfem::xGeomElem geo_appro(e);
               xfem::xGeomElem geo_integ(se);
               // Loop on nodes
               for (int n = 0; n < se->size(0); ++n)
               {
                  mVertex *v = (mVertex *)se->get(0, n);
                  v = (mVertex *)xfem::getSource(v);
                  geo_appro.setUVWForXYZ(v->point());
                  geo_integ.setUVWForXYZ(v->point());
                  eval_lsn(&geo_appro, &geo_integ, lsnv);
                  eval_lst(&geo_appro, &geo_integ, lstv);
                  double precision = 1e-12;
                  bool lsnNull = (lsnv > -precision) && (lsnv < precision);
                  bool lstNeg = (lstv < -precision);
                  if ((!lsnNull) || (!lstNeg))
                  {
                     ////////////////// Node not on crack  //////////////////
                     if (!v->getAttachedInt(new_id_tag))
                     {
                        nodesList.push_back(v);
                        v->attachInt(new_id_tag, nodesList.size());
                        approNodeList.push_back(e);
                        integNodeList.push_back(se);
                     }
                     elementConnectivity.push_back(v->getAttachedInt(new_id_tag) - 1);
                  }
                  else
                  {
                     //////////////////// Node on crack  ////////////////////
                     int side = lsn->side_of(&geo_appro, &geo_integ);
                     if (!v->getAttachedInt(new_id_tag))
                     {  // nodes not in nodesList
                        nodesList.push_back(v);
                        nodesList.push_back(v);
                        approNodeList.push_back(e);
                        approNodeList.push_back(e);
                        integNodeList.push_back(se);
                        integNodeList.push_back(se);
                        int id = 0;
                        if (side == -1)
                        {
                           id = nodesList.size() - 1;
                        }
                        else if (side == 1)
                        {
                           id = nodesList.size();
                        }
                        v->attachInt(new_id_tag, id);
                        v->attachInt(side_tag, side);
                        elementConnectivity.push_back(v->getAttachedInt(new_id_tag) - 1);
                     }
                     else
                     {  // root node in list
                        int id;
                        if (side == v->getAttachedInt(side_tag))
                        {  // actual node and root node on the same side
                           id = v->getAttachedInt(new_id_tag);
                        }
                        else
                        {  // actual node and root node on opposite side
                           id = v->getAttachedInt(new_id_tag);
                           id = v->getAttachedInt(new_id_tag);
                           id += side;
                           approNodeList[id - 1] = e;
                           integNodeList[id - 1] = se;
                        }
                        elementConnectivity.push_back(id - 1);
                     }
                  }
               }
               element_type.push_back(e->getType());
               connectivity.push_back(elementConnectivity);
               approElemList.push_back(e);
               integElemList.push_back(se);
            }
         }
      }
   }
   for (vector<mVertex *>::iterator it = nodesList.begin(); it != nodesList.end(); ++it)
   {
      (*it)->deleteData(new_id_tag);
      if ((*it)->getAttachedInt(side_tag)) (*it)->deleteData(side_tag);
   }
}

}  // namespace xexport
