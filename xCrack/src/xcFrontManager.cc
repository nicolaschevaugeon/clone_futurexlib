/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
// std
#include <iostream>
// xcrack
#include "xcFrontManager.h"
// AOMD
#include "AOMD.h"
#include "mVertex.h"
// xfem
#include "ParUtil.h"
#include "xGatherMesh.h"
#include "xMesh.h"
#include "xParseData.h"
#include "xSubMesh.h"

#ifdef PARALLEL
#include "autopack.h"
#endif

using std::cout;
using std::endl;

using namespace AOMD;
using namespace xfem;

xcFrontManager::xcFrontManager(const xMesh &_front_mesh, const string &_frontname, const xMesh &_mesh,
                               std::function<double(mVertex *)> f, const xParseData &_parameters)
    : front_mesh(_front_mesh), front_name(_frontname), mesh(_mesh), front_distance(f), parameters(_parameters)
{
   int verbose = parameters.getInt("front_verbosity");
   if (verbose)
   {
      std::string filename = "frontmanagerlog.txt";
      std::ofstream out(filename.c_str(), std::ios_base::app);
      out << "Entering xcFrontManager constructor" << endl;
      out.close();
   }

   createFronts();

   if (verbose)
   {
      std::string filename = "frontmanagerlog.txt";
      std::ofstream out(filename.c_str(), std::ios_base::app);
      out << "------- nb front parts = " << fronts.size() << endl;
      int i = 1;
      for (xcFrontPartBase *fp : *this)
      {
         out << " part " << i++ << " : nb nodes = " << fp->local_vertices.size() << endl;
      }
      out.close();
   }
}

xcFrontManager::~xcFrontManager()
{
   iterator it = fronts.begin();
   iterator itend = fronts.end();
   while (it != itend)
   {
      ((*front_mesh_global)).deleteSubMesh((*it)->getFrontName());
      delete (*it);
      (*it) = nullptr;
      ++it;
   }

   fronts.clear();
   delete front_mesh_global;
   front_mesh_global = nullptr;
}

double xcFrontManager::getTotalFrontLength() const
{
   double length = 0.;
   for (xcFrontPartBase *fp : *this) length += fp->getFrontLenght();
   return length;
}

void xcFrontManager::createFronts()
{
   int verbose = parameters.getInt("front_verbosity");
   xMesh &front_mesh_mod = (const_cast<xMesh &>(front_mesh));
   // front_mesh_mod.getMesh()->modifyAllState();
   AOMD::mMesh &mmesh = front_mesh_mod.getMesh();
   mmesh.modifyState(3, 2, true);
   mmesh.modifyState(3, 1, true);
   mmesh.modifyState(3, 0, true);
   mmesh.modifyState(2, 1, true);
   mmesh.modifyState(2, 0, true);
   mmesh.modifyState(1, 0, true);
   mmesh.modifyState(0, 1, true);
   mmesh.modifyState(0, 2, true);
   mmesh.modifyState(0, 3, true);
   mmesh.modifyState(1, 2, true);
   mmesh.modifyState(1, 3, true);
   mmesh.modifyState(2, 3, true);
   front_mesh_global = new xMesh;
   if (verbose) cout << "broadcasting_front_mesh" << endl;
   xEntityCopyCallBackAttachEntity<xMesh::datamanager_t> callback(xMesh::get_const_was_created_by(), xMesh::get_was_created_by());
   BroadCastMesh(&front_mesh_mod, front_mesh_global, &callback);

   if (verbose) cout << "separting front_mesh in connected parts" << endl;
   std::map<std::string, std::pair<mVertex *, mVertex *>> frontParts = Separate1dMeshBranch();
   std::map<std::string, std::pair<mVertex *, mVertex *>>::iterator frontPartsIt = frontParts.begin();
   std::map<std::string, std::pair<mVertex *, mVertex *>>::iterator frontPartsItEnd = frontParts.end();
   if (verbose) cout << " xcInteractionIntegralsOnCrack::createFront : nb front " << frontParts.size() << endl;
   if (verbose) cout << "creating xcFronts" << endl;
   while (frontPartsIt != frontPartsItEnd)
   {
      xcFrontPartBase::frontType fronttype;
      std::string frontsubname = frontPartsIt->first;
      // int dim = front_mesh_global->dim();// should be the size of the submesh describing the frontpart... to separate points
      // from lineclosed
      xSubMesh &temp = front_mesh_global->getSubMesh(frontsubname);
      int dim = temp.dim();
      if (dim == 0)
         fronttype = xcFrontPartBase::Point;
      else if (dim == 1)
      {
         if (frontPartsIt->second.first == frontPartsIt->second.second)
            fronttype = xcFrontPartBase::LineClose;
         else
            fronttype = xcFrontPartBase::LineOpen;
      }
      else
         abort();
      mVertex *start = frontPartsIt->second.first;
      mVertex *end = frontPartsIt->second.second;
      switch (fronttype)
      {
         case xcFrontPartBase::Point:
         {
            fronts.push_back(new xcFrontPartPoint(*front_mesh_global, frontsubname, mesh, front_distance, parameters));
            break;
         }
         case xcFrontPartBase::LineOpen:
         {
            fronts.push_back(
                new xcFrontPartLineOpen(*front_mesh_global, frontsubname, mesh, front_distance, parameters, start, end));
            break;
         }
         case xcFrontPartBase::LineClose:
         {
            fronts.push_back(new xcFrontPartLineClose(*front_mesh_global, frontsubname, mesh, front_distance, parameters, start));
            break;
         }
         case xcFrontPartBase::None:
         {
            cout << "warning can't create  front part with type xcFrontPartBase::None" << endl;
            break;
         }
      }
      ++frontPartsIt;
   }
   if (verbose) cout << "done creating front parts" << endl;
};

std::map<std::string, std::pair<mVertex *, mVertex *>> xcFrontManager::Separate1dMeshBranch()
{
   // NB: when considering multiple frontparts merging, a node might be connected to more than 2 edges
   // infinte loops are then observed with the algorithm below.
   // Quick Fix used: slightly shift the level set before calling xcFront...
   const int verbose = parameters.getInt("front_verbosity");
   std::string subsetnamebase(front_name + "_part_");
   std::map<std::string, std::pair<mVertex *, mVertex *>> parts;
   if (verbose) cout << front_mesh_global->dim() << endl;
   switch (front_mesh_global->dim())
   {
      case 0:
      {
         int partnumber = 1;
         for (mEntity *pv : front_mesh_global->range(0))
         {
            std::stringstream subsetname;
            subsetname << subsetnamebase << partnumber;
            xSubMesh &subset = front_mesh_global->createSubMesh(subsetname.str());
            subset.add(pv);
            parts.insert(std::make_pair(subsetname.str(), std::make_pair((AOMD::mVertex *)(pv), (AOMD::mVertex *)nullptr)));
            ++partnumber;
         }
         break;
      }
      case 1:
      {
         std::list<mVertex *> totreatv;
         front_mesh_global->getMesh().modifyState(0, 1, true);
         for (mEntity *pv : front_mesh_global->range(0)) totreatv.push_back(static_cast<mVertex *>(pv));
         int partnumber = 0;
         xinterface::aomd::xAttachedDataManagerAOMD<int> visited;
         if (verbose) cout << "Separate 1d mesh Branch start " << endl;
         while (totreatv.size() != 0)
         {  // continue untill every node has been stored in a front part
            ++partnumber;
            std::stringstream subsetname;
            subsetname << subsetnamebase << partnumber;
            xSubMesh &subset = front_mesh_global->createSubMesh(subsetname.str());
            // cout << "totreat " <<totreatv.size() << endl;
            std::list<mVertex *>::iterator itv = totreatv.begin();
            mVertex *starttmp = (*itv);
            mVertex *start = nullptr;
            mVertex *theend = nullptr;
            mVertex *current = starttmp;
            AOMD::mEdge *previousedge = nullptr;
            bool isloop = false;
            subset.add(starttmp);
            visited.setData(*starttmp) = 1;
            if (starttmp->size(1) == 1)
            {
               if (verbose) cout << "Separate 1d mesh Branch start on boundary" << endl;
               start = starttmp;
               AOMD::mEdge *currentedge = (AOMD::mEdge *)current->get(1, 0);
               if (currentedge == previousedge)
               {
                  currentedge = (AOMD::mEdge *)current->get(1, 1);
               }
               current = (mVertex *)E_otherVertex((mEntity *)currentedge, current);
               subset.add(current);
               subset.add((mEntity *)currentedge);
               visited.setData(*current) = 1;
               previousedge = currentedge;
               while (current->size(1) > 1)
               {
                  AOMD::mEdge *currentedge = (AOMD::mEdge *)current->get(1, 0);
                  if (currentedge == previousedge)
                  {
                     currentedge = (AOMD::mEdge *)current->get(1, 1);
                  }
                  current = (mVertex *)E_otherVertex((mEntity *)currentedge, current);

                  subset.add(current);
                  subset.add((mEntity *)currentedge);
                  visited.setData(*current) = 1;
                  previousedge = currentedge;
               }
               theend = current;
            }
            else
            {
               if (verbose) cout << "Separate 1d mesh Branch start on a middle" << endl;
               while (current->size(1) > 1)
               {
                  AOMD::mEdge *currentedge = (AOMD::mEdge *)current->get(1, 0);
                  if (currentedge == previousedge)
                  {
                     currentedge = (AOMD::mEdge *)current->get(1, 1);
                  }
                  subset.add((mEntity *)currentedge);
                  current = (mVertex *)E_otherVertex((mEntity *)currentedge, current);
                  if (current == starttmp) break;
                  subset.add(current);
                  visited.setData(*current) = 1;
                  previousedge = currentedge;
               }
               if (current == starttmp)
               {
                  start = starttmp;
                  theend = starttmp;
               }
               else
               {
                  if (verbose) cout << "Separate 1d mesh Branch going the other way" << endl;
                  theend = current;
                  previousedge = nullptr;
                  current = starttmp;
                  while (current->size(1) > 1)
                  {
                     AOMD::mEdge *currentedge = (AOMD::mEdge *)current->get(1, 1);
                     if (currentedge == previousedge)
                     {
                        currentedge = (AOMD::mEdge *)current->get(1, 0);
                     }
                     current = (mVertex *)E_otherVertex((mEntity *)currentedge, current);
                     visited.setData(*current) = 1;
                     subset.add(current);
                     subset.add((mEntity *)currentedge);
                     previousedge = currentedge;
                  }
                  start = current;
               }
            }
            isloop = (theend == start);
            if (verbose) cout << "Separate 1d mesh Branch find branch " << start << " " << theend << " " << isloop << endl;
            parts.insert(std::make_pair(subsetname.str(), std::make_pair(start, theend)));
            totreatv.erase(std::remove_if(totreatv.begin(), totreatv.end(),
                                          [&visited](AOMD::mVertex *pv) { return bool(visited.getData(*pv)); }),
                           totreatv.end());
            subset.modifyState(0, 1, true);
         }
         for (mEntity *pe : front_mesh_global->range(0)) visited.deleteData(*pe);
         break;
      }
      default:
      {
         cout << "Error : front_mesh_global should be of dimension 0 or 1 " << endl;
         throw;
      }
   }
   return parts;
};
