/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xcFrontDomainManager.h"
#include "xMesh.h"
#include "xSubMeshManager.h"

using namespace AOMD;
using namespace std;
using namespace xfem;

xcFrontDomainManager::xcFrontDomainManager(xcFrontManager *_fmanager, int physical_process)
    : fmanager(_fmanager), globalSubMesh(nullptr)
{
   // creationg a domain for each front part
   for (xcFrontPartBase *fp : *fmanager)
   {
      switch (fp->getFrontType())
      {
         case xcFrontPartBase::Point:
            domains.push_back(new xcFrontDomainPoint(fp, physical_process));
            break;
         case xcFrontPartBase::LineOpen:
            domains.push_back(new xcFrontDomainLineOpen(fp, physical_process));
            break;
         case xcFrontPartBase::LineClose:
            domains.push_back(new xcFrontDomainLineClosed(fp, physical_process));
            break;
         case xcFrontPartBase::None:
            cout << "warning can't create domain for front part with type xcFrontPartBase::None" << endl;
            break;
      }
   }
};

xSubMesh &xcFrontDomainManager::createGlobalSubMesh()
{
   if (globalSubMesh) return *globalSubMesh;
   const xMesh &mesh(fmanager->getMesh());
   globalSubMesh = &(mesh.createSubMesh("globalsubmesh"));
   iterator it = begin_domains_iter();
   iterator itend = end_domains_iter();
   while (it != itend)
   {
      // const xcFrontPartBase *currentPart = (*it)->getFrontPart();
      const string subset_name = (*it)->subset_for_int_label;
      xSubMesh *mysubmesh = (const_cast<xSubMesh *>(&mesh.getSubMesh(subset_name)));
      for (xIter itf = mysubmesh->begin(2); itf != mysubmesh->end(2); ++itf)
      {
         mEntity *e = (mEntity *)*itf;
         globalSubMesh->add(e);
         // seules les faces sont ajoutées par defaut... manque les noeuds...
         int nb_vertices_on_face = e->size(0);
         for (int inode = 0; inode < nb_vertices_on_face; inode++)
         {
            globalSubMesh->add(e->get(0, inode));
         }
      }
      it++;
   }
   return *globalSubMesh;
}

void xcFrontDomainManager::extendParametrization()
{
   createGlobalSubMesh();
   iterator it = begin_domains_iter();
   iterator itend = end_domains_iter();
   while (it != itend)
   {
      (*it)->extendParametrization();
      it++;
   }
}
