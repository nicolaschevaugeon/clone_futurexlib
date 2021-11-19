/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xcFrontDomainManager_
#define _xcFrontDomainManager_

#include "xcFrontDomain.h"
#include "xcFrontManager.h"

namespace xfem
{
class xIntegrationRule;
class xSpaceRegular;
class xSubMesh;
}  // namespace xfem

class xcFrontDomainManager
{
  public:
   xcFrontDomainManager(xcFrontManager *_fmanager, int physical_process);
   ~xcFrontDomainManager() { domains.clear(); }

   xfem::xSubMesh &createGlobalSubMesh();
   void extendParametrization();

   typedef std::list<xcFrontDomain *>::iterator iterator;
   typedef std::list<xcFrontDomain *>::const_iterator const_iterator;
   size_t size_domains_parts() const { return domains.size(); }
   iterator begin_domains_iter() { return domains.begin(); }
   iterator end_domains_iter() { return domains.end(); }
   const_iterator begin_domains_iter() const { return domains.begin(); }
   const_iterator end_domains_iter() const { return domains.end(); }
   xcFrontManager *fmanager;
   /// exports the modes in gmsh format with a filename (basefilename)_part_(frontpart)_mod_(modenumber)_iter_(iter).pos

  private:
   std::list<xcFrontDomain *> domains;
   xfem::xSubMesh *globalSubMesh;
};

#endif
