/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xcFrontManager_
#define _xcFrontManager_

#include <functional>
#include <string>
#include "xcFrontPart.h"

class xParseData;

namespace xfem
{
class xMesh;
}

class xcFrontManager
{
  public:
   xcFrontManager(const xfem::xMesh& front_mesh, const std::string& _frontname, const xfem::xMesh& mesh,
                  std::function<double(AOMD::mVertex*)> front_distance, const xParseData& _parameters);
   ~xcFrontManager();
   typedef std::list<xcFrontPartBase*>::iterator iterator;
   typedef std::list<xcFrontPartBase*>::const_iterator const_iterator;
   int size() const { return fronts.size(); };
   iterator begin() { return fronts.begin(); };
   iterator end() { return fronts.end(); };
   const_iterator begin() const { return fronts.begin(); };
   const_iterator end() const { return fronts.end(); };

   const xfem::xMesh& getFrontMesh() const { return *front_mesh_global; };
   const xfem::xMesh& getMesh() const { return mesh; };
   double getTotalFrontLength() const;

  private:
   void createFronts();
   std::map<std::string, std::pair<AOMD::mVertex*, AOMD::mVertex*>> Separate1dMeshBranch();

   const xfem::xMesh& front_mesh;
   xfem::xMesh* front_mesh_global;
   const std::string front_name;
   const xfem::xMesh& mesh;
   std::function<double(AOMD::mVertex*)> front_distance;
   const xParseData& parameters;
   std::list<xcFrontPartBase*> fronts;
};

#endif
