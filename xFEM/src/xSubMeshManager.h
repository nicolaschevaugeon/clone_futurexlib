/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xSubMeshManager_H_
#define _xSubMeshManager_H_

#include <map>
#include <memory>
#include <string>

#include "mAOMD.h"
#include "mEntity.h"
#include "mIterator.h"
#include "mMesh.h"
#include "mVertex.h"
#include "xMesh.h"

namespace xfem
{
// this class should in the near future replace
// all the sub capabilities contained in xMesh.

class xSubMeshManager
{
  public:
   ~xSubMeshManager();
   void allocateSubsetEntities(const string& sub);
   // to_code void createSubsetEntities(const string& subset, xSubsetCreator& creator);
   void removeSubsetEntities(const string& subset);
   void printSubsetEntities() const;
   void removeAllSubsets();

   AOMD::mEntity* find_sub(AOMD::mEntity*, const string& sub) const;
   void add_sub(AOMD::mEntity*, const string& subset);
   void del_sub(AOMD::mEntity*, const string& subset);
   xIter begin_sub(int what, const string& sub) const;
   xIter end_sub(int what, const string& sub) const;
   xtool::xRange<xIter> range_sub(int what, const string& sub) const
   {
      return xtool::make_range<xIter>(begin_sub(what, sub), end_sub(what, sub));
   }
   xIterall beginall_sub(int what, const string& sub) const;
   xIterall endall_sub(int what, const string& sub) const;
   int size_sub(int what, const string& sub) const;
   int dim_sub(const string& sub) const;

   // to_code (if needed) void modifyState_sub(int from, int to, bool state, const string& name = "all",int with = 0);

  private:
   std::map<string, AOMD::mMeshEntityContainer*> subsetEntities;
   typedef std::map<string, AOMD::mMeshEntityContainer*>::const_iterator const_iter_subsets;
   typedef std::map<string, AOMD::mMeshEntityContainer*>::iterator iter_subsets;
   AOMD::mMeshEntityContainer* getEntities(const string& sub);
   const AOMD::mMeshEntityContainer* getEntities(const string& sub) const;
};

}  // namespace xfem

#endif
