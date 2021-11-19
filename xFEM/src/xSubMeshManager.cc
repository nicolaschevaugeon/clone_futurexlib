/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xSubMeshManager.h"
#include <cstdio>

using std::cerr;

namespace xfem
{
using AOMD::mEntity;
using AOMD::mMeshEntityContainer;

xSubMeshManager::~xSubMeshManager()
{
   for (iter_subsets it = subsetEntities.begin(); it != subsetEntities.end(); ++it)
   {
      delete it->second;
   }
}

void xSubMeshManager::del_sub(mEntity* e, const string& subset)
{
   mMeshEntityContainer* m = getEntities(subset);
   m->del(e);
}

void xSubMeshManager::add_sub(mEntity* e, const string& subset)
{
   mMeshEntityContainer* m = getEntities(subset);
   m->add(e);
}

// mEntity* xSubMeshManager::find(mEntity *e) { return mesh->find(e);}

mMeshEntityContainer* xSubMeshManager::getEntities(const string& sub)
{
   iter_subsets it = subsetEntities.find(sub);
   if (it == subsetEntities.end())
   {
      cerr << "The entities set with name " << sub << " does not exist\n";
   }
   return it->second;
}
const mMeshEntityContainer* xSubMeshManager::getEntities(const string& sub) const
{
   const_iter_subsets it = subsetEntities.find(sub);
   if (it == subsetEntities.end())
   {
      cerr << "The entities set with name " << sub << " does not exist\n";
   }
   return it->second;
}

void xSubMeshManager::allocateSubsetEntities(const string& sub)
{
   removeSubsetEntities(sub);
   subsetEntities[sub] = new mMeshEntityContainer;
}

void xSubMeshManager::printSubsetEntities() const
{
   printf("The subset entities are:\n");
   for (const auto& name_pentitycont : subsetEntities)
   {
      std::string name = name_pentitycont.first;
      AOMD::mMeshEntityContainer* pcont = name_pentitycont.second;
      printf("%s : %d %d %d %d\n", name.c_str(), pcont->size(0), pcont->size(1), pcont->size(2), pcont->size(3));
      printf("ids of the vertices\n");
      for (mEntity* pe : range_sub(0, name)) pe->print();
      printf("ids of the edges\n");
      for (mEntity* pe : range_sub(1, name)) pe->print();
      printf("ids of the faces\n");
      for (mEntity* pe : range_sub(2, name)) pe->print();
      printf("ids of the volumes\n");
      for (mEntity* pe : range_sub(3, name)) pe->print();
   }
   return;
}

void xSubMeshManager::removeSubsetEntities(const string& sub)
{
   iter_subsets it = subsetEntities.find(sub);
   if (it != subsetEntities.end())
   {
      delete it->second;
      subsetEntities.erase(it);
   }
}

void xSubMeshManager::removeAllSubsets()
{
   for (iter_subsets it = subsetEntities.begin(); it != subsetEntities.end(); ++it)
   {
      delete it->second;
   }
   subsetEntities.clear();
}

mEntity* xSubMeshManager::find_sub(mEntity* e, const string& sub) const
{
   const mMeshEntityContainer* m = getEntities(sub);
   return m->find(e);
}

int xSubMeshManager::dim_sub(const string& sub) const
{
   const mMeshEntityContainer* m = getEntities(sub);
   return (m->size(3)) ? 3 : ((m->size(2)) ? 2 : ((m->size(1)) ? 1 : 0));
}

// xIter
xIter xSubMeshManager::begin_sub(int what, const string& sub) const
{
   const mMeshEntityContainer* m = getEntities(sub);
   return AOMD::mLeavesIterator(m->begin(what), m->end(what));
}
// xIter
xIter xSubMeshManager::end_sub(int what, const string& sub) const
{
   const mMeshEntityContainer* m = getEntities(sub);
   return AOMD::mLeavesIterator(m->end(what), m->end(what));
}

// xIter
xIterall xSubMeshManager::beginall_sub(int what, const string& sub) const
{
   const mMeshEntityContainer* m = getEntities(sub);
   return m->begin(what);
   //  return beginall(what);
}
// xIter
xIterall xSubMeshManager::endall_sub(int what, const string& sub) const
{
   const mMeshEntityContainer* m = getEntities(sub);
   return m->end(what);
   //  return endall(what);
}

int xSubMeshManager::size_sub(int what, const string& sub) const
{
   const mMeshEntityContainer* m = getEntities(sub);
   return m->size(what);
}

// void xSubMeshManager::modifyState_sub (int i , int j , bool state, const string& subset, int with)
// {
//   if(!state)
//   {
//     std::for_each (beginall_sub(i,subset), endall_sub(i,subset), deleteAdjFunctor (j));
//   }
//   else
//   if(i>j)
//   {
//     std::for_each (beginall_sub(i,subset), endall_sub(i,subset), xCreateDownwardFunctor_sub (j,with,this,subset));
//   }
//   else
//   if (i<j)
//   {
//     std::for_each (begin_sub(j,subset), end_sub(j,subset), createUpwardFunctor (i));
//   }
// }

}  // namespace xfem
