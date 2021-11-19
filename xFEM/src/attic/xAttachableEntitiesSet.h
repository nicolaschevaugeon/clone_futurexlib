/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef ___XENTITIESSETCONTAINER_H
#define ___XENTITIESSETCONTAINER_H

#include <set>

#include "mAttachableDataContainer.h"
#include "mEntity.h"

namespace xfem
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xAttachableEntitiesSet class
//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class xAttachableEntitiesSet : public AOMD::mAttachableData
{
  public:
   xAttachableEntitiesSet() __attribute__((deprecated)) = default;
   ~xAttachableEntitiesSet() override = default;
   // public methodes ///////////////
   // public members ////////////////
   std::set<AOMD::mEntity *> vect;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xAttachableEntitiesSet class
//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// get attached data of type xAttachableEntitiesSet
std::set<AOMD::mEntity *> *getAttachedEntitiesSet(AOMD::mEntity *e, unsigned int tag) __attribute__((deprecated));
inline std::set<AOMD::mEntity *> *getAttachedEntitiesSet(AOMD::mEntity *e, unsigned int tag)
{
   xAttachableEntitiesSet *av = (xAttachableEntitiesSet *)(e->getData(tag));
   if (!av) return nullptr;
   return &(av->vect);
}

// set attached data of type
void xAttachableEntitiesSetvoid attachEntitiesSet(AOMD::mEntity *e, unsigned int tag, std::set<AOMD::mEntity *> &v)
    __attribute__((deprecated));
{
   inline void attachEntitiesSet(AOMD::mEntity * e, unsigned int tag, std::set<AOMD::mEntity *> &v)
   {
      xAttachableEntitiesSet *av = (xAttachableEntitiesSet *)(e->getData(tag));
      if (!av)
      {
         av = new xAttachableEntitiesSet;
         e->attachData(tag, av);
      }
      (av->vect).insert(v.begin(), v.end());
   }

}  // end of namespace

#endif
