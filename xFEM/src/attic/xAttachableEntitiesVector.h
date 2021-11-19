/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef ___XENTITIESVECTORCONTAINER_H
#define ___XENTITIESVECTORCONTAINER_H

#include <vector>

#include "mAttachableDataContainer.h"
#include "mEntity.h"

namespace xfem
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xAttachableEntitiesVector class
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class xAttachableEntitiesVector : public AOMD::mAttachableData
{
  public:
   xAttachableEntitiesVector() __attribute__((deprecated)) = default;
   ~xAttachableEntitiesVector() override = default;
   // public methodes ///////////////
   // public members ////////////////
   std::vector<AOMD::mEntity *> vect;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xAttachableEntitiesVector class
///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// get attached data of type xAttachableEntitiesVector
std::vector<AOMD::mEntity *> *getAttachedEntitiesVector(AOMD::mEntity *e, unsigned int tag) __attribute__((deprecated));
inline std::vector<AOMD::mEntity *> *getAttachedEntitiesVector(AOMD::mEntity *e, unsigned int tag)
{
   xAttachableEntitiesVector *av = (xAttachableEntitiesVector *)(e->getData(tag));
   if (!av) return nullptr;
   return &(av->vect);
}

// set attached data of type  xAttachableEntitiesVector
void attachEntitiesVector(AOMD::mEntity *e, unsigned int tag, std::vector<AOMD::mEntity *> &v) __attribute__((deprecated));

inline void attachEntitiesVector(AOMD::mEntity *e, unsigned int tag, std::vector<AOMD::mEntity *> &v)
{
   xAttachableEntitiesVector *av = (xAttachableEntitiesVector *)(e->getData(tag));
   if (!av)
   {
      av = new xAttachableEntitiesVector;
      e->attachData(tag, av);
   }
   (av->vect).resize(v.size());
   std::copy(v.begin(), v.end(), (av->vect).begin());
}

}  // namespace xfem

#endif
