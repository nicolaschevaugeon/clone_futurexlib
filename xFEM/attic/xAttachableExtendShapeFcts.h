/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef ___XATTACHABLEEXTENDSHAPEFCTS_H
#define ___XATTACHABLEEXTENDSHAPEFCTS_H

#include "mAttachableDataContainer.h"
#include "mEntity.h"

namespace xfem
{
class xExtendShapeFcts;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xAttachableExtendShapeFcts class
///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class xAttachableExtendShapeFcts : public AOMD::mAttachableData
{
  public:
   ~xAttachableExtendShapeFcts() override = default;
   // public methodes ///////////////
   // public members ////////////////
   xExtendShapeFcts *m;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xAttachableExtendShapeFcts class
//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO ---
// TODO --- TODO --- TODO see if it's the right place for these functions

// get attached data of type  xAttachableExtendShapeFcts
inline xExtendShapeFcts *getAttachedExtendShapeFcts(AOMD::mEntity *e, unsigned int tag)
{
   xAttachableExtendShapeFcts *ac = (xAttachableExtendShapeFcts *)(e->getData(tag));
   if (!ac) return nullptr;
   return ac->m;
};

// set attached data of type xAttachableExtendShapeFcts
inline void attachExtendShapeFcts(AOMD::mEntity *e, unsigned int tag, xExtendShapeFcts *m)
{
   xAttachableExtendShapeFcts *ac = (xAttachableExtendShapeFcts *)(e->getData(tag));
   if (!ac)
   {
      ac = new xAttachableExtendShapeFcts;
      e->attachData(tag, ac);
   }
   ac->m = m;
};

// nullify attached data pointer of type xAttachableExtendShapeFcts if any
inline void nullifyAttachedExtendShapeFcts(AOMD::mEntity *e, unsigned int tag)
{
   xAttachableExtendShapeFcts *ac = (xAttachableExtendShapeFcts *)(e->getData(tag));
   if (ac) ac->m = nullptr;
};

// TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO ---
// TODO --- TODO --- TODO

}  // namespace xfem

#endif
