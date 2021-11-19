/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef ___XCHARCONTAINER_H
#define ___XCHARCONTAINER_H

#include "mAttachableDataContainer.h"
#include "mEntity.h"

namespace xinterface
{
namespace aomd
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xAttachableChar class
//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class xAttachableChar : public AOMD::mAttachableData
{
  public:
   ~xAttachableChar() override = default;
   // public methodes ///////////////
   // public members ////////////////
   char c;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xAttachableChar class
//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO ---
// TODO --- TODO --- TODO see if it's the right place for these functions

// get attached data of type xAttachableChar
char getAttachedChar(const AOMD::mEntity* e, unsigned int tag);
inline char getAttachedChar(const AOMD::mEntity* e, unsigned int tag)
{
   xAttachableChar* ac = (xAttachableChar*)(e->getData(tag));
   if (!ac) return 0;
   return ac->c;
}
// set attached data of type xAttachableChar
void attachChar(AOMD::mEntity* e, unsigned int tag, char c);
inline void attachChar(AOMD::mEntity* e, unsigned int tag, char c)
{
   xAttachableChar* ac = (xAttachableChar*)(e->getData(tag));
   if (!ac)
   {
      ac = new xAttachableChar;
      e->attachData(tag, ac);
   }
   ac->c = c;
}
// TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO ---
// TODO --- TODO --- TODO

}  // namespace aomd
}  // namespace xinterface

#endif
