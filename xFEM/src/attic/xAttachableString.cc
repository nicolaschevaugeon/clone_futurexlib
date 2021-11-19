/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xAttachableString.h"

// AOMD
#include <mEntity.h>

namespace xfem
{
// get attached data of type xAttachableString
std::string &getAttachedString(AOMD::mEntity *e, unsigned int tag)
{
   xAttachableString *a = (xAttachableString *)(e->getData(tag));
   if (!a)
   {
      a = new xAttachableString;
      e->attachData(tag, a);
   }
   return a->s;
}
// set attached data of type xAttachableString
void attachString(AOMD::mEntity *e, unsigned int tag, const std::string &s)
{
   xAttachableString *a = (xAttachableString *)(e->getData(tag));
   if (!a)
   {
      a = new xAttachableString;
      e->attachData(tag, a);
   }
   a->s = s;
   return;
}

}  // namespace xfem
