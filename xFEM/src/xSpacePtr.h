/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _SPACE_POINTER_H
#define _SPACE_POINTER_H

#include <memory>

namespace xfem
{

class xSpace;
class xSpaceBase;

using spacePtr_t = std::shared_ptr < xSpace > ;
using spaceBasePtr_t = std::shared_ptr < xSpaceBase > ;

}

#endif
