/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _XTENSORS_POINTER_H
#define _XTENSORS_POINTER_H

#include <memory>

namespace xfem
{

class xTensors;

typedef std::shared_ptr < xTensors > tensorsPtr_t;

}

#endif
