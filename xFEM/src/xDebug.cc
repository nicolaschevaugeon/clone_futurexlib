/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xDebug.h"

namespace xfem {
void SetXfemDebugFlag(bool b) {xDebugSingleton::instance().setFlag(b); }
}

