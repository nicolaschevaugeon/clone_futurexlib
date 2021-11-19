/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __VARIAB__MANAGER__H
#define __VARIAB__MANAGER__H

#include "xValManager.h"
#include "xTensors.h"
#include "xTensorsPtr.h"

namespace xfem
{
template <class T> class xValue;
class xValKey;
struct xHashValKey;
struct xEqualValKey;
typedef xValManager<xValKey, xValue< tensorsPtr_t >, xHashValKeyGauss, xEqualValKey> xVariabManager;
} // end of namespace

#endif
