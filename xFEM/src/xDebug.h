/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __xdebug_h
#define __xdebug_h

#include "xSingleton.h"

namespace xfem {

class xDebug 
{
public:
  xDebug() : flag(false) {}
  const bool& getFlag() const {return flag;}
private:
  void setFlag(bool b)  {flag = b;}  
  bool flag;
  friend void SetXfemDebugFlag(bool b);
};

typedef xtool::xSingleton<xDebug> xDebugSingleton;

static const bool& xdebug_flag = xDebugSingleton::instance().getFlag();
void SetXfemDebugFlag(bool b);

} // end namespace



#endif
