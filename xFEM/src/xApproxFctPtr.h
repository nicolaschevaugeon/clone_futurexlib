/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _APPROX_FUNCTION_POINTER_H
#define _APPROX_FUNCTION_POINTER_H

#include <vector>
#include <memory>


namespace xfem
{
  class xApproxFunction;

typedef std::shared_ptr < xApproxFunction > approxFctPtr_t;
typedef std::vector<approxFctPtr_t> femFcts_t;
}

#endif
