/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XATTACHABLE_GP_H
#define _XATTACHABLE_GP_H

// Trellis
#include "mAttachableDataContainer.h"

// xtensor
#include "xPoint.h"

namespace xinterface
{
namespace aomd
{
class xAttachableGaussPoints : public AOMD::mAttachableData
{
  public:
   typedef std::vector<std::pair<xtensor::xPoint, double>> Container;
   typedef Container::iterator Iter;
   xAttachableGaussPoints() __attribute__((deprecated)) = default;

   Container gauss_points;
};

}  // namespace aomd
}  // namespace xinterface

#endif
