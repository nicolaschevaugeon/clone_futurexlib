/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef POINT_TO_DOUBLE__H
#define POINT_TO_DOUBLE__H

#include <functional>

// xtensor
#include "xPoint.h"

namespace xfem
{
typedef std::function<double(const xtensor::xPoint&)> xPointToDouble;
}
#endif
