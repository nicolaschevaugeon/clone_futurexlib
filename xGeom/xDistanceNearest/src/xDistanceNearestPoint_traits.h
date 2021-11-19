/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _DISTANCE_CLOSEST_POINT_TRAITS_H
#define _DISTANCE_CLOSEST_POINT_TRAITS_H

namespace xgeom
{
/// \brief Traits structures for specialisations of template
/// class xDistanceNearestPoint_internal implemented in
/// xDistanceNearestPoint_internal_imp.h. These structures are
/// to be declared as a typedef inside generator classes
/// implemented in xDistanceNearestPoint.h. Generator classes
/// are then used in combination with xDistanceNearestPoint
/// template class. xDistanceNearestPoint_internal must not be
/// used in applications.
struct xCGALTrait
{
};
struct xANNTrait
{
};
struct xBruteForceTrait
{
};
}  // namespace xgeom

#endif
