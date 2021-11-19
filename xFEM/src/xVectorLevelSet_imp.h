/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include "xVectorLevelSet.h"

namespace xfem
{
/// Implentation of template classes
template <class VERTEXITERATOR>
xVectorLevelSet::xVectorLevelSet(const VERTEXITERATOR &_vbegin, const VERTEXITERATOR &_vend, xPointToVectorLevelSetData vlseval)
{
   for (VERTEXITERATOR itv = _vbegin; itv != _vend; ++itv)
   {
      AOMD::mVertex *v = static_cast<AOMD::mVertex *>(*itv);
      data.emplaceSetData(*v, vlseval(v->point()));
   }
}

template <typename VERTEXITERATOR, typename GETCLOSEST>
xVectorLevelSet::xVectorLevelSet(const VERTEXITERATOR &_vbegin, const VERTEXITERATOR &_vend, const xLevelSet &ls,
                                 GETCLOSEST &getclosestpoint, const double &shift)
{
   AOMD::mEntity *nearest_e;
   bool exist;

   // loop on mesh vertex
   for (VERTEXITERATOR itv = _vbegin; itv != _vend; ++itv)
   {
      exist = true;
      AOMD::mVertex *v = static_cast<AOMD::mVertex *>(*itv);
      double lsx = ls(v);
      int inout;
      double dist;
      xtensor::xPoint current_point = v->point();
      xtensor::xPoint closest_point;
      xtensor::xPoint closest_point_from_isolc;
      getclosestpoint.nearestEntityDistance((*itv), nearest_e, closest_point, dist);
      xtensor::xVector<> ulc(current_point, closest_point);
      double delta = dist - shift;
      if (shift == 0.)
      {
         if (lsx < 0.)
            inout = -1;
         else if (lsx > 0.)
            inout = 1;
         else
            inout = 0;
      }
      else
      {
         xtensor::xVector<> u(ulc);
         u.norm();
         const double eps = 1.e-7;
         if (lsx < (-eps))
         {
            if (delta > eps)
               inout = -1;
            else if (delta < (-eps))
               inout = 1;
            else
               inout = 0;
            ulc = u * delta;
            // check not in sckeleton vicinity
            xtensor::xPoint on_iso_lc(ulc[0] + current_point(0), ulc[1] + current_point(1), ulc[2] + current_point(2));
            getclosestpoint.nearestEntityDistance(on_iso_lc, nearest_e, closest_point_from_isolc, dist);
            if (dist < shift && inout == 1)
            {
               // check that the new closest point is locate in the opposite half space
               xtensor::xVector<> dvect(closest_point_from_isolc, closest_point);
               if (dvect * ulc > 0. && 0)
               {
                  ulc = u * eps;
                  exist = false;
               }
            }
         }
         else if (lsx > eps)
         {
            ulc = u * (dist + shift);
            inout = 1;
         }
         else
         {
            ulc = ls.getGrad(v);
            ulc.norm();
            ulc *= -shift;
            inout = 1;
         }
      }

      data.emplaceSetData(*v, inout, ulc, exist);
   }
   return;
}

}  // end namespace xfem
