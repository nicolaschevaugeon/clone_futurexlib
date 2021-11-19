/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xNearestNeighborInterface.h"

#include "xMesh.h"

xFindClosestTo1dMesh_BruteForce::xFindClosestTo1dMesh_BruteForce(const xfem::xMesh *_mesh) : mesh(_mesh){};
void xFindClosestTo1dMesh_BruteForce::nearestpointto(const double *P, double *Cp)
{
   double A[2], B[2], Cpe[2];
   double dist = std::numeric_limits<double>::max();
   for (auto e : mesh->range(1))
   {
      const xtensor::xPoint p0 = static_cast<const AOMD::mVertex *>(e->get(0, 0))->point();
      const xtensor::xPoint p1 = static_cast<const AOMD::mVertex *>(e->get(0, 1))->point();
      A[0] = p0(0);
      A[1] = p0(1);
      B[0] = p1(0);
      B[1] = p1(1);
      xClosestPointToEdge2d(A, B, P, Cpe);
      const double diste = (Cpe[0] - P[0]) * (Cpe[0] - P[0]) + (Cpe[1] - P[1]) * (Cpe[1] - P[1]);
      if (diste < dist)
      {
         Cp[0] = Cpe[0];
         Cp[1] = Cpe[1];
         dist = diste;
      }
   }
}
