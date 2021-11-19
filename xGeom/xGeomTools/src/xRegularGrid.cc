/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// xgeom
#include "xRegularGrid.h"

namespace xgeom
{
using xgeom::xBoundingBox;
using xtensor::xPoint;

xRegularGrid::xRegularGrid(const xgeom::xBoundingBox &bb, int nx, int ny, int nz) : N{nx, ny, nz}, numbox{nx * ny * nz}
{
   initTable();
   const xtensor::xPoint PTOL{TOL, TOL, TOL};
   Pmin = bb.min - PTOL;
   Pmax = bb.max + PTOL;
}

void xRegularGrid::print(std::ostream &out) const
{
   out << "Here are the infos in the regular grid\n";
   out << "#####################################\n";
   for (int i = 0; i < N[0]; ++i)
      for (int j = 0; j < N[1]; ++j)
         for (int k = 0; k < N[2]; ++k)
         {
            int index = getIndex(i, j, k);
            const xBrick &brick = getBrick(index);
            int nb = brick.Objects.size();
            out << "The brick at location " << i << " " << j << " " << k << " with index " << index << " contains " << nb
                << "elts\n";
            /*            if (nb)
                        {
                           out << "The list is ";
                           for (const void *e : brick.Objects) out << e << " ";
                           out << "\n";
                        }
            */
         }
   out << "End of printing the regular grid infos\n";
   out << "#####################################\n";
   return;
}

bool xRegularGrid::addObject(const void *obj, const xgeom::xBoundingBox &bb)
{
   for (auto index : getIndicesForBox(bb))
   {
      bool known = false;
      xBrick &brick = getBrick(index);
      for (const void *pobj : brick.Objects)
         if (pobj == obj)
         {
            known = true;
            break;
         }
      if (!known) brick.Objects.push_back(obj);
   }
   return true;
}

// returns the indices of the little boxes covered by the box given in
// argument
std::vector<int> xRegularGrid::getIndicesForBox(const xgeom::xBoundingBox &bb) const
{
   int Imin[3], Imax[3];
   std::vector<int> indices;
   for (int i = 0; i < 3; ++i)
   {
      Imin[i] = (int)((double)N[i] * (bb.min(i) - Pmin(i)) / (Pmax(i) - Pmin(i)));
      Imin[i] = std::max(Imin[i], 0);
      Imax[i] = (int)((double)N[i] * (bb.max(i) - Pmin(i)) / (Pmax(i) - Pmin(i)));
      Imax[i] = std::min(Imax[i], N[i] - 1);
   }
   for (int i = Imin[0]; i <= Imax[0]; ++i)
      for (int j = Imin[1]; j <= Imax[1]; ++j)
         for (int k = Imin[2]; k <= Imax[2]; ++k) indices.push_back(getIndex(i, j, k));
   return indices;
}

void xRegularGrid::initTable()
{
   Table.clear();
   Table.reserve(numbox);
   for (int i = 0; i < numbox; ++i) Table.push_back(xBrick(i));
}

int xRegularGrid::getBrickId(const xtensor::xPoint &P) const
{
   int I[3];
   for (int i = 0; i < 3; ++i)
   {
      int k = (int)((double)N[i] * (P(i) - Pmin(i)) / (Pmax(i) - Pmin(i)));
      k = std::min(k, N[i] - 1);
      if (k < 0) k = 0;
      I[i] = k;
   }
   return getIndex(I[0], I[1], I[2]);
}

}  // namespace xgeom
