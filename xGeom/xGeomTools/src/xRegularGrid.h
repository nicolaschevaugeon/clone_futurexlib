/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XREGULARGRID_H
#define _XREGULARGRID_H

// std
#include <vector>
// xtensor
#include "xPoint.h"
// xgeom
#include "xBoundingBox.h"

namespace xgeom
{
class xBrick
{
  public:
   std::vector<const void *> Objects;
   int operator<(const xBrick &other) { return iNumBrick < other.iNumBrick; }
   int operator>(const xBrick &other) { return iNumBrick > other.iNumBrick; }
   int operator==(const xBrick &other) { return iNumBrick == other.iNumBrick; }
   xBrick(int i) : iNumBrick(i) {}
   const void *operator[](int i)
   {
      if (i < 0 || (unsigned)i >= Objects.size()) throw i;
      return Objects[i];
   }
   int size() { return Objects.size(); }

  private:
   int iNumBrick;
};

/// a grid to store pointer to elements that cover part of each grid cell. USed for Fast retrivial of element that cover a point.
class xRegularGrid
{
   // a vector of bricks is created.
   // The grid is regular (Nx,Ny,Nz) so that finding
   // in what brick you are takes cosntant time.
  public:
   /// constructor from a bounding box. nx, ny, nz are the number of regular division in each direction.
   //! This constructor set up the data structures, but is empty (no element is associated to the xBricks)
   xRegularGrid(const xgeom::xBoundingBox &bb, int nx = 40, int ny = 40, int nz = 40);
   /// output the contents of the datastructure. usefull for debuging
   void print(std::ostream &out = std::cout) const;
   /// return the brick containing the point if any. should throw if no brick contains the point. const version
   const xBrick &getBrick(const xtensor::xPoint &P) const { return (getBrick(getBrickId(P))); }
   /// add the entity to each bricks that cover the given bb
   bool addObject(const void *, const xgeom::xBoundingBox &bb);
   /// return a vector containing the indices of xBricks covering the given bounding box.
   std::vector<int> getIndicesForBox(const xgeom::xBoundingBox &bb) const;

  private:
   constexpr static double TOL = 1.e-06;
   std::vector<xBrick> Table;
   const std::array<int, 3> N;
   const int numbox;
   xtensor::xPoint Pmin, Pmax;
   void initTable();
   int getIndex(int i, int j, int k) const { return i + j * N[0] + k * N[0] * N[1]; }
   int getBrickId(const xtensor::xPoint &P) const;
   /// return the brick containing the point if any. should throw if no brick contains the point
   xBrick &getBrick(const xtensor::xPoint &P) { return (getBrick(getBrickId(P))); }
   const xBrick &getBrick(int index) const
   {
      if (index < 0 || index >= numbox) throw index;
      return Table[index];
   }
   xBrick &getBrick(int index)
   {
      if (index < 0 || index >= numbox) throw index;
      return Table[index];
   }
};

}  // namespace xgeom

#endif
