/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xBoundary.h"

#include <cstdio>
#include <string>

namespace xfem
{
bool xBoundary::AddToBoundary(const int &pToInsert) { return (boundary.insert(pToInsert).second); }

void xBoundary::EmptyBoundary() { boundary.erase(boundary.begin(), boundary.end()); }
bool xBoundary::IsInBoundary(const int &pToFind) const
{
   if (boundary.find(pToFind) == boundary.end())
      return false;
   else
      return true;
}

void xBoundary::RemoveFromBoundary(const int &pToErase) { boundary.erase(pToErase); }

xBoundary *xBoundaryContainer::getBoundaryWithId(int id)
{
   BoundaryCont_t::iterator found;

   found = boundaries.find(id);

   if (found == boundaries.end()) return nullptr;

   return &((*found).second);
}

void xBoundaryContainer::AddToBoundaryWithId(int &id, const int &pToInsert)
{
   xBoundary boundary;
   xBoundary *boundary_ptr = getBoundaryWithId(id);
   if (boundary_ptr != nullptr)
      boundary_ptr->AddToBoundary(pToInsert);
   else
   {
      boundary.AddToBoundary(pToInsert);
      addBoundary(id, boundary);
   }
   return;
}

void xBoundaryContainer::addBoundary(int &id, xBoundary &Boundary)
{
   std::pair<const int, xBoundary> to_add(id, Boundary);

   boundaries.insert(to_add);
}

xBoundary *xBoundaryContainer::getBoundaryAtPos(int p_pos)
{
   BoundaryCont_t::iterator it = boundaries.begin();
   BoundaryCont_t::iterator ite = boundaries.end();

   xBoundary *ret = nullptr;

   int pos = 0;

   while (it != ite)
   {
      if (pos == p_pos)
      {
         ret = &((*it).second);
         break;
      }

      pos++;
      it++;
   }

   return ret;
}

int xBoundaryContainer::getNbBoundaries() { return boundaries.size(); }

}  // namespace xfem
