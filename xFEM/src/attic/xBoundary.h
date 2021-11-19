/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __GROUP_H
#define __GROUP_H

#include <algorithm>
#include <map>
#include <set>

namespace xfem
{
class xBoundary
{
  public:
   xBoundary(){};
   bool AddToBoundary(const int &pToAdd);
   void EmptyBoundary();
   void RemoveFromBoundary(const int &pToErase);
   bool IsInBoundary(const int &pToFind) const;
   int size() const { return boundary.size(); };
   std::set<int>::const_iterator begin() const { return boundary.begin(); };
   std::set<int>::const_iterator end() const { return boundary.end(); };
   typedef std::set<int>::const_iterator const_iterator;

  private:
   std::set<int> boundary;
};

class xBoundaryContainer
{
  public:
   int getNbBoundaries();
   void addBoundary(int &id, xBoundary &Boundary);
   xBoundary *getBoundaryWithId(int id);
   xBoundary *getBoundaryAtPos(int pos);
   void AddToBoundaryWithId(int &id, const int &pToInsert);

  private:
   typedef std::map<int, xBoundary> BoundaryCont_t;

   BoundaryCont_t boundaries;
};

}  // namespace xfem

#endif
