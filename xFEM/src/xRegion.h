/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef REGION__H
#define REGION__H

#include <boost/iterator/filter_iterator.hpp>
#include <cassert>
#include <string>

#include "xIteratorTools.h"
#include "xMesh.h"

namespace AOMD
{
class mEntity;
}

namespace xfem
{
class xSubMesh;

class xRegion
{
  public:
   xRegion();
   xRegion(const xMesh* m);
   xRegion(const xMesh* m, const std::string& subset);
   xRegion(const xSubMesh* sm);
   // mandatory interface of region
   xIter begin() const;
   xIter end() const;
   xtool::xRange<xIter> range() const { return xtool::make_range(begin(), end()); }
   int size() const;
   xIter begin(int what) const;
   xIter end(int what) const;
   xtool::xRange<xIter> range(int what) const { return xtool::make_range(begin(what), end(what)); }
   int size(int what) const;
   // end mandatory
   int dim() const { return dimension; }
   bool IsInRegion(AOMD::mEntity* e) const;
   xMesh::partman_t& getPartitionManager() const;
   AOMD::mEntity* find(AOMD::mEntity*) const;
   xMesh* getMesh() const;
   const xSubMesh* getSubMesh() const { return sub; }
   void setMesh(xMesh* m);
   std::string getSubName() const;

  private:
   xMesh* mesh;
   const xSubMesh* sub;
   int dimension;
};

template <class ITER, class F>
class xFilteredRegion
{  //
  public:
   typedef boost::filter_iterator<F, ITER> FilterIter;
   // typedef FilterIter::policies_type Policy;
  public:
   xFilteredRegion(ITER b, ITER e, const F& f) : iter_beg(f, b, e), iter_end(f, e, e) {}
   xFilteredRegion(ITER b, ITER e) : iter_beg(F(), b, e), iter_end(F(), e, e) {}
   // mandatory interface  of region
   FilterIter begin() const { return iter_beg; }
   FilterIter end() const { return iter_end; }

  private:
   //  typename FilterIter::policies_type policy;
   // F filter;
   FilterIter iter_beg, iter_end;
};

class xClassRegion
{
  private:
   xMesh* mesh;
   int entity_id, entity_dim, dimension;

  public:
   xClassRegion(xMesh* m, int e_id, int e_dim);
   xClassIter begin();
   xClassIter end();
   xClassIter begin(int what);
   xClassIter end(int what);
   const xMesh* getMesh() {
       if (mesh) return mesh;
       else throw 75864;
   }
   size_t size() const;
   size_t size(int what) const;
};

}  // namespace xfem

#endif
