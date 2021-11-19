/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef REGION__H
#define REGION__H

#include <cassert>
#include <string>
#include <boost/iterator/filter_iterator.hpp>
#include "xMesh.h"


#ifdef PARALLEL 
#include "xParallel.h"
#endif
namespace AOMD{
  class mEntity;
}

namespace xfem
{
  class xMesh;
  class xSubMesh;
  //la classe xRegion est à virer à terme.
  
 
class xRegion {
	
public:
	xRegion();
	xRegion(const xMesh* m);
	xRegion(const xMesh* m, const std::string& subset);
	xRegion(const xSubMesh *sm);
        //mandatory interface of region
        xIter begin() const ;
        xIter   end() const ;
        int    size() const ;
        xIter begin(int what) const ;
        xIter   end(int what) const ;
        int    size(int what) const ;
        //end mandatory
        int dim() const { return dimension; }
	bool IsInRegion(AOMD::mEntity* e) const;
#ifdef PARALLEL 
	const xPartitionBoundaryInfo & getPartitionBoundaryInfo(void) const;
#endif
        AOMD::mEntity* find(AOMD::mEntity*) const ;
        xMesh* getMesh() const;
        const xSubMesh *getSubMesh() const{ return sub;}
        void setMesh(xMesh *m);
	
	std::string getSubName() const;
private:	
	xMesh* mesh;
        const xSubMesh *sub;
	int dimension;
      };

template <class ITER, class F>
class xFilteredRegion {		// 
public:
  typedef boost::filter_iterator<F, ITER> FilterIter;
  //typedef FilterIter::policies_type Policy;
public:
  xFilteredRegion(ITER b,  ITER e, const F& f) 
     : iter_beg(f, b, e), iter_end(f, e, e) {}
  xFilteredRegion(ITER b,  ITER e) 
     : iter_beg(F(),b,e), iter_end(F(), e, e) {}
  //mandatory interface  of region 
  FilterIter begin() const { return iter_beg; }
  FilterIter   end() const { return iter_end; }
  
private :
//  typename FilterIter::policies_type policy;
  //F filter;
  FilterIter  iter_beg, iter_end;
};

class xClassRegion {
	
private:	
	xMesh* mesh;
        int    entity_id, entity_dim, dimension;
public:
	xClassRegion(xMesh* m, int e_id, int e_dim);
        xClassIter begin();
        xClassIter   end();
        xClassIter begin(int what);
        xClassIter   end(int what);
	const xMesh* getMesh(){return mesh;};
	//	size_t size();
	//size_t size(in what);
      };


} // end of namespace

#endif







