/*
   xfem : C++ Finite Element Library
   developed under the GNU Lesser General Public License
   See the NOTICE & LICENSE files for conditions.
*/

#ifndef Interface_OCTREEDB_TO_AOMD__H
#define Interface_OCTREEDB_TO_AOMD__H

#include <string>

#include "oLevelSet.h"
#include "oOctree.h"
#include "xAttachedDataManagerAOMD.h"
#include "xLevelSet.h"
#include "xMesh.h"

// InterfaceOctreeToAOMD
namespace xinterface
{
namespace xoctree
{
using ::xoctree::oField;
using ::xoctree::oKey;
using ::xoctree::oKeyManager;
using ::xoctree::oMapping;
using ::xoctree::oOctree;
using ::xoctree::oTopo;

const int init_physical_index = 100;
// base class to custumize how xMesh entities will be classifyed when created from an oOctree
class iClassifyCriteria
{
  public:
   iClassifyCriteria() = default;
   virtual ~iClassifyCriteria() = default;
   virtual int operator()(oOctree::const_iterator cell, int level, const int* ijk, oOctree::const_iterator children_beg,
                          oOctree::const_iterator children_end) const = 0;
};

/// from an oOctree and a oField, create an xMesh and a levelSet on the xMesh
/*!note         alternate : if alternate set to true, quad spliting in triangle in 2d wil alternate from one cel to the next.
                the parameter replace the previous preprocessing directive ALTERNATE which was set to 1 by default in
   InterfaceOctreeToAOMD It seems, at the time of factorizing this ALTERNATE directive into a parameter, that the code can only
   work with alternate =true and that alternate =false would give bad answer... I think that we need to either rethink the
   alternate/notalternate or just remove the possibility.
*/
template <template <class> class DATAMANAGER = xfem::xMesh::datamanager_t>
void InterfaceOctreeToAOMD(const oOctree& octree, const oField& lsoct, xfem::xMesh& mesh, xfem::xLevelSet& ls,
                           const bool simplex = true, iClassifyCriteria* ptr_clfCriteria = nullptr, bool alternate = true,
                           DATAMANAGER<int>& octreeLevel = xfem::xMesh::get_octree_level(),
                           DATAMANAGER<AOMD::mEntity*>& isHangingOn = xfem::xMesh::get_is_hanging_on(),
                           DATAMANAGER<AOMD::mEntity*>& isHangingBy = xfem::xMesh::get_is_hanging_by(),
                           DATAMANAGER<std::vector<AOMD::mEntity*>>& downGroup = xfem::xMesh::get_down_group(),
                           DATAMANAGER<std::vector<AOMD::mEntity*>>& bndGroup = xfem::xMesh::get_bnd_group());
// same as above, but without hanging nodes (more cuts)
template <template <class> class DATAMANAGER = xfem::xMesh::datamanager_t>
void InterfaceOctreeToAOMDNoHanging(const oOctree& octree, const oField& lsoct, xfem::xMesh& mesh, xfem::xLevelSet& ls,
                                    const bool coarser = false, const bool simplex = true,
                                    iClassifyCriteria* ptr_clfCriteria = nullptr,
                                    DATAMANAGER<int>& octreeLevel = xfem::xMesh::get_octree_level()

);
// push the lsoct value to the ls values.
void OctreeToAOMDProjection(const oOctree& octree, const oField& lsoct, xfem::xMesh* mesh, xfem::xLevelSet& ls);
void OctreeToAOMDNoHangingProjection(const oOctree& octree, const oField& lsoct, xfem::xMesh* mesh, xfem::xLevelSet& ls);
// push the ls value to the lsoct value.... of course lsoct must be on an octree compatible with the mesh
void InterfaceAOMDToOctree(const xfem::xMesh& mesh, const xfem::xLevelSet& ls, oField& lsoct);

}  // namespace xoctree

}  // namespace xinterface

#include "InterfaceOctreeToAOMD_imp.h"

#endif
