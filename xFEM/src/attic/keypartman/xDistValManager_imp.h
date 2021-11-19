/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _DIST_VAL_MANAGER_IMP_H
#define _DIST_VAL_MANAGER_IMP_H

namespace xfem
{


template < typename K, typename V, typename HASH, typename EQUAL, typename MPM>
xValManagerKeyPartitionManager < K,V,HASH,EQUAL,MPM >::xValManagerKeyPartitionManager(const xValManager < K,V,HASH,EQUAL > &vm_,const MPM& mesh_partman_) :
    vm(vm_)
    ,mesh_partman(mesh_partman_)
    ,partition_manager(mesh_partman->getComm())
{}

/// function that create distributed informations
template < typename K, typename V, typename HASH, typename EQUAL, typename MPM>
void xValManagerKeyPartitionManager < K,V,HASH,EQUAL,MPM >::genPartitionManager(void)
{
    partition_manager.clear();
    throw -1;
}
template < typename K, typename V, typename HASH, typename EQUAL, typename MPM>
auto xValManagerKeyPartitionManager < K,V,HASH,EQUAL,MPM >::getPartitionManager(void) const->const partman_t &
{
    return partition_manager;
}

} // end of namespace

#endif









