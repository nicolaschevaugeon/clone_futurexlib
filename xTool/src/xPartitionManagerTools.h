/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _xPartitionManagerTools_
#define _xPartitionManagerTools_
#include <cassert>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include "mpi.h"
#include "xDataExchanger.h"
#include "xPartitionManager.h"

namespace xtool
{


/*! createPartitionManagerForSubGroup shrink an initial partition manager to reflect behavior of a subset of element
    The subset of element is viewed as a mesh with its own partition manager which is constructed by this function.
    Connectivity of this subset rely on the underlying mesh where sub group of elements have been picked out.
    \param [in] sub_group_boundary_range range of entity of TopLevel-1 describing  subset of element boundary common to
    many proc
    \param [in] general_partman actual partition manager managing underlying mesh partition used as source of element subset
    \param [in,out] sub_group_partman partition manager filled by this function. It is user responsibility to create an empty one
                    acting on the MPI communicator that he wants. This function then add all connections related to the given subset
                    and communicator.
                
 */
template <  template < typename > class DM, typename T, typename RANGE, typename GPM, typename SGPM >
void createPartitionManagerForSubGroup(RANGE sub_group_boundary_range, const GPM & general_partman, SGPM & sub_group_partman);

/*! Generate a graphviz dof file from given partition manager
 *  One file is create for each dimension (i.e. for node, edge,face)
 *  All process send there graph part
 *  Graph node are node, edge or face described by there address, there id, there type 
 *  Directed edge show the counterpart of an entity on remote process
 */
template <  typename PM >
void exportInGraphvizDof(const PM & partman, std::string name, size_t filter=0);

} // end namespace

#include "xPartitionManagerTools_imp.h"

#endif
