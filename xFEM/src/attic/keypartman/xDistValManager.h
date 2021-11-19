/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _DIST_VAL_MANAGER_H
#define _DIST_VAL_MANAGER_H

#include "xValManager.h"
#include "xUnorderedMapDataManager.h"
#include "xPartitionManager.h"


namespace xfem
{


template < typename K, typename V, typename HASH, typename EQUAL, typename MPM>
class xValManagerKeyPartitionManager
{

    private:
        template < typename T >
        using data_manager_t = xtool::xUnorderedMapDataManager < K, T >;


    public:
        typedef xtool::xPartitionManager < data_manager_t > partman_t;

        /// Constructor that store a reference to the value manager for which we
        //! will need distributed information
        xValManagerKeyPartitionManager(const xValManager < K,V,HASH,EQUAL > &vm_, const MPM& mesh_partman_);

        /// function that create distributed information
        void genPartitionManager(void);

        /// function to obtain the const version of the distributed information
        auto getPartitionManager(void) const->const partman_t &;

    private:
        const xValManager < K,V,HASH,EQUAL > &vm;
        const MPM &mesh_partman;
        partman_t partition_manager;
};

} // end of namespace

#include "xDistValManager_imp.h"

#endif
