/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __DOUBLE__MANAGER__H
#define __DOUBLE__MANAGER__H

#include "xValKey.h"
#include "xValManager.h"
#include "xDistValManager.h"
#include "xValue.h"
#include "xValKeyDataManager.h"
#include "xPartitionManager.h"

namespace xfem 
{
typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
struct xCloneDoubleValue
{
    xValue<double>* operator () (const  xValue<double>* val, const std::map<xValue<double> *,xValue<double> *> &coresp) const
    {
        return xCloneValue(val,coresp);
    }
};

/// xValManagerKeyPartitionManager specialization for double manager with AOMD mesh partition manager
template < >
class xValManagerKeyPartitionManager < xValKey,  xValue<double>, xHashValKey, xEqualValKey, xMesh::partman_t>
{

    private:
        template < typename T >
        using data_manager_t = xValKeyDataManager < T >;


    public:
        typedef xtool::xPartitionManager < data_manager_t > partman_t;


        /// Constructor that store a reference to the value manager for which we
        //! will need distributed information
        xValManagerKeyPartitionManager(const xValueManagerDist<double> &dm_, const xMesh::partman_t& mesh_partman_);

        /// function that create distributed information
        void genPartitionManager(void);

        /// function to obtain the const version of the distributed information
        auto getPartitionManager(void) const->const partman_t &;

    private:
        const xValueManagerDist<double> &dm;
        const xMesh::partman_t &mesh_partman;
        partman_t partition_manager;
};


typedef xValManagerKeyPartitionManager < xValKey,  xValue<double>, xHashValKey, xEqualValKey, xMesh::partman_t> xValueManagerDist<double>KeyPartitionManager;

} // end of namespace

#endif
