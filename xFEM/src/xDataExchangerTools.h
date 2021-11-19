/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XDATAEXCHANGERTOOLS_H
#define _XDATAEXCHANGERTOOLS_H

#include <iostream>

// xtools
#include "xPartitionManager.h"

// xinterface aomd
#include "xAttachedDataManagerAOMD.h"


// trellis
#include "mEntity.h"

namespace xfem
{
/// For many user this partition manager type is well suited for there need. Most (but not all) class of xfem library rely on it.
  typedef xtool::xPartitionManager < xinterface::aomd::xAttachedDataManagerAOMD > partmanAOMD_t;

/// This class is to be consider as a default concrete implementation of a key manager for accumulate method of xKeyContainerSendAndRecv class
class keyManagerSendAndReceive
{
    public:
        typedef const AOMD::mEntity * information_key_t;

        keyManagerSendAndReceive ( const partmanAOMD_t & _partman );

        information_key_t localObjectKey( const AOMD::mEntity & o) const;

        information_key_t remoteObjectKey(const xtool::xRemoteObject < AOMD::mEntity > & ro, const AOMD::mEntity  &lo ) const;

        xtool::xConstPartitionObject < AOMD::mEntity > getConstPartitionObject(const AOMD::mEntity &e) const;
    protected:
        const partmanAOMD_t &partman;
};

/// This class is to be consider as a default concrete implementation of a key manager for accumulate method of xKeyContainerSendOrRecv class
class keyManagerSendOrReceive
{
    public:
        typedef const AOMD::mEntity * information_key_t;

        keyManagerSendOrReceive ( const partmanAOMD_t & _partman );

        information_key_t localObjectKey( const AOMD::mEntity  & o) const;

        std::set < int > getMessageRanks( const information_key_t & lk) const;
    private:
        const partmanAOMD_t &partman;
};

} // end of namespace

#endif
