/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xDataExchangerTools.h"



namespace xfem
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// keyManagerSendAndReceive implementation ////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
keyManagerSendAndReceive::keyManagerSendAndReceive ( const partmanAOMD_t & _partman ) : partman(_partman){}
auto keyManagerSendAndReceive::localObjectKey( const AOMD::mEntity & o) const->information_key_t
{
    return &o;
}
auto keyManagerSendAndReceive::remoteObjectKey(const xtool::xRemoteObject < AOMD::mEntity > &ro, const AOMD::mEntity  &lo ) const->information_key_t
{
    return ro.getObjectAddress();
}
xtool::xConstPartitionObject < AOMD::mEntity > keyManagerSendAndReceive::getConstPartitionObject(const AOMD::mEntity &e) const
{
    return partman.getConstPartitionObject(e);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// keyManagerSendOrReceive implementation /////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
keyManagerSendOrReceive::keyManagerSendOrReceive ( const partmanAOMD_t & _partman ) : partman(_partman){}
auto keyManagerSendOrReceive::localObjectKey( const AOMD::mEntity & o) const->information_key_t
{
    return &o;
}
std::set < int > keyManagerSendOrReceive::getMessageRanks( const information_key_t &lk) const
{
    // lk is the pointer to local object
    AOMD::mEntity &lo = *( const_cast < AOMD::mEntity * >( lk ));
    // So we can used it to access to remote object
    xtool::xConstPartitionObject < AOMD::mEntity > po = partman.getConstPartitionObject(lo);

    std::set < int > res;

    // if object have remote copy
    if (po.hasRemoteObject())
    {
        // grabe remote proc id for object and store them in res
        for (const auto ro : po.getRemoteObjectsCollectionRange( ))
        {
            res.insert(ro.getProcId());
        }
    }

    return res;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// function implementation ////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
} // end of namespace