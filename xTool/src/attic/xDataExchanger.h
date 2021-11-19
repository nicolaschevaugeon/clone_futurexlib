/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef XDATAEXCHANGER_H
#define XDATAEXCHANGER_H

//std
#include <cassert>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <string>
#include <exception>
#include <sstream>
#include <type_traits>

// parrallel
#include <mpi.h>

//eXlibris_tools
#include "xDataExchangerTraits.h"
#include "xPartitionManager.h"
#include "xMPITag.h"


// notes and Ideas of semantics ...
//  we have iterator to objects
//    for each object o we get from Info (usally throught a xPArtitionManager inside info) a xPartitionObject (po)
//      from the xPArtitionObject, we can get the remotes object associated to object throught iterator itro
//          for each ro we can get a key from Infomanager (supposed to be unique and consistant accroos procs) ...
//               in most cases the key is just the address inside the remote object ....
//           Then the communication process start
//               For each key Infomanager will be prompted for some data to exchange.
//               For each key marqued to be receive for, the infomanager will get some data.
//  Proposition :  an object and it's remote object could give rise to more than one key ... That could fit neatly
//     With the way dof are created in Xfem  object are element, form element and fields keyssss are created, for //     each of which one has data (the dofs status, value and number )



namespace xtool
{

class xMpiInputBuffer
{
    public:
        inline xMpiInputBuffer(const MPI_Comm & _comm, int reserve = 0);
        inline int   pack(const void *inbuff, int incount, MPI_Datatype datatype);
        inline int   size() const;
        inline void *data();
    private:
        std::vector < char > dat;
        const MPI_Comm comm;
        mutable int position = 0;
};

class xMpiOutputBuffer
{
    public:
        inline xMpiOutputBuffer(const MPI_Comm & _comm, int size);
        inline int         unPack(void *outbuff, int outcount, MPI_Datatype datatype) const;
        inline const void *data() const;
        inline void       *data();
        inline bool exhausted() const {return position >= (int) ( dat.size()); }
    private:
        const MPI_Comm comm;
        std::vector < char > dat;
        mutable int position = 0;
};


template < class K, class KEYDATAASSOCIATION >
class exchangeInformationPolicy
{
    public:
        template < typename  I >
        void exchangeInformation(const int tag_info, void ( *work )(void), I & info_manager, MPI_Comm mpi_comm) {throw; }
};

template < typename K >
class xDataExchanger
{
    public:
        typedef typename K::information_key_t information_key_t;
        typedef exchangeInformationPolicy < K, typename K::data_key_association_trait  > exchange_information_policy_t;
        // both info_key_recv_container_t and  info_key_send_container_t are map first indexed with the processorId
        //	typedef std::map < int, std::set < information_key_t > >                    info_key_recv_container_t;
        // typedef std::map < int, std::map < information_key_t, information_key_t > > info_key_send_container_t;

        // typedef typename std::set < information_key_t >::const_iterator info_key_iter_t;
        //typedef typename std::map < information_key_t, information_key_t >::const_iterator info_key_pair_iter_t;
        xDataExchanger(K & _key_info_manager, MPI_Comm _mpi_comm);
        typename exchange_information_policy_t::info_key_send_container_t getSendKeys()
        {
            static_assert(!std::is_same < typename K::data_key_association_trait, recv_only_keys_association_trait >::value,
                          " getSendKey not valid  for recv_only_keys_association_trait ");
            static_assert(!std::is_same < typename K::data_key_association_trait, no_keys_association_trait >::value,
                          " getSendKey not valid for no_keys_association_trait ");
            return expolicy.getSendKeys();
        }

        void setRecvKeys(const typename exchange_information_policy_t::info_key_recv_container_t &  _recv_keys_from )
        {
            static_assert(!std::is_same < typename K::data_key_association_trait, send_only_keys_association_trait >::value,
                          " setRecvKeys not valid  for send_only_keys_association_trait ");
            static_assert(!std::is_same < typename K::data_key_association_trait, no_keys_association_trait >::value,
                          " setRecvKeys not valid for no_keys_association_trait ");
            expolicy.recv_keys_from = _recv_keys_from;
        }

        typename exchange_information_policy_t::info_key_recv_container_t getRecvKeys()
        {
            static_assert(!std::is_same < typename K::data_key_association_trait, send_only_keys_association_trait >::value,
                          " getRecvKey not valid  for send_only_keys_association_trait ");
            static_assert(!std::is_same < typename K::data_key_association_trait, no_keys_association_trait >::value,
                          " getRecvKey not valid for no_keys_association_trait ");
            return expolicy.getRecvKeys();
        };

        void clearKeys( )
        {
            expolicy.clearKeys();
        }

        /// each local object in the iterator range send to all it's remote object if the object is owned by the current proc. the data send to the remote object is not necesseraly the same.
        /// Collective on the MPI_Comm of the exchanger
        /// Only avalaible for exchangeInformationPolicy with trait send_recv_keys_association_trait
        // We might need to consider the  BroadCast case, where the same data is send to all the remote of the object
        template < typename ITER >
        void accumulateKeysOwnerScatter( ITER itb, ITER ite );
        // Question : should Owner send to himself ... not really send but everything done as if it had send ...
        //            This would insure same semantic as an MPI_BCast and perhaps simplify some algo

        /// each key in the iterator range send to all it's owners's remote key if the key is not owned
        /// Collective on the MPI_Comm of the exchanger
        /// Only avalaible for exchangeInformationPolicy with trait send_recv_keys_association_trait
        template < typename ITER >
        void accumulateKeysOwnerGather( ITER itb, ITER ite );
        // Question : should Owner send to himself ... not really send but everything done as if it had send ...
        //            This would insure same semantic as an MPI_BCast and perhaps simplify some algo

        /// each key in the iterator range send to all it's remote key
        /// Collective on the MPI_Comm of the exchanger
        template < typename ITER >
        void accumulateKeysAllGather( ITER itb, ITER ite );
        // Experimental : not ready yet .....
        /// each local object in the iterator range send to all it's remote object if the object is owned by the current proc. the data send to the remote object is not necesseraly the same.
        /// Collective on the MPI_Comm of the exchanger
        /// Only avalaible for exchangeInformationPolicy with trait send_only_keys_association_trait
        template < typename ITER >
        void accumulateKeysSenderOnly(ITER itb, ITER ite);
        template < typename OBJECTYPE >
        void accumulateKeysSenderOnly(const OBJECTYPE & lo);

        template < typename ITER >
        void accumulateKeysReceiverOnly(ITER itb, ITER ite);

        void exchangeInformation( void ( *work )(void) = NULL)
        {
            exchangeInformation(key_info_manager,work);
        }
        template < typename  I >
        void exchangeInformation( I & info_manager,  void ( *work )(void) = NULL)
        {
            // if there is less then 2 proc we don't need to exchange anything
            if (mpi_size < 2) return;
            // if not part of a communication world
            if (mpi_comm == MPI_COMM_NULL) return;
            // arbitrary tag used for this sequence
            const int tag_info = xtool::xMPITag::getNewTag(mpi_comm);
            // exchange sequence
            expolicy.exchangeInformation( tag_info, work, info_manager, mpi_comm );
        }

    private:
        exchange_information_policy_t expolicy;
        K &key_info_manager;
        const MPI_Comm mpi_comm;
        int mpi_rank, mpi_size;

        template < typename OBJECTYPE >
        void accumulateKeysSenderOnlyPrivate(const OBJECTYPE & lo);

};


class xDataExchangerException : public std::exception
{
    public:
        xDataExchangerException(std::string message, std::string file, int line, std::string date, std::string time );
        ~xDataExchangerException() throw( );
        const char * what() const throw( );

    private:
        std::string msg;
};

/// Function to get an object out of an object. By itself it doesn't look very helpful. But with overload version with pointer it helps
/// getting object out iterator used in accumulate methods
template < typename OBJECTYPE >
inline OBJECTYPE & getObject(OBJECTYPE &o) { return o; }
/// Function to get object out of pointer argument. By itself it doesn't look very helpful. But with overload version with reference it helps
/// getting object out iterator used in accumulate methods
template < typename OBJECTYPE >
inline OBJECTYPE & getObject(OBJECTYPE *o) { return *o; }


} // end namspace


//***********************************************************
//     All Implementation start here !
//***********************************************************


namespace xtool
{
//***********************************************************
// xMpiInputBuffer Implementation
//***********************************************************
xMpiInputBuffer::xMpiInputBuffer(const MPI_Comm & _comm, int reserve ) : dat(), comm(_comm)
{
    dat.reserve(reserve);
}

int xMpiInputBuffer::pack(const void *inbuff, int incount, MPI_Datatype datatype)
{
    int packsize;
    MPI_Pack_size(incount, datatype, comm, &packsize);
    const int datsize = dat.size();
    if (position+packsize > datsize ) dat.resize(2*( position+packsize ));
    int tmp = MPI_Pack(const_cast < void * >( inbuff ), incount, datatype, data(), dat.size(), &position, comm);
    /*if (datatype == MPI_AINT){
       std::cout << "IN pack "<< inbuff << std::endl;
       }*/
    return tmp;
}

int xMpiInputBuffer::size() const {return position; }

void * xMpiInputBuffer::data(){return dat.data(); }

//***********************************************************
// xMpiOutputBuffer Implementation
//***********************************************************
xMpiOutputBuffer::xMpiOutputBuffer(const MPI_Comm & _comm, int size) : comm(_comm), dat(size){}
int xMpiOutputBuffer::unPack(void *outbuff, int outcount, MPI_Datatype datatype) const
{
    int tmp = MPI_Unpack( const_cast < void * >( data()),  dat.size(), &position, outbuff, outcount, datatype,  MPI_COMM_WORLD  );
    /*if (datatype == MPI_AINT){
       std::cout << "IN un pack "<< outbuff << std::endl;
       }*/
    return tmp;


}
const void *xMpiOutputBuffer::data() const {return dat.data(); }
void *xMpiOutputBuffer::data() {return dat.data(); }


//***********************************************************
// xDataExchanger Implementation
//***********************************************************
template < typename K >
xDataExchanger < K >::xDataExchanger(  K & _key_info_manager, MPI_Comm _mpi_comm) : key_info_manager(_key_info_manager), mpi_comm(_mpi_comm)
{
    // getting mpi_rank and mpi_size in mpi_comm
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm,&mpi_size);
}


template < typename K >
template < typename ITER >
void xDataExchanger < K >::accumulateKeysOwnerScatter( ITER itb, ITER ite )
{
    static_assert(std::is_same < typename K::data_key_association_trait, send_recv_keys_association_trait >::value,
                  " accumulateKeysOwnerScatter only valid for send_recv_keys_association_trait ");

    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;


    // loop
    for (ITER it = itb; it != ite; ++it)
    {
        const auto &lo = getObject(*it);       // lo is the reference to the local object
        const information_key_t lk = key_info_manager.localObjectKey(lo);     // lk is the key to the local object
        const auto po = key_info_manager.getConstPartitionObject(lo);     // po is the partition object associated to lo


        if (po.hasRemoteObject())
        {
            if (!po.isOwner())
            {
                const auto ro_owner = po.getOwner(); //ro_owner is the remote object on the process that own the object.
                expolicy.recv_keys_from[ro_owner.getProcId()].insert(lk);
            }
            else
            {
                for (const auto ro : po.getRemoteObjectsCollectionRange( ))
                {
                    // ro is a remoteObject associated to lo.
                    const information_key_t rk = key_info_manager.remoteObjectKey(ro, lo);
                    expolicy.send_keys_to[ro.getProcId()].insert( std::make_pair( rk,  lk ));
                }
            }
        }
    }
    return;
}

template < typename K >
template < typename ITER >
void xDataExchanger < K >::accumulateKeysOwnerGather( ITER itb, ITER ite )
{
    static_assert(std::is_same < typename K::data_key_association_trait, send_recv_keys_association_trait >::value,
                  " accumulateKeysOwnerGather only valid for send_recv_keys_association_trait ");
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;

    // loop on iterator to collect information key and put them in correct order
    for (ITER it = itb; it != ite; ++it)
    {
        const auto &lo = getObject(*it);     // lo is the reference to the local object
        const information_key_t lk = key_info_manager.localObjectKey(lo);   // lk is the key to the local object
        const auto po = key_info_manager.getConstPartitionObject( lo);               // po is the partition object associated to po
        if (po.hasRemoteObject())
        {
            if (po.isOwner())
            {
                // ro is a remoteObject associated to lo.
                for (auto ro : po.getRemoteObjectsCollectionRange( )) expolicy.recv_keys_from[ro.getProcId()].insert(lk);
            }
            else
            {
                const auto ro = po.getOwner(); //ro is the remote object on the process that own the object.
                const information_key_t rk = key_info_manager.remoteObjectKey(ro, lo);
                expolicy.send_keys_to[ro.getProcId()].insert ( std::make_pair( rk,  lk ) );
            }
        }
    }
    return;
}

template < typename K >
template < typename ITER >
void xDataExchanger < K >::accumulateKeysAllGather( ITER itb, ITER ite )
{
    static_assert(std::is_same < typename K::data_key_association_trait, send_recv_keys_association_trait >::value,
                  " accumulateKeysOwnerGather only valid for send_recv_keys_association_trait ");
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;

    // loop on iterator to collect information key and put them in correct order
    for (ITER it = itb; it != ite; ++it)
    {
        const auto &lo = getObject(*it);           // lo is the reference to the local object
        const information_key_t lk = key_info_manager.localObjectKey(lo);         // lk is the key to the local object
        const auto po = key_info_manager.getConstPartitionObject( lo);                     // po is the partition object associated to po
        if (po.hasRemoteObject())
        {
            // ro is a remoteObject associated to lo.
            for (auto ro : po.getRemoteObjectsCollectionRange( ))
            {
                expolicy.recv_keys_from[ro.getProcId()].insert(lk);
                const information_key_t rk = key_info_manager.remoteObjectKey(ro, lo);
                expolicy.send_keys_to[ro.getProcId()].insert ( std::make_pair( rk,  lk ) );
            }
        }
    }
    return;
}

template < typename K >
template < typename ITER >
void xDataExchanger < K >::accumulateKeysSenderOnly(ITER itb, ITER ite)
{
    static_assert(std::is_same < typename K::data_key_association_trait, send_only_keys_association_trait >::value,
                  " accumulateKeySenderOnly only valid for send_only_keys_association_trait ");
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;
    for (ITER it = itb; it != ite; ++it)
    {
        accumulateKeysSenderOnlyPrivate(getObject(*it));
    }
}
template < typename K >
template < typename OBJECTYPE >
void xDataExchanger < K >::accumulateKeysSenderOnly( const OBJECTYPE & lo)
{
    static_assert(std::is_same < typename K::data_key_association_trait, send_only_keys_association_trait >::value,
                  " accumulateKeySenderOnly only valid for send_only_keys_association_trait ");
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;

    // do work for lo
    accumulateKeysSenderOnlyPrivate(lo);
}

template < typename K >
template < typename OBJECTYPE >
void xDataExchanger < K >::accumulateKeysSenderOnlyPrivate( const OBJECTYPE & lo)
{
    const information_key_t lk = key_info_manager.localObjectKey(lo);         // lk is the key to the local object
    std::set < int > targets = key_info_manager.getMessageTargetRanks(lk);
    for (auto targetproc : targets)
    {
        expolicy.send_keys_to[targetproc].insert(lk);
    }
}


template < typename K >
template < typename ITER >
void xDataExchanger < K >::accumulateKeysReceiverOnly(ITER itb, ITER ite)
{
    static_assert(std::is_same < typename K::data_key_association_trait, recv_only_keys_association_trait >::value,
                  " accumulateKeyReceiverOnly only valid for recv_only_keys_association_trait ");
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;
    for (ITER it = itb; it != ite; ++it)
    {
        const auto &lo = getObject(*it);       // lo is the reference to the local object
        const information_key_t lk = key_info_manager.localObjectKey(lo);     // lk is the key to the local object
        std::set < int > sources = key_info_manager.getMessageSourceRanks(lk);
        for (auto sourceproc : sources)
        {
            expolicy.recv_keys_from[sourceproc].insert(lk);
        }
    }
}

template < typename K, typename I, typename DATASTYLE >
class sendRecvKeysExchangeInformationPolicy
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void), I & info_manager, K & key_info_manager, MPI_Comm mpi_comm) {throw; }
};
//***********************************************************
//     exchangeInformationPolicy for send_recv_keys_association_trait
//***********************************************************
template < class K >
class exchangeInformationPolicy < K, send_recv_keys_association_trait >
{
    public:
        typedef typename K::information_key_t information_key_t;
        typedef std::map < int, std::set < information_key_t > >                    info_key_recv_container_t;
        typedef std::map < int, std::map < information_key_t, information_key_t > > info_key_send_container_t;

        info_key_send_container_t send_keys_to;
        info_key_recv_container_t recv_keys_from;
        info_key_send_container_t getSendKeys()
        {
            return send_keys_to;
        }
        info_key_send_container_t getRecvKeys()
        {
            return recv_keys_from;
        }
        void clearKeys( )
        {
            send_keys_to.clear();
            recv_keys_from.clear();
        }
        template < typename  I >
        void exchangeInformation(const int tag_info, void ( *work )(void), I & info_manager, MPI_Comm mpi_comm)
        {
            sendRecvKeysExchangeInformationPolicy < exchangeInformationPolicy < K, send_recv_keys_association_trait >, I, typename I::data_style_trait >::exchangeInformation(tag_info,work,info_manager,*this,mpi_comm);
        }
};
//***********************************************************
//      sendRecvKeysExchangeInformationPolicy for homogeneous_data_style_trait
//***********************************************************
template < typename K, typename I  >
class sendRecvKeysExchangeInformationPolicy < K,I, homogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void),I & info_manager,K &key_setting,MPI_Comm mpi_comm)
        {
            typedef typename I::information_t information_t;
            int mpi_rank;
            MPI_Comm_rank(mpi_comm, &mpi_rank);
            MPI_Status status;

            const int nb_send_message = key_setting.send_keys_to.size();
            const int nb_recv_message = key_setting.recv_keys_from.size();
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);
            std::vector < MPI_Request > request_from (nb_recv_message, MPI_REQUEST_NULL);
            const size_t size_info = sizeof( information_t );
            // storage space to pack information to send. indexed by target processor.
            std::map < int, std::vector < information_t >  > send_buffers_to;
            // for a priori known pattern
            // storage space to receive information. indexed by target processor.
            std::map < int, std::vector < information_t > > recv_buffers_from;
            int recv_req_count = 0;
            // posting receive for each proc we are expecting message from
            for (const auto &recv_keys : key_setting.recv_keys_from)
            {
                const int recv_from = recv_keys.first;
                const size_t recv_size = recv_keys.second.size();

                //  std::cout << "PROC " << mpi_rank << " Prepare to receive " << recv_size << " Data " << " from proc "  << recv_from << std::endl;

                auto & recv_buff = recv_buffers_from[recv_from];
                recv_buff.resize( recv_size );
                MPI_Irecv( recv_buff.data(), recv_size*size_info, MPI_BYTE, recv_from, tag_info, mpi_comm, &request_from[recv_req_count++]);
            }
            int send_req_count = 0;
            for (const auto &send_keys : key_setting.send_keys_to)
            {
                const int send_to = send_keys.first;
                const size_t send_size = send_keys.second.size();
                auto & send_buff = send_buffers_to[send_to];
                send_buff.reserve(send_size);
                //  std::cout << "PROC " << mpi_rank << " Prepare to send " << send_size << " Data " << " to  proc "  << send_to << std::endl;
                for ( const auto &key_lr : send_keys.second) send_buff.push_back(info_manager.getInfo(key_lr.second, send_to) );
                MPI_Isend(send_buff.data(), send_size*size_info, MPI_BYTE, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);


            }
            // overlaping computation and communication
            if (work) work();
            int index;
            MPI_Waitany(nb_recv_message, request_from.data(), &index, &status);

            while (index != MPI_UNDEFINED)
            {
                assert(tag_info == status.MPI_TAG);
                const int from = status.MPI_SOURCE;
                auto & recv_keys = key_setting.recv_keys_from.at(from);
                auto & recv_buff = recv_buffers_from.at(from);
                int ikey = 0;
                for (const auto &rkey : recv_keys  ) info_manager.setInfo( rkey, recv_buff[ikey++], from);
                MPI_Waitany(nb_recv_message, request_from.data(), &index, &status);


            }
            // check every thing been sent (you don't want to exit the function before that ... You might destroy data that where not already in send buffer)
            MPI_Waitall(nb_send_message, request_to.data(), MPI_STATUS_IGNORE);

        }
};
//***********************************************************
//      sendRecvKeysExchangeInformationPolicy for nonhomogeneous_data_style_trait
//***********************************************************
template < typename K, typename I  >
class sendRecvKeysExchangeInformationPolicy < K,I, nonhomogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void),I & info_manager,K &key_setting,MPI_Comm mpi_comm)
        {
            int mpi_rank;
            MPI_Comm_rank(mpi_comm, &mpi_rank);
            MPI_Status status;
            const int nb_send_message = key_setting.send_keys_to.size();
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);
            const size_t size_info = info_manager.getApproxDataSize();
            // tag : arbitrary.
            // storage space to pack information to send. indexed by target processor.
            std::map < int, xMpiInputBuffer  > send_buffers_to;
            int send_req_count = 0;
            for (const auto &send_keys : key_setting.send_keys_to)
            {
                const int send_to = send_keys.first;
                const size_t send_nb = send_keys.second.size();
                const size_t send_size = send_nb*size_info;
                auto sendbuffit = send_buffers_to.emplace( send_to, xMpiInputBuffer(mpi_comm, send_size) );
                xMpiInputBuffer & send_buff = ( *( sendbuffit.first )).second;

                for ( const auto &key_lr : send_keys.second) info_manager.getInfo(key_lr.second,  send_buff,  send_to);
                MPI_Isend(send_buff.data(),send_buff.size(), MPI_PACKED, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);
            }
            std::set < int > recved_from;
            std::list < std::pair < int, xMpiOutputBuffer > > recv_buffers_to_treat;
            //Receiving Loop  : probing for  message, if message ready and conform, receive it in outbuff and add the out buff into the queue for setInfo treatment.
            while ( ( key_setting.recv_keys_from.size() != recved_from.size()) || ( !recv_buffers_to_treat.empty()) )
            {
                int flag;
                MPI_Iprobe(MPI_ANY_SOURCE, tag_info, mpi_comm, &flag, &status);
                if (flag)
                {
                    int from = status.MPI_SOURCE;
                    if (key_setting.recv_keys_from.find(from) == key_setting.recv_keys_from.end())
                        throw xDataExchangerException("Error receiving from unexpected proc "
                                                      + std::to_string(mpi_rank),
                                                      __FILE__,__LINE__,__DATE__,__TIME__);

                    if (recved_from.find(from) != recved_from.end() )
                        throw xDataExchangerException("Error receiving a second time from proc "
                                                      + std::to_string(mpi_rank),
                                                      __FILE__,__LINE__,__DATE__,__TIME__);

                    recved_from.insert(from);
                    int count;
                    MPI_Get_count( &status, MPI_PACKED, &count );
                    // not sure I really need to store the receive buffers ...
                    // I choose to just push them in a list, but store the fact that I received it, to know when I'm done.
                    recv_buffers_to_treat.emplace_back(from, xMpiOutputBuffer(mpi_comm, count));
                    MPI_Recv(recv_buffers_to_treat.back().second.data(), count, MPI_PACKED, from, tag_info,mpi_comm,&status);
                }
                else if (!recv_buffers_to_treat.empty())
                {
                    const xMpiOutputBuffer & buff = recv_buffers_to_treat.front().second;
                    const int from = recv_buffers_to_treat.front().first;
                    for ( auto rkey : key_setting.recv_keys_from.at(from) )
                    {
                        info_manager.setInfo( rkey, buff, from);
                    }
                    recv_buffers_to_treat.pop_front();
                }
                else if (work)
                    work();
            }                                                    // at the end of the loop, all the messages are received and the recv_buffers_to_treat is empty.
                                                                 // Wait for all assynchronus send to be finished.
            MPI_Waitall(nb_send_message, request_to.data(), MPI_STATUS_IGNORE);
        }
};

template < typename K, typename I, typename DATASTYLE >
class sendOnlyKeysExchangeInformationPolicy
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void), I & info_manager, K & key_info_manager, MPI_Comm mpi_comm) {throw; }
};
//***********************************************************
//     exchangeInformationPolicy for send_only_keys_association_trait
//***********************************************************
template < class K >
class exchangeInformationPolicy < K,  send_only_keys_association_trait >
{
    public:
        typedef typename K::information_key_t information_key_t;
        typedef std::map < int, std::unordered_set < information_key_t > >                    info_key_send_container_t;
        typedef bool info_key_recv_container_t;
        info_key_send_container_t send_keys_to;
        info_key_send_container_t getSendKeys()
        {
            return send_keys_to;
        }
        info_key_send_container_t getRecvKeys()
        {
            static_assert(std::true_type::value, "exchangeInformationPolicy with send_only_keys_association_trait has no recv keys; ");
        }
        void clearKeys( )
        {
            send_keys_to.clear();
        }
        template < typename  I >
        void exchangeInformation(const int tag_info, void ( *work )(void), I & info_manager, MPI_Comm mpi_comm)
        {
            sendOnlyKeysExchangeInformationPolicy <  exchangeInformationPolicy < K,  send_only_keys_association_trait >, I, typename I::data_style_trait >::exchangeInformation(tag_info,work,info_manager,*this,mpi_comm);
        }
};
//***********************************************************
//      sendOnlyKeysExchangeInformationPolicy for homogeneous_data_style_trait
//***********************************************************
template < typename K, typename I  >
class sendOnlyKeysExchangeInformationPolicy < K,I, homogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void),I & info_manager,K &key_setting,MPI_Comm mpi_comm)
        {
            typedef typename I::information_t information_t;
            MPI_Status status;
            const int nb_send_message = key_setting.send_keys_to.size();
            const size_t size_info = sizeof( information_t );
            int mpi_size;
            MPI_Comm_size(mpi_comm, &mpi_size);

            // storage space to pack information to send. indexed by target processor.
            std::map < int, std::vector < information_t >  > send_buffers_to;

            //vector of mpi requests to check for completion of each send
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);

#ifdef FORCE_MPI2
            std::vector < int > receive_sizes_sendbuf(mpi_size, 0);
            for (auto i : key_setting.send_keys_to) receive_sizes_sendbuf[i.first] = i.second.size( );
            std::vector < int > receive_sizes_receivebuf(mpi_size, 0);
            MPI_Alltoall(receive_sizes_sendbuf.data(), 1, MPI_INT,  receive_sizes_receivebuf.data(), 1, MPI_INT, mpi_comm );
            int nb_recv_message = 0;
            for (int i = 0; i < mpi_size; ++i)
                if (receive_sizes_receivebuf[i]) ++nb_recv_message;

            // storage space to receive information. indexed by target processor.
            std::map < int, std::vector < information_t > > recv_buffers_from;
            //vector of mpi requests to check for completion of each send and receive
            std::vector < MPI_Request > request_from (nb_recv_message, MPI_REQUEST_NULL);

            int recv_req_count = 0;
            for (int i = 0; i < mpi_size; ++i)
            {
                int recv_size = receive_sizes_receivebuf[i];
                if (recv_size)
                {
                    auto & recv_buff = recv_buffers_from[i];
                    recv_buff.resize( recv_size );
                    MPI_Irecv( recv_buff.data(), recv_size*size_info, MPI_BYTE, i, tag_info, mpi_comm, &request_from[recv_req_count++]);
                }
            }
            int send_req_count = 0;
            for (const auto &send_keys : key_setting.send_keys_to)
            {
                const int send_to = send_keys.first;
                const size_t send_size = send_keys.second.size();
                auto & send_buff = send_buffers_to[send_to];
                send_buff.reserve(send_size);
                for ( const auto &lkey_send : send_keys.second) send_buff.push_back(info_manager.getInfo(lkey_send, send_to) );
                MPI_Isend(send_buff.data(), send_size*size_info, MPI_BYTE, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);
            } // overlaping computation and communication
            if (work) work();
            int index;
            MPI_Waitany(nb_recv_message, request_from.data(), &index, &status);


            while (index != MPI_UNDEFINED)
            {
                assert(tag_info == status.MPI_TAG);
                const int from = status.MPI_SOURCE;
                //auto & recv_k =  recv_keys_from.at(from);
                auto & infos = recv_buffers_from.at(from);
                info_manager.setInfo( infos, from);
                MPI_Waitany(nb_recv_message, request_from.data(), &index, &status);
            }
            // check every thing been sent (you don't want to exit the function before that ... You might destroy data that where not already in send buffer)
            MPI_Waitall(nb_send_message, request_to.data(), MPI_STATUS_IGNORE);
#else

            // receive buffers
            std::vector < information_t > recv_buffer;

            // Loop to send information (synchro send)
            int send_req_count = 0;
            for (const auto &send_keys : key_setting.send_keys_to)
            {
                const int send_to = send_keys.first;
                const size_t send_size = send_keys.second.size();
                auto & send_buff = send_buffers_to[send_to];
                send_buff.reserve(send_size);
                for ( const auto &lkey_send : send_keys.second) send_buff.push_back(info_manager.getInfo(lkey_send, send_to) );
                MPI_Issend((void *) send_buff.data(), send_size*size_info, MPI_BYTE, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);
            }
            assert(send_req_count == nb_send_message);

            // overlaping computation and communication
            if (work) work();

            bool do_test = false;
            int flag,from;
            // Receiving loop
            MPI_Request barrier_request(MPI_REQUEST_NULL);
            do
            {
                // probing pending message without blocking to be able to check barrier achievement (end of loop)
                MPI_Iprobe(MPI_ANY_SOURCE,tag_info,mpi_comm,&flag,&status);

                // there is a pending message : get it and unpack
                if (flag)
                {
                    // get size and then message
                    assert(tag_info == status.MPI_TAG);
                    from = status.MPI_SOURCE;
                    int count_byte = 0;
                    MPI_Get_count( &status, MPI_BYTE, &count_byte );
                    if (count_byte%size_info)
                    {
                        int mpi_rank;
                        MPI_Comm_rank(mpi_comm, &mpi_rank);
                        throw xDataExchangerException("Problem with incoming message size, it is not a multiple of the information size! proc "
                                                      + std::to_string(mpi_rank),
                                                      __FILE__,__LINE__,__DATE__,__TIME__);
                    }
                    size_t recv_message_size = count_byte/size_info;
                    recv_buffer.resize(recv_message_size);
                    MPI_Recv((void *) recv_buffer.data(),count_byte,MPI_BYTE,from,tag_info,mpi_comm,&status);

                    // unpack
                    info_manager.setInfo( recv_buffer, from);
                }

                // Does all sended mesages are arrived ? we test ... we don't want to wait for that as we want to loop
                // so it is a test i.e. non blocking
                // We do this test only if we haven't already entered barrier
                if (!do_test)
                {
#ifdef TEST_ALL
                    MPI_Testall(nb_send_message,&request_to[0],&flag,MPI_STATUS_IGNORE);
                    if (flag)
#else
                    MPI_Testany(nb_send_message,&request_to[0],&from,&flag,MPI_STATUS_IGNORE);
                    if (flag && from == MPI_UNDEFINED )
#endif
                    {
                        // This process has sent all its messages and all receiving process have treated them
                        // so it enters the barrier. We know that it is finished in receiving process
                        // because we use a synchronous send. When all process of the communicator will
                        // have enter the barrier we will know that all message sent have been received
                        // by all target process and we will not need to loop anymore to wait for nothing.
                        MPI_Ibarrier(mpi_comm,&barrier_request);

                        // we can now test if barrier request is achieved and this avoid also to reenter barrier many time
                        // and test if all senders have been tested
                        do_test = true;
                    }
                }

                // reset flag if no test is done to loop again
                flag = 0;

                // This test is done only if at least current proc have entered the barrier. Otherwise it's useless (N.C. idea)
                // If all proc enter the barrier flag will be set to a non zero which will finish the loop
                if (do_test)
                    MPI_Test(&barrier_request,&flag,MPI_STATUS_IGNORE);

            } while (!flag);


#endif
        }

};

//***********************************************************
//      sendOnlyKeysExchangeInformationPolicy for nonhomogeneous_data_style_trait
//***********************************************************
template < typename K, typename I  >
class sendOnlyKeysExchangeInformationPolicy < K,I, nonhomogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void),I & info_manager,K &key_setting,MPI_Comm mpi_comm)
        {
            int mpi_rank, mpi_size;
            MPI_Comm_rank(mpi_comm, &mpi_rank);
            MPI_Comm_size(mpi_comm, &mpi_size);
            MPI_Status status;
            const int nb_send_message = key_setting.send_keys_to.size();
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);
            const size_t size_info = info_manager.getApproxDataSize();
            // each send to remote if it has a message .
            std::vector < int > nbkeys_send_to(mpi_size, 0);
            // storage space to pack information to send. indexed by target processor.
            std::map < int, xMpiInputBuffer  > send_buffers_to;

#ifdef FORCE_MPI2
            std::vector < int > nbkeys_recv_from(mpi_size, 0);
            for (auto i : key_setting.send_keys_to) nbkeys_send_to[i.first] = i.second.size( );
            MPI_Alltoall(nbkeys_send_to.data(), 1, MPI_INT,  nbkeys_recv_from.data(), 1, MPI_INT, mpi_comm );
            int nb_messages_to_recv = 0;
            for (int i = 0; i < mpi_size; ++i)
                if (nbkeys_recv_from[i]) ++nb_messages_to_recv;
#endif

            int send_req_count = 0;
            for (const auto &to_send_keys_pair : key_setting.send_keys_to)
            {
                const int send_to = to_send_keys_pair.first;
                const auto & send_keys = to_send_keys_pair.second;
                const size_t send_nb = send_keys.size();
                const size_t send_size = send_nb*size_info;
                auto sendbuffit = send_buffers_to.emplace( send_to, xMpiInputBuffer(mpi_comm, send_size) );
                xMpiInputBuffer & send_buff = ( *( sendbuffit.first )).second;

                for ( const auto &key : send_keys) info_manager.getInfo(key,  send_buff,  send_to);
#ifdef FORCE_MPI2
                MPI_Isend(send_buff.data(),send_buff.size(), MPI_PACKED, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);
#else
                MPI_Issend((void *) send_buff.data(),send_buff.size(),MPI_PACKED, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);
#endif
            }
            assert(send_req_count == nb_send_message);

            // overlapping computation and communication
            if (work) work();

#ifdef FORCE_MPI2
            //std::set<int > recved_from;
            //	std::list< std::pair <int, xMpiOutputBuffer > > recv_buffers_to_treat;
            //Receiving Loop  : probing for  message,  receive it in outbuff and add the out buff into the queue for setInfo treatment.
            int nb_received_message = 0;
            //std::cout << "Proc " << mpi_rank << " will receive "<<  nb_messages_to_recv << std::endl;
            while (nb_received_message < nb_messages_to_recv  )
            {
                //  int flag;
                MPI_Probe(MPI_ANY_SOURCE, tag_info, mpi_comm, &status);
                int from = status.MPI_SOURCE;
                int count;
                MPI_Get_count( &status, MPI_PACKED, &count );
                xMpiOutputBuffer recv_buff(mpi_comm, count);
                MPI_Recv(recv_buff.data(), count, MPI_PACKED, from, tag_info, mpi_comm, &status);
                info_manager.setInfo( recv_buff, from);
                ++nb_received_message;
                // std::cout  << "Proc " << mpi_rank << " has received one message from "<< from   << std::endl;
            }
            //std::cout <<  "Proc " << mpi_rank << "has received all messages" << std::endl;
            MPI_Waitall(nb_send_message, request_to.data(), MPI_STATUS_IGNORE);
            //std::cout <<  "Proc " << mpi_rank << "All messages have been send" << std::endl;
#else
            bool do_test = false;
            int flag,from;
            // Receiving loop
            MPI_Request barrier_request(MPI_REQUEST_NULL);
            do
            {
                // probing pending message without blocking to be able to check barrier achievement (end of loop)
                MPI_Iprobe(MPI_ANY_SOURCE,tag_info,mpi_comm,&flag,&status);

                // there is a pending message : get it and unpack
                if (flag)
                {
                    // get size and then message
                    assert(tag_info == status.MPI_TAG);
                    from = status.MPI_SOURCE;
                    int count = 0;
                    MPI_Get_count( &status, MPI_PACKED, &count);
                    xMpiOutputBuffer recv_buff(mpi_comm, count);
                    MPI_Recv((void *) recv_buff.data(), count, MPI_PACKED, from, tag_info, mpi_comm, &status);

                    // unpack
                    info_manager.setInfo( recv_buff, from);
                }

                // Does all sended mesages are arrived ? we test ... we don't want to wait for that as we want to loop
                // so it is a test i.e. non blocking
                // We do this test only if we haven't already entered barrier
                if (!do_test)
                {
#ifdef TEST_ALL
                    MPI_Testall(nb_send_message,&request_to[0],&flag,MPI_STATUS_IGNORE);
                    if (flag)
#else
                    MPI_Testany(nb_send_message,&request_to[0],&from,&flag,MPI_STATUS_IGNORE);
                    if (flag && from == MPI_UNDEFINED)
#endif
                    {
                        // This process has sent all its messages and all receiving process have treated them
                        // so it enters the barrier. We know that it is finished in receiving process
                        // because we use a synchronous send. When all process of the communicator will
                        // have enter the barrier we will know that all message sent have been received
                        // by all target process and we will not need to loop anymore to wait for nothing.
                        MPI_Ibarrier(mpi_comm,&barrier_request);

                        // we can now test if barrier request is achieved and this avoid also to reenter barrier many time
                        do_test = true;
                    }
                }

                // reset flag if no test is done to loop again
                flag = 0;

                // This test is done only if at least current proc have entered the barrier. Otherwise it's useless (N.C. idea)
                // If all proc enter the barrier flag will be set to a non zero which will finish the loop
                if (do_test)
                    MPI_Test(&barrier_request,&flag,MPI_STATUS_IGNORE);

            } while (!flag);
#endif
        }
};

template  < typename K, typename I, typename DATASTYLE >
class recvOnlyKeysExchangeInformationPolicy
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void), I & info_manager, K & key_info_manager, MPI_Comm mpi_comm) {throw; }
};
//***********************************************************
//     exchangeInformationPolicy for recv_only_keys_association_trait
//***********************************************************
template < class K >
class exchangeInformationPolicy < K,  recv_only_keys_association_trait >
{
    public:
        typedef typename K::information_key_t information_key_t;
        typedef std::map < int, std::unordered_set < information_key_t > >                    info_key_recv_container_t;
        typedef bool info_key_send_container_t;

        info_key_recv_container_t recv_keys_from;

        void clearKeys( )
        {
            recv_keys_from.clear();
        }
        template < typename  I >
        void exchangeInformation(const int tag_info, void ( *work )(void), I & info_manager, MPI_Comm mpi_comm)
        {
            recvOnlyKeysExchangeInformationPolicy < exchangeInformationPolicy < K,  recv_only_keys_association_trait >, I, typename I::data_style_trait >::exchangeInformation(tag_info,work,info_manager,*this,mpi_comm);
        }
};
//***********************************************************
//      recvOnlyKeysExchangeInformationPolicy for homogeneous_data_style_trait
//***********************************************************
template < typename K, typename I  >
class recvOnlyKeysExchangeInformationPolicy < K,I, homogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(void),I & info_manager,K &key_setting,MPI_Comm mpi_comm)
        {
            typedef typename I::information_t information_t;
            MPI_Status status;
            const int nb_recv_message = key_setting.recv_keys_from.size();
            const size_t size_info = sizeof( information_t );
            int mpi_size,  mpi_rank;
            MPI_Comm_size(mpi_comm, &mpi_size);
            MPI_Comm_rank(mpi_comm, &mpi_rank);

            std::map < int, std::vector < information_t > > recv_buffers_from;
            std::vector < MPI_Request > request_from (nb_recv_message, MPI_REQUEST_NULL);

            int recv_req_count = 0;
            for (auto recv_keys : key_setting.recv_keys_from)
            {
                const int from = recv_keys.first;
                const int recv_size = recv_keys.second.size();
                auto & recv_buff = recv_buffers_from[from];
                recv_buff.resize( recv_size );
                //  std::cout << "P" << mpi_rank <<  "Post a receive from " << from << " of size " << recv_size << std::endl;
                MPI_Irecv( recv_buff.data(), recv_size*size_info, MPI_BYTE, from, tag_info, mpi_comm, &request_from[recv_req_count++]);
            }

            std::map < int, std::vector < information_t > > send_buffers_to;
            for (int send_to = 0; send_to < mpi_size; ++send_to)
            {
                std::vector < information_t >  send_buff = info_manager.getInfo( send_to);
                if (!send_buff.empty()) send_buffers_to.emplace(send_to, std::move(send_buff) );
                //	  if (!send_buff.empty()) send_buffers_to[send_to] = send_buff;
            }
            const int nb_send_message = send_buffers_to.size();
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);
            int send_req_count = 0;
            //	std::cout << "P" << mpi_rank <<  " Will Post a send to "<< send_buffers_to.size() << " Remotes procs" << std::endl;
            for (const auto &to_send_buff_pair : send_buffers_to )
            {
                const int send_to = to_send_buff_pair.first;
                const auto & send_buff = to_send_buff_pair.second;
                const int send_size = send_buff.size();
                //  std::cout << "P" << mpi_rank <<  " Post a send to " << send_to << " of size " << send_buff.size() << " : " << std::endl;
                // for (information_t data : send_buff) std::cout << " " << data ;
                //std::cout << std::endl;
                MPI_Isend( const_cast <  information_t  * >( send_buff.data()), send_size*size_info, MPI_BYTE, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);
            }

            int index;
            MPI_Waitany(nb_recv_message, request_from.data(), &index, &status);
            while (index != MPI_UNDEFINED)
            {
                assert(tag_info == status.MPI_TAG);
                const int from = status.MPI_SOURCE;

                auto & recv_keys = key_setting.recv_keys_from.at(from);
                auto & recv_buff = recv_buffers_from.at(from);
                //std::cout << " P" << mpi_rank <<  " Received message from " << from << " Of size " <<recv_buff.size() << " " << recv_keys.size() << " And start to set info from it" <<  std::endl;
                int ikey = 0;
                for (const auto &rkey : recv_keys  ) info_manager.setInfo( rkey, recv_buff[ikey++], from);

                MPI_Waitany(nb_recv_message, request_from.data(), &index, &status);
            }
            // check every thing has been sent (you don't want to exit the function before that ... You might destroy data that where not already in send buffer)
            MPI_Waitall(nb_send_message, request_to.data(), MPI_STATUS_IGNORE);
        }
};


}                                                    // end namespace

namespace xtool
{
//***********************************************************
// xDataExchangerException Implementation
//***********************************************************

inline xDataExchangerException::xDataExchangerException(std::string message, std::string file,  int line, std::string date, std::string time)
{
    std::ostringstream oss;
    oss << " In file "<< file << " line " << line << " compiled "<<date<<" at "<<time<<std::endl;
    oss << "xDataExchangerException : "<< message << std::endl;
    msg = oss.str();
}
/// general exception object : destructor
inline xDataExchangerException :: ~xDataExchangerException() throw( ) {}

/// mandatory what method
inline const char * xDataExchangerException::what() const throw( )
{
    return this->msg.c_str();
}

}                                                    // end namespace

#endif
