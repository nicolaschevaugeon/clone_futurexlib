/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef XDATAEXCHANGER_IMP_H
#define XDATAEXCHANGER_IMP_H

namespace xtool
{

//***********************************************************
// helping function
//***********************************************************
/// Function to get an object out of an object. By itself it doesn't look very helpful. But with overload version with pointer it helps
/// getting object out iterator used in accumulate methods
template < typename OBJECTYPE >
inline OBJECTYPE & getObject(OBJECTYPE &o) { return o; }
/// Function to get object out of pointer argument. By itself it doesn't look very helpful. But with overload version with reference it helps
/// getting object out iterator used in accumulate methods
template < typename OBJECTYPE >
inline OBJECTYPE & getObject(OBJECTYPE *o) { return *o; }

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
    int tmp = MPI_Unpack( const_cast < void * >( data()),  dat.size(), &position, outbuff, outcount, datatype,  comm  );
    /*if (datatype == MPI_AINT){
       std::cout << "IN un pack "<< outbuff << std::endl;
       }*/
    return tmp;


}
const void *xMpiOutputBuffer::data() const {return dat.data(); }
void *xMpiOutputBuffer::data() {return dat.data(); }


//***********************************************************
// xKeyContainerSendAndRecv Implementation
//***********************************************************
template < typename K >
xKeyContainerSendAndRecv < K >::xKeyContainerSendAndRecv(  MPI_Comm _mpi_comm) : mpi_comm(_mpi_comm)
{
    // getting mpi_rank and mpi_size in mpi_comm
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
}

template < typename K >
void xKeyContainerSendAndRecv < K >::clearKeys( )
{
    send_keys_to.clear();
    recv_keys_from.clear();
}

template < typename K >
MPI_Comm xKeyContainerSendAndRecv < K >::getComm() { return mpi_comm; }

template < typename K >
template < typename ITER, typename KM >
void xKeyContainerSendAndRecv < K >::accumulateKeysOwnerScatter( ITER itb, ITER ite, KM & key_manager )
{
    static_assert(std::is_same < K, typename KM::information_key_t >::value, " The key type of your key_container is not the same as the one of the key_manager" );
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;


    // loop
    for (ITER it = itb; it != ite; ++it)
    {
        const auto &lo = getObject(*it);       // lo is the reference to the local object
        const information_key_t lk = key_manager.localObjectKey(lo);     // lk is the key to the local object
        const auto po = key_manager.getConstPartitionObject(lo);     // po is the partition object associated to lo


        if (po.hasRemoteObject())
        {
            if (!po.isOwner())
            {
                const auto ro_owner = po.getOwner(); //ro_owner is the remote object on the process that own the object.
                recv_keys_from[ro_owner.getProcId()].insert(lk);
            }
            else
            {
                for (const auto ro : po.getRemoteObjectsCollectionRange( ))
                {
                    // ro is a remoteObject associated to lo.
                    const information_key_t rk = key_manager.remoteObjectKey(ro, lo);
                    send_keys_to[ro.getProcId()].insert( std::make_pair( rk,  lk ));
                }
            }
        }
    }
    return;
}

template < typename K >
template < typename ITER, typename KM >
void xKeyContainerSendAndRecv < K >::accumulateKeysOwnerGather( ITER itb, ITER ite, KM & key_manager )
{
    static_assert(std::is_same < K, typename KM::information_key_t >::value, " The key type of your key_container is not the same as the one of the key_manager" );
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;

    // loop on iterator to collect information key and put them in correct order
    for (ITER it = itb; it != ite; ++it)
    {
        const auto &lo = getObject(*it);     // lo is the reference to the local object
        const information_key_t lk = key_manager.localObjectKey(lo);   // lk is the key to the local object
        const auto po = key_manager.getConstPartitionObject( lo);               // po is the partition object associated to po
        if (po.hasRemoteObject())
        {
            if (po.isOwner())
            {
                // ro is a remoteObject associated to lo.
                for (auto ro : po.getRemoteObjectsCollectionRange( )) recv_keys_from[ro.getProcId()].insert(lk);
            }
            else
            {
                const auto ro = po.getOwner(); //ro is the remote object on the process that own the object.
                const information_key_t rk = key_manager.remoteObjectKey(ro, lo);
                send_keys_to[ro.getProcId()].insert ( std::make_pair( rk,  lk ) );
            }
        }
    }
    return;
}

template < typename K >
template < typename ITER, typename KM >
void xKeyContainerSendAndRecv < K >::accumulateKeysAllGather( ITER itb, ITER ite, KM & key_manager)
{
    static_assert(std::is_same < K, typename KM::information_key_t >::value, " The key type of your key_container is not the same as the one of the key_manager" );
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;

    // loop on iterator to collect information key and put them in correct order
    for (ITER it = itb; it != ite; ++it)
    {
        const auto &lo = getObject(*it);           // lo is the reference to the local object
        const information_key_t lk = key_manager.localObjectKey(lo);         // lk is the key to the local object
        const auto po = key_manager.getConstPartitionObject( lo);                     // po is the partition object associated to po
        if (po.hasRemoteObject())
        {
            // ro is a remoteObject associated to lo.
            for (auto ro : po.getRemoteObjectsCollectionRange( ))
            {
                recv_keys_from[ro.getProcId()].insert(lk);
                const information_key_t rk = key_manager.remoteObjectKey(ro, lo);
                send_keys_to[ro.getProcId()].insert ( std::make_pair( rk,  lk ) );
            }
        }
    }
    return;
}

//***********************************************************
// xKeyContainerSendAndRecvFromTo Implementation
//***********************************************************
template < typename K >
xKeyContainerSendAndRecvFromTo < K >::xKeyContainerSendAndRecvFromTo( int from_, int to_, MPI_Comm _mpi_comm) : mpi_comm(_mpi_comm),from(from_),to(to_)
{
    assert(from != to);
    // getting mpi_rank and mpi_size in mpi_comm
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
}

template < typename K >
void xKeyContainerSendAndRecvFromTo < K >::clearKeys( )
{
    send_keys_to.clear();
    recv_keys_from.clear();
}

template < typename K >
MPI_Comm xKeyContainerSendAndRecvFromTo < K >::getComm() { return mpi_comm; }



template < typename K >
template < typename ITER, typename KM >
void xKeyContainerSendAndRecvFromTo < K >::accumulateKeys( ITER itb, ITER ite, KM & key_manager)
{
    static_assert(std::is_same < K, typename KM::information_key_t >::value, " The key type of your key_container is not the same as the one of the key_manager" );
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;
    // if not part of the "from to" condition
    if (mpi_rank != from && mpi_rank != to) return;
    // if it is the from of "from to" condition
    else if (mpi_rank == from )
    {
        // loop on iterator to collect information key and put them in correct order
        for (ITER it = itb; it != ite; ++it)
        {
            const auto &lo = getObject(*it);       // lo is the reference to the local object
            const information_key_t lk = key_manager.localObjectKey(lo);     // lk is the key to the local object
            const auto po = key_manager.getConstPartitionObject( lo);                 // po is the partition object associated to po
            if (po.hasRemoteObject())
            {
                // ro is a remoteObject associated to lo.
                for (auto ro : po.getRemoteObjectsCollectionRange( ))
                {
                    if (ro.getProcId() == to)
                    {
                        const information_key_t rk = key_manager.remoteObjectKey(ro, lo);
                        send_keys_to[to].insert ( std::make_pair( rk,  lk ) );
                    }
                }
            }
        }
    }
    // if it is the to of "from to" condition
    // Note: no test as from#to and pid#from so pid==to
    else
    {
        // loop on iterator to collect information key and put them in correct order
        for (ITER it = itb; it != ite; ++it)
        {
            const auto &lo = getObject(*it);       // lo is the reference to the local object
            const information_key_t lk = key_manager.localObjectKey(lo);     // lk is the key to the local object
            const auto po = key_manager.getConstPartitionObject( lo);                 // po is the partition object associated to po
            if (po.hasRemoteObject())
            {
                // ro is a remoteObject associated to lo.
                for (auto ro : po.getRemoteObjectsCollectionRange( ))
                {
                    if (ro.getProcId() == from)
                        recv_keys_from[from].insert(lk);
                }
            }
        }
    }
    return;
}
//***********************************************************
// xKeyContainerSendOrRecv Implementation
//***********************************************************
template < typename K, typename HASHFCN, typename EQUALKEY >
xKeyContainerSendOrRecv < K, HASHFCN, EQUALKEY >::xKeyContainerSendOrRecv(  MPI_Comm _mpi_comm) : mpi_comm(_mpi_comm)
{
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
}

template < typename K, typename HASHFCN, typename EQUALKEY >
void xKeyContainerSendOrRecv < K, HASHFCN, EQUALKEY >::clearKeys( ) { keys_from_to.clear(); }

template < typename K, typename HASHFCN, typename EQUALKEY >
MPI_Comm xKeyContainerSendOrRecv < K, HASHFCN, EQUALKEY >::getComm() { return mpi_comm; }

template < typename K, typename HASHFCN, typename EQUALKEY >
template < typename ITER, typename KM >
void xKeyContainerSendOrRecv < K, HASHFCN, EQUALKEY >::accumulateKeys(ITER itb, ITER ite, KM & key_manager)
{
    static_assert(std::is_same < K, typename KM::information_key_t >::value, " The key type of your key_container is not the same as the one of the key_manager" );
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;
    for (ITER it = itb; it != ite; ++it)
    {
        accumulateKeysPrivate(getObject(*it), key_manager);
    }
}
template < typename K, typename HASHFCN, typename EQUALKEY >
template < typename OBJECTYPE, typename KM >
void xKeyContainerSendOrRecv < K, HASHFCN, EQUALKEY >::accumulateKeys( const OBJECTYPE & lo, KM & key_manager)
{
    static_assert(std::is_same < K, typename KM::information_key_t >::value, " The key type of your key_container is not the same as the one of the key_manager" );
    // if there is less then 2 proc we don't need to exchange anything
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;

    // do work for lo
    accumulateKeysPrivate(lo,key_manager);
}

template < typename K, typename HASHFCN, typename EQUALKEY >
template < typename OBJECTYPE, typename KM >
void xKeyContainerSendOrRecv < K, HASHFCN, EQUALKEY >::accumulateKeysPrivate( const OBJECTYPE & lo,  KM  &key_manager)
{
    static_assert(std::is_same < K, typename KM::information_key_t >::value, " The key type of your key_container is not the same as the one of the key_manager" );
    const information_key_t lk = key_manager.localObjectKey(lo);         // lk is the key to the local object key
    std::set < int > procs = key_manager.getMessageRanks(lk);
    for (auto proc : procs)
    {
        keys_from_to[proc].insert(lk);
    }
}

//***********************************************************
//     exchangeInformationPolicy non specialized class
//***********************************************************
template < typename C, typename IM, typename COMMPATTERN, typename DATASTYLE >
class exchangeInformationPolicy
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(), C & key_info_container, IM & info_manager, MPI_Comm mpi_comm) {throw; }
};
//***********************************************************
//     exchangeInformation function
//***********************************************************
template < typename C, typename  IM >
void exchangeInformation( C & key_info_container, IM & info_manager,  void ( *work )(void))
{
    // usage checking
    static_assert(
        ( std::is_same < typename IM::communication_trait, send_and_recv_keys_communication_trait >::value ) ?
        ( std::is_same < typename C::keys_association_trait, send_and_recv_keys_association_trait >::value )
        :
        !( std::is_same < typename C::keys_association_trait, send_and_recv_keys_association_trait >::value )
        , "  mixing usage of traits !\n You must use  xKeyContainerSendOrRecv container with send_only_keys_communication_trait or recv_only_keys_communication_trait\nYou must use  xKeyContainerSendAndRecv container with send_and_recv_keys_communication_trait");

    static_assert(std::is_same < typename C::information_key_t, typename IM::information_key_t >::value,
                  " key must be the same in betwen your 'information manager' and your  'key information container' (i.e. your key information manager) ");

    // if there is less then 2 proc we don't need to exchange anything
    MPI_Comm mpi_comm = key_info_container.getComm();
    int mpi_size;
    MPI_Comm_size(mpi_comm,&mpi_size);
    if (mpi_size < 2) return;
    // if not part of a communication world
    if (mpi_comm == MPI_COMM_NULL) return;
    // arbitrary tag used for this sequence
    const int tag_info = xtool::xMPITag::getNewTag(mpi_comm);
    // exchange sequence
    exchangeInformationPolicy < C,IM,typename IM::communication_trait,typename IM::data_style_trait >::exchangeInformation( tag_info, work, key_info_container, info_manager, mpi_comm );
}

//***********************************************************
//     exchangeInformationPolicy specialized for send_and_recv_keys_communication_trait and homogeneous_data_style_trait
//***********************************************************
template < typename C, typename IM >
class exchangeInformationPolicy < C,IM, send_and_recv_keys_communication_trait, homogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(), C & key_info_container, IM & info_manager, MPI_Comm mpi_comm)
        {
            typedef typename IM::information_t information_t;
            int mpi_rank;
            MPI_Comm_rank(mpi_comm, &mpi_rank);
            MPI_Status status;

            const int nb_send_message = key_info_container.send_keys_to.size();
            const int nb_recv_message = key_info_container.recv_keys_from.size();
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);
            std::vector < MPI_Request > request_from (nb_recv_message, MPI_REQUEST_NULL);
            const size_t size_info = sizeof( information_t );
            // storage space to pack information to send. indexed by target processor.
            std::map < int, std::vector < information_t >  > send_buffers_to;
            // for a priori known pattern
            // storage space to receive information. indexed by source processor.
            std::map < int, std::vector < information_t > > recv_buffers_from;
            int recv_req_count = 0;
            // posting receive for each proc we are expecting message from
            for (const auto &recv_keys : key_info_container.recv_keys_from)
            {
                const int recv_from = recv_keys.first;
                const size_t recv_size = recv_keys.second.size();

                //  std::cout << "PROC " << mpi_rank << " Prepare to receive " << recv_size << " Data " << " from proc "  << recv_from << std::endl;

                auto & recv_buff = recv_buffers_from[recv_from];
                recv_buff.resize( recv_size );
                MPI_Irecv( recv_buff.data(), recv_size*size_info, MPI_BYTE, recv_from, tag_info, mpi_comm, &request_from[recv_req_count++]);
            }
            int send_req_count = 0;
            for (const auto &send_keys : key_info_container.send_keys_to)
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
                auto & recv_keys = key_info_container.recv_keys_from.at(from);
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
//     exchangeInformationPolicy specialized for send_and_recv_keys_communication_trait and nonhomogeneous_data_style_trait
//***********************************************************
template < typename C, typename IM >
class exchangeInformationPolicy < C,IM, send_and_recv_keys_communication_trait, nonhomogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(), C & key_info_container, IM & info_manager, MPI_Comm mpi_comm)
        {
            int mpi_rank;
            MPI_Comm_rank(mpi_comm, &mpi_rank);
            MPI_Status status;
            const int nb_send_message = key_info_container.send_keys_to.size();
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);
            const size_t size_info = info_manager.getApproxDataSize();
            // tag : arbitrary.
            // storage space to pack information to send. indexed by target processor.
            std::map < int, xMpiInputBuffer  > send_buffers_to;
            int send_req_count = 0;
            for (const auto &send_keys : key_info_container.send_keys_to)
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
            while ( ( key_info_container.recv_keys_from.size() != recved_from.size()) || ( !recv_buffers_to_treat.empty()) )
            {
                int flag;
                MPI_Iprobe(MPI_ANY_SOURCE, tag_info, mpi_comm, &flag, &status);
                if (flag)
                {
                    int from = status.MPI_SOURCE;
                    if (key_info_container.recv_keys_from.find(from) == key_info_container.recv_keys_from.end())
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
                    for ( auto rkey : key_info_container.recv_keys_from.at(from) )
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

//***********************************************************
//     exchangeInformationPolicy specialized for send_only_keys_communication_trait and homogeneous_data_style_trait
//***********************************************************
template < typename C, typename IM >
class exchangeInformationPolicy < C,IM, send_only_keys_communication_trait, homogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(), C & key_info_container, IM & info_manager, MPI_Comm mpi_comm)
        {
            typedef typename IM::information_t information_t;
            MPI_Status status;
            const int nb_send_message = key_info_container.keys_from_to.size();
            const size_t size_info = sizeof( information_t );
            int mpi_size;
            MPI_Comm_size(mpi_comm, &mpi_size);

            // storage space to pack information to send. indexed by target processor.
            std::map < int, std::vector < information_t >  > send_buffers_to;

            //vector of mpi requests to check for completion of each send
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);

#ifdef FORCE_MPI2
            std::vector < int > receive_sizes_sendbuf(mpi_size, 0);
            for (auto i : key_info_container.keys_from_to) receive_sizes_sendbuf[i.first] = i.second.size( );
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
            for (const auto &send_keys : key_info_container.keys_from_to)
            {
                const int send_to = send_keys.first;
                const size_t send_size = send_keys.second.size();
                auto & send_buff = send_buffers_to[send_to];
                send_buff.reserve(send_size);
                for ( const auto &lkey_send : send_keys.second) send_buff.push_back(info_manager.getInfo(lkey_send, send_to) );
                MPI_Isend(send_buff.data(), send_size*size_info, MPI_BYTE, send_to, tag_info, mpi_comm, &request_to[send_req_count++]);
            }     // overlaping computation and communication
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
            for (const auto &send_keys : key_info_container.keys_from_to)
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
//     exchangeInformationPolicy specialized for send_only_keys_communication_trait and nonhomogeneous_data_style_trait
//***********************************************************
template < typename C, typename IM >
class exchangeInformationPolicy < C,IM, send_only_keys_communication_trait, nonhomogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(), C & key_info_container, IM & info_manager, MPI_Comm mpi_comm)
        {
            int mpi_rank, mpi_size;
            MPI_Comm_rank(mpi_comm, &mpi_rank);
            MPI_Comm_size(mpi_comm, &mpi_size);
            MPI_Status status;
            const int nb_send_message = key_info_container.keys_from_to.size();
            std::vector < MPI_Request > request_to   (nb_send_message, MPI_REQUEST_NULL);
            const size_t size_info = info_manager.getApproxDataSize();
            // each send to remote if it has a message .
            std::vector < int > nbkeys_send_to(mpi_size, 0);
            // storage space to pack information to send. indexed by target processor.
            std::map < int, xMpiInputBuffer  > send_buffers_to;

#ifdef FORCE_MPI2
            std::vector < int > nbkeys_recv_from(mpi_size, 0);
            for (auto i : key_info_container.keys_from_to) nbkeys_send_to[i.first] = i.second.size( );
            MPI_Alltoall(nbkeys_send_to.data(), 1, MPI_INT,  nbkeys_recv_from.data(), 1, MPI_INT, mpi_comm );
            int nb_messages_to_recv = 0;
            for (int i = 0; i < mpi_size; ++i)
                if (nbkeys_recv_from[i]) ++nb_messages_to_recv;
#endif

            int send_req_count = 0;
            for (const auto &to_send_keys_pair : key_info_container.keys_from_to)
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

//***********************************************************
//     exchangeInformationPolicy specialized for recv_only_keys_communication_trait and homogeneous_data_style_trait
//***********************************************************
template < typename C, typename IM >
class exchangeInformationPolicy < C,IM, recv_only_keys_communication_trait, homogeneous_data_style_trait >
{
    public:
        static void exchangeInformation(const int tag_info, void ( *work )(), C & key_info_container, IM & info_manager, MPI_Comm mpi_comm)
        {
            typedef typename IM::information_t information_t;
            MPI_Status status;
            const int nb_recv_message = key_info_container.keys_from_to.size();
            const size_t size_info = sizeof( information_t );
            int mpi_size,  mpi_rank;
            MPI_Comm_size(mpi_comm, &mpi_size);
            MPI_Comm_rank(mpi_comm, &mpi_rank);

            std::map < int, std::vector < information_t > > recv_buffers_from;
            std::vector < MPI_Request > request_from (nb_recv_message, MPI_REQUEST_NULL);

            int recv_req_count = 0;
            for (auto recv_keys : key_info_container.keys_from_to)
            {
                const int from = recv_keys.first;
                const int recv_size = recv_keys.second.size();
                auto & recv_buff = recv_buffers_from[from];
                recv_buff.resize( recv_size );
                // std::cout << "P" << mpi_rank <<  "Post a receive from " << from << " of size " << recv_size << std::endl;
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
                // std::cout << "P" << mpi_rank <<  " Post a send to " << send_to << " of size " << send_buff.size() << " : " << std::endl;
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

                auto & recv_keys = key_info_container.keys_from_to.at(from);
                auto & recv_buff = recv_buffers_from.at(from);
                // std::cout << " P" << mpi_rank <<  " Received message from " << from << " Of size " <<recv_buff.size() << " " << recv_keys.size() << " And start to set info from it" <<  std::endl;
                int ikey = 0;
                for (const auto &rkey : recv_keys  ) info_manager.setInfo( rkey, recv_buff[ikey++], from);

                MPI_Waitany(nb_recv_message, request_from.data(), &index, &status);
            }
            // check every thing has been sent (you don't want to exit the function before that ... You might destroy data that where not already in send buffer)
            MPI_Waitall(nb_send_message, request_to.data(), MPI_STATUS_IGNORE);
        }
};


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
inline xDataExchangerException :: ~xDataExchangerException() throw( ) = default;

/// mandatory what method
inline const char * xDataExchangerException::what() const throw( )
{
    return this->msg.c_str();
}

}                                                    // end namespace

#endif
