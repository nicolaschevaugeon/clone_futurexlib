/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef XCOMMEXCHANGERIMP_H
#define XCOMMEXCHANGERIMP_H

namespace distmesh
{

template < typename I >
xCommExchanger < I >::xCommExchanger(void) : world(MPI_COMM_NULL),proc_id(-1),nb_proc(-1),exchange_type(0),info_manager(NULL) {}

template < typename I >
void xCommExchanger < I >::initIndex( I & info_manager_, const MPI_Comm &world_)
{
    info_manager = &info_manager_;
    send.clear();
    receive.clear();
    world = world_;
    exchange_type = 0;
    int err_mpi;
    err_mpi = MPI_Comm_rank(world,&proc_id);
    // getting rank in world
    if ( err_mpi != MPI_SUCCESS)
    {
        std::ostringstream oss;
        oss << "Problem with MPI_Comm_rank !\n error="<<err_mpi<<std::endl;
        throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }
    // getting size of world
    err_mpi = MPI_Comm_size(world,&nb_proc);
    if ( err_mpi != MPI_SUCCESS)
    {
        std::ostringstream oss;
        oss << "Problem with MPI_Comm_size !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
        throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }

    // if there is less then 2 proc we don't need to exchange anything
    if (nb_proc < 2) return;

    // resizing index container
    send.resize(nb_proc);
    receive.resize(nb_proc);
}

template < typename I >
template < typename ITER >
void xCommExchanger < I >::acumulateIndexOwnerOnly( ITER itb, ITER ite )
{
    // if there is less then 2 proc we don't need to exchange anything
    if (nb_proc < 2) return;
    // if not part of a communication world
    if (world == MPI_COMM_NULL) return;
    assert(info_manager);
    // set exchange type and check if not miss used
    if (exchange_type != 1 && exchange_type != 0)
    {
        std::ostringstream oss;
        oss << "You can not accumulate index information of different kind: proc "<<proc_id<<std::endl;
        throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }
    exchange_type = 1;

    // loop on iterator to collect information key and put them in correct order
    xdistmesh::xCommEntity* ce;
    size_t target_pid;
    for (ITER it = itb; it != ite; ++it)
    {
        typename I::information_key_t ik = info_manager->infKey(it);
        ce = info_manager->getComm(ik);
        if (ce)
        {
            for (int i = 0, m = ce->size(); i < m; ++i)
            {
                void *adress = ( *ce )[i].get(target_pid);
                if (target_pid < proc_id && !i)
                {
                    receive[target_pid].insert (ik);
                    break;
                }
                else if (target_pid > proc_id)
                {
                    send[target_pid].insert ( std::make_pair(info_manager->transformInfKey(ik,adress),ik) );
                }
            }
        }
    }

    return;

}

template < typename I >
template < typename ITER >
void xCommExchanger < I >::acumulateIndexAllToAll( ITER itb, ITER ite )
{
    // if there is less then 2 proc we don't need to exchange anything
    if (nb_proc < 2) return;
    // if not part of a communication world
    if (world == MPI_COMM_NULL) return;
    assert(info_manager);
    // set exchange type and check if not miss used
    if (exchange_type != 2 && exchange_type != 0)
    {
        std::ostringstream oss;
        oss << "You can not accumulate index information of different kind: proc "<<proc_id<<std::endl;
        throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }
    exchange_type = 2;

    // loop on iterator to collect information key and put them in correct order
    xdistmesh::xCommEntity* ce;
    size_t target_pid;
    for (ITER it = itb; it != ite; ++it)
    {
        typename I::information_key_t ik = info_manager->infKey(it);
        ce = info_manager->getComm(ik);
        if (ce)
        {
            for (int i = 0, m = ce->size(); i < m; ++i)
            {
                void *adress = ( *ce )[i].get(target_pid);
                receive[target_pid].insert (ik);
                send[target_pid].insert ( std::make_pair(info_manager->transformInfKey(ik,adress),ik) );
            }
        }
    }

    return;

}

template < typename I >
template < typename ITER >
void xCommExchanger < I >::acumulateIndexOneSide( ITER itb, ITER ite )
{
    // if there is less then 2 proc we don't need to exchange anything
    if (nb_proc < 2) return;
    // if not part of a communication world
    if (world == MPI_COMM_NULL) return;
    assert(info_manager);
    // set exchange type and check if not miss used
    if (exchange_type != 3 && exchange_type != 0)
    {
        std::ostringstream oss;
        oss << "You can not accumulate index information of different kind: proc "<<proc_id<<std::endl;
        throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }
    exchange_type = 3;

    // loop on iterator to collect information key and put them in some order
    xdistmesh::xCommEntity* ce;
    size_t target_pid;
    for (ITER it = itb; it != ite; ++it)
    {
        typename I::information_key_t ik = info_manager->infKey(it);
        ce = info_manager->getComm(ik);
        if (ce)
        {
            for (int i = 0, m = ce->size(); i < m; ++i)
            {
                void *adress = ( *ce )[i].get(target_pid);
                send[target_pid].insert ( std::make_pair(info_manager->transformInfKey(ik,adress),ik) );
            }
        }
    }

    return;

}
template < typename I >
void xCommExchanger < I >::exchangeInformation( void ( *work )(void) )
{
    // if there is less then 2 proc we don't need to exchange anything
    if (nb_proc < 2) return;
    // if not part of a communication world
    if (world == MPI_COMM_NULL) return;

    assert(info_manager);
    assert(exchange_type);

    int err_mpi;
    MPI_Status status;
    int index;
    int flag;
    int count;
    // 'request'
    std::vector < MPI_Request > requestfrom(nb_proc,MPI_REQUEST_NULL);
    std::vector < MPI_Request > requestto(nb_proc,MPI_REQUEST_NULL);
    // size
    size_t size_info = sizeof( typename I::information_t );
    // tag : arbitrary. Ideally it should be unique for each exchange TODO
    int tag_info = size_info/nb_proc;
    // storage vectors to pack information
    std::vector < std::vector < typename I::information_t > > info_to(nb_proc);

    // for a priori known pattern
    if (exchange_type < 3)
    {
        // storage vectors to pack information
        std::vector < std::vector < typename I::information_t > > info_from(nb_proc);

        // loop to start receiving information
        for (int i = 0; i < nb_proc; ++i)
        {
            size_t r = receive[i].size();
            if (r)
            {
                info_from[i].resize(r);

                // start receiving information set in other proc
                err_mpi = MPI_Irecv((void *) &info_from[i][0],r*size_info,MPI_BYTE,i,tag_info,world,&requestfrom[i]);
                if ( err_mpi != MPI_SUCCESS)
                {
                    std::ostringstream oss;
                    oss << "Problem with MPI_Irecv when posting receive from !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }
            }
        }
        // loop to pack and send information where it have to
        for (size_t i = 0; i < nb_proc; ++i)
        {
            size_t s = send[i].size();
            if (s)
            {
                info_to.reserve(s);
                for (auto iti = send[i].begin(),itie = send[i].end(); iti != itie; ++iti)
                {
                    info_to[i].push_back(info_manager->getInfo(iti->second,i));
                }
                err_mpi = MPI_Isend((void *) &info_to[i][0],s*size_info,MPI_BYTE,i,tag_info,world,&requestto[i]);
                if ( err_mpi != MPI_SUCCESS)
                {
                    std::ostringstream oss;
                    oss << "Problem with MPI_Isend when posting send to !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }
            }
        }

        // overlaping computation and comunication
        if (work) work();


        // loop to treat received information
        err_mpi = MPI_Waitany(nb_proc,&requestfrom[0],&index,&status);
        if ( err_mpi != MPI_SUCCESS)
        {
            std::ostringstream oss;
            oss << "Problem with MPI_Waitany with from request !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
            throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
        }

        if (index != MPI_UNDEFINED)
            do
            {

                size_t r = receive[index].size();
                assert(index == status.MPI_SOURCE);
                assert(tag_info == status.MPI_TAG);

                int k = 0;
                for (auto iti = receive[index].begin(),itie = receive[index].end(); iti != itie; ++iti,++k)
                {
                    info_manager->setInfo(*iti,info_from[index][k],index);
                }

                err_mpi = MPI_Waitany(nb_proc,&requestfrom[0],&index,&status);
                if ( err_mpi != MPI_SUCCESS)
                {
                    std::ostringstream oss;
                    oss <<"Problem with MPI_Waitany with from request !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }

            } while (index != MPI_UNDEFINED);
        // check every thing been sent
        err_mpi = MPI_Waitall(nb_proc,&requestto[0],MPI_STATUS_IGNORE);
        if ( err_mpi != MPI_SUCCESS )
        {
            std::ostringstream oss;
            oss <<"Problem with MPI_Waitall with to request !\n error = "<<err_mpi<<" proc "<<proc_id<<std::endl;
            throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
        }
    }
    // for unknown receive pattern
    else if (exchange_type < 4)
    {
        std::vector < typename I::information_t >  info_from;
        // loop to pack and send information where it have to
        for (size_t i = 0; i < nb_proc; ++i)
        {
            size_t s = send[i].size();
            if (s)
            {
                info_to.reserve(s);
                for (auto iti = send[i].begin(),itie = send[i].end(); iti != itie; ++iti)
                {
                    info_to[i].push_back(info_manager->getInfo(iti->second,i));
                }
                // use of non blocking synchronize send to insure that when barrier is entered there is no packet send and not yet
                // receive somewhere (avoid finishing loop somewhere an loosing a entering packet).
                // This is done by sender because only sender know that it is sending
                err_mpi = MPI_Issend((void *) &info_to[i][0],s*size_info,MPI_BYTE,i,tag_info,world,&requestto[i]);
                if ( err_mpi != MPI_SUCCESS)
                {
                    std::ostringstream oss;
                    oss << "Problem with MPI_Issend when posting send to !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }
            }
        }

        // overlaping computation and comunication
        if (work) work();

        // Receiving loop
        MPI_Request barrier_request(MPI_REQUEST_NULL);
        do
        {
            // probing pending message without blocking to be able to check barrier achievement (end of loop)
            err_mpi = MPI_Iprobe(MPI_ANY_SOURCE,tag_info,world,&flag,&status);
            if ( err_mpi != MPI_SUCCESS)
            {
                std::ostringstream oss;
                oss << "Problem with MPI_Iprob when testing incomming message !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
            }

            // there is a pending message : get it and unpack
            if (flag)
            {
                index = status.MPI_SOURCE;
                err_mpi = MPI_Get_count( &status, MPI_BYTE, &count );
                if ( ( err_mpi != MPI_SUCCESS ) && ( count != MPI_UNDEFINED ) )
                {
                    std::ostringstream oss;
                    oss << "Problem with MPI_Get_count when computing incomming message size !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }
                if (count%size_info)
                {
                    std::ostringstream oss;
                    oss << "Problem with incomming message size, it is not a multiple of the information size!\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }
                size_t r = count/size_info;
                info_from.resize(r);
                err_mpi = MPI_Recv((void *) &info_from[0],count,MPI_BYTE,index,tag_info,world,&status);
                if ( err_mpi != MPI_SUCCESS)
                {
                    std::ostringstream oss;
                    oss << "Problem with MPI_Recv when receiving incoming message from !\n error="<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }
                for (int k = 0; k < r; ++k)
                {
                    info_manager->setInfo(info_from[k],index);
                }
            }

            // Does all sended mesages are arrived ? we test ... we don't want to wait for that as we want to loop
            err_mpi = MPI_Testany(nb_proc,&requestto[0],&index,&flag,MPI_STATUS_IGNORE);
            if ( err_mpi != MPI_SUCCESS )
            {
                std::ostringstream oss;
                oss << "Problem with MPI_Testany with to request !\n error = "<<err_mpi<<" proc "<<proc_id<<std::endl;
                throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
            }
            // what we want to know is if all ssend request are finished or not
            // if yes this proc enter the barrier
            if (flag && index == MPI_UNDEFINED)
            {
                err_mpi = MPI_Ibarrier(world,&barrier_request);
                if ( err_mpi != MPI_SUCCESS )
                {
                    std::ostringstream oss;
                    oss << "Problem with MPI_Testany with to request !\n error = "<<err_mpi<<" proc "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                }
            }

            err_mpi = MPI_Test(&barrier_request,&flag,MPI_STATUS_IGNORE);
            if ( err_mpi != MPI_SUCCESS )
            {
                std::ostringstream oss;
                oss << "Problem with MPI_Test with barrier !\n error = "<<err_mpi<<" proc "<<proc_id<<std::endl;
                throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
            }

        } while (!flag);
    }
    else
    {
        std::ostringstream oss;
        oss << "Unknown exchange type proc: "<<proc_id<<std::endl;
                    throw xCommExchangerException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }

}

} // end namspace
#endif
