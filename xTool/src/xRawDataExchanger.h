/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/


/*
 *   This file define some simple functions to exchange data between processors using MPI.
 *   It differ from xDataExchanger in the sense that it only exchange data, their is no key - information association mechanism. It is usefull when no concept of key are needed.
 */
#ifndef XRAWDATAEXCHANGER_H
#define XRAWDATAEXCHANGER_H
// std includes
#include <vector>
#include <map>
#include <functional>
// mpi includes
#include "mpi.h"
// xtool includes
#include "xMPITag.h"
namespace xtool{
/*! Collective on all the procs of mpi_comm
 *  Each proc send data to some of the procs of mpi_comm, and receive the same amount of data
 *  from the proc it send to. When data is received it is treated using the recvfunction
 *  IN : data_to_send, std map of data to send. data are indexed with target proc id.
 *       The data is an std::vector of DATA
 *  IN : recvfunction : a function that take as input an int that represent the source proc of       * the message and an std::vector of DATA tha reprensent the received DATA. The user has to define in this function what he wants to do with the data.
!*/
template < typename DATA >
void exchangeRawData(const std::map< int, std::vector <DATA > > &data_to_send, std::function< void (int, const std::vector< DATA> &) > recvfunction , MPI_Comm mpi_comm ){

    const int tag_info = xtool::xMPITag::getNewTag(mpi_comm);
    const size_t size_info = sizeof( DATA );
    std::map< int, std::vector <DATA > > received_data;
    const int nb_send_recv_message = data_to_send.size();
    std::vector < MPI_Request > request_from (nb_send_recv_message, MPI_REQUEST_NULL);
    std::vector < MPI_Request > request_to   (nb_send_recv_message, MPI_REQUEST_NULL);
    size_t send_recv_count = 0;

    for ( const auto & key_data : data_to_send){
        const int send_recv_id = key_data.first;
        const size_t send_recv_size = key_data.second.size();
        received_data.emplace(send_recv_id, std::vector<DATA>(send_recv_size) );
        MPI_Irecv( received_data[send_recv_id].data(), send_recv_size*size_info, MPI_BYTE, send_recv_id, tag_info, mpi_comm, &request_from[send_recv_count]);
        MPI_Isend( const_cast< DATA *>(data_to_send.at(send_recv_id).data()), send_recv_size*size_info, MPI_BYTE, send_recv_id, tag_info, mpi_comm, &request_to[send_recv_count]);
        ++send_recv_count;
    }
    int index;
    MPI_Status status;
    MPI_Waitany(nb_send_recv_message, request_from.data(), &index, &status);

    while (index != MPI_UNDEFINED){
        if (tag_info != status.MPI_TAG) throw;
        const int from = status.MPI_SOURCE;
        recvfunction(from, received_data[from] );
        received_data.erase(from);
        MPI_Waitany(nb_send_recv_message, request_from.data(), &index, &status);
    }
    MPI_Waitall(nb_send_recv_message, request_to.data(),   MPI_STATUS_IGNORE);
    return;
}

/*! Collective on all the procs of mpi_comm
 *  Each proc send data to some of the procs of mpi_comm, and receive any  amount of data from any proc. When data is received it is treated using the recvfunction
 *  IN : data_to_send, std map of data to send. data are indexed with target proc id.
 *       The data is an std::vector of DATA
 *  IN : recvfunction : a function that take as input an int that represent the source proc of       * the message and an std::vector of DATA tha reprensent the received DATA. The user has to define in this function what he wants to do with the data.
!*/
template < typename DATA >
void sendRecvRawData(const std::map< int, std::vector <DATA > > &data_to_send, std::function< void (int, const std::vector< DATA> &) > recvfunction , MPI_Comm mpi_comm ){
    int mpi_rank, mpi_size;
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
    const int tag_info = xtool::xMPITag::getNewTag(mpi_comm);
    const size_t size_info = sizeof( DATA );
    int nb_to_send = data_to_send.size();
    std::vector < MPI_Request > request_to   (nb_to_send, MPI_REQUEST_NULL);
    MPI_Request request_allsenddone = MPI_REQUEST_NULL;
    int nb_send = 0;
    for (const auto & to_data :data_to_send){
        const int  to      = to_data.first;
        const auto & data  = to_data.second;
        MPI_Issend( const_cast< DATA *>(data.data()), size_info*data.size(), MPI_BYTE, to, tag_info, mpi_comm, &request_to[nb_send] );
        ++nb_send;
    }
    MPI_Status status;
    int messageready = 0;
    bool pass_barrier = false;
    int all_send_done = 0;
    while(!all_send_done){
        MPI_Iprobe(MPI_ANY_SOURCE, tag_info, mpi_comm, &messageready, &status );
        if (messageready){
            int from =   status.MPI_SOURCE;
            int size;
            MPI_Get_count(&status, MPI_BYTE, &size);
            std::vector<DATA> recv_buff(size/size_info);
            MPI_Recv(recv_buff.data(), size, MPI_BYTE, from, tag_info, mpi_comm, MPI_STATUSES_IGNORE);
            recvfunction(from, recv_buff);
        }
        if(!pass_barrier){
            int all_send_received = 0;
            MPI_Testall(nb_to_send, request_to.data(), &all_send_received,  MPI_STATUSES_IGNORE);
            if (all_send_received){
                MPI_Ibarrier(mpi_comm, &request_allsenddone);
                pass_barrier =true;
            }
        }
        if (pass_barrier)
            MPI_Test(&request_allsenddone, &all_send_done, MPI_STATUS_IGNORE);
    }
}

} // end namespace xtool
#endif 
