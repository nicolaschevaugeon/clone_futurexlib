/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef XCOMMEXCHANGER_H
#define XCOMMEXCHANGER_H

//std
#include <cassert>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <exception>
#include <sstream>

// eXlibris

// parrallel
#include <mpi.h>

// distmesh
#include "xComm.h"
#include "xStoreCommEntity.h"


namespace distmesh
{

template < typename I >
class xCommExchanger
{
    public:
        typedef std::vector < std::set < typename I::information_key_t > > info_key_container_t;
        typedef std::vector < std::map < typename I::information_key_t, typename I::information_key_t > > info_key_pair_container_t;
        typedef typename std::set < typename I::information_key_t >::const_iterator info_key_iter_t;
        typedef typename std::map < typename I::information_key_t, typename I::information_key_t >::const_iterator info_key_pair_iter_t;

        xCommExchanger(void);

        void initIndex( I & info_manager_, const MPI_Comm &world_);
        template < typename ITER >
        void acumulateIndexOwnerOnly( ITER itb, ITER ite );
        template < typename ITER >
        void acumulateIndexAllToAll( ITER itb, ITER ite );
        template < typename ITER >
        void acumulateIndexOneSide( ITER itb, ITER ite );
        void exchangeInformation( void (*work)(void) = NULL);

    private:
        MPI_Comm world;
        int proc_id,nb_proc,exchange_type;
        info_key_pair_container_t send;
        info_key_container_t receive;
        I *info_manager;
};

class xCommExchangerException : public std::exception
{
    public:
        xCommExchangerException(std::string,std::string,int,std::string,std::string);
        ~xCommExchangerException() throw( );
        const char * what() const throw( );

    private:
        std::string msg;
};


} // end namspace

#include "xCommExchanger_imp.h"

#endif
