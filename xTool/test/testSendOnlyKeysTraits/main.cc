//#define FORCE_MPI2
#include "xMPIEnv.h"
#include "xDataExchanger.h"
#include "mpi.h"
#include <iostream>
#include <array>
#include <set>
#include <vector>

// message size
#define MSGSIZE 50

using namespace std;



typedef array < int,MSGSIZE > info_t;

class dummyKeyAndInfoManager
{
    public:
        // mandatory types for information class familly
        typedef xtool::homogeneous_data_style_trait data_style_trait;
        typedef xtool::send_only_keys_communication_trait communication_trait;
        typedef const info_t * information_key_t;
        typedef  info_t information_t;

        /// Constructor
        dummyKeyAndInfoManager(int proc_id_) : proc_id(proc_id_)
        {
            if (!proc_id)
                one_instance.fill(9);

        }

        // mandatory method for key manager class familly
        information_key_t localObjectKey( const info_t & o)
        {
            return &o;
        }
        std::set < int > getMessageRanks( const information_key_t & lk)
        {
            set < int > tmp;
            tmp.insert(1);
            return tmp;
        }
        // mandatory method for information manager class familly
        information_t getInfo(information_key_t key, int sendto)
        {
            if (!proc_id)
                return one_instance;
            else
                throw -1;

        }
        void setInfo(const std::vector < information_t > &info, int receivedfrom )
        {
            // we just do nothing to be able to measure only communication and implementation overhead
    /*        if (receivedfrom == 0 && proc_id == 1)
                cout<<"Ok me proc 1 i got something from 0 : beg="<<info[0][0]<<" end="<<info[0][MSGSIZE-1]<<endl;
                */
        }
        information_t one_instance;
    private:
        int proc_id;
};

int main(int argc, char *argv[])
{
    // initialize mpi universe
    xtool::xMPIEnv::init(argc,argv);

    int proc_id, nb_proc;
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size( MPI_COMM_WORLD, &nb_proc);

    // manager creation : it contains the information that will be send
    dummyKeyAndInfoManager ikm(proc_id);

    // exchanger : we are going to exchange in a send only manner so we use approriate key container
    xtool::xKeyContainerSendOrRecv < dummyKeyAndInfoManager::information_key_t > key_container(MPI_COMM_WORLD );

    // set keys for information to send
    if (proc_id == 0)
    {
        key_container.accumulateKeys(ikm.one_instance,ikm);
    }

    for (int i = 0; i < 5000; ++i)
    {
        // synchronisation to profile only exchange time, not waiting period
        MPI_Barrier(MPI_COMM_WORLD);

        // exchange
        xtool::exchangeInformation(key_container,ikm);
    }


    return  xtool::xMPIEnv::finalize();

}

