/*
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <pthread.h>
#include <cmath>
#include <cstdio>
#include "xMPIEnv.h"


using namespace xtool;
using namespace std;
class data_
{
    public:
        MPI_Comm world;
        int proc_id,nb_proc,thread_id,nb;
        pthread_t id;
};
void * ring_com(void *arg)
{
    data_ &d = *( static_cast < data_ * >( arg ));
    MPI_Comm_rank(d.world, &d.proc_id);
    MPI_Comm_size(d.world,&d.nb_proc);
    int tag = 489,err_mpi;
    int nb = ( d.proc_id+1 )*d.nb_proc*( d.thread_id+1 );
    int offset;
    cout<<"Start ring_com proc "<<d.proc_id<<" thread "<<d.thread_id<<" with nb "<<nb<<endl;
    if (d.proc_id > 0)
    {
        err_mpi = MPI_Recv((void *) &offset,1,MPI_INT,d.proc_id-1,tag,d.world,MPI_STATUS_IGNORE);
        if ( err_mpi != MPI_SUCCESS)
        {
            cout<<"Prb with MPI_RECV"<<endl;
            throw d.proc_id;
        }
        nb += offset;
    }
    /* to slow down an see somethin on top for example
       double y=1.2;
       for (int i=0;i<1000000;++i)
       {
        y=pow(i*1.,2.);
        y=pow(i*1.,3.);
        y=pow(i*1.,4.);
       }
       cout<<"y = "<<y<<endl;
     */
    if (d.proc_id < d.nb_proc-1)
    {
        err_mpi = MPI_Send((void *) &nb,1,MPI_INT,d.proc_id+1,tag,d.world);
        if ( err_mpi != MPI_SUCCESS)
        {
            cout<<"Prb with MPI_SEND"<<endl;
            throw d.proc_id;
        }
    }

    d.nb = nb;

    return nullptr;
}


int main(int argc, char *argv[])
{
    // initialize mpi univers
    char ini_type = xMPIEnv::THREAD_MULTIPLE;
    //char ini_type = xMPIEnv::THREAD_FUNNELED;
    xMPIEnv::init(argc,argv,ini_type);


    int nb_thread = 4;
    int proc_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    // all outputs are now recorded in individual files
    std::string pid = std::to_string(proc_id);
    std::string fo = "proc_"+pid+"_output.txt";
    auto ok = freopen (fo.c_str(),"w",stdout);
    if(!ok){std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << std::endl; throw;}
  
    


    // thread
    pthread_attr_t attr;
    if (pthread_attr_init(&attr)) throw -1;

    vector < data_ > vdata(nb_thread);


    // Prepare  data
    for (int i = 0; i < nb_thread; ++i)
    {
        MPI_Comm_dup(MPI_COMM_WORLD, &vdata[i].world);
        vdata[i].thread_id = i;
    }
    // starting threads
    for (int i = 0; i < nb_thread; ++i)
    {
        cout<<"Start thread "<<i<<endl;
        pthread_create(&vdata[i].id,&attr,ring_com,(void *) &vdata[i]);
    }

    // main process open file during thread works
    ofstream ofs("proc_"+pid+"_ref.txt",ios::out | ios::trunc);


    // waiting threads to join
    for (int i = 0; i < nb_thread; ++i)
    {
        void *res;
        pthread_join(vdata[i].id, &res);
        ofs<<"proc "<<proc_id<<" ( proc dump "<<vdata[i].proc_id<<" )thread "<<vdata[i].thread_id<<" got "<<vdata[i].nb<<endl;
    }

    ofs.close();

    // waiting all proc reach this point (all pthread finished)
    MPI_Barrier(MPI_COMM_WORLD);

    // cleanning 
    for (int i = 0; i < nb_thread; ++i)
    {
        MPI_Comm_free(&vdata[i].world);
    }

    // don't forget this one
    return xMPIEnv::finalize();

}
/*
 * Below is given a version that can be used for testing a more shuffled pattern of communication.
 * Says with this version all thread send messages that can be received by any thread of target MPI processor.
 * By commenting out a special tag value you can have similar  behavior compare to official version above.
 *
#include <iostream>
#include <sstream>
#include <vector>
#include <pthread.h>
#include <cmath>
#include "xMPIEnv.h"


using namespace xfem;
using namespace std;
class data
{
    public:
        int proc_id,nb_proc,thread_id,nb;
        pthread_t id;
};
void * ring_com(void *arg)
{
    data &d = *( static_cast < data * >( arg ));
    int err_mpi;
    // full unordered :  threat 0 of proc k send info that can be receive by any thread of proc k+1
     int tag = 489 ;
    // to order communication : threat 0 of proc k send info and only thread 0 of proc k+1 can receive it
    // fall back to something rather similar to official version with dumped communicator
    // int tag = 489+d.thread_id ;
    int nb = ( d.proc_id+1 )*d.nb_proc*( d.thread_id+1 );
    int offset;
    cout<<"Start ring_com proc "<<d.proc_id<<" thread "<<d.thread_id<<" with nb "<<nb<<endl;
    if (d.proc_id > 0)
    {
        err_mpi = MPI_Recv((void *) &offset,1,MPI_INT,d.proc_id-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if ( err_mpi != MPI_SUCCESS)
        {
            cout<<"Prb with MPI_RECV"<<endl;
            throw d.proc_id;
        }
        nb += offset;
    }
    // to slow down an see something in top for example
    //   double y=1.2;
    //   for (int i=0;i<1000000;++i)
    //   {
    //    y=pow(i*1.,2.);
    //    y=pow(i*1.,3.);
    //    y=pow(i*1.,4.);
    //   }
    //   cout<<"y = "<<y<<endl;
    //
    if (d.proc_id < d.nb_proc-1)
    {
        err_mpi = MPI_Send((void *) &nb,1,MPI_INT,d.proc_id+1,tag,MPI_COMM_WORLD);
        if ( err_mpi != MPI_SUCCESS)
        {
            cout<<"Prb with MPI_SEND"<<endl;
            throw d.proc_id;
        }
    }

    d.nb = nb;

    return NULL;
}


int main(int argc, char *argv[])
{
    // initialize mpi univers
    //char ini_type = xMPIEnv::THREAD_MULTIPLE;
    char ini_type = xMPIEnv::THREAD_FUNNELED;
    xMPIEnv::init(argc,argv,ini_type);


    int nb_thread = 8;
    int proc_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    // all outputs are now recorded in individual files
    std::ostringstream oss;
    oss<<"proc_"<<proc_id<<"_"<<"output.txt";
    freopen (oss.str().c_str(),"w",stdout);


    int nb_proc;
    MPI_Comm_size(MPI_COMM_WORLD,&nb_proc);

    // thread
    pthread_attr_t attr;
    if (pthread_attr_init(&attr)) throw -1;

    vector < data > vdata(nb_thread);


    for (int i = 0; i < nb_thread; ++i)
    {
        vdata[i].proc_id = proc_id;
        vdata[i].nb_proc = nb_proc;
        vdata[i].thread_id = i;
        cout<<"Start thread "<<i<<endl;
        pthread_create(&vdata[i].id,&attr,ring_com,(void *) &vdata[i]);
    }


    for (int i = 0; i < nb_thread; ++i)
    {
        void *res;
        pthread_join(vdata[i].id, &res);
        cout<<"proc "<<proc_id<<"thread "<<vdata[i].thread_id<<" got "<<vdata[i].nb<<endl;
    }

    // don't forget this one
    xMPIEnv::finalize();

    return 0;
}
*/
