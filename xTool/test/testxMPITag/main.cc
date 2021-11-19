/*
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include "xMPIEnv.h"
#include "xMPITag.h"


using namespace xtool;
using namespace std;


int main(int argc, char *argv[])
{
    // initialize mpi universe
    xMPIEnv::init(argc,argv);

    // local variable
    MPI_Comm world = MPI_COMM_WORLD;
    int proc_id, proc_idc;

    // get rank for world
    MPI_Comm_rank(world, &proc_id);

    // first call for world comm
    int tag = xMPITag::getNewTag(world);

    // open reference
    ofstream ofs("proc_"+std::to_string(proc_id)+"_ref.txt",ios::out | ios::trunc);

    // output
    ofs<<"Tag for world proc "<<proc_id<<" is "<<tag<<endl;

    // duplcation of world
    MPI_Comm dupworld;
    MPI_Comm_dup(MPI_COMM_WORLD, &dupworld);

    // get rank for dupworld
    MPI_Comm_rank(dupworld, &proc_id);

    // first call for dupworld comm
    tag = xMPITag::getNewTag(dupworld);
    ofs<<"Tag for dubworld proc "<<proc_id<<" is "<<tag<<endl;

    // new com based on color
    MPI_Comm colworld;
    MPI_Comm_split ( dupworld,proc_id%2,0, &colworld );

    // get rank for colworld
    MPI_Comm_rank(colworld, &proc_idc);

    // first call for colworld comm
    tag = xMPITag::getNewTag(colworld);
    ofs<<"Tag for colworld proc "<<proc_idc<<" is "<<tag<<endl;

    // second call
    ofs<<"new serie"<<endl;
    tag = xMPITag::getNewTag(dupworld);
    ofs<<"Tag for dubworld proc "<<proc_id<<" is "<<tag<<endl;
    tag = xMPITag::getNewTag(colworld);
    ofs<<"Tag for colworld proc "<<proc_idc<<" is "<<tag<<endl;

    // first call for SELF
    tag = xMPITag::getNewTag(MPI_COMM_SELF);
    ofs<<"Tag for selfworld proc "<<proc_idc<<" is "<<tag<<endl;

    // test that xMPITag work also with a null comm
    tag = xMPITag::getNewTag(MPI_COMM_NULL);
    ofs<<"Tag for nullworld proc "<<proc_idc<<" is "<<tag<<endl;


    // up to internal limit of xMPITag
    // at time of creation it is 32767.
    for(int i=0; i< 32767+3; ++i)
        tag = xMPITag::getNewTag(colworld);
    ofs<<"Tag for colworld after loop  proc "<<proc_idc<<" is "<<tag<<endl;


    tag = xMPITag::getNewTag(dupworld);
    ofs<<"Tag for dubworld proc "<<proc_id<<" is "<<tag<<endl;


    ofs.close();

    return xMPIEnv::finalize();
}
