/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include <sstream>
#include <fstream>

// eXlibris_tools
#include "xMPIEnv.h"
// xlinalg
#include "xDistIndex.h"

/*
 * Testing xDistIndex class with following example :
 *  On 3 procs each proc have the following set of Fortran index
 *
 *       p0        p1     p2
 *                        16
 *       1                1
 *       2
 *                 3      3
 *                 5
 *                 13
 *                 15
 *       8         8      8
 *       10
 *
 *
 *  Tested : all  public method and created embedded partition manager
 *
 */

using namespace std;

int main(int argc, char *argv[])
{
    // initialize mpi universe
    xtool::xMPIEnv::init(argc,argv);

    // local variable
    MPI_Comm world = MPI_COMM_WORLD;
    int proc_id,nb_proc;

    // get rank for world
    MPI_Comm_rank(world, &proc_id);
    MPI_Comm_size(world, &nb_proc);


    if (nb_proc != 3)
    {
        cout << argv[0] << " Need 3 mpi process to run "<< endl;
        MPI_Abort(world,-1);
    }

    // open reference
    ofstream ofs("proc_"+std::to_string(proc_id)+"_ref.txt",ios::out | ios::trunc);

    // output
    ofs<<"In proc "<<proc_id<<" Starting"<<endl;

    // constructor test
    xlinalg::xDistIndex idx(world);

    // insertIndex,  insertToFrom, finalize test
    switch (proc_id)
    {
        case 0 :
         {
             idx.insertIndex(1);
             idx.insertIndex(2);
             idx.insertIndex(8);
             idx.insertToFrom(1,2); // to test use of this method in the middle of insertion
             idx.insertIndex(10);
             idx.insertToFrom(8,1);
             idx.insertToFrom(8,2);
             idx.finalize(4,10);
             break;
         }
        case 1 :
         {
             idx.insertIndex(3);
             idx.insertIndex(5);
             idx.insertIndex(13);
             idx.insertIndex(15);
             idx.insertIndex(8);
             idx.insertToFrom(3,2);
             idx.insertToFrom(8,0);
             idx.insertToFrom(8,2);
             idx.finalize(4,15);
             break;
         }
        case 2 :
         {
             idx.insertIndex(16);
             idx.insertIndex(1);
             idx.insertIndex(8);
             idx.insertIndex(3); // here we intentionally don't respect order of p1 to check that it does no make problem
             idx.insertToFrom(1,0);
             idx.insertToFrom(3,1);
             idx.insertToFrom(8,0);
             idx.insertToFrom(8,1);
             idx.finalize(1,16);
             break;
         }
    }

    ofs <<"getGlobalIndexSize : "<<idx.getGlobalIndexSize()<<endl;
    ofs <<"getPackedIndexSize : "<<idx.getPackedIndexSize()<<endl;
    ofs <<"getPackedLocalIndexSize : "<<idx.getPackedLocalIndexSize()<<endl;

    ofs <<"Does 3 index present in this proc : "<<idx.testIndexPresent(3)<<endl;

    int check;
    MPI_Comm_compare(idx.getComm(),world,&check);
    ofs <<"check getComm : "<<( check == MPI_IDENT )<<endl;

    xlinalg::xDistIndex::idx_t const *all_local = idx.getAllPackedLocalIndexSize();
    ofs <<"All local size : ";
    for (int i = 0; i < nb_proc; ++i) ofs <<" "<<i<<","<<all_local[i];
    ofs <<endl;

    xlinalg::xDistIndex::idx_t const *all_dist = idx.getAllPackedLocalIndexDisp();
    ofs <<"All dist size : ";
    for (int i = 0; i < nb_proc; ++i) ofs <<" "<<i<<","<<all_dist[i];
    ofs <<endl;

    ofs <<"All index : ";
    for (auto i : idx) ofs <<" "<<i;
    ofs <<endl;

    ofs <<"All index location in packed index container: ";
    for (auto i : idx) ofs <<" "<<idx.getPackedIndex(i);
    ofs <<endl;

    // choosing 8 has it is common to all proc
    auto b = idx.getPackedAdress(8); // test with user parameter value (not like in loop with iterators)
    ofs <<"All offset from adress of 8 index: ";
    for (auto i : idx) ofs <<" "<<idx.getPackedAdress(i)-b;
    ofs <<endl;

    ofs <<"All owner  : ";
    for (auto i : idx)
    {
        xtool::xConstPartitionObject < xlinalg::xDistIndex::idx_t > po = idx.getConstPartitionObject(i);
        ofs <<" "<<i<<","<<po.isOwner();
    }
    ofs <<endl;

    ofs <<"All id remote  : "<<endl;
    // Choosing 8 has it is common to all proc give a remote address reference which can be used to
    // compute offset. This is to avoid an instable reference file which show address. Using 8 avoid
    // any extra communication
    xtool::xConstPartitionObject < xlinalg::xDistIndex::idx_t > po8 = idx.getConstPartitionObject(8);
    for (auto i : idx)
    {
        xtool::xConstPartitionObject < xlinalg::xDistIndex::idx_t > po = idx.getConstPartitionObject(i);
        ofs <<"for index "<<i<<endl;
        if (po.hasRemoteObject())
            for (auto ro : po.getRemoteObjectsCollectionRange() )
            {
                const xlinalg::xDistIndex::idx_t * p8 = po8.getRemoteObjectOn(ro.getProcId());
                ofs <<" in remote proc "<<ro.getProcId()<<" offset from address of 8 index is "<<ro.getObjectAddress()-p8<<endl;
            }
        else
            ofs <<" no remote equivalent index"<<endl;
    }

    // test copy constructor
    xlinalg::xDistIndex cidx(idx);
    ofs <<"All index location in packed index container for copy: ";
    for (auto i : cidx) ofs <<" "<<cidx.getPackedIndex(i);
    ofs <<endl;

    ofs.close();


    return xtool::xMPIEnv::finalize();
}
