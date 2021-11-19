/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include <iostream>
#include <fstream>
// eXlibris_tools
#include "xMPIEnv.h"
// xlinalg
#include "xGenericSparseMatrix.h"
#include "xGraphMatrix.h"
#include "xTraitsMatrix.h"

using namespace std;
using namespace xlinalg;

// column wise functor
struct cfunctor {

    int nbr;
    int operator () (const int &i, const int &j) {return i+j*nbr; }
};




int main(int argc, char *argv[])
{
    // initialize mpi universe
    xtool::xMPIEnv::init(argc,argv);
    int nb;
    MPI_Comm_size(MPI_COMM_WORLD,&nb);
    if (nb > 1)
    {
        std::cout<<"This test is sequential !"<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // graph for a unsymetric matrix
    int n = 5;
    int m = 5;
    xGraphMatrix g(n,m,0);
    g.add(0,0);
    g.add(0,3);
    g.add(0,4);
    g.add(1,1);
    g.add(1,3);
    g.add(1,4);
    g.add(2,1);
    g.add(3,2);
    g.add(4,2);
    g.add(4,4);
    g.countNNZ();

    // result file
    ofstream oss("out.txt",ios::out | ios::trunc);
/*

   |1|2|3|4|5|
   |0|1|2|3|4|
     ===========
   |1|0||1| | |2|3|
   |2|1|| |4| |3|5|
   |3|2|| |6| | | |
   |4|3|| | |7| | |
   |5|4|| | |8| |9|
     ===========

 */

    // starting unsymetric square matrix creation
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat(g);
    mat(0,0) = 1.;
    mat(0,3) = 2;
    mat(0,4) = 3;
    mat(1,1) = 4;
    mat(1,3) = 3;
    mat(1,4) = 5;
    mat(2,1) = 6;
    mat(3,2) = 7;
    mat(4,2) = 8;
    mat(4,4) = 9;
    oss<<"Starting unsymetrique square matrix (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    mat.printMatrixMarket(oss);

    // maximun storage for dense matrix
    double *dense = new double[25];

    std::fill(dense,dense+25,0.);

    // select patern
    int idx_r[] = {-1,2,-1,0,1};
    int idx_c[] = {-1,-1,-1,1,0};


    oss<<"sub bloc : selr+selc"<<endl;
    cfunctor cf;
    cf.nbr = 3;
    mat.toDense(3,2,dense,idx_r,idx_c,cf);
    oss<<"full storage"<<endl;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            oss << dense[i+5*j] << ' ';
        }
        oss << endl;
    }
    oss<<"packed storage"<<endl;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            oss << dense[cf(i,j)] << ' ';
        }
        oss << endl;
    }

    oss<<"sub bloc : selr+full"<<endl;
    mat.toDense(3,5,dense,idx_r,nullptr,cf);
    oss<<"full storage"<<endl;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            oss << dense[i+5*j] << ' ';
        }
        oss << endl;
    }
    oss<<"packed storage"<<endl;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            oss << dense[cf(i,j)] << ' ';
        }
        oss << endl;
    }

    oss<<"sub bloc : full+selc"<<endl;
    cf.nbr = 5;
    mat.toDense(5,2,dense,nullptr,idx_c,cf);
    oss<<"full storage"<<endl;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            oss << dense[i+5*j] << ' ';
        }
        oss << endl;
    }
    oss<<"packed storage"<<endl;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            oss << dense[cf(i,j)] << ' ';
        }
        oss << endl;
    }

    oss<<"full bloc : full+full"<<endl;
    mat.toDense(5,5,dense,nullptr,nullptr,cf);
    oss<<"full storage"<<endl;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            oss << dense[i+5*j] << ' ';
        }
        oss << endl;
    }

    // closing results
    oss.close();

    return xtool::xMPIEnv::finalize();
}
