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


    oss<<"Simple copy"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat1(mat,false);
    mat1.printMatrixMarket(oss);

    oss<<"Copy and transpose"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat2(mat,true);
    mat2.printMatrixMarket(oss);

    oss<<"Copy and transpose again"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat3(mat2,true);
    mat3.printMatrixMarket(oss);

    int sel_n[5] = {0,1,1,0,0};
    oss<<"Extracting sub 2x1 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat4(mat,sel_n,nullptr);
    mat4.printMatrixMarket(oss);


    oss<<"Copy and transpose retangular"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat5(mat4,true);
    mat5.printMatrixMarket(oss);

    // closing results
    oss.close();

    return xtool::xMPIEnv::finalize();
}
