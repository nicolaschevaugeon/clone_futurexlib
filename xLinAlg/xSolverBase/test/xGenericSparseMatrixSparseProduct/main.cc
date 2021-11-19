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

    // result file
    ofstream oss("out.txt",ios::out | ios::trunc);

    // =================================================================================================================================================
    // unsymetric matrix
    // =================================================================================================================================================
    //
/*

    matrix A

        |1|2|3|4|5|
        |0|1|2|3|4|
        ===========
   |1|0||1| | |2|3|
   |2|1|| |4| |3|5|
   |3|2|| |6| | | |
   |4|3|| | |7| | |
        ===========

 */
    int n = 4;
    int m = 5;
    xGraphMatrix gA(n,m,0);
    gA.add(0,0);
    gA.add(0,3);
    gA.add(0,4);
    gA.add(1,1);
    gA.add(1,3);
    gA.add(1,4);
    gA.add(2,1);
    gA.add(3,2);
    gA.countNNZ();
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > A(gA);
    gA.clear();
    A(0,0) = 1.;
    A(0,3) = 2;
    A(0,4) = 3;
    A(1,1) = 4;
    A(1,3) = 3;
    A(1,4) = 5;
    A(2,1) = 6;
    A(3,2) = 7;
    oss<<"A matrix"<<endl;
    A.printMatrixMarket(oss);

/*

    matrix B

        |1|2|3|4|5|6|7|8|9|
        |0|1|2|3|4|5|6|7|8|
        ===================
   |1|0||1| | |2|3| | |5| |
   |2|1|| |4| |3|5| | | |7|
   |3|2|| |6| | | |8| | | |
   |4|3|| | |7| | | |9| | |
   |5|4|| | |7| | | |6| |3|
        ===================

 */
    n = 5;
    m = 9;
    xGraphMatrix gB(n,m,0);
    gB.add(0,0);
    gB.add(0,3);
    gB.add(0,4);
    gB.add(0,7);
    gB.add(1,1);
    gB.add(1,3);
    gB.add(1,4);
    gB.add(1,8);
    gB.add(2,1);
    gB.add(2,5);
    gB.add(3,2);
    gB.add(3,6);
    gB.add(4,2);
    gB.add(4,6);
    gB.add(4,8);
    gB.countNNZ();
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > B(gB);
    gB.clear();
    B(1,1) = 1.;
    B(1,4) = 2.;
    B(1,5) = 3.;
    B(1,8) = 5.;
    B(2,2) = 4.;
    B(2,4) = 3.;
    B(2,5) = 5.;
    B(2,9) = 7.;
    B(3,2) = 6.;
    B(3,6) = 8.;
    B(4,3) = 7.;
    B(4,7) = 9.;
    B(5,3) = 7.;
    B(5,7) = 6.;
    B(5,9) = 3.;
    oss<<"B matrix"<<endl;
    B.printMatrixMarket(oss);

/*

    matrix C

        |1|2|3|4|5|6|7|8|9|
        |0|1|2|3|4|5|6|7|8|
        ===================
   |1|0|| | | | | | | | | |
   |2|1|| | | | |2| | | | |
   |3|2|| | | | | | | | | |
   |4|3|| | | | | | | | | |
        ===================

 */
    n = 4;
    m = 9;
    xGraphMatrix gC(n,m,0);
    gC.add(1,4);
    gC.countNNZ();
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > C(gC);
    gC.clear();
    C(1,4) = 2.;
    oss<<"initial C matrix"<<endl;
    C.printMatrixMarket(oss);

    double alpha=1.;
    double beta=0.;

    // test alpha=1 beta=0
    A.gemm(alpha,beta,B,C);
    oss<<"C matrix test 1 0"<<endl;
    C.printMatrixMarket(oss);

    // test alpha=2 beta=0
    alpha=2.;
    A.gemm(alpha,beta,B,C);
    oss<<"C matrix test 2 0"<<endl;
    C.printMatrixMarket(oss);

    // test alpha=1 beta=-1
    alpha=1.;
    beta=-1.;
    A.gemm(alpha,beta,B,C);
    oss<<"C matrix test 1 -1"<<endl;
    C.printMatrixMarket(oss);

    // test alpha=2 beta=1
    alpha=2.;
    beta=1.;
    A.gemm(alpha,beta,B,C);
    oss<<"C matrix test 2 1"<<endl;
    C.printMatrixMarket(oss);

    // test alpha=0 beta=1
    alpha=0.;
    beta=1.;
    A.gemm(alpha,beta,B,C);
    oss<<"C matrix test 0 1"<<endl;
    C.printMatrixMarket(oss);

    oss<<"Copy and transpose A"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > At(A,true);
    At.printMatrixMarket(oss);


    oss<<"Copy and transpose B"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > B2(B,nullptr,nullptr);
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > Bt(B2,true);
    Bt.printMatrixMarket(oss);

    oss<<"Copy and transpose C"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > Ct(C,true);
    Ct.printMatrixMarket(oss);

    // test alpha=1 beta=1
    alpha=1.;
    beta=-1.;
    Bt.gemm(alpha,beta,At,Ct);
    oss<<"C matrix test transpose"<<endl;
    Ct.printMatrixMarket(oss);


    // closing results
    oss.close();

    return xtool::xMPIEnv::finalize();
}
