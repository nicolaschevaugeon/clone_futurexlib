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
    int n = 5;
    int m = 9;
    xGraphMatrix gA(n,m,0);
    gA.add(0,0);
    gA.add(0,3);
    gA.add(0,4);
    gA.add(0,7);
    gA.add(1,1);
    gA.add(1,3);
    gA.add(1,4);
    gA.add(1,8);
    gA.add(2,1);
    gA.add(2,5);
    gA.add(3,2);
    gA.add(3,6);
    gA.add(4,2);
    gA.add(4,6);
    gA.add(4,8);
    gA.countNNZ();
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > A(gA);
    gA.clear();
    A(1,1) = 1.;
    A(1,4) = 2.;
    A(1,5) = 3.;
    A(1,8) = 5.;
    A(2,2) = 4.;
    A(2,4) = 3.;
    A(2,5) = 5.;
    A(2,9) = 7.;
    A(3,2) = 6.;
    A(3,6) = 8.;
    A(4,3) = 7.;
    A(4,7) = 9.;
    A(5,3) = 7.;
    A(5,7) = 6.;
    A(5,9) = 3.;
    oss<<"A matrix"<<endl;
    A.printMatrixMarket(oss);
/*

    matrix B

        |1|2|3|4|5|
        |0|1|2|3|4|
        ===========
   |1|0||1| | |2|3|
   |2|1|| |4| |3|5|
   |3|2|| |6| | | |
   |4|3|| | |7| | |
        ===========

 */
    n = 4;
    m = 5;
    xGraphMatrix gB(n,m,0);
    gB.add(0,0);
    gB.add(0,3);
    gB.add(0,4);
    gB.add(1,1);
    gB.add(1,3);
    gB.add(1,4);
    gB.add(2,1);
    gB.add(3,2);
    gB.countNNZ();
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > B(gB);
    gB.clear();
    B(0,0) = 1.;
    B(0,3) = 2;
    B(0,4) = 3;
    B(1,1) = 4;
    B(1,3) = 3;
    B(1,4) = 5;
    B(2,1) = 6;
    B(3,2) = 7;
    oss<<"B matrix"<<endl;
    B.printMatrixMarket(oss);


/*

    matrix C : empty of size A

 */
    n = 5;
    m = 9;
    xGraphMatrix gC(n,m,0);
    gC.countNNZ();
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixCindex > C(gC);
    gC.clear();
    oss<<"initial C matrix"<<endl;
    C.printMatrixMarket(oss);

    // test offset_r=0 offset_c=0 addition
    A.AoperB(0,0,B,C,[](double a, double b){return a+b;});
    oss<<"C matrix test offset_r=0 offset_c=0 addition"<<endl;
    C.printMatrixMarket(oss);

    // test offset_r=1 offset_c=0 substraction
    A.AoperB(1,0,B,C,[](double a, double b){return a-b;});
    oss<<"C matrix test offset_r=1 offset_c=0 substraction"<<endl;
    C.printMatrixMarket(oss);

    // test offset_r=1 offset_c=2 replacement
    A.AoperB(1,2,B,C,[](double a, double b){return b;});
    oss<<"C matrix test offset_r=1 offset_c=2 replacement"<<endl;
    C.printMatrixMarket(oss);

    // test offset_r=1 offset_c=4 replacement if not null
    A.AoperB(1,4,B,C,[](double a, double b){return (b)?b:a;});
    oss<<"C matrix test offset_r=1 offset_c=4 replacement if not null"<<endl;
    C.printMatrixMarket(oss);

    // test same matrix multiply
    A.AoperB(0,0,A,C,[](double a, double b){return a*b;});
    oss<<"C matrix test same matrix multiply"<<endl;
    C.printMatrixMarket(oss);
    // closing results
    oss.close();

    return xtool::xMPIEnv::finalize();
}
