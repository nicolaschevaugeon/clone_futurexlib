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

    // starting unsymetric matrix creation
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


    oss<<"Simple copy with template change (xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat1(mat,nullptr,nullptr);
    mat1.printMatrixMarket(oss);

    // starting selecting vector
    int sel_n[5] = {0,0,1,1,0};
    int sel_m[5] = {0,1,0,0,0};

    oss<<"Extracting sub 2x1 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat2(mat,sel_n,sel_m);
    mat2.printMatrixMarket(oss);

    oss<<"Filter sub 2x1 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat2b(mat,sel_n,sel_m,false);
    mat2b.printMatrixMarket(oss);

    oss<<"Extracting sub 2x1 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat3(mat,sel_n,sel_m);
    mat3.printMatrixMarket(oss);


    sel_n[2] = 0;
    sel_n[4] = 1;

    oss<<"Extracting sub 2x5 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat4(mat,sel_n,nullptr);
    mat4.printMatrixMarket(oss);

    oss<<"Extracting sub 5x2 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat5(mat,nullptr,sel_n);
    mat5.printMatrixMarket(oss);

    oss<<"Extracting sub 2x2 matrix (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat6(mat,sel_n,sel_n);
    mat6.printMatrixMarket(oss);

    oss<<"Filtering sub 2x2 matrix (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat6b(mat,sel_n,sel_n,false);
    mat6b.printMatrixMarket(oss);
    sel_m[0] = 1;

    oss<<"Extracting sub 2x2 matrix (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat7(mat,sel_m,sel_n);
    mat7.printMatrixMarket(oss);

    oss<<"Extracting sub 2x2 matrix from 5x2 (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat8(mat5,sel_m,nullptr);
    mat8.printMatrixMarket(oss);

    oss<<"Extracting sub 2x2 matrix from 2x5 (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat9(mat4,nullptr,sel_n);
    mat9.printMatrixMarket(oss);

    oss<<"Extracting with merging 2 to 3 an 5 to 4 in sub 2x5 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    sel_n[0] = 0;
    sel_n[1] = -3;
    sel_n[2] = 1;
    sel_n[3] = 1;
    sel_n[4] = -4;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat10(mat,sel_n,nullptr);
    mat10.printMatrixMarket(oss);

    oss<<"Extracting with merging 2 to 3 an 5 to 4 in sub 2x2 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat11(mat,sel_n,sel_n);
    mat11.printMatrixMarket(oss);

    oss<<"Filtering with merging 2 to 3 an 5 to 4 in sub 2x5 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat12(mat,sel_n,nullptr,false);
    mat12.printMatrixMarket(oss);

    oss<<"Filtering with merging 2 to 3 an 5 to 4 in sub 2x2 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat13(mat,sel_n,sel_n,false);
    mat13.printMatrixMarket(oss);

    oss<<"Extracting with merging 2 to 3 an 5 to 4 in sub 2x2 matrix (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat14(mat,sel_n,sel_n);
    mat14.printMatrixMarket(oss);
/*

   |1|2|3|4|5|
   |0|1|2|3|4|
     ===========
   |1|0||1| | |2|3|
   |2|1|| |4| |3|5|
   |3|2|| | |6| | |
   |4|3||2|3| |7|8|
   |5|4||3|5| |8|9|
     ===========

 */

    // graph for a symetric matrix
    xGraphMatrix gs(n,m,1);
    gs.add(0,0);
    gs.add(3,0);
    gs.add(4,0);
    gs.add(1,1);
    gs.add(3,1);
    gs.add(4,1);
    gs.add(2,2);
    gs.add(3,3);
    gs.add(4,3);
    gs.add(4,4);
    gs.countNNZ();

    // starting symetric matrix creation
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mats(gs);
    mats(0,0) = 1;
    mats(3,0) = 2;
    mats(4,0) = 3;
    mats(1,1) = 4;
    mats(3,1) = 3;
    mats(4,1) = 5;
    mats(2,2) = 6;
    mats(3,3) = 7;
    mats(4,3) = 8;
    mats(4,4) = 9;
    oss<<"Starting symetrique square matrix (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    mats.printMatrixMarket(oss);

    sel_m[0] = 0;
    sel_m[2] = 1;
    sel_n[0] = 0;
    sel_n[1] = 1;
    sel_n[2] = 1;
    sel_n[3] = 1;
    sel_n[4] = 0;

    oss<<"Extracting sub 3x2 matrix  (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat100(mats,sel_n,sel_m);
    mat100.printMatrixMarket(oss);

    oss<<"Extracting sub 2x3 matrix  (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat101(mats,sel_m,sel_n);
    mat101.printMatrixMarket(oss);

    oss<<"Extracting sub 2x5 matrix  (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat102(mats,sel_m,nullptr);
    mat102.printMatrixMarket(oss);

    oss<<"Extracting sub 5x2 matrix  (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat103(mats,nullptr,sel_m);
    mat103.printMatrixMarket(oss);

    oss<<"Filtering sub 5x2 matrix  (xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat103b(mats,nullptr,sel_m,false);
    mat103b.printMatrixMarket(oss);

    oss<<"Extracting sub 3x3 matrix  (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat104(mats,sel_n,sel_n);
    mat104.printMatrixMarket(oss);

    oss<<"Filtering sub 3x3 matrix  (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat104b(mats,sel_n,sel_n,false);
    mat104b.printMatrixMarket(oss);

    oss<<"Extracting sub 2x2 matrix in diag block (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat105(mats,sel_m,sel_m);
    mat105.printMatrixMarket(oss);

    oss<<"Extracting sub 2x2 matrix in diag block (xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat106(mats,sel_m,sel_m);
    mat106.printMatrixMarket(oss);

    sel_m[0] = 0;
    sel_m[1] = 0;
    sel_m[2] = 0;
    sel_m[3] = 1;
    sel_m[4] = 1;
    sel_n[0] = 1;
    sel_n[1] = 1;
    sel_n[2] = 0;
    sel_n[3] = 0;
    sel_n[4] = 0;
    oss<<"Extracting sub 2x2 matrix from upper block  (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat107(mats,sel_n,sel_m);
    mat107.printMatrixMarket(oss);
    oss<<"Extracting sub 2x2 matrix from lower block  (xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > mat108(mats,sel_m,sel_n);
    mat108.printMatrixMarket(oss);


    oss<<"Extracting with merging 2 to 3 an 5 to 4 in sub 2x5 matrix (xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    sel_n[0] = 0;
    sel_n[1] = -3;
    sel_n[2] = 1;
    sel_n[3] = 1;
    sel_n[4] = -4;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat109(mats,sel_n,nullptr);
    mat109.printMatrixMarket(oss);

    oss<<"Extracting with merging 2 to 3 an 5 to 4 in sub 2x2 sym matrix (xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat110(mats,sel_n,sel_n);
    mat110.printMatrixMarket(oss);

    oss<<"Extracting with merging 2 to 3 an 5 to 1 in sub 3x3 sym matrix (xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex)"<<endl;
    sel_n[0] = 1;
    sel_n[4] = -1;
    xGenericSparseMatrix < double,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > mat111(mats,sel_n,sel_n);
    mat111.printMatrixMarket(oss);

    // to be continued
    // somme will not work
    // .... from CSC to CSR
    // .... from CSR to CSC
    // .... from COO to CSC
    // etc




    // closing results
    oss.close();
    return xtool::xMPIEnv::finalize();
}
