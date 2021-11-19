/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

// eXlibris_tools
#include "xMPIEnv.h"
// xlinalg
#include "xDenseMatrix.h"
#include "xCSRVector.h"
#include <iostream>
#include <fstream>

using namespace xlinalg;

int main(int argc, char *argv[]){
  // initialize mpi universe
  xtool::xMPIEnv::init(argc,argv);
  int nb;
  MPI_Comm_size(MPI_COMM_WORLD,&nb);
  if (nb>1) {
      std::cout<<"This test is sequential !"<<std::endl; 
      MPI_Abort(MPI_COMM_WORLD, -1);
  }
  std::ofstream out("res.txt");
  xCSRVector x (2);
  x(0) = 1;
  x(1) = 1;

  xCSRVector y(2);
    
  xDenseMatrix A(2,2);
  A(0,0) = 1.;
  A(0,1) = 1.;
  A(1,0) = 3.;
  A(1,1) = 1.;

  
  //testing level 2 blas using gemv
  gemv(false, 1. ,  A , x, 0., y );
  out << "y = " << y(0) << " " << y(1) << std::endl; // res is y =  2 4
  
  //testing level 2 blas using operator *
  y = A*x;
  out << "y = " << y(0) << " " << y(1) << std::endl; // res is y =  2 4

  xDenseMatrix B(2,2);
  B(0,0) = 1.;
  B(0,1) = 2.;
  B(1,0) = 3.;
  B(1,1) = 4.;
  xDenseMatrix C(2,2);

  //testing level 3 blas using gemm call
  gemm(false, false, 1., A, B, 0., C);
  out << "C = " << std::endl;
  out << C(0,0) << " " << C(0,1) << std::endl;
  out << C(1,0) << " " << C(1,1) << std::endl; //res is C = 4 6; 6 10
  
  // test level 2 blas using operator *
  C = A*B;
  out << "C = " << std::endl;
  out << C(0,0) << " " << C(0,1) << std::endl;
  out << C(1,0) << " " << C(1,1) << std::endl; //res is C = 4 6; 6 10
  
  return xtool::xMPIEnv::finalize();
 
}
