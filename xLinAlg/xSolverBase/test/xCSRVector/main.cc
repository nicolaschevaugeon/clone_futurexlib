/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include <fstream>

// eXlibris_tools
#include "xMPIEnv.h"

// xlinalg
#include "xCSRVector.h"


using std::ofstream;
using std::endl;
using namespace xlinalg;


// This test BLAS1 call on xCSRVector
int main(int argc, char *argv[]){
  // initialize mpi universe
  xtool::xMPIEnv::init(argc,argv);
  int nb;
  MPI_Comm_size(MPI_COMM_WORLD,&nb);
  if (nb>1) {
      std::cout<<"This test is sequential !"<<std::endl; 
      MPI_Abort(MPI_COMM_WORLD, -1);
  }

  ofstream out("res.txt");
  xCSRVector test(3);
  test(0) = 1.;
  test(1) = 1.;
  test(2) = 1.;
  xCSRVector test2(3);
  test2(0) = 2.;
  test2(1) = 3.;
  test2(2) = 4.;
  xCSRVector test3(test2);
  out << nrm2(test) << " " << "1.73205" <<  endl;
  out << dot(test, test) << " " << "3 " << endl;
  out << dot(test, test2) << " " << "9"<<   endl;  
  axpy(2, test, test3);
  out << test3[0] <<" "  << test3[1]<< " " << test3[2] << " " << "4 5 6  " <<   endl;

  return xtool::xMPIEnv::finalize();

}
