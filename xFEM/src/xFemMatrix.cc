#include "xFemMatrix.h"
#include "xBlasDef.h"


namespace xfem {

void axpy(const double &alpha, const xFemMatrix<double> &A, xFemMatrix<double> &B)
{
  int N =  A.getNbRow()* A.getNbCol();
  int M =  B.getNbRow()* B.getNbCol();
  if (M!=N) {
    std::cout << "Dimensions not compatible in " << __FILE__<< ":"<<  __LINE__ << std::endl;
    throw ;
  }
  xlinalg::xCPPBlasDef<double>::axpy(N,alpha,A.ReturnValPointer(),B.ReturnValPointer());
}

}
