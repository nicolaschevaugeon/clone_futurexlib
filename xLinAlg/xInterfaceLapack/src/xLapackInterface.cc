/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

  
#include "xLapackInterface.h"
#include <algorithm> 
#include <iostream>

namespace xlinalg {

 void solve(const xDenseMatrix &A, const xCSRVector &b, xCSRVector &x){
   int M = A.nline();
   int N= A.ncolumn();
   int LDA = M;
   double *a = new double [M*N];
   const double *ain = A.ReturnValPointer ();
   std::copy(&ain[0],&ain[M*N], &a[0]  );
   int *IPIV = new int[M>N? N:M];
   int INFO; 
   dgetrf_(&M, &N, a, &LDA, IPIV, &INFO ); 
   //std::cout << INFO << std::endl;
   if (INFO !=0){
     if (INFO > 0){
       std::cout << "Can't factorise : singular matrice nul pivot at line  "<< INFO << std::endl; 
     }
     else{
       std::cout << "Problem in the interface with lapack "<< __FILE__ << " " << __LINE__ << std::endl; 
     }
     throw;  
   }
   /*
   // the following compute the reciproc Condition number ... Not needed in most case.
   double matnorm =1.;
   double rcond;
   std::vector <double > WORK(4*M);
   std::vector< int>     IWORK(M);
   dgecon_("1", &M, a, &LDA, &matnorm , &rcond, &WORK[0], &IWORK[0], &INFO );
   std::cout << "Reciproc Condition Number (rcond) : " << rcond<<" " << INFO<<  std::endl; 
   */
   int nrhs =1;
   std::copy(&b[0],&b[M], &x[0]  );
   dgetrs_("N", &M, &nrhs, a, &LDA,  IPIV, &x[0], &M, &INFO);  
   delete[] a;
   delete[] IPIV;
 }
  
  void solve (const xDenseLUFactor &LU, const xCSRVector & b, xCSRVector &x){
    std::copy(b.begin(),b.end(), x.begin()  );
    int INFO;
    int nrhs =1;

    dgetrs_("N", const_cast<int *>(&LU.nl), &nrhs, const_cast<double *>(LU.LU), const_cast<int *>(&LU.nl),  const_cast<int *>(LU.IPIV), &x[0], const_cast<int *>(&LU.nl), &INFO);  
  }

  void solve (const xDenseLUFactor &LU, double *x, int nrhs, int ldx){
    int INFO;
    dgetrs_("N", const_cast<int *>(&LU.nl), &nrhs, const_cast<double *>(LU.LU), const_cast<int *>(&LU.nl),  const_cast<int *>(LU.IPIV), &x[0], const_cast<int *>(&LU.nl), &INFO);  
  }
  
  void solve(xDenseMatrix &A, xCSRVector &b){
   int M = A.nline();
   int N=  A.ncolumn();
   int LDA = M;
   double *a = A.ReturnValPointer ();
   int *IPIV = new int[M>N? N:M];
   int INFO;
   dgetrf_(&M, &N, a, &LDA, IPIV, &INFO );
   int nrhs =1;
   dgetrs_("N", &M, &nrhs, a, &LDA,  IPIV, &b[0], &M, &INFO);  
   delete[] IPIV;
  }

 void solve(xDenseMatrix &A, xDenseMatrix &b){
   int M = A.nline();
   int N=  A.ncolumn();
   int LDA = M;
   double *a = A.ReturnValPointer ();
   int *IPIV = new int[M>N? N:M];
   int INFO;
   dgetrf_(&M, &N, a, &LDA, IPIV, &INFO );
   int nrhs = b.ncolumn();
   dgetrs_("N", &M, &nrhs, a, &LDA,  IPIV, b.ReturnValPointer (), &M, &INFO);  
   delete[] IPIV;
  }

  xDenseLUFactor::xDenseLUFactor(const xDenseMatrix & A): nl(A.nline()), nc(A.ncolumn()), LU(nullptr), IPIV(nullptr){
    int M = A.nline();
    int N = A.ncolumn();
    int LDA = nl;
    LU = new double [M*N];
    const double *ain = A.ReturnValPointer ();
    std::copy(&ain[0],&ain[M*N], &LU[0]  );
    IPIV = new int[M>N? N:M];
    dgetrf_(&M, &N, LU, &LDA, IPIV, &INFO ); 
    // std::cout << INFO << std::endl;
    if (INFO !=0){
      if (INFO > 0){
	std::cout << "Can't factorise : singular matrice nul pivot at line  "<< INFO << std::endl; 
      }
      else{
	std::cout << "Problem in the interface with lapack "<< __FILE__ << " " << __LINE__ << std::endl; 
      }
      throw;  
    }
  }
  
  xDenseLUFactor::~xDenseLUFactor(){
    delete [] LU;
    delete [] IPIV;
  } 

  void invert(xDenseMatrix &A){
    int M = A.nline();
    int N=  A.ncolumn();
    int LDA = M;
    double *a = A.ReturnValPointer ();
    int *IPIV = new int[M>N? N:M];
    int INFO;
    int NB = 4;
    int LWORK = std::max(1,N*NB);
    double *WORK = new double[LWORK];
    dgetrf_(&M, &N, a, &LDA, IPIV, &INFO );
    dgetri_(&N, a, &LDA,  IPIV, WORK, &LWORK , &INFO);  
    delete[] IPIV;
    delete[] WORK;
  }

} // end of namespace
  
