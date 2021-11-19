/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include <limits>
#include <algorithm>
#include <functional>

#include "xDistBlockDenseMatrix.h"

#include "xBlasDef.h"
namespace xlinalg
{

/*
   xDistBlockDenseMatrix::xDistBlockDenseMatrix(const xDistBlockDenseMatrix &other) :
   dat(other.dat)
   ,dist_index(other.dist_index)
   ,reduce_keys(dist_index.getComm() )
   ,set_keys(dist_index.getComm() )
   ,status(other.status)
   ,proc_id(other.proc_id)
   {
   //create key manager
   xKeyManagerDistBlockDenseMatrix key_manager(dist_index,dat);
   //create keys : here keys have changed from other has dat is not the same
   reduce_keys.accumulateKeysOwnerGather(dist_index.begin(),dist_index.end(),key_manager);
   set_keys.accumulateKeysOwnerScatter(dist_index.begin(),dist_index.end(),key_manager);

   }
 */

//===== xPolicyxDistBlockDenseMatrixGemv ================================================================================
void xPolicyxDistBlockDenseMatrixGemv < double,xTraitMatrixLowerSym >::mv( const int & n, const double & alpha, const double*A, const int & LDA,const double *X, const double & beta, double *Y)
{
    xlinalg::xCPPBlasDef<double>::symv(LOWERTRIANGULAR, n, alpha, A, LDA, X, beta, Y);
}
void xPolicyxDistBlockDenseMatrixGemv < double,xTraitMatrixUpperSym >::mv( const int & n, const double & alpha, const double*A, const int & LDA,const double *X, const double & beta, double *Y)
{
    xlinalg::xCPPBlasDef<double>::symv(UPPERTRIANGULAR, n, alpha, A, LDA, X, beta, Y);

}
void xPolicyxDistBlockDenseMatrixGemv < double,xTraitMatrixUnSym >::mv( const int &n, const double & alpha, const double*A, const int & LDA,const double *X, const double & beta, double *Y)
{
    xlinalg::xCPPBlasDef<double>::gemv(NOTRANSPOSE, n, n, alpha, A, LDA, X,  beta,  Y);
}
}             // end of namespace

