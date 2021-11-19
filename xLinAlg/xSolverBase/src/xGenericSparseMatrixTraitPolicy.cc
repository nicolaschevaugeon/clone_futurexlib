/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xGenericSparseMatrix.h"

namespace xlinalg
{

int xGenericSparseMatrixOpIncrease (int i) { return ++i; }

int * xPolicyGenericSparseMatrix<xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex>::rowidx_stat= nullptr;
int * xPolicyGenericSparseMatrix<xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex>::colidx_stat= nullptr;
int   xPolicyGenericSparseMatrix<xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex>::nnz_stat= 0;

int * xPolicyGenericSparseMatrix<xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex>::rowidx_stat= nullptr;
int * xPolicyGenericSparseMatrix<xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex>::colidx_stat= nullptr;
int   xPolicyGenericSparseMatrix<xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex>::nnz_stat= 0;

int * xPolicyGenericSparseMatrix<xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixFindex>::rowidx_stat= nullptr;
int * xPolicyGenericSparseMatrix<xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixFindex>::colptr_stat= nullptr;

int * xPolicyGenericSparseMatrix<xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex>::rowidx_stat= nullptr;
int * xPolicyGenericSparseMatrix<xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex>::colptr_stat= nullptr;

int * xPolicyGenericSparseMatrix<xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex>::rowidx_stat= nullptr;
int * xPolicyGenericSparseMatrix<xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex>::colptr_stat= nullptr;

int * xPolicyGenericSparseMatrix<xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixCindex>::rowptr_stat= nullptr;
int * xPolicyGenericSparseMatrix<xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixCindex>::colidx_stat= nullptr;

} // end of namespace
