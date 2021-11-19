/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef _XCSRMATRIXIMP_H
#define _XCSRMATRIXIMP_H
//
template < typename STORAGE_TYPE_REQUESTED , typename INDEXING_REQUESTED >
void xCSRMatrix::getMemoryAccess(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, STORAGE_TYPE_REQUESTED &matrix_storage_requested,INDEXING_REQUESTED &indexing_requested)
{ 

   transfertToConnectionParameter<matrix_pattern,matrix_storage,STORAGE_TYPE_REQUESTED,INDEXING_REQUESTED>(ns,ms,nnzs,idx1,idx2,dat,matrix_storage_requested,indexing_requested);

   return;

}

//general template function to transfert from XCSRMatrix pointer suitable for solver to do the connection
// see specialization below to see what is possible (implemented)
template <typename PATTERN, typename STORAGE_TYPE , typename STORAGE_TYPE_REQUESTED ,typename INDEXING_REQUESTED >
void xCSRMatrix::transfertToConnectionParameter(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, STORAGE_TYPE_REQUESTED &matrix_storage_requested,INDEXING_REQUESTED &indexing_requested)
{
    std::cout<<"look's like you are trying to connect to a xCSRMatrix with traits defining pattern,storage and index not covered by the actual implementation !\n";
    throw;

}
// specialisation template of transfertToConnectionParameter :
// xCSRMatrix considered as xTraitMatrixUnSym and xTraitMatrixSparceCSR
// connecting object : xTraitMatrixSparceCOO, xTraitMatrixFindex (Mumps for example)
template <>
void xCSRMatrix::transfertToConnectionParameter<xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixSparceCOO,xTraitMatrixFindex>(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, xTraitMatrixSparceCOO &matrix_storage_requested,xTraitMatrixFindex &indexing_requested);

// specialisation template of transfertToConnectionParameter :
// xCSRMatrix considered as xTraitMatrixUnSym and xTraitMatrixSparceCSR
// connecting object : xTraitMatrixSparceCSR, xTraitMatrixCindex (SuperLu for example)
template <>
void xCSRMatrix::transfertToConnectionParameter<xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixSparceCSR,xTraitMatrixCindex>(int &ns,int &ms,int &nnzs, int **idx1,int **idx2,double **dat, xTraitMatrixSparceCSR &matrix_storage_requested,xTraitMatrixCindex &indexing_requested);
      
#endif
