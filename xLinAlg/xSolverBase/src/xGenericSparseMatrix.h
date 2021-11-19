/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */
#ifndef _XGENERICSPARSEMATRIX_H
#define _XGENERICSPARSEMATRIX_H

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "xCSRVector.h"
#include "xDataType.h"
#include "xGenericSparseMatrixTraitPolicy.h"
#include "xGraphMatrix.h"
#include "xTraitsMatrix.h"

namespace xlinalg
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xGenericSparseMatrix class ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// This container is a kind of general container for sparse matrix. It support many format via traits template parameter.
//! It represent a matrix A which may be square or not.
//
//! <br/> The connection to the assembler is done via a method given by this class : AddMatrix
//!       In Xfem this is used by xAssembler/xAssembler_dispatcher class which call at low level AddMatrix to add at term in a
//!       container. For now only atomic term are added but sooner or latter a coarser chunk of data may be added at once (block
//!       strategy of Nicolas Chevaugeon).
//
//! <br/> The connection to the solver is done via a method given by this class :  getMemoryAccess
//!       In solver interface this method is called by the connectMatrix method. getMemoryAccess give pointer that will be connect
//!       to the specific place the solver is using to access input matrix. Some how this can be considered as a violation of the
//!       private space of the class. But it is the best way to avoid data duplication to give access to the matrix for a solver
//!       by pointer.
//!
//! <br/>template parameter of the class are of 2 kinds :
//! <br/>
//! <br/>     type characterization :
//! <br/>                  T for double,float,complex, ...
//! <br/>     matrix traits characterization :
//! <br/>                  DEFINED for  defined, defined positive not defined (use by the solver)
//! <br/>                  PATTERN for  lower sym triangle, unsym, .... (use by this class, by some solver interface and by
//! assembler) <br/>                  STORAGE_TYPE Compressed Row storage, Compressed Column Storage, Coordinate storage,... (use
//! by this class and by some matrix graph class) <br/>                  INDEXING for  C or Fortran indexation in this container
//! (starting from 0 or 1) (use by this class, by some solver interface and by some matrix graph class)
//!
//! <br/>
//! Bellow is given the list of actual implemented method against template parameter and some useful information.
//! <br/> This table concern T=double.
//! <table border="1" cellpadding="2" cellspacing="2" width="100%" style="text-align:center">
//! <tr ><td>  ST </td><td> I </td><td>   P     </td><td>      index1      </td><td>    index2       </td><td>     extra1
//! </td><td>      solver     </td><td> AddMatrix </td><td> gemv </td><td> printMatrixMarket </td>><td> fillZeroDiag </td></tr>
//! <tr ><td> COO </td><td> F </td><td>  UNSYM   </td><td> row index (nnz) </td><td> col index (nnz) </td><td> col ptr (m+1)
//! </td><td>      mumps      </td><td>     X     </td><td>      </td><td>         X         </td><td>         ?     </td></tr>
//! <tr ><td> COO </td><td> F </td><td> LOWERSYM </td><td> row index (nnz) </td><td> col index (nnz) </td><td> col ptr (m+1)
//! </td><td>      mumps      </td><td>     X     </td><td>  X   </td><td>         X         </td><td>         X     </td></tr>
//! <tr ><td> CSC </td><td> F </td><td> LOWERSYM </td><td> col ptr (m+1)   </td><td> row index (nnz) </td><td> </td><td> pastix
//! </td><td>     X     </td><td>  X   </td><td>         X         </td><td>         X     </td></tr> <tr ><td> CSC </td><td> C
//! </td><td> LOWERSYM </td><td> col ptr (m+1)   </td><td> row index (nnz) </td><td>               </td><td>      taucs </td><td>
//! X     </td><td>  X   </td><td>         X         </td><td>         X     </td></tr> <tr ><td> CSC </td><td> C </td><td>  UNSYM
//! </td><td> col ptr (m+1)   </td><td> row index (nnz) </td><td>               </td><td> superLU/umfpack </td><td>     X
//! </td><td>  X   </td><td>         X         </td><td>         ?     </td></tr> <tr ><td> CSR </td><td> C </td><td>  UNSYM
//! </td><td> row ptr (m+1)   </td><td> col index (nnz) </td><td>               </td><td>     superLU     </td><td>     X
//! </td><td>  X   </td><td>         X         </td><td>         X     </td></tr>
//! </table>
//! <br/> where :
//! <ul><li>   ST = STORAGE_TYPE template parameter where :
//! <ul><li>                COO = coordinate storage ; index and value are stored
//! </li><li>               CSC = Compressed Sparse Column
//! </li><li>               CSR = Compressed Sparse Row
//! </li></ul>
//! </li><li>  I = INDEXING template parameter where :
//! <ul><li>                F = Fortran indexing (start at 1)
//! </li><li>               C = C indexing (start at 0)
//! </li></ul>
//! </li><li>  P = PATTERN template parameter where :
//! <ul><li>                UNSYM = all non zero term are stored
//! <ul><li>                LOWERSYM = all non zero term of the lower triangular are stored considering matrix as symmetric
//! (aij=aji)
//! </li></ul>
//! </li><li>  X = if present this method is implemented for storage,pattern and index of the line
//! </li><li>  ? = if present this method is implemented for storage,pattern and index of the line but not verified with a test
//! </li><li>  index1,index2,extra1 = member of this class which store structure of the matrix
//! </li></ul>
//
//
template <typename T, typename DEFINED, typename PATTERN, typename STORAGE_TYPE, typename INDEXING>
class xGenericSparseMatrix
{
  public:
   //
   // Traits :
   //
   //   defined directly by the user
   typedef T matrix_data_type;
   typedef PATTERN matrix_pattern;
   typedef DEFINED matrix_defined;  // for now used only with Mumps and Pastix
   //
   //   defined directly by the user which must be coherent with future connection to real
   //   storage space
   typedef STORAGE_TYPE matrix_storage;
   typedef INDEXING matrix_indexing;
   //
   //   this one is static as logically we want a fixed structure with this kind of matrix
   //   as separation of symbolic and numerical phase impose constant structure
   typedef xTraitMatrixForcedAssemblyOnZero matrix_assembly_on_zero;

   // constructor/destructor

   /// Constructor using a precomputed graph to set internal structure of the matrix
   //! The construction of the matrix graph structure is done outside this container and given to it via a argument of this
   //! constructor All graph information needed by this containers are used at construction time. It is user responsibility to
   //! clean graph made outside this container after it's use in this constructor. For the container a transformed copy suitable
   //! for the solver involved, have been made from the graph
   //!
   //! <br/> It is important to note that it is the graph which give the sparse pattern of A. No assumption is made to see if this
   //! graph represent or not only none null value.
   //!
   //! <br/>  Requirement for the graph structure of type G is that it must provide the following methods :
   //! <ul><li>  getN() : give the number of lines of the matrix described by the graph
   //! </li><li>  getM() : give the number of columns of the matrix described by the graph
   //! </li><li>  getNNZ() : give the number of non zero terms of the matrix described by the graph
   //! </li><li>  getCol(k) : give indices of lines for column k where non zero term is present as a ordered packed array of int
   //! (return a pointer to this array)
   //! </li><li>  getSizeCol(k) : give number of non zero termes of column k (return a int corresponding to the size of the array
   //! given by getCol)
   //! </li><li>  isSym() : respond true if only half part of the graph is stored, otherwise false if full graph is stored
   //! </li></ul>
   //! <br/>   The graph itself is considered as a sorted structure stored by column accessed by method above. Index of lines are
   //! given in C indexing and argument of "Col" function are
   //!         also in C indexing. If the graph represent a symmetric matrix only lower half of it have to be given (in this case
   //!         line index start from diagonal or below)
   //!
   template <typename G>
   xGenericSparseMatrix(G &graph);

   /// Constructor using an "other" matrix and extract a subset of it to construct a new matrix.
   //! Two vector are used to filter row and column with the following encoding :
   //! <ul><li>  value strictly greater then 0  => keep row/column
   //! <ul><li>  null value => remove row/column
   //! <ul><li>  value strictly less then 0 => merge row/column with "-value-1" row/column. In this case row/column does not exist
   //! anymore
   //! </li></ul>
   //!
   //! <br/>
   //! If sel_n is a null pointer then all row are selected
   //! If sel_m is a null pointer then all column are selected
   //! "other" matrix is considered of any kind with use of XXXO template parameter
   //! If both sel vector are null this constructor transcript "other"  from it's traits to the one of the constructor
   //! If both sel vector are null and "other"  got same template parameter has the constructor this is just equivalent to a copy
   //! constructor <br/> The resulting matrix is renumbered according to new dimension if filtration is involved and pack is true
   //! (default). It follow original order of "other" matrix numeration starting at 1/0. For example in fortran indexing if term
   //! (3,4) is selected as the first term of submatrix then it's index become (1,1). And if term (3,6) is the next to be selected
   //! when inspecting row then in submatrix it's index is (1,2). This imply that it's user responsibility to keep numbering
   //! tracking in between original and extracted matrix, if needed. If pack is false extracted submatrix terms are kept embedded
   //! in a matrix of the same size as the given matrix "other". It correspond more to a filtration of terms in a matrix then
   //! generating a new submatrix. This permit to keep matrix vector operation identical between "other" and genarated matrix if
   //! needed. <br/> When merging is involved it is user responsibility to insure that all given target position given in negated
   //! Fortran indexing correspond to retained row/column. Says sel_n[-sel_n[i]-1]==1 for all i with sel_n[i]<0. Same holds for
   //! sel_m. Anyway in debug mode assert are checking this condition. <br/> Be careful with traits transcription as not all
   //! combination have been tested and some may be bugged or impossible. See test case to see how it works and what have been
   //! tested. From those test you'll see that you can extract : <ul><li>       a non symmetric rectangle or square block from any
   //! part of lower sym or unsym matrix
   //! </li><li>      a symmetric square block around diagonal of lower sym matrix
   //! </li><li>      a symmetric square block from above or below diagonal of lower sym or unsym matrix (be careful to what you
   //! do here)
   //! </li></ul>
   //!
   //! Be careful also to pack set to false. It have almost not been tested.
   //!
   template <typename TO, typename DEFINEDO, typename PATTERNO, typename STORAGE_TYPEO, typename INDEXINGO>
   xGenericSparseMatrix(const xGenericSparseMatrix<TO, DEFINEDO, PATTERNO, STORAGE_TYPEO, INDEXINGO> &other, int *sel_n,
                        int *sel_m, bool pack = true);

   /// Constructor using an "other" matrix to construct a new equal or transposed matrix.
   //! Compare to generic copy constructor with filtering this constructor is giving a new matrix of the
   //! same type as "other" matrix. This version offers also a transpose functionality for non symmetric matrix.
   xGenericSparseMatrix(const xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING> &other, bool trans);

   ~xGenericSparseMatrix();

   // public methods ///////////////

   /// method needed to connect this matrix to a solver or whatever
   void getMemoryAccess(int &ns, int &ms, int &nnzs, int **idx1, int **idx2, T **dat, STORAGE_TYPE &storage, INDEXING &indexing);

   /// reseting values to zero
   void resetMatrixToZero();

   /// storing intermediate state of matrix values
   void saveState();

   /// restoring last saved state of matrix values
   void restorState();

   /// do Y = alpha.A.X + beta.Y  where A is this matrix
   //! For now this use a algorithm rather correct (mainly for CSC) but not necessary optimal as
   //! it  is coming from a sparse product of a full term matrix. Here we adapt
   //! it to symmetric case when only half matrix is given.
   //! comment the macro "#define PREFETCH 1" if you encounter problems while compiling
   //! this will remove __builtin_prefetch which is non portable (built-in functions provided by GCC and maybe other compiler)
   //! for unsym matrix it should be ok
   //! For CSR storage it have to be reviewed (using for example openMP implementation of N. Chevaugeon), for now it's a basic
   //! ugly implementation size_X  must be at least of size m (this mean that alpha.A.X is only done with X sub-vector 1:m )
   //! size_Y  must be at least of size n (this mean that alpha.A.X is only added to Y sub-vector of size 1:n; beta scaling is
   //! done on full Y) sparse_X is a boolean indicating if X is sparse (having few non null term) or not. By default X is
   //! considered as dense. If set to true user tell gemv method to test value of X before doing operation (i.e. if Xi is null
   //! nothing is done) Only few implementation is taking into account this parameter transpose is a boolean indicating if X has
   //! to be transposed in gemv.
   void gemv(int size_X, int size_Y, T alpha, T beta, T *X, T *Y, bool sparse_X = false, bool transpose = false);

   ///  do Y = alpha.A.X + beta.Y  where A is this matrix and X,Y xDistVector instance
   //! The set of local non null column of A must be included in the set of local X index
   //! The set of local non null row of A must be included in the set of local Y index
   void gemv(T alpha, T beta, const xDistVector<T> &X, xDistVector<T> &Y, bool sparse_X = false, bool transpose = false);

   /// do C = alpha.A.op(B) + beta.C  where A is this matrix and  B and C are dense matrix given by column (column 1 followed by
   /// column 2 ,....)
   //! op(B) = Bt or B depending on parameter trans_B = T or N (in fact only op(B)=B is possible/implemented for now)
   //! nb_col_OPB is the number of column of B and C
   //! nb_row_OPB must be at least of size m (this mean that product A.op(B) is only done with sub-bloc  m x nb_col_OPB of op(B))
   //! nb_row_C must be at least of size n (this imply adding to sub-bloc n x nb_col_OPB of C the product A.op(B) ; beta scaling
   //! is done one full matrix C ) see gemm blas level3 documentation for ldb,ldc parameter definition experimental : use at your
   //! own risk
   void gemm(char *trans_B, int nb_row_OPB, int nb_col_OPB, int nb_row_C, int ldb, int ldc, T alpha, T beta, T *B, T *C);

   /// do C = alpha.A.B + beta.C  where A is this matrix and  B and C are instance of this class (with their own traits)
   //! C is reshaped on output and have a pattern coresponding to pattern of product of A by B plus original pattern of C if beta
   //! not null Use algorithm : Efficient Sparse Matrix-Matrix Products Using Colorings, 2013, M. McCourt and B. F. Smith and H.
   //! Zhang (Argonne)
   template <typename DEFINEDB, typename PATTERNB, typename STORAGE_TYPEB, typename INDEXINGB, typename DEFINEDC,
             typename PATTERNC, typename STORAGE_TYPEC, typename INDEXINGC>
   void gemm(T alpha, T beta, const xGenericSparseMatrix<T, DEFINEDB, PATTERNB, STORAGE_TYPEB, INDEXINGB> &B,
             xGenericSparseMatrix<T, DEFINEDC, PATTERNC, STORAGE_TYPEC, INDEXINGC> &C);

   /// Do C = A oper B where :
   //! A must be able to contain (dimension point of view) B some how
   //! C must have the same dimension as A and is replaced by the result the operation
   //! oper is a given function that take 2 T and return a T. The first
   //! T argument is the value of A at a (i,j) location. The second double argument
   //! correspond the value of B at the same (i,j) location. The return value correspond to
   //! the operation that you want in between those two terms and it is stored in C at the (i,j) location.
   //! It could be a addition, a multiplication, a replacement , a tensorial
   //! product (if T is of  tensor type and you add appropriate xTraitsGenericSparseMatrixType specialization) ..... whatever you
   //! want. Note that oper is call for any position where at least one of A or B term is not null. You cannot create new term not
   //! already present in A or B. A symbolic phase create non nul term in C. This is not taking into account the fact that oper
   //! may create zero value.
   //!
   //! offset_r,offset_c are offset from location (0,0) of A to place B and do the operation. If dimensions of A are nA=10 x
   //! mA=20. If dimensions of B are nB=5 x nB=5 and offset_r=5,offset_c=15. Then B will be viewed by A as a matrix 10x20 with all
   //! it's term in the bloc [0 4]x[0 19] and [5 9]x[0 14] null and terms in bloc [5 9]x[15 19] corresponding to given B terms.

   template <typename DEFINEDB, typename PATTERNB, typename STORAGE_TYPEB, typename INDEXINGB, typename DEFINEDC,
             typename PATTERNC, typename STORAGE_TYPEC, typename INDEXINGC>
   void AoperB(int offset_r, int offset_c, const xGenericSparseMatrix<T, DEFINEDB, PATTERNB, STORAGE_TYPEB, INDEXINGB> &B,
               xGenericSparseMatrix<T, DEFINEDC, PATTERNC, STORAGE_TYPEC, INDEXINGC> &C,
               std::function<T(const T &, const T &)> oper);

   //
   /// do Y = alpha.K-1.X + beta.Y  where K-1 is the inverse of this matrix
   //! WARNING : this is only possible if numerical factorization of A have been made before
   //! WARNING : this is a quite restrictive implementation as xCSRVector are of double type !!!!!!!!
   // !  No check made ..... ugly and dangerous .... to do ....
   template <typename S>
   void gemvfac(S &, int, double, double, xCSRVector &, xCSRVector &, xCSRVector &);

   // here by convention i,j is in Fortran indexing (comes from xfem Double manager which is in Fortran indexing)
   void AddMatrix(const int &, const int &, const T &val);

   // here i,j is in traits indexing
   T &operator()(const int &i, const int &j);
   T operator()(const int &i, const int &j) const;

   //! print in Matrix Market format data of this container
   void printMatrixMarket(std::ostream &os);

   /// When matrix have explicitly null equation this method set diagonal zero term to a arbitrary value using "val"
   //! The  null equation is then replace by a arbitrary isolated equation which will give arbitrary solution
   //! depending on the way matrix is used
   //! "norm" indicate what value is set :
   //! <ul><li>   1 : Aii = val*norm1(A) = val * maximum absolute column sum of the matrix (to be done)
   //! </li><li>  2 : Aii = val*normInf(A) = val * maximum absolute row sum of the matrix (to be done)
   //! </li><li>  3 : Aii = val*normMax(A) = val * maximum absolute value of the matrix
   //! </li></ul>
   //! A automatic threshold (depending on chosen norm above) is calculated to consider diag term as null.
   //! If you have almost zero diag termes ( <threshold ) but this is not a null equation don't use this method
   //! or modify threshold
   void fillZeroDiag(const T &val, const int norm = 3);

   /// norm max of A :
   //! ||A|| =  maximum absolute value of the matrix
   T normMax();

   //! Scale the matrix values
   void scaleMatrix(T scale);

   //! Operator +=
   //! Warning this operator is unsafe but quick. Unsafe has it presuppose that pattern of calling
   //! instance and rhs is the same. Quick as this presupposition let this operation reduce to
   //! the addition of 2 vectors.
   void operator+=(const xGenericSparseMatrix<T, DEFINED, PATTERN, STORAGE_TYPE, INDEXING> &rhs);

   /// Store selected row and/or lines in packed form in a dense matrix of size nbrxnbc
   //! nbr number of row of the packed dense matrix
   //! nbc number of column of the packed dense matrix
   //! ids_r and ids_c contains information to filter and sort informations from sparse matrix to dense matrix :
   //! <ul><li>   There contains are -1 if row/column has to be skipped and index where to store in dense_matrix otherwise.
   //! <ul><li>   index given by those vector correspond  to the order chosen by the user. Says, sparse matrix columns 12 and 13
   //! may be the 5th and 3th ones respectively in dense matrix <ul><li>   if ids_r or ids_c are null corresponding dimension is
   //! fully taken (in sparse matrix order) and zero are stored in dense matrix <ul><li>   The only security test available for
   //! those vector is to check that given index is less then nbr or nbc
   //! </li></ul>
   //! <br/>
   //!  denseIndex is a function that, giving i and j argument coming from ids_r or ids_c, return the index of corresponding term
   //!  in array dense matrix. User is providing this function to sort dense matrix in row or column wise form.
   //
   void toDense(int nbr, int nbc, T *dense_matrix, int *ids_r, int *ids_c, std::function<int(int, int)> denseIndex);

   /// Pass in review all matrix terms in traits order. For each it apply functor given in argument to (i,j,x) triplet
   //! where i, j correspond to term indexes (in FORTRAN indexing) and x to term value.
   //! functor possible signature: void (int,int,const T&);
   template <typename F>
   void forEach(const F &func);

   //! Return the number of rows, columns and non-zero terms
   int getN() const;
   int getM() const;
   int getNNZ() const;

   //! Tells by inspecting template parameter PATTERN if graph for constructor have to be sym or not.
   static bool graphHaveToBeSym() { return xTraitsGenericSparseMatrixPattern<PATTERN>::symGraph(); }

   // public members ////////////////
   // none

   // Force all matrix type to be friend of each other
   // This is simplifying some methodes
   template <typename TO, typename DEFINEDO, typename PATTERNO, typename STORAGE_TYPEO, typename INDEXINGO>
   friend class xGenericSparseMatrix;

  protected:
   // private members ////////////////
   /// dimensions of the matrix (n for lines, m for columns) and number of non zero in the graph
   int n, m, nnz;
   /// Generic index container. May store index (i or j index of A(i,j)) or index pointer (starting/ending index of row/column
   /// i/j)
   //! interpretation of those generic index is subject to traits used
   int *index1;
   int *index2;
   /// In case of COO storage extra1 is used as a index pointer to keep a quick access to data
   std::vector<int> extra1;

   /// Container of non zero (of the graph) values of A
   T *data;

   /// extra storage used to keep data value when user ask for it (see storeMatrixValues and setMatrixValuesFromLastStore above)
   std::vector<T> stored_data;

   /// indicator of connection status
   bool connected_status;

   // private methods ///////////////
   template <typename G>
   void initMemoryWithGraph(G &g);
   template <typename... VARG>
   int gemvDispatch(const T &alpha, const T *X, T *Y, bool sparse_X, bool transpose, VARG &... args);
   // none
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xGenericSparseMatrix class
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xGenericSparseMatrixException class
///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// interface derived class of standard exception for xGenericSparseMatrix class
class xGenericSparseMatrixException : public std::exception
{
  public:
   xGenericSparseMatrixException(std::string, std::string, int, std::string, std::string);
   ~xGenericSparseMatrixException() throw() override;
   const char *what() const throw() override;

  private:
   std::string msg;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xGenericSparseMatrixException  class
//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace xlinalg

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xGenericSparseMatrix implementation
///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "xGenericSparseMatrix_imp.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xGenericSparseMatrix implementation
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
