/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
/*  -- ScaLAPACK routine (version 1.x) (x range from 1 to 8) --
 *     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
 *     and University of California, Berkeley.
 *     Jan 30, 2006
 */

#ifndef _SCALAPACKDEF_
#define _SCALAPACKDEF_

// ==================================================================================================================================
// Initialize the Process Grid
// ===========================
/*
   ICTXT
    (global output) INTEGER
    ICTXT specifies the BLACS context identifying the created process grid.
   NPROW
    (global input) INTEGER
    NPROW specifies the number of process rows in the process grid to be created.
   NPCOL
    (global input) INTEGER
    NPCOL specifies the number of process columns in the process grid to be created.
 */
extern "C" void sl_init_(int* ICTXT, int *NPROW, int *NPCOL);

// ==================================================================================================================================
// BLACS 
// ======

extern "C" void  Cblacs_get(int ConTxt, int what, int *val);
extern "C" int   Cblacs_gridinit(int *ConTxt, char *order, int nprow, int npcol);
extern "C" void  Cblacs_gridexit(int ConTxt);
extern "C" void  Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);


// ==================================================================================================================================
// Array descriptor initialization
// ===============================
/*
 *  DESC    (output) INTEGER array of dimension DLEN_.
 *          The array descriptor of a distributed matrix to be set.
 *
 *  M       (global input) INTEGER
 *          The number of rows in the distributed matrix. M >= 0.
 *
 *  N       (global input) INTEGER
 *          The number of columns in the distributed matrix. N >= 0.
 *
 *  MB      (global input) INTEGER
 *          The blocking factor used to distribute the rows of the
 *          matrix. MB >= 1.
 *
 *  NB      (global input) INTEGER
 *          The blocking factor used to distribute the columns of the
 *          matrix. NB >= 1.
 *
 *  IRSRC   (global input) INTEGER
 *          The process row over which the first row of the matrix is
 *          distributed. 0 <= IRSRC < NPROW.
 *
 *  ICSRC   (global input) INTEGER
 *          The process column over which the first column of the
 *          matrix is distributed. 0 <= ICSRC < NPCOL.
 *
 *  ICTXT   (global input) INTEGER
 *          The BLACS context handle, indicating the global context of
 *          the operation on the matrix. The context itself is global.
 *
 *  LLD     (local input)  INTEGER
 *          The leading dimension of the local array storing the local
 *          blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)).
 *
 *  INFO    (output) INTEGER
 *          = 0: successful exit
 *          < 0: if INFO = -i, the i-th argument had an illegal value
 *
 */
extern "C" void descinit_( int*DESCA, int*M, int*N, int*MB, int*NB, int*IRSRC, int *ICSRC, int *ICTXT, int *LLD, int *INFO );

// ==================================================================================================================================
// numerical
// =========

// BLAS LEVEL 1
// ============
/*
 *  Purpose
 *  =======
 *
 *  PDAXPY  adds one subvector to another,
 *
 *     sub( Y ) := sub( Y ) + alpha * sub( X ),
 *
 *  where
 *
 *     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
 *                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X, and,
 *
 *     sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
 *                      Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
 *
 *  Notes
 *  =====
 *
 *  A description  vector  is associated with each 2D block-cyclicly dis-
 *  tributed matrix.  This  vector  stores  the  information  required to
 *  establish the  mapping  between a  matrix entry and its corresponding
 *  process and memory location.
 *
 *  In  the  following  comments,   the character _  should  be  read  as
 *  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
 *  block cyclicly distributed matrix.  Its description vector is DESC_A:
 *
 *  NOTATION         STORED IN       EXPLANATION
 *  ---------------- --------------- ------------------------------------
 *  DTYPE_A (global) DESCA[ DTYPE_ ] The descriptor type.
 *  CTXT_A  (global) DESCA[ CTXT_  ] The BLACS context handle, indicating
 *                                   the NPROW x NPCOL BLACS process grid
 *                                   A  is  distributed over. The context
 *                                   itself  is  global,  but  the handle
 *                                   (the integer value) may vary.
 *  M_A     (global) DESCA[ M_     ] The  number of rows in the distribu-
 *                                   ted matrix A, M_A >= 0.
 *  N_A     (global) DESCA[ N_     ] The number of columns in the distri-
 *                                   buted matrix A, N_A >= 0.
 *  IMB_A   (global) DESCA[ IMB_   ] The number of rows of the upper left
 *                                   block of the matrix A, IMB_A > 0.
 *  INB_A   (global) DESCA[ INB_   ] The  number  of columns of the upper
 *                                   left   block   of   the  matrix   A,
 *                                   INB_A > 0.
 *  MB_A    (global) DESCA[ MB_    ] The blocking factor used to  distri-
 *                                   bute the last  M_A-IMB_A  rows of A,
 *                                   MB_A > 0.
 *  NB_A    (global) DESCA[ NB_    ] The blocking factor used to  distri-
 *                                   bute the last  N_A-INB_A  columns of
 *                                   A, NB_A > 0.
 *  RSRC_A  (global) DESCA[ RSRC_  ] The process row over which the first
 *                                   row of the matrix  A is distributed,
 *                                   NPROW > RSRC_A >= 0.
 *  CSRC_A  (global) DESCA[ CSRC_  ] The  process column  over  which the
 *                                   first column of  A  is  distributed.
 *                                   NPCOL > CSRC_A >= 0.
 *  LLD_A   (local)  DESCA[ LLD_   ] The  leading dimension  of the local
 *                                   array  storing  the  local blocks of
 *                                   the distributed matrix A,
 *                                   IF( Lc( 1, N_A ) > 0 )
 *                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
 *                                   ELSE
 *                                      LLD_A >= 1.
 *
 *  Let K be the number of  rows of a matrix A starting at the global in-
 *  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
 *  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
 *  receive if these K rows were distributed over NPROW processes.  If  K
 *  is the number of columns of a matrix  A  starting at the global index
 *  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
 *  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
 *  these K columns were distributed over NPCOL processes.
 *
 *  The values of Lr() and Lc() may be determined via a call to the func-
 *  tion PB_Cnumroc:
 *  Lr( IA, K ) = PB_Cnumroc( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
 *  Lc( JA, K ) = PB_Cnumroc( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
 *
 *  Arguments
 *  =========
 *
 *  N       (global input) INTEGER.
 *          On entry,  N  specifies the  length of the  subvectors to  be
 *          added. N must be at least zero.
 *
 *  ALPHA   (global input) DOUBLE PRECISION
 *          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *          supplied as zero then the local entries of the array  X  cor-
 *          responding to the entries of the subvector sub( X ) need  not
 *          be set on input.
 *
 *  X       (local input) DOUBLE PRECISION array
 *          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
 *          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
 *          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
 *          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
 *          Before  entry,  this array  contains the local entries of the
 *          matrix X.
 *
 *  IX      (global input) INTEGER
 *          On entry, IX  specifies X's global row index, which points to
 *          the beginning of the submatrix sub( X ).
 *
 *  JX      (global input) INTEGER
 *          On entry, JX  specifies X's global column index, which points
 *          to the beginning of the submatrix sub( X ).
 *
 *  DESCX   (global and local input) INTEGER array
 *          On entry, DESCX  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix X.
 *
 *  INCX    (global input) INTEGER
 *          On entry,  INCX   specifies  the  global  increment  for  the
 *          elements of  X.  Only two values of  INCX   are  supported in
 *          this version, namely 1 and M_X. INCX  must not be zero.
 *
 *  Y       (local input/local output) DOUBLE PRECISION array
 *          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
 *          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
 *          MAX( 1, Lr( 1, IY+N-1 ) )  otherwise,  and,  Ky  is  at least
 *          Lc( 1, JY+N-1 )  when  INCY = M_Y  and Lc( 1, JY ) otherwise.
 *          Before  entry,  this  array contains the local entries of the
 *          matrix Y.  On  exit, sub( Y ) is overwritten with the updated
 *          subvector.
 *
 *  IY      (global input) INTEGER
 *          On entry, IY  specifies Y's global row index, which points to
 *          the beginning of the submatrix sub( Y ).
 *
 *  JY      (global input) INTEGER
 *          On entry, JY  specifies Y's global column index, which points
 *          to the beginning of the submatrix sub( Y ).
 *
 *  DESCY   (global and local input) INTEGER array
 *          On entry, DESCY  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix Y.
 *
 *  INCY    (global input) INTEGER
 *          On entry,  INCY   specifies  the  global  increment  for  the
 *          elements of  Y.  Only two values of  INCY   are  supported in
 *          this version, namely 1 and M_Y. INCY  must not be zero.
 */
extern "C" void  pdaxpy_ (int *N, double *ALPHA, double *X, int *IX, int *JX, int *DESCX, int *INCX, double *Y, int *IY, int *JY, int *DESCY, int *INCY);

// BLAS LEVEL 2
// ============
/*
 *  Purpose
 *  =======
 *
 *  PDGEMV  performs one of the matrix-vector operations
 *
 *     sub( Y ) := alpha*sub( A ) *sub( X )  + beta*sub( Y ),  or
 *     sub( Y ) := alpha*sub( A )'*sub( X )  + beta*sub( Y ),
 *
 *  where
 *
 *     sub( A ) denotes A(IA:IA+M-1,JA:JA+N-1).
 *
 *  When TRANS = 'N',
 *
 *     sub( X ) denotes X(IX:IX,JX:JX+N-1), if INCX = M_X,
 *                      X(IX:IX+N-1,JX:JX), if INCX = 1 and INCX <> M_X,
 *     and,
 *
 *     sub( Y ) denotes Y(IY:IY,JY:JY+M-1), if INCY = M_Y,
 *                      Y(IY:IY+M-1,JY:JY), if INCY = 1 and INCY <> M_Y,
 *  and, otherwise
 *
 *     sub( X ) denotes X(IX:IX,JX:JX+M-1), if INCX = M_X,
 *                      X(IX:IX+M-1,JX:JX), if INCX = 1 and INCX <> M_X,
 *     and,
 *
 *     sub( Y ) denotes Y(IY:IY,JY:JY+N-1), if INCY = M_Y,
 *                      Y(IY:IY+N-1,JY:JY), if INCY = 1 and INCY <> M_Y.
 *
 *  Alpha and beta are scalars, and sub( X ) and sub( Y ) are  subvectors
 *  and sub( A ) is an m by n submatrix.
 *
 *  Notes
 *  =====
 *
 *  A description  vector  is associated with each 2D block-cyclicly dis-
 *  tributed matrix.  This  vector  stores  the  information  required to
 *  establish the  mapping  between a  matrix entry and its corresponding
 *  process and memory location.
 *
 *  In  the  following  comments,   the character _  should  be  read  as
 *  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
 *  block cyclicly distributed matrix.  Its description vector is DESC_A:
 *
 *  NOTATION         STORED IN       EXPLANATION
 *  ---------------- --------------- ------------------------------------
 *  DTYPE_A (global) DESCA[ DTYPE_ ] The descriptor type.
 *  CTXT_A  (global) DESCA[ CTXT_  ] The BLACS context handle, indicating
 *                                   the NPROW x NPCOL BLACS process grid
 *                                   A  is  distributed over. The context
 *                                   itself  is  global,  but  the handle
 *                                   (the integer value) may vary.
 *  M_A     (global) DESCA[ M_     ] The  number of rows in the distribu-
 *                                   ted matrix A, M_A >= 0.
 *  N_A     (global) DESCA[ N_     ] The number of columns in the distri-
 *                                   buted matrix A, N_A >= 0.
 *  IMB_A   (global) DESCA[ IMB_   ] The number of rows of the upper left
 *                                   block of the matrix A, IMB_A > 0.
 *  INB_A   (global) DESCA[ INB_   ] The  number  of columns of the upper
 *                                   left   block   of   the  matrix   A,
 *                                   INB_A > 0.
 *  MB_A    (global) DESCA[ MB_    ] The blocking factor used to  distri-
 *                                   bute the last  M_A-IMB_A  rows of A,
 *                                   MB_A > 0.
 *  NB_A    (global) DESCA[ NB_    ] The blocking factor used to  distri-
 *                                   bute the last  N_A-INB_A  columns of
 *                                   A, NB_A > 0.
 *  RSRC_A  (global) DESCA[ RSRC_  ] The process row over which the first
 *                                   row of the matrix  A is distributed,
 *                                   NPROW > RSRC_A >= 0.
 *  CSRC_A  (global) DESCA[ CSRC_  ] The  process column  over  which the
 *                                   first column of  A  is  distributed.
 *                                   NPCOL > CSRC_A >= 0.
 *  LLD_A   (local)  DESCA[ LLD_   ] The  leading dimension  of the local
 *                                   array  storing  the  local blocks of
 *                                   the distributed matrix A,
 *                                   IF( Lc( 1, N_A ) > 0 )
 *                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
 *                                   ELSE
 *                                      LLD_A >= 1.
 *
 *  Let K be the number of  rows of a matrix A starting at the global in-
 *  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
 *  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
 *  receive if these K rows were distributed over NPROW processes.  If  K
 *  is the number of columns of a matrix  A  starting at the global index
 *  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
 *  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
 *  these K columns were distributed over NPCOL processes.
 *
 *  The values of Lr() and Lc() may be determined via a call to the func-
 *  tion PB_Cnumroc:
 *  Lr( IA, K ) = PB_Cnumroc( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
 *  Lc( JA, K ) = PB_Cnumroc( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
 *
 *  Arguments
 *  =========
 *
 *  TRANS   (global input) CHARACTER*1
 *          On entry,  TRANS  specifies the  operation to be performed as
 *          follows:
 *
 *             TRANS = 'N' or 'n'
 *                sub( Y ) := alpha*sub( A )  * sub( X ) + beta*sub( Y ),
 *
 *             TRANS = 'T' or 't',
 *                sub( Y ) := alpha*sub( A )' * sub( X ) + beta*sub( Y ),
 *
 *             TRANS = 'C' or 'c',
 *                sub( Y ) := alpha*sub( A )' * sub( X ) + beta*sub( Y ).
 *
 *  M       (global input) INTEGER
 *          On entry,  M  specifies the number of rows of  the  submatrix
 *          sub( A ). M  must be at least zero.
 *
 *  N       (global input) INTEGER
 *          On entry, N  specifies the number of columns of the submatrix
 *          sub( A ). N  must be at least zero.
 *
 *  ALPHA   (global input) DOUBLE PRECISION
 *          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *          supplied  as  zero  then  the  local entries of the arrays  A
 *          and X corresponding to the entries of the submatrix  sub( A )
 *          and the subvector sub( X ) need not be set on input.
 *
 *  A       (local input) DOUBLE PRECISION array
 *          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
 *          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
 *          the local entries of the matrix A.
 *
 *  IA      (global input) INTEGER
 *          On entry, IA  specifies A's global row index, which points to
 *          the beginning of the submatrix sub( A ).
 *
 *  JA      (global input) INTEGER
 *          On entry, JA  specifies A's global column index, which points
 *          to the beginning of the submatrix sub( A ).
 *
 *  DESCA   (global and local input) INTEGER array
 *          On entry, DESCA  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix A.
 *
 *  X       (local input) DOUBLE PRECISION array
 *          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
 *          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
 *          MAX( 1, Lr( 1, IX+Lx-1 ) )  otherwise, and,  Kx  is  at least
 *          Lc( 1, JX+Lx-1 )  when INCX = M_X  and Lc( 1, JX ) otherwise.
 *          Lx is N when TRANS = 'N' or 'n' and  M  otherwise. Before en-
 *          try, this array  contains the local entries of the  matrix X.
 *
 *  IX      (global input) INTEGER
 *          On entry, IX  specifies X's global row index, which points to
 *          the beginning of the submatrix sub( X ).
 *
 *  JX      (global input) INTEGER
 *          On entry, JX  specifies X's global column index, which points
 *          to the beginning of the submatrix sub( X ).
 *
 *  DESCX   (global and local input) INTEGER array
 *          On entry, DESCX  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix X.
 *
 *  INCX    (global input) INTEGER
 *          On entry,  INCX   specifies  the  global  increment  for  the
 *          elements of  X.  Only two values of  INCX   are  supported in
 *          this version, namely 1 and M_X. INCX  must not be zero.
 *
 *  BETA    (global input) DOUBLE PRECISION
 *          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
 *          supplied  as  zero  then  the  local entries of  the array  Y
 *          corresponding to the entries of the subvector  sub( Y )  need
 *          not be set on input.
 *
 *  Y       (local input/local output) DOUBLE PRECISION array
 *          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
 *          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
 *          MAX( 1, Lr( 1, IY+Ly-1 ) )  otherwise, and,  Ky  is at  least
 *          Lc( 1, JY+Ly-1 )  when  INCY = M_Y and Lc( 1, JY ) otherwise.
 *          Ly is  M  when TRANS = 'N' or 'n' and N otherwise. Before en-
 *          try, this  array  contains the local entries of the matrix Y.
 *          On exit, sub( Y ) is overwritten by the updated subvector.
 *
 *  IY      (global input) INTEGER
 *          On entry, IY  specifies Y's global row index, which points to
 *          the beginning of the submatrix sub( Y ).
 *
 *  JY      (global input) INTEGER
 *          On entry, JY  specifies Y's global column index, which points
 *          to the beginning of the submatrix sub( Y ).
 *
 *  DESCY   (global and local input) INTEGER array
 *          On entry, DESCY  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix Y.
 *
 *  INCY    (global input) INTEGER
 *          On entry,  INCY   specifies  the  global  increment  for  the
 *          elements of  Y.  Only two values of  INCY   are  supported in
 *          this version, namely 1 and M_Y. INCY  must not be zero.
 */
extern "C" void  pdgemv_ (char* TRANS, int *M, int *N, double *ALPHA, double *A, int *IA, int *JA, int *DESCA, double *X, int *IX, int *JX, int *DESCX, int *INCX, double *BETA, double *Y, int *IY, int *JY, int *DESCY, int *INCY);


// Parallel lapack
// ===============
/*
 * Purpose
 *  =======
 *
 *  PDGESV computes the solution to a real system of linear equations
 *
 *                        sub( A ) * X = sub( B ),
 *
 *  where sub( A ) = A(IA:IA+N-1,JA:JA+N-1) is an N-by-N distributed
 *  matrix and X and sub( B ) = B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS
 *  distributed matrices.
 *
 *  The LU decomposition with partial pivoting and row interchanges is
 *  used to factor sub( A ) as sub( A ) = P * L * U, where P is a permu-
 *  tation matrix, L is unit lower triangular, and U is upper triangular.
 *  L and U are stored in sub( A ). The factored form of sub( A ) is then
 *  used to solve the system of equations sub( A ) * X = sub( B ).
 *
 *  Notes
 *  =====
 *
 *  Each global data object is described by an associated description
 *  vector.  This vector stores the information required to establish
 *  the mapping between an object element and its corresponding process
 *  and memory location.
 *
 *  Let A be a generic term for any 2D block cyclicly distributed array.
 *  Such a global array has an associated description vector DESCA.
 *  In the following comments, the character _ should be read as
 *  "of the global array".
 *
 *  NOTATION        STORED IN      EXPLANATION
 *  --------------- -------------- --------------------------------------
 *  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
 *                                 DTYPE_A = 1.
 *  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
 *                                 the BLACS process grid A is distribu-
 *                                 ted over. The context itself is glo-
 *                                 bal, but the handle (the integer
 *                                 value) may vary.
 *  M_A    (global) DESCA( M_ )    The number of rows in the global
 *                                 array A.
 *  N_A    (global) DESCA( N_ )    The number of columns in the global
 *                                 array A.
 *  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
 *                                 the rows of the array.
 *  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
 *                                 the columns of the array.
 *  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
 *                                 row of the array A is distributed.
 *  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
 *                                 first column of the array A is
 *                                 distributed.
 *  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
 *                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
 *
 *  Let K be the number of rows or columns of a distributed matrix,
 *  and assume that its process grid has dimension p x q.
 *  LOCr( K ) denotes the number of elements of K that a process
 *  would receive if K were distributed over the p processes of its
 *  process column.
 *  Similarly, LOCc( K ) denotes the number of elements of K that a
 *  process would receive if K were distributed over the q processes of
 *  its process row.
 *  The values of LOCr() and LOCc() may be determined via a call to the
 *  ScaLAPACK tool function, NUMROC:
 *          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
 *          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
 *  An upper bound for these quantities may be computed by:
 *          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
 *          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
 *
 *  This routine requires square block decomposition ( MB_A = NB_A ).
 *  Arguments
 *  =========
 *
 *  N       (global input) INTEGER
 *          The number of rows and columns to be operated on, i.e. the
 *          order of the distributed submatrix sub( A ). N >= 0.
 *
 *  NRHS    (global input) INTEGER
 *          The number of right hand sides, i.e., the number of columns
 *          of the distributed submatrix sub( B ). NRHS >= 0.
 *
 *  A       (local input/local output) DOUBLE PRECISION pointer into the
 *          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
 *          On entry, the local pieces of the N-by-N distributed matrix
 *          sub( A ) to be factored. On exit, this array contains the
 *          local pieces of the factors L and U from the factorization
 *          sub( A ) = P*L*U; the unit diagonal elements of L are not
 *          stored.
 *
 *  IA      (global input) INTEGER
 *          The row index in the global array A indicating the first
 *          row of sub( A ).
 *
 *  JA      (global input) INTEGER
 *          The column index in the global array A indicating the
 *          first column of sub( A ).
 *
 *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
 *          The array descriptor for the distributed matrix A.
 *
 *  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
 *          This array contains the pivoting information.
 *          IPIV(i) -> The global row local row i was swapped with.
 *          This array is tied to the distributed matrix A.
 *
 *  B       (local input/local output) DOUBLE PRECISION pointer into the
 *          local memory to an array of dimension
 *          (LLD_B,LOCc(JB+NRHS-1)).  On entry, the right hand side
 *          distributed matrix sub( B ). On exit, if INFO = 0, sub( B )
 *          is overwritten by the solution distributed matrix X.
 *
 *  IB      (global input) INTEGER
 *          The row index in the global array B indicating the first
 *          row of sub( B ).
 *
 *  JB      (global input) INTEGER
 *          The column index in the global array B indicating the
 *          first column of sub( B ).
 *
 *  DESCB   (global and local input) INTEGER array of dimension DLEN_.
 *          The array descriptor for the distributed matrix B.
 *
 *  INFO    (global output) INTEGER
 *          = 0:  successful exit
 *          < 0:  If the i-th argument is an array and the j-entry had
 *                an illegal value, then INFO = -(i*100+j), if the i-th
 *                argument is a scalar and had an illegal value, then
 *                INFO = -i.
 *          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
 *                The factorization has been completed, but the factor U
 *                is exactly singular, so the solution could not be
 *                computed.
 *
 */
extern "C" void pdgesv_( int *N, int *NRHS, double *A, int *IA, int *JA, int *DESCA, int *IPIV, double *B, int *IB, int *JB, int *DESCB, int *INFO );

/*
 *  Purpose
 *  =======
 *
 *  PDPOSV computes the solution to a real system of linear equations
 *
 *                        sub( A ) * X = sub( B ),
 *
 *  where sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1) and is an N-by-N
 *  symmetric distributed positive definite matrix and X and sub( B )
 *  denoting B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS distributed
 *  matrices.
 *
 *  The Cholesky decomposition is used to factor sub( A ) as
 *
 *                     sub( A ) = U**T * U,  if UPLO = 'U', or
 *
 *                     sub( A ) = L * L**T,  if UPLO = 'L',
 *
 *  where U is an upper triangular matrix and L is a lower triangular
 *  matrix.  The factored form of sub( A ) is then used to solve the
 *  system of equations.
 *
 *  Notes
 *  =====
 *
 *  Each global data object is described by an associated description
 *  vector.  This vector stores the information required to establish
 *  the mapping between an object element and its corresponding process
 *  and memory location.
 *
 *  Let A be a generic term for any 2D block cyclicly distributed array.
 *  Such a global array has an associated description vector DESCA.
 *  In the following comments, the character _ should be read as
 *  "of the global array".
 *
 *  NOTATION        STORED IN      EXPLANATION
 *  --------------- -------------- --------------------------------------
 *  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
 *                                 DTYPE_A = 1.
 *  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
 *                                 the BLACS process grid A is distribu-
 *                                 ted over. The context itself is glo-
 *                                 bal, but the handle (the integer
 *                                 value) may vary.
 *  M_A    (global) DESCA( M_ )    The number of rows in the global
 *                                 array A.
 *  N_A    (global) DESCA( N_ )    The number of columns in the global
 *                                 array A.
 *  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
 *                                 the rows of the array.
 *  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
 *                                 the columns of the array.
 *  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
 *                                 row of the array A is distributed.
 *  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
 *                                 first column of the array A is
 *                                 distributed.
 *  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
 *                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
 *
 *  Let K be the number of rows or columns of a distributed matrix,
 *  and assume that its process grid has dimension p x q.
 *  LOCr( K ) denotes the number of elements of K that a process
 *  would receive if K were distributed over the p processes of its
 *  process column.
 *  Similarly, LOCc( K ) denotes the number of elements of K that a
 *  process would receive if K were distributed over the q processes of
 *  its process row.
 *  The values of LOCr() and LOCc() may be determined via a call to the
 *  ScaLAPACK tool function, NUMROC:
 *          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
 *          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
 *  An upper bound for these quantities may be computed by:
 *          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
 *          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
 *
 *  This routine requires square block decomposition ( MB_A = NB_A ).
 *
 *  Arguments
 *  =========
 *
 *  UPLO    (global input) CHARACTER
 *          = 'U':  Upper triangle of sub( A ) is stored;
 *          = 'L':  Lower triangle of sub( A ) is stored.
 *
 *  N       (global input) INTEGER
 *          The number of rows and columns to be operated on, i.e. the
 *          order of the distributed submatrix sub( A ). N >= 0.
 *
 *  NRHS    (global input) INTEGER
 *          The number of right hand sides, i.e., the number of columns
 *          of the distributed submatrix sub( B ). NRHS >= 0.
 *
 *  A       (local input/local output) DOUBLE PRECISION pointer into the
 *          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
 *          On entry, this array contains the local pieces of the
 *          N-by-N symmetric distributed matrix sub( A ) to be factored.
 *          If UPLO = 'U', the leading N-by-N upper triangular part of
 *          sub( A ) contains the upper triangular part of the matrix,
 *          and its strictly lower triangular part is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          sub( A ) contains the lower triangular part of the distribu-
 *          ted matrix, and its strictly upper triangular part is not
 *          referenced. On exit, if INFO = 0, this array contains the
 *          local pieces of the factor U or L from the Cholesky factori-
 *          zation sub( A ) = U**T*U or L*L**T.
 *
 *  IA      (global input) INTEGER
 *          The row index in the global array A indicating the first
 *          row of sub( A ).
 *
 *  JA      (global input) INTEGER
 *          The column index in the global array A indicating the
 *          first column of sub( A ).
 *
 *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
 *          The array descriptor for the distributed matrix A.
 *
 *  B       (local input/local output) DOUBLE PRECISION pointer into the
 *          local memory to an array of dimension (LLD_B,LOC(JB+NRHS-1)).
 *          On entry, the local pieces of the right hand sides distribu-
 *          ted matrix sub( B ). On exit, if INFO = 0, sub( B ) is over-
 *          written with the solution distributed matrix X.
 *
 *  IB      (global input) INTEGER
 *          The row index in the global array B indicating the first
 *          row of sub( B ).
 *
 *  JB      (global input) INTEGER
 *          The column index in the global array B indicating the
 *          first column of sub( B ).
 *
 *  DESCB   (global and local input) INTEGER array of dimension DLEN_.
 *          The array descriptor for the distributed matrix B.
 *
 *  INFO    (global output) INTEGER
 *          = 0:  successful exit
 *          < 0:  If the i-th argument is an array and the j-entry had
 *                an illegal value, then INFO = -(i*100+j), if the i-th
 *                argument is a scalar and had an illegal value, then
 *                INFO = -i.
 *          > 0:  If INFO = K, the leading minor of order K,
 *                A(IA:IA+K-1,JA:JA+K-1) is not positive definite, and
 *                the factorization could not be completed, and the
 *                solution has not been computed.
 */
extern "C" void pdposv_( char *UPLO, int *N, int *NRHS, double *A, int *IA, int *JA, int *DESCA, double *B, int *IB, int *JB, int *DESCB, int *INFO);


// ==================================================================================================================================
// tools
// =====

extern "C" void Cpdgemr2d( int m, int n, double *ptrmyblock, int ia, int ja, int *ma, double *ptrmynewblock, int ib, int jb, int * mb, int globcontext );

#endif
