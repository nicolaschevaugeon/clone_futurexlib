/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#ifndef _XTRAITPOLICYSOLVERMATRIX_H
#define _XTRAITPOLICYSOLVERMATRIX_H

#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include <memory>
#include "xTraitsMatrix.h"
#include "xBlasDef.h"
#include "xGraphMatrix.h"
#include "xDistVector.h"

// traits and policies for GEMV class
#include "xGenericSparseMatrixGemvPolicy.h"


namespace xlinalg
{

// local operator
int xGenericSparseMatrixOpIncrease (int i);

//------------------------------------------------------ Pattern dependante
template < typename PATTERN >
class xTraitsGenericSparseMatrixPattern
{
    public:
        static std::string spatern() { return "general"; }
        static bool symGraph() { return false; }
};

template < >
class xTraitsGenericSparseMatrixPattern < xTraitMatrixLowerSym >
{
    public:
        static std::string spatern() { return "symmetric"; }
        static bool symGraph() { return true; }
};
template < >
class xTraitsGenericSparseMatrixPattern < xTraitMatrixUpperSym >
{
    public:
        static std::string spatern() { return "symmetric"; }
        static bool symGraph() { return true; }
};
//------------------------------------------------------ End Pattern dependante

//------------------------------------------------------ Pattern,storage and indexing dependante
template < typename PATTERN, typename STORAGE_TYPE, typename INDEXING >
class xPolicyGenericSparseMatrix;


// CSC Compressed Sparce Column
// index1 = colone pointeur , size m+1
// index2 = row index , size nnz
// extra1 = unused
//
// Lower sym : only lower terms and diagonal are stored (i>=j)
//
// C indexing
//
template < >
class xPolicyGenericSparseMatrix < xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex >
{
    public:
        static void  allocMemory(int n,int m,int nnz,int **index1,int **index2,std::vector < int > &extra1)
        {
            ( *index1 ) = new int[m+1];
            ( *index2 ) = new int[nnz];
        }
        static int  getIndexFromFortran(int i, int j, int *ptr, int *rowidx, int *extra)
        {
            const int istart = ptr[j-1];
            const int iend = ptr[j];
            int *pos = std::lower_bound(&rowidx[istart], &rowidx[iend],i-1);
            return std::distance(rowidx, pos );
        }
        static int  getIndex(const int i, const int j, const int *ptr, const int *rowidx, const int *extra)
        {
            const int istart = ptr[j];
            const int iend = ptr[j+1];
            const int *pos = std::lower_bound(&rowidx[istart], &rowidx[iend],i);
            return std::distance(rowidx, pos );

        }
        static int  firstIndexEnd(int n_, int m_, int nnz_, int *colptr_, int *rowidx_)
        {
            rowidx_stat = rowidx_;
            colptr_stat = colptr_;
            return ( m_ );
        }
        static int  secondIndexStart(int k)
        {
            return( colptr_stat[k] );
        }
        static int  secondIndexEnd(int k)
        {
            return( colptr_stat[k+1] );
        }
        static int  indexI(int k, int l)
        {
            return( rowidx_stat[l]+1 );
        }
        static int  indexJ(int k, int l)
        {
            return( k+1 );
        }
        static int  indexVal(int k, int l)
        {
            return( l );
        }
        static void  initDiagSearch(int n_, int m_, int nnz_, int *colptr_, int *rowidx_, int *ptr)
        {
            rowidx_stat = rowidx_;
            colptr_stat = colptr_;
        }
        static int  diagIndex(int k)
        {
            const int idx = colptr_stat[k];
            if (rowidx_stat[idx] != k) return -1;
            else return( idx );
        }
        static const bool twoPass = false;
        static inline int setStorage(int k, int n, int *idx, int *ptr, int *ridx, int * extra1, int &offset)
        {
            // if null column
            if (!n)
            {
                // set column pointeur
                ptr[k] = offset;
                return( 0 );
            }

            // check consitancy betewn given graph and storage.
            // As graph is in lower column and here storage is of the same type the only consitancy test
            // possible is to verifie that i>=j assuming that index are correctely ordered (ascending)
            if (idx[0] < k) return( 1 );

            // set column pointeur
            ptr[k] = offset;

            // copy row index directely as graph is in C index
            memcpy ( (void *) ( &ridx[offset] ), (void *) ( &idx[0] ), sizeof( int )*n );

            // update offset
            offset += n;

            return( 0 );
        }
        static inline int setStoragePassTwo(int k, int n, int *idx, int *ptr, int *ridx, int * extra1, int &offset) {return( 0 ); }
        static inline void endSetStorage(int n, int m, int *ptr, int *ridx, int * extra1, int &offset)
        {
            // set column pointeur
            ptr[m] = offset;
        }
        static inline bool checkGraph(bool sym) {return ( !sym ); }

    private:
        static int *rowidx_stat;
        static int *colptr_stat;
};

// CSC Compressed Sparce Column
// index1 = colone pointeur , size m+1
// index2 = row index , size nnz
// extra1 = unused
//
// unsym : all term are stored
//
// C indexing
//
template < >
class xPolicyGenericSparseMatrix < xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex >
{
    public:
        static void  allocMemory(int n,int m,int nnz,int **index1,int **index2,std::vector < int > &extra1)
        {
            ( *index1 ) = new int[m+1];
            ( *index2 ) = new int[nnz];
        }
        static int  getIndexFromFortran(int i, int j, int *ptr, int *rowidx, int *extra)
        {
            const int istart = ptr[j-1];
            const int iend = ptr[j];
            int *pos = std::lower_bound(&rowidx[istart], &rowidx[iend],i-1);
            return std::distance(rowidx, pos );
        }
        static int  getIndex(const int i, const int j, const int *ptr, const int *rowidx, const int *extra)
        {
            const int istart = ptr[j];
            const int iend = ptr[j+1];
            const int *pos = std::lower_bound(&rowidx[istart], &rowidx[iend],i);
            return std::distance(rowidx, pos );

        }
        static int  firstIndexEnd(int n_, int m_, int nnz_, int *colptr_, int *rowidx_)
        {
            rowidx_stat = rowidx_;
            colptr_stat = colptr_;
            return ( m_ );
        }
        static int  secondIndexStart(int k)
        {
            return( colptr_stat[k] );
        }
        static int  secondIndexEnd(int k)
        {
            return( colptr_stat[k+1] );
        }
        static int  indexI(int k, int l)
        {
            return( rowidx_stat[l]+1 );
        }
        static int  indexJ(int k, int l)
        {
            return( k+1 );
        }
        static int  indexVal(int k, int l)
        {
            return( l );
        }
        static void  initDiagSearch(int n_, int m_, int nnz_, int *colptr_, int *rowidx_, int *ptr)
        {
            rowidx_stat = rowidx_;
            colptr_stat = colptr_;
        }
        static int  diagIndex(int k)
        {
            const int istart = colptr_stat[k];
            const int iend = colptr_stat[k+1];
            int *pos = std::lower_bound(&rowidx_stat[istart], &rowidx_stat[iend],k);
            const int d = std::distance(rowidx_stat, pos );
            if (rowidx_stat[d] != k) return -1;
            else return( d );
        }
        static const bool twoPass = false;
        static inline int setStorage(int k, int n, int *idx, int *ptr, int *ridx, int * extra1, int &offset)
        {
            // check consitancy betewn given graph and storage.
            // here no check possible


            // set column pointeur
            ptr[k] = offset;

            // if null column
            if (!n) return( 0 );

            // copy row index directely as graph is in C index
            memcpy ( (void *) ( &ridx[offset] ), (void *) ( &idx[0] ), sizeof( int )*n );

            // update offset
            offset += n;

            return( 0 );
        }
        static inline int setStoragePassTwo(int k, int n, int *idx, int *ptr, int *ridx, int * extra1, int &offset) {return( 0 ); }
        static inline void endSetStorage(int n,int m, int *ptr, int *ridx, int * extra1, int &offset)
        {
            // set column pointeur
            ptr[m] = offset;
        }
        static inline bool checkGraph(bool sym) {return ( sym ); }

    private:
        static int *rowidx_stat;
        static int *colptr_stat;
};

// CSC Compressed Sparce Column
// index1 = colone pointeur , size m+1
// index2 = row index , size nnz
// extra1 = unused
//
// Lower sym : only lower terms and diagonal are stored (i>=j)
//
// fortran indexing
//
template < >
class xPolicyGenericSparseMatrix < xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixFindex >
{
    public:
        static void  allocMemory(int n,int m,int nnz,int **index1,int **index2,std::vector < int > &extra1)
        {
            ( *index1 ) = new int[m+1];
            ( *index2 ) = new int[nnz];
        }
        static int  getIndexFromFortran(int i, int j, int *ptr, int *rowidx, int *extra)
        {
            return getIndex( i, j, ptr,rowidx,extra);
        }
        static int  getIndex(const int i, const int j, const int *ptr, const int *rowidx, const int *extra)
        {
            const int istart = ptr[j-1]-1;
            const int iend = ptr[j]-1;
            const int *pos = std::lower_bound(&rowidx[istart], &rowidx[iend],i);
            return std::distance(rowidx, pos );

        }
        static int  firstIndexEnd(int n_, int m_, int nnz_, int *colptr_, int *rowidx_)
        {
            rowidx_stat = rowidx_;
            colptr_stat = colptr_;
            return ( m_ );
        }
        static int  secondIndexStart(int k)
        {
            return( colptr_stat[k]-1 );
        }
        static int  secondIndexEnd(int k)
        {
            return( colptr_stat[k+1]-1 );
        }
        static int  indexI(int k, int l)
        {
            return( rowidx_stat[l] );
        }
        static int  indexJ(int k, int l)
        {
            return( k+1 );
        }
        static int  indexVal(int k, int l)
        {
            return( l );
        }
        static void  initDiagSearch(int n_, int m_, int nnz_, int *colptr_, int *rowidx_, int *ptr)
        {
            rowidx_stat = rowidx_;
            colptr_stat = colptr_;
        }
        static int  diagIndex(int k)
        {
            const int idx = colptr_stat[k]-1;
            if (rowidx_stat[idx] != k+1) return -1;
            else return( idx );
        }
        static const bool twoPass = false;
        static inline int setStorage(int k, int n, int *idx, int *ptr, int *ridx, int * extra1, int &offset)
        {
            // if null column
            if (!n)
            {
                // set column pointeur
                ptr[k] = offset+1;
                return( 0 );
            }
            // check consitancy betewn given graph and storage.
            // As graph is in lower column and here storage is of the same type the only consitancy test
            // possible is to verifie that i>=j assuming that index are correctely ordered (ascending)
            if (idx[0] < k) return( 1 );

            // set column pointeur (fortran)
            ptr[k] = offset+1;

            // copy row index and shift index as graph is in C index and here we want fortran indexing
            std::transform(&idx[0],&idx[n],&ridx[offset],xGenericSparseMatrixOpIncrease);

            // update offset
            offset += n;

            return( 0 );
        }
        static inline int setStoragePassTwo(int k, int n, int *idx, int *ptr, int *ridx, int * extra1, int &offset) {return( 0 ); }
        static inline void endSetStorage(int n, int m, int *ptr, int *ridx, int * extra1, int &offset)
        {
            // set column pointeur
            ptr[m] = offset+1;
        }
        static inline bool checkGraph(bool sym) {return ( !sym ); }
    private:
        static int *rowidx_stat;
        static int *colptr_stat;
};

// CSR Compressed Sparce Row
// index1 = row pointeur , size n+1
// index2 = colone index , size nnz
// extra1 = unused
//
// unsym : all term are stored
//
// C indexing
//
template < >
class xPolicyGenericSparseMatrix < xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixCindex >
{
    public:
        static void  allocMemory(int n,int m,int nnz,int **index1,int **index2,std::vector < int > &extra1)
        {
            int *ptr = ( *index1 ) = new int[n+1];
            ( *index2 ) = new int[nnz];

            // in this context as ptr will be used to count ligne size (number of nz per ligne)
            // it have to be set to zero before the first storage pass
            int i;
            for (i = 0; i < n; ++i) ptr[i] = 0;
            ptr[i] = 0;

        }
        static int  getIndexFromFortran(int i, int j, int *ptr, int *colidx, int *extra)
        {
            const int jstart = ptr[i-1];
            const int jend = ptr[i];
            int *pos = std::lower_bound(&colidx[jstart], &colidx[jend],j-1);
            return std::distance(colidx, pos );
        }
        static int  getIndex(const int i, const int j, const int *ptr, const int *colidx, const int *extra)
        {
            const int jstart = ptr[i];
            const int jend = ptr[i+1];
            const int *pos = std::lower_bound(&colidx[jstart], &colidx[jend],j);
            return std::distance(colidx, pos );

        }
        static int  firstIndexEnd(int n_, int m_, int nnz_, int *rowptr_, int *colidx_)
        {
            colidx_stat = colidx_;
            rowptr_stat = rowptr_;
            return ( n_ );
        }
        static int  secondIndexStart(int k)
        {
            return( rowptr_stat[k] );
        }
        static int  secondIndexEnd(int k)
        {
            return( rowptr_stat[k+1] );
        }
        static int  indexI(int k, int l)
        {
            return( k+1 );
        }
        static int  indexJ(int k, int l)
        {
            return( colidx_stat[l]+1 );
        }
        static int  indexVal(int k, int l)
        {
            return( l );
        }
        static void  initDiagSearch(int n_, int m_, int nnz_, int *rowptr_, int *colidx_, int *ptr)
        {
            colidx_stat = colidx_;
            rowptr_stat = rowptr_;
        }
        static int  diagIndex(int k)
        {
            const int jstart = rowptr_stat[k];
            const int jend = rowptr_stat[k+1];
            int *pos = std::lower_bound(&colidx_stat[jstart], &colidx_stat[jend],k);
            const int d = std::distance(colidx_stat, pos );
            if (colidx_stat[d] != k) return -1;
            else return( d );
        }
        static const bool twoPass = true;
        static inline int setStorage(int k, int n, int *idx, int *ptr, int *cidx, int * extra1, int &offset)
        {
            // graph and storage are in different order but with same indexing
            // swapping from column-wise to  row-wise
            // on this first pass we count ligne size (number of nz per ligne)

            int l;
            int *ptr_p_1 = ptr+1;

            for (l = 0; l < n; ++l)
            {
                ++ptr_p_1[idx[l]];
            }

            return( 0 );
        }
        static inline int setStoragePassTwo(int k, int n, int *idx, int *ptr, int *cidx, int * extra1, int &offset)
        {
            // local
            int i,l;

            // graph and storage are in different order but with same indexing
            // swapping from column-wise to  row-wise
            // on this second pass we insert in cidx correct value at correct place using ligne size previously
            // computed

            // before any use of ptr it have to be summed (once) to obtaine real offset
            if (offset)
            {
                for (l = 1; l < offset; ++l) ptr[l] += ptr[l-1];
                offset = 0;
            }

            for (l = 0; l < n; ++l)
            {

                i = idx[l];
                cidx[ptr[i]] = k;
                ++ptr[i];
            }

            return( 0 );
        }
        static inline void endSetStorage(int n, int m, int *ptr, int *cidx, int * extra1, int &offset)
        {
            // set row pointeur
            for (int l = n; l > 0; --l) ptr[l] = ptr[l-1];
            ptr[0] = 0;
        }
        static inline bool checkGraph(bool sym) {return ( sym ); }

    private:
        static int *colidx_stat;
        static int *rowptr_stat;
};

// COO Coordinate list
// index1 = row index , size nnz
// index2 = colone index , size nnz
// extra1 = colone pointeur (in c indexing whatever traits indexing is) for underlying CSC storage   , size m+1
//
// lower sym : only lower terms and diagonal are stored (i>=j)
//
// fortran indexing
//
template < >
class xPolicyGenericSparseMatrix < xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex >
{
    public:
        static void  allocMemory(int n,int m,int nnz,int **index1,int **index2,std::vector < int > &extra1)
        {
            ( *index1 ) = new int[nnz];
            ( *index2 ) = new int[nnz];
            extra1.resize(m+1);
        }
        static int  getIndexFromFortran(int i, int j, int *rowidx, int *colidx, int *ptr)
        {
            return getIndex( i, j, rowidx,colidx,ptr);
        }
        static int  getIndex(const int i, const int j, const int *rowidx, const int *colidx, const int *ptr)
        {
            const int istart = ptr[j-1];
            const int iend = ptr[j];
            const int *pos = std::lower_bound(&rowidx[istart], &rowidx[iend],i);
            return std::distance(rowidx, pos );

        }
        static int  firstIndexEnd(int n_, int m_, int nnz_, int *rowidx_, int *colidx_)
        {
            rowidx_stat = rowidx_;
            colidx_stat = colidx_;
            nnz_stat = nnz_;
            return ( 1 );
        }
        static int  secondIndexStart(int k)
        {
            return( 0 );
        }
        static int  secondIndexEnd(int k)
        {
            return( nnz_stat );
        }
        static int  indexI(int b, int k)
        {
            return( rowidx_stat[k] );
        }
        static int  indexJ(int b, int k)
        {
            return( colidx_stat[k] );
        }
        static int  indexVal(int b, int k)
        {
            return( k );
        }
        static void  initDiagSearch(int n_, int m_, int nnz_, int *rowidx_, int *colidx_, int *ptr)
        {
            rowidx_stat = rowidx_;
            colidx_stat = ptr; // here we use the underlying CSC storage and use colidx_stat for conveniancy. colptr_stat would have been a better name
        }
        static int  diagIndex(int k)
        {
            const int idx = colidx_stat[k];
            if (rowidx_stat[idx] != k+1) return -1;
            else return( idx );
        }
        static const bool twoPass = false;
        static inline int setStorage(int k, int n, int *idx, int *ridx, int *cidx, int * ptr, int &offset)
        {
            // if null column
            if (!n)
            {
                // set column pointeur (C indexing)
                ptr[k] = offset;
                return( 0 );
            }
            // check consitancy betewn given graph and storage.
            // As graph is in lower column and here storage is made the same type the only consitancy test
            // possible is to verifie that i>=j assuming that index are correctely ordered (ascending)
            if (idx[0] < k) return( 1 );

            // set column pointeur (C indexing)
            ptr[k] = offset;

            // copy row index and shift index as graph is in C index and here we want fortran indexing
            std::transform(&idx[0],&idx[n],&ridx[offset],xGenericSparseMatrixOpIncrease);

            // set increment for colone index generation
            int i = offset;

            // update offset
            offset += n;

            // generate colone index
            ++k;
            for (; i < offset; ++i) cidx[i] = k;

            return( 0 );

        }
        static inline int setStoragePassTwo(int k, int n, int *idx, int *ridx, int *cidx, int * ptr, int &offset) {return( 0 ); }
        static void endSetStorage(int n, int m, int *ridx, int *cidx, int * ptr, int &offset)
        {
            // set column pointeur
            ptr[m] = offset;
        }
        static inline bool checkGraph(bool sym) {return ( !sym ); }

    private:
        static int *rowidx_stat;
        static int *colidx_stat;
        static int nnz_stat;
};

// COO Coordinate list
// index1 = row index , size nnz
// index2 = colone index , size nnz
// extra1 = colone pointeur (in c indexing whatever traits indexing is) for underlying CSC storage   , size m+1
//
// unsym : all term are stored
//
// fortran indexing
//
template < >
class xPolicyGenericSparseMatrix < xTraitMatrixUnSym,xTraitMatrixSparceCOO,xTraitMatrixFindex >
{
    public:
        static void  allocMemory(int n,int m,int nnz,int **index1,int **index2,std::vector < int > &extra1)
        {
            ( *index1 ) = new int[nnz];
            ( *index2 ) = new int[nnz];
            extra1.resize(m+1);
        }
        static int  getIndexFromFortran(int i, int j, int *rowidx, int *colidx, int *ptr)
        {
            return getIndex( i, j, rowidx,colidx,ptr);
        }
        static int  getIndex(const int i, const int j, const int *rowidx, const int *colidx, const int *ptr)
        {
            const int istart = ptr[j-1];
            const int iend = ptr[j];
            const int *pos = std::lower_bound(&rowidx[istart], &rowidx[iend],i);
            return std::distance(rowidx, pos );
        }
        static int  firstIndexEnd(int n_, int m_, int nnz_, int *rowidx_, int *colidx_)
        {
            rowidx_stat = rowidx_;
            colidx_stat = colidx_;
            nnz_stat = nnz_;
            return ( 1 );
        }
        static int  secondIndexStart(int k)
        {
            return( 0 );
        }
        static int  secondIndexEnd(int k)
        {
            return( nnz_stat );
        }
        static int  indexI(int b, int k)
        {
            return( rowidx_stat[k] );
        }
        static int  indexJ(int b, int k)
        {
            return( colidx_stat[k] );
        }
        static int  indexVal(int b, int k)
        {
            return( k );
        }
        static void  initDiagSearch(int n_, int m_, int nnz_, int *rowidx_, int *colidx_, int *ptr)
        {
            rowidx_stat = rowidx_;
            colidx_stat = ptr; // here we use the underlying CSC storage and use colidx_stat for conveniancy. colptr_stat would have been a better name
        }
        static int  diagIndex(int k)
        {
            const int istart = colidx_stat[k];
            const int iend = colidx_stat[k+1];
            int *pos = std::lower_bound(&rowidx_stat[istart], &rowidx_stat[iend],k+1);
            const int d = std::distance(rowidx_stat, pos );
            if (rowidx_stat[d] != k+1) return -1;
            else return( d );
        }
        static const bool twoPass = false;
        static inline int setStorage(int k, int n, int *idx, int *ridx, int *cidx, int * ptr, int &offset)
        {
            // check consitancy betewn given graph and storage.
            // here no check possible

            // set column pointeur (C indexing)
            ptr[k] = offset;

            // if null column
            if (!n) return( 0 );

            // copy row index and shift index as graph is in C index and here we want fortran indexing
            std::transform(&idx[0],&idx[n],&ridx[offset],xGenericSparseMatrixOpIncrease);

            // set increment for colone index generation
            int i = offset;

            // update offset
            offset += n;

            // generate colone index
            ++k;
            for (; i < offset; ++i) cidx[i] = k;

            return( 0 );

        }
        static inline int setStoragePassTwo(int k, int n, int *idx, int *ridx, int *cidx, int * ptr, int &offset) {return( 0 ); }
        static void endSetStorage(int n, int m, int *ridx, int *cidx, int * ptr, int &offset)
        {
            // set column pointeur
            ptr[m] = offset;
        }
        static inline bool checkGraph(bool sym) {return ( sym ); }

    private:
        static int *rowidx_stat;
        static int *colidx_stat;
        static int nnz_stat;
};
//------------------------------------------------------ End Pattern,storage and indexing dependante

//------------------------------------------------------ Type,Pattern,storage and indexing dependante
template <  typename T, typename PATTERN, typename STORAGE_TYPE, typename INDEXING >
class xPolicyGenericSparseMatrixGemv;

// COO Coordinate list
// index1 = row index , size nnz
// index2 = colone index , size nnz
// extra1 = colone pointeur (in c indexing whatever traits indexing is) for underlying CSC storage   , size m+1
//
// lower sym : only lower terms and diagonal are stored (i>=j)
//
// T type product policy
//
// fortran indexing
//
template <  typename T >
class xPolicyGenericSparseMatrixGemv < T,xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex >
{
    public:
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *rowidx, int *colidx, int *ptr, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixFindex,xTraitRHSFull > X_index(rowidx);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCOO,xTraitMatrixFindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,xTraitMatrixFindex,xTraitRHSFull >::CscLowerSym(n,m,A,X,Y, X_index, ptr_idx, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *rowidx, int *colidx, int *ptr, const PROD &prod_alpha, const SPARSEX &test_v_zero )
        {
            // symetrique matrix => transposed==original => same function
            std::cout <<"xGenericSparseMatrix::gemv: Why using the transposed option with a symetric matrix ??"<<std::endl;
            return gemv1(n,m,A,X,Y,rowidx,colidx,ptr,prod_alpha,test_v_zero);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *rowidx, int *colidx, int *ptr, const PROD &prod_alpha, const SPARSEX & test_v_zero, const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixFindex,xTraitRHSDistIndex > X_index(rowidx,idx_X);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCOO,xTraitMatrixFindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,xTraitMatrixFindex,xTraitRHSDistIndex >::CscLowerSym(n,m,A,X,Y, X_index, ptr_idx,idx_X,idx_Y, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *rowidx, int *colidx, int *ptr, const PROD &prod_alpha, const SPARSEX & test_v_zero, const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            // symetrique matrix => transposed==original => same function
            std::cout <<"xGenericSparseMatrix::gemv: Why using the transposed option with a symetric matrix ??"<<std::endl;
            return gemv1(n,m,A,X,Y,rowidx,colidx,ptr,prod_alpha,test_v_zero,idx_X,idx_Y);
        }


};

// CSC Compressed Sparce Column
// index1 = colone pointeur , size n+1
// index2 = row index , size nnz
// extra1 = unused
//
// Lower sym : only lower terms and diagonal are stored (i>=j)
//
// T type, PROD product policy
//
// C indexing
//
template <  typename T >
class xPolicyGenericSparseMatrixGemv < T,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex >
{
    public:
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSFull > X_index(rowidx);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixCindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSFull >::CscLowerSym(n,m,A,X,Y, X_index, ptr_idx, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            // symetrique matrix => transposed==original => same function
            std::cout <<"xGenericSparseMatrix::gemv: Why using the transposed option with a symetric matrix ??"<<std::endl;
            return gemv1(n,m,A,X,Y,ptr,rowidx,extra1,prod_alpha,test_v_zero);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1,  const PROD &prod_alpha, const SPARSEX & test_v_zero,const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSDistIndex > X_index(rowidx,idx_X);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixCindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSDistIndex >::CscLowerSym(n,m,A,X,Y, X_index, ptr_idx,idx_X,idx_Y, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1,  const PROD &prod_alpha, const SPARSEX & test_v_zero,const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            // symetrique matrix => transposed==original => same function
            std::cout <<"xGenericSparseMatrix::gemv: Why using the transposed option with a symetric matrix ??"<<std::endl;
            return gemv1(n,m,A,X,Y,ptr,rowidx,extra1,prod_alpha,test_v_zero,idx_X,idx_Y);
        }
};

// CSR Compressed Sparce Row
// index1 = row pointeur , size n+1
// index2 = colone index , size nnz
// extra1 = unused
//
// unsym : all term are stored
//
// T type, PROD product policy
//
// C indexing
//
template <  typename T >
class xPolicyGenericSparseMatrixGemv < T,xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixCindex >
{
    public:
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *colidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSFull > X_index(colidx);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSR,xTraitMatrixCindex > ptr_idx(ptr,colidx);
            return xPolicyGenericSparseMatrixGemvCsrUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSFull >::CsrUnSym(n,m,A,X,Y, X_index, ptr_idx, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *colidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSFull > Y_index(colidx);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixCindex > ptr_idx(ptr,colidx);
            return xPolicyGenericSparseMatrixGemvCscUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSFull >::CscUnSym(m,n,A,X,Y, Y_index, ptr_idx, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *colidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero, const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSDistIndex > X_index(colidx,idx_X);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSR,xTraitMatrixCindex > ptr_idx(ptr,colidx);
            return xPolicyGenericSparseMatrixGemvCsrUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSDistIndex >::CsrUnSym(n,m,A,X,Y, X_index, ptr_idx,idx_X,idx_Y, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *colidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero, const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSDistIndex > Y_index(colidx,idx_Y);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSR,xTraitMatrixCindex > ptr_idx(ptr,colidx);
            return xPolicyGenericSparseMatrixGemvCscUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSDistIndex >::CscUnSym(m,n,A,X,Y, Y_index, ptr_idx,idx_X,idx_Y, test_v_zero,prod_alpha);
        }
};

// CSC Compressed Sparce Column
// index1 = colone pointeur , size n+1
// index2 = row index , size nnz
// extra1 = unused
//
// unsym : all term are stored
//
// T type, PROD product policy
//
// C indexing
//
template <  typename T >
class xPolicyGenericSparseMatrixGemv < T,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex >
{
    public:
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSFull > Y_index(rowidx);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixCindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSFull >::CscUnSym(n,m,A,X,Y, Y_index, ptr_idx, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSFull > X_index(rowidx);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSR,xTraitMatrixCindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCsrUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSFull >::CsrUnSym(m,n,A,X,Y, X_index, ptr_idx, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1,  const PROD &prod_alpha, const SPARSEX & test_v_zero,const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSDistIndex > Y_index(rowidx,idx_Y);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixCindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSDistIndex >::CscUnSym(n,m,A,X,Y, Y_index, ptr_idx,idx_X,idx_Y, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1,  const PROD &prod_alpha, const SPARSEX & test_v_zero,const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex,xTraitRHSDistIndex > X_index(rowidx,idx_X);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixCindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCsrUnSym < PROD,SPARSEX,T,xTraitMatrixCindex,xTraitRHSDistIndex >::CsrUnSym(m,n,A,X,Y, X_index, ptr_idx,idx_X,idx_Y, test_v_zero,prod_alpha);
        }
};

// CSC Compressed Sparce Column
// index1 = colone pointeur , size n+1
// index2 = row index , size nnz
// extra1 = unused
//
// Lower sym : only lower terms and diagonal are stored (i>=j)
//
// T type, PROD product policy
//
// F indexing
//
template <  typename T >
class xPolicyGenericSparseMatrixGemv < T,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixFindex >
{
    public:
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixFindex,xTraitRHSFull > X_index(rowidx);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixFindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,xTraitMatrixFindex,xTraitRHSFull >::CscLowerSym(n,m,A,X,Y, X_index, ptr_idx, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1, const PROD &prod_alpha, const SPARSEX & test_v_zero)
        {
            // symetrique matrix => transposed==original => same function
            std::cout <<"xGenericSparseMatrix::gemv: Why using the transposed option with a symetric matrix ??"<<std::endl;
            return gemv1(n,m,A,X,Y,ptr,rowidx,extra1,prod_alpha,test_v_zero);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1,  const PROD &prod_alpha, const SPARSEX & test_v_zero,const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixFindex,xTraitRHSDistIndex > X_index(rowidx,idx_X);
            xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixFindex > ptr_idx(ptr,rowidx);
            return xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,xTraitMatrixFindex,xTraitRHSDistIndex >::CscLowerSym(n,m,A,X,Y, X_index, ptr_idx,idx_X,idx_Y, test_v_zero,prod_alpha);
        }
        template <  typename PROD, typename SPARSEX >
        static int  gemv1Transposed(int n, int m, T *A, const T *X, T *Y, int *ptr, int *rowidx, int *extra1,  const PROD &prod_alpha, const SPARSEX & test_v_zero,const xDistIndex &idx_X, const xDistIndex &idx_Y)
        {
            // symetrique matrix => transposed==original => same function
            std::cout <<"xGenericSparseMatrix::gemv: Why using the transposed option with a symetric matrix ??"<<std::endl;
            return gemv1(n,m,A,X,Y,ptr,rowidx,extra1,prod_alpha,test_v_zero,idx_X,idx_Y);
        }
};
//------------------------------------------------------ End Prod,Type,Pattern,storage and indexing dependante

//------------------------------------------------------ Pattern dependante
class xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopyBase(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            sel_n(sel_n_),
            sel_m(sel_m_),
            id_pack_n(id_pack_n_),
            id_pack_m(id_pack_m_)
        {}
        virtual void addGraph(xGraphMatrix & g, int i, int j) const = 0;
        virtual void addMatrix(int i, int j, std::function < void (int,int) > add_m) const = 0;

        virtual ~xPolicyGenericSparseMatrixCopyBase() = default;
    protected:
        int * sel_n;
        int * sel_m;
        int * id_pack_n;
        int * id_pack_m;
};

template < typename PATTERN, typename PATTERNO, bool f_n, bool f_m >
class xPolicyGenericSparseMatrixCopy : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_){}
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            g.add(i,j);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            add_m(i+1,j+1);
            return;
        }

};
template < typename PATTERN, typename PATTERNO >
class xPolicyGenericSparseMatrixCopy < PATTERN,PATTERNO,true,true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_){}
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            i = sel_n[i];
            j = sel_m[j];
            if (i && j)
                g.add(id_pack_n[i],id_pack_m[j]);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            i = sel_n[i];
            j = sel_m[j];
            if (i && j)
                add_m(id_pack_n[i]+1,id_pack_m[j]+1);
            return;
        }
};
template < typename PATTERN, typename PATTERNO >
class xPolicyGenericSparseMatrixCopy < PATTERN,PATTERNO,true,false > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_){}
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            i = sel_n[i];
            if (i)
                g.add(id_pack_n[i],j);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            i = sel_n[i];
            if (i)
                add_m(id_pack_n[i]+1,j+1);
            return;
        }
};
template < typename PATTERN, typename PATTERNO >
class xPolicyGenericSparseMatrixCopy < PATTERN,PATTERNO,false,true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_){}
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            j = sel_m[j];
            if (j)
                g.add(i,id_pack_m[j]);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            j = sel_m[j];
            if (j)
                add_m(i+1,id_pack_m[j]+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixLowerSym,xTraitMatrixUnSym,true,false > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        // USE WITH CARE !
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            assert(sel_n);
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            i = sel_n[i];
            if ( i && ( id_pack_n[i] >= j ) )
                g.add(id_pack_n[i],j);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            i = sel_n[i];
            if ( i && ( id_pack_n[i] >= j ) )
                add_m(id_pack_n[i]+1,j+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixLowerSym,xTraitMatrixUnSym,false,true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        // USE WITH CARE !
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            assert(sel_m);
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            j = sel_m[j];
            if ( j && ( i >= id_pack_m[j] ) )
                g.add(i,id_pack_m[j]);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            j = sel_m[j];
            if ( j && ( i >= id_pack_m[j] ) )
                add_m(i+1,id_pack_m[j]+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixLowerSym,xTraitMatrixUnSym, true, true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        // USE WITH CARE !
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            assert(sel_n);
            assert(sel_m);
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            i = sel_n[i];
            j = sel_m[j];
            if ( i && j && ( id_pack_n[i] >= id_pack_m[j] ) )
                g.add(id_pack_n[i],id_pack_m[j]);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            i = sel_n[i];
            j = sel_m[j];
            if ( i && j && ( id_pack_n[i] >= id_pack_m[j] ) )
                add_m(id_pack_n[i]+1,id_pack_m[j]+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixUnSym,xTraitMatrixLowerSym, true, true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            assert(sel_n);
            assert(sel_m);
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            const int k = sel_m[i];
            const int l = sel_n[j];
            if ( ( i != j ) && k && l)
                g.add(id_pack_n[l],id_pack_m[k]);
            i = sel_n[i];
            j = sel_m[j];
            if (i && j)
                g.add(id_pack_n[i],id_pack_m[j]);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            const int k = sel_m[i];
            const int l = sel_n[j];
            if ( ( i != j ) && k && l)
                add_m(id_pack_n[l]+1,id_pack_m[k]+1);
            i = sel_n[i];
            j = sel_m[j];
            if (i && j)
                add_m(id_pack_n[i]+1,id_pack_m[j]+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixUnSym,xTraitMatrixLowerSym, true, false > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            assert(sel_n);
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            const int l = sel_n[j];
            if ( ( i != j ) && l)
                g.add(id_pack_n[l],i);
            i = sel_n[i];
            if (i)
                g.add(id_pack_n[i],j);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            const int l = sel_n[j];
            if ( ( i != j ) && l)
                add_m(id_pack_n[l]+1,i+1);
            i = sel_n[i];
            if (i)
                add_m(id_pack_n[i]+1,j+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixUnSym,xTraitMatrixLowerSym, false, true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            assert(sel_m);
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            const int k = sel_m[i];
            if ( ( i != j ) && k )
                g.add(j,id_pack_m[k]);
            j = sel_m[j];
            if (j)
                g.add(i,id_pack_m[j]);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            const int k = sel_m[i];
            if ( ( i != j ) && k )
                add_m(j+1,id_pack_m[k]+1);
            j = sel_m[j];
            if (j)
                add_m(i+1,id_pack_m[j]+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixLowerSym,xTraitMatrixLowerSym, true, true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            assert(sel_n);
            assert(sel_m);
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            const int k = sel_m[i];
            const int l = sel_n[j];
            if ( ( i != j ) && k && l && ( id_pack_n[l] >= id_pack_m[k] ) )
                g.add(id_pack_n[l],id_pack_m[k]);
            i = sel_n[i];
            j = sel_m[j];
            if ( i && j && ( id_pack_n[i] >= id_pack_m[j] ) )
                g.add(id_pack_n[i],id_pack_m[j]);
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            const int k = sel_m[i];
            const int l = sel_n[j];
            if ( ( i != j ) && k && l && ( id_pack_n[l] >= id_pack_m[k] ) )
                add_m(id_pack_n[l]+1,id_pack_m[k]+1);
            i = sel_n[i];
            j = sel_m[j];
            if ( i && j && ( id_pack_n[i] >= id_pack_m[j] ) )
                add_m(id_pack_n[i]+1,id_pack_m[j]+1);
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixLowerSym,xTraitMatrixLowerSym, true, false > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            // imposible ??
            throw -1;
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            // imposible ??
            throw -1;
            return;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            // imposible ??
            throw -1;
            return;
        }
};
template < >
class xPolicyGenericSparseMatrixCopy < xTraitMatrixLowerSym,xTraitMatrixLowerSym, false, true > : public xPolicyGenericSparseMatrixCopyBase
{
    public:
        xPolicyGenericSparseMatrixCopy(int *sel_n_,int * sel_m_,int * id_pack_n_,int * id_pack_m_) :
            xPolicyGenericSparseMatrixCopyBase(sel_n_,sel_m_,id_pack_n_,id_pack_m_)
        {
            // imposible ??
            throw -1;
        }
        void addGraph(xGraphMatrix & g, int i, int j) const override
        {
            // imposible ??
            throw -1;
        }
        void addMatrix(int i, int j,std::function < void (int,int) > add_m) const override
        {
            // imposible ??
            throw -1;
        }
};
class xPolicyDenseMatrixCopyBase
{
    public:
        xPolicyDenseMatrixCopyBase(int nbr_,int nbc_,int *ids_r_,int * ids_c_,std::function < int (int,int) > denseIndex_) :
            ids_r(ids_r_-1),
            ids_c(ids_c_-1),
            nbr(nbr_),nbc(nbc_),
            denseIndex(denseIndex_)
        {}
        virtual int denseSelectedIndex(int i, int j) const = 0;
        virtual ~xPolicyDenseMatrixCopyBase() = default;
    protected:
        int * ids_r;
        int * ids_c;
        int nbr,nbc;
        std::function < int (int,int) > denseIndex;
};

template < typename PATTERN, bool f_n, bool f_m >
class xPolicyDenseMatrixCopy : public xPolicyDenseMatrixCopyBase
{
    public:
        xPolicyDenseMatrixCopy(int nbr,int nbc,int *ids_r_,int * ids_c_,std::function < int (int,int) > denseIndex_) :
            xPolicyDenseMatrixCopyBase(nbr,nbc,ids_r_,ids_c_,denseIndex_){}
        // By not defining denseSelectedIndex method we keep the generic derided class
        // abstract so that user may only compile with approved specialize class.
        // May change in future

};
template < >
class xPolicyDenseMatrixCopy < xTraitMatrixUnSym,true,true > : public xPolicyDenseMatrixCopyBase
{
    public:
        xPolicyDenseMatrixCopy(int nbr,int nbc,int *ids_r_,int * ids_c_,std::function < int (int,int) > denseIndex_) :
            xPolicyDenseMatrixCopyBase(nbr,nbc,ids_r_,ids_c_,denseIndex_){}
        int denseSelectedIndex( int i, int j) const override
        {
            const int k = ids_r[i];
            const int l = ids_c[j];
            assert(k < nbr);
            assert(l < nbc);
            if (k > -1 && l > -1)
                return denseIndex(k,l);

            return -1;
        }

};

template < >
class xPolicyDenseMatrixCopy < xTraitMatrixUnSym,false,true > : public xPolicyDenseMatrixCopyBase
{
    public:
        xPolicyDenseMatrixCopy(int nbr,int nbc,int *ids_r_,int * ids_c_,std::function < int (int,int) > denseIndex_) :
            xPolicyDenseMatrixCopyBase(nbr,nbc,ids_r_,ids_c_,denseIndex_){}
        int denseSelectedIndex( int i, int j) const override
        {
            const int l = ids_c[j];
            assert(l < nbc);
            if (l > -1)
                return denseIndex(i-1,l);

            return -1;
        }

};
template < >
class xPolicyDenseMatrixCopy < xTraitMatrixUnSym,true,false > : public xPolicyDenseMatrixCopyBase
{
    public:
        xPolicyDenseMatrixCopy(int nbr,int nbc,int *ids_r_,int * ids_c_,std::function < int (int,int) > denseIndex_) :
            xPolicyDenseMatrixCopyBase(nbr,nbc,ids_r_,ids_c_,denseIndex_){}
        int denseSelectedIndex( int i, int j) const override
        {
            const int k = ids_r[i];
            assert(k < nbr);
            if (k > -1)
                return denseIndex(k,j-1);

            return -1;
        }

};
template < >
class xPolicyDenseMatrixCopy < xTraitMatrixUnSym,false,false > : public xPolicyDenseMatrixCopyBase
{
    public:
        xPolicyDenseMatrixCopy(int nbr,int nbc,int *ids_r_,int * ids_c_,std::function < int (int,int) > denseIndex_) :
            xPolicyDenseMatrixCopyBase(nbr,nbc,ids_r_,ids_c_,denseIndex_){}
        int denseSelectedIndex( int i, int j) const override
        {
            return denseIndex(i-1,j-1);
        }

};
//------------------------------------------------------ End Pattern dependante
} // end of namespace
#endif
