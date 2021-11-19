/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#ifndef _GENERICSPARSEMATRIXTRAITPOLICY_H
#ifdef _XTRAITPOLICYSOLVERMATRIX_H
#define _GENERICSPARSEMATRIXTRAITPOLICY_H


#define PREFETCH 1

namespace xlinalg
{

// RHST traits
struct xTraitRHSFull { };
struct xTraitRHSDistIndex {};

//------------------------------------------------------ Indexing dependent
template < typename INDEXING >
class xPolicyGenericSparseMatrixIndexing;

template < >
class xPolicyGenericSparseMatrixIndexing < xTraitMatrixCindex >
{
    public:
        static int toC(const int & i) {return i; }
        static int CtoI(const int & i) {return i; }
};
template < >
class xPolicyGenericSparseMatrixIndexing < xTraitMatrixFindex >
{
    public:
        static int toC(const int & i) {return i-1; }
        static int CtoI(const int & i) {return i+1; }
};
//------------------------------------------------------ End Indexing dependent
//------------------------------------------------------ STORAGE_TYPE,INDEXING Dependant
template < typename STORAGE_TYPE, typename INDEXING >
class xPolicyGenericSparseMatrixGemvPtrIdx;

template < >
class xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCOO,xTraitMatrixFindex >
{
    public:
        xPolicyGenericSparseMatrixGemvPtrIdx (int*ptri_,int*indx_) : ptri(ptri_),indx(indx_) {}
        int ptr(int j) const
        {
            return ptri[j];
        }
        int idx(int i) const
        {
            return indx[i]-1;
        }
    private:
        int *ptri,*indx;
};
template < >
class xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixFindex >
{
    public:
        xPolicyGenericSparseMatrixGemvPtrIdx (int*ptri_,int*indx_) : ptri(ptri_),indx(indx_) {}
        int ptr(int j) const
        {
            return ptri[j]-1;
        }
        int idx(int i) const
        {
            return indx[i]-1;
        }
    private:
        int *ptri,*indx;
};
template < >
class xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSC,xTraitMatrixCindex >
{
    public:
        xPolicyGenericSparseMatrixGemvPtrIdx (int*ptri_,int*indx_) : ptri(ptri_),indx(indx_) {}
        int ptr(int j) const
        {
            return ptri[j];
        }
        int idx(int i) const
        {
            return indx[i];
        }
    private:
        int *ptri,*indx;
};
template < >
class xPolicyGenericSparseMatrixGemvPtrIdx < xTraitMatrixSparceCSR,xTraitMatrixCindex >
{
    public:
        xPolicyGenericSparseMatrixGemvPtrIdx (int*ptri_,int*indx_) : ptri(ptri_),indx(indx_) {}
        int ptr(int j) const
        {
            return ptri[j];
        }
        int idx(int i) const
        {
            return indx[i];
        }
    private:
        int *ptri,*indx;
};

//------------------------------------------------------ End STORAGE_TYPE,INDEXING Dependant
//------------------------------------------------------ INDEXING,RHST Dependant
template < typename INDEXING, typename RHST >
class xPolicyGenericSparseMatrixGemvVindx;

//--------RHS Full
template < >
class xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixFindex, xTraitRHSFull  >
{
    public:
        xPolicyGenericSparseMatrixGemvVindx (const int*idx_) : idx(idx_) {}
        int index(int i) const
        {
            return idx[i]-1;
        }
    private:
        const int *idx;
};
template < >
class xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex, xTraitRHSFull  >
{
    public:
        xPolicyGenericSparseMatrixGemvVindx (const int*idx_) : idx(idx_) {}
        int index(int i) const
        {
            return idx[i];
        }
    private:
        const int *idx;
};

//--------RHS xDistIndex
template < >
class xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixFindex, xTraitRHSDistIndex  >
{
    public:
        xPolicyGenericSparseMatrixGemvVindx (const int*idx_,const xDistIndex &didx_) : idx(idx_),didx(didx_) {}
        int index(int i) const
        {
            return didx.getPackedIndex(idx[i]);
        }
    private:
        const int *idx;
        const xDistIndex &didx;
};
template < >
class xPolicyGenericSparseMatrixGemvVindx < xTraitMatrixCindex, xTraitRHSDistIndex  >
{
    public:
        xPolicyGenericSparseMatrixGemvVindx (const int*idx_,const xDistIndex &didx_) : idx(idx_),didx(didx_) {}
        int index(int i) const
        {
            return didx.getPackedIndex(idx[i]+1);
        }
    private:
        const int *idx;
        const xDistIndex &didx;
};
//------------------------------------------------------ End INDEXING,RHST Dependant
//----------------------------------------- template helping function
template < typename SPARSEX, typename PROD, typename T, typename INDEXING, typename RHST, typename ID >
int  CscLowerSymColumn(int j
                       , int q
                       , T *A
                       , const T *X
                       , T *Y
                       , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,RHST > &XY_index
                       , const ID &ptr_idx
                       , const SPARSEX  &test_v_zero
                       , const PROD  &prod_alpha)
{

    // j,q in C indexing


    // local
    int i,k,i1,i2,i3,i4,l,kf,id0,id1;
    const int delta = 4;
    T vxj,symj,v1,v2,v3,v4;

    // initial row index
    i = ptr_idx.ptr(j);
    // final row index
    k = ptr_idx.ptr(j+1);

    // no term in column => nothing to do
    // in the xDistVector context if it append it is rather strange
    // It means that selected column from X index is null which mean either that
    // it is treated on another process or it is null in absolute. In both case
    // it means that this process doesn't contribute. In first case if an other
    // process got a non null part of this column it forcedly imply that
    // X index is not owned by at least one process (this one or the one having
    // non null term) which is not so easy to check.
    if (k == i) return 0;

    // value of X corresponding to this column
    vxj = X[q];

    // skip null terms if sparceX or nothing otherwise
    if (test_v_zero(vxj)) return 0;

    // size calculation
    l = ( k-i-1 )%delta;
    kf = k-l;

    // first term always treated apart as it is on diagonal
    i1 = ptr_idx.idx(i);
    if (i1 != j)             //for now error if not diag term
    {
        std::cout <<" i1="<<i1<<" j="<<j<<std::endl;
        return -1;
    }
    symj = A[i]*vxj;

    // eventualy scale vxj by alpha
    prod_alpha(vxj);

    // loop on rows by delta
    for (++i; i < kf; i += delta)
    {
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( A+i ),0,3);
#endif
        i1 = XY_index.index(i);
        id1 = i+1;
        i2 = XY_index.index(id1);
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( X+i1 ),0,1);
        __builtin_prefetch((void * ) ( X+i2 ),0,1);
#endif
        v1 = vxj*A[i];
        symj += X[i1]*A[i];
        v2 = vxj*A[id1];
        symj += X[i2]*A[id1];
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( Y+i1 ),1,2);
        __builtin_prefetch((void * ) ( Y+i2 ),1,2);
#endif
        id0 = i+2;
        id1 = i+3;
        i3 = XY_index.index(id0);
        i4 = XY_index.index(id1);
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( X+i3 ),0,1);
#endif
        v3 = vxj*A[id0];
        symj += X[i3]*A[id0];
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( X+i4 ),0,1);
        __builtin_prefetch((void * ) ( Y+i3 ),1,1);
        __builtin_prefetch((void * ) ( Y+i4 ),1,1);
#endif
        v4 = vxj*A[id1];
        symj += X[i4]*A[id1];
        Y[i1] += v1;
        Y[i2] += v2;
        Y[i3] += v3;
        Y[i4] += v4;
    }

    // finish loop
    for (; i < k; ++i)
    {
        i1 = XY_index.index(i);
        v1 = vxj*A[i];
        symj += X[i1]*A[i];
        Y[i1] += v1;
    }

    // eventualy scale symj by alpha
    prod_alpha(symj);

    // Add sym contribution
    Y[q] += symj;

    return 0;
}
template < typename SPARSEX, typename PROD, typename T, typename INDEXING, typename RHST, typename ID >
int  CscUnSymColumn(int j
                    , int q
                    , T *A
                    , const T *X
                    , T *Y
                    , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,RHST > &Y_index
                    , const ID &ptr_idx
                    , const SPARSEX  &test_v_zero
                    , const PROD &prod_alpha)
{

    // j,q in C indexing


    // local
    int i,k,i1,i2,i3,i4,l,kf,id0,id1;
    const int delta = 4;
    T vxj,v1,v2,v3,v4;

    // initial row index
    i = ptr_idx.ptr(j);
    // final row index
    k = ptr_idx.ptr(j+1);

    // no term in column => nothing to do
    // in the xDistVector context if it append it is rather strange
    // It means that selected column from X index is null which mean either that
    // it is treated on another process or it is null in absolute. In both case
    // it means that this process doesn't contribute. In first case if an other
    // process got a non null part of this column it forcedly imply that
    // X index is not owned by at least one process (this one or the one having
    // non null term) which is not so easy to check.
    if (k == i) return 0;

    // value of X corresponding to this column
    vxj = X[q];

    // skip null terms if sparceX or nothing otherwise
    if (test_v_zero(vxj)) return 0;

    // size calculation
    l = ( k-i )%delta;
    kf = k-l;

    // eventualy scale vxj by alpha
    prod_alpha(vxj);

    // loop on rows by delta
    for (; i < kf; i += delta)
    {
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( A+i ),0,3);
#endif
        i1 = Y_index.index(i);
        id1 = i+1;
        i2 = Y_index.index(id1);
        v1 = vxj*A[i];
        v2 = vxj*A[id1];
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( Y+i1 ),1,2);
        __builtin_prefetch((void * ) ( Y+i2 ),1,2);
#endif
        id0 = i+2;
        id1 = i+3;
        i3 = Y_index.index(id0);
        i4 = Y_index.index(id1);
        v3 = vxj*A[id0];
#ifdef PREFETCH
        __builtin_prefetch((void * ) ( Y+i3 ),1,1);
        __builtin_prefetch((void * ) ( Y+i4 ),1,1);
#endif
        v4 = vxj*A[id1];
        Y[i1] += v1;
        Y[i2] += v2;
        Y[i3] += v3;
        Y[i4] += v4;
    }

    // finish loop
    for (; i < k; ++i)
    {
        i1 = Y_index.index(i);
        v1 = vxj*A[i];
        Y[i1] += v1;
    }

    return 0;
}
template < typename SPARSEX, typename PROD, typename T, typename INDEXING, typename RHST, typename ID >
int  CsrUnSymRow(int i
                 , int q
                 , T *A
                 , const T *X
                 , T *Y
                 , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,RHST > &X_index
                 , const ID &ptr_idx
                 , const SPARSEX  &test_v_zero
                 , const PROD  &prod_alpha)
{
    // i,q in C indexing
    // test_v_zero not used here but might be introduced in future

    // local
    int j,k;
    T dot_val;

    // initial column index
    j = ptr_idx.ptr(i);
    // final column index
    k = ptr_idx.ptr(i+1);

    // no term in row => nothing to do
    if (k == j) return 0;

    dot_val = xtool::xDataType < T >::zero();

    // loop on columns
    // NOTE: horrible, enrolling might help, packing+blas ?
    for (; j < k; ++j)
    {
        dot_val += A[j]*X[X_index.index(j)];
    }

    // eventually scale dot_val by alpha
    prod_alpha(dot_val);

    // add term to Y
    Y[q] += dot_val;

    return 0;
}
//------------------------------------------------------ Type,RHST Dependant
// ----- CSC LowerSym
template < typename PROD, typename SPARSEX, typename T, typename INDEXING, typename RHST >
class xPolicyGenericSparseMatrixGemvCscLowerSym
{};

template < typename PROD, typename SPARSEX, typename T, typename INDEXING >
class xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,INDEXING,xTraitRHSFull  >
{
    public:
        template < typename ID >
        static int  CscLowerSym(int n
                                , int m
                                , T *A
                                , const T *X
                                , T *Y
                                , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,xTraitRHSFull > &XY_index
                                , const ID &ptr_idx
                                , const SPARSEX  &test_v_zero
                                , const PROD  &prod_alpha)
        {

            // local
            int j,q;

            // for now implementation impose that sym imply n=m
            if (n != m)
            {
                std::cout <<" Symetric pattern is possible only with n=m"<<std::endl;
                return -1;
            }
            // loop on columns
            for (j = 0; j < m; ++j)
            {
                // j location in X or Y
                q = j;

                // compute columns
                if (  CscLowerSymColumn(j,q,A,X,Y,XY_index,ptr_idx,test_v_zero,prod_alpha) < 0 ) return -1;
            }
            return 0;
        }
};
template < typename PROD, typename SPARSEX, typename T, typename INDEXING >
class xPolicyGenericSparseMatrixGemvCscLowerSym < PROD,SPARSEX,T,INDEXING,xTraitRHSDistIndex >
{
    public:
        template < typename ID >
        static int  CscLowerSym(int n
                                , int m
                                , T *A
                                , const T *X
                                , T *Y
                                , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,xTraitRHSDistIndex > &XY_index
                                , const ID &ptr_idx
                                , const xDistIndex &idx_X
                                , const xDistIndex &idx_Y
                                , const SPARSEX  &test_v_zero
                                , const PROD  &prod_alpha)
        {

            // local
            int q;

            // for now implementation impose that sym imply n=m
            if (n != m)
            {
                std::cout <<" Symetric pattern is possible only with n=m"<<std::endl;
                return -1;
            }
            // for now implementation impose that sym imply idx_Y = idx_X
            // for now xDistIndex do not have != operator so test is done directly on adresses
            if (&idx_Y != &idx_X)
            {
                std::cout <<" Symetric pattern is possible only with  idx_Y = idx_X"<<std::endl;
                return -1;
            }
            // loop on columns
            for (auto j : idx_X)
            {
                // j in fortran indexing

                // j location in X or Y
                q = idx_X.getPackedIndex(j);

                // compute columns
                if (  CscLowerSymColumn(j-1,q,A,X,Y,XY_index,ptr_idx,test_v_zero,prod_alpha) < 0 ) return -1;
            }
            return 0;
        }
};
// ----- CSC unsym
template < typename PROD, typename SPARSEX, typename T, typename INDEXING, typename RHST >
class xPolicyGenericSparseMatrixGemvCscUnSym
{};

template < typename PROD, typename SPARSEX, typename T, typename INDEXING >
class xPolicyGenericSparseMatrixGemvCscUnSym < PROD,SPARSEX,T,INDEXING,xTraitRHSFull  >
{
    public:
        template < typename ID >
        static int  CscUnSym(int n
                             , int m
                             , T *A
                             , const T *X
                             , T *Y
                             , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,xTraitRHSFull > &Y_index
                             , const ID &ptr_idx
                             , const SPARSEX  &test_v_zero
                             , const PROD  &prod_alpha)
        {

            // local
            int j,q;

            // loop on columns
            for (j = 0; j < m; ++j)
            {
                // j location in X
                q = j;

                // compute columns
                if (  CscUnSymColumn(j,q,A,X,Y,Y_index,ptr_idx,test_v_zero,prod_alpha) < 0 ) return -1;
            }
            return 0;
        }
};
template < typename PROD, typename SPARSEX, typename T, typename INDEXING >
class xPolicyGenericSparseMatrixGemvCscUnSym < PROD,SPARSEX,T,INDEXING,xTraitRHSDistIndex >
{
    public:
        template < typename ID >
        static int  CscUnSym(int n
                             , int m
                             , T *A
                             , const T *X
                             , T *Y
                             , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,xTraitRHSDistIndex > &Y_index
                             , const ID &ptr_idx
                             , const xDistIndex &idx_X
                             , const xDistIndex &idx_Y
                             , const SPARSEX  &test_v_zero
                             , const PROD &prod_alpha)
        {

            // local
            int q;

            // loop on columns
            for (auto j : idx_X)
            {
                // j in fortran indexing

                // j location in X or Y
                q = idx_X.getPackedIndex(j);

                // compute columns
                if (  CscUnSymColumn(j-1,q,A,X,Y,Y_index,ptr_idx,test_v_zero,prod_alpha) < 0 ) return -1;
            }
            return 0;
        }
};
// ----- CSR unsym
template < typename PROD, typename SPARSEX, typename T, typename INDEXING, typename RHST >
class xPolicyGenericSparseMatrixGemvCsrUnSym
{};

template < typename PROD, typename SPARSEX, typename T, typename INDEXING >
class xPolicyGenericSparseMatrixGemvCsrUnSym < PROD,SPARSEX,T,INDEXING,xTraitRHSFull  >
{
    public:
        template < typename ID >
        static int  CsrUnSym(int n
                             , int m
                             , T *A
                             , const T *X
                             , T *Y
                             , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,xTraitRHSFull > &X_index
                             , const ID &ptr_idx
                             , const SPARSEX  &test_v_zero
                             , const PROD  &prod_alpha)
        {

            // local
            int i,q;

            // loop on rows
            for (i = 0; i < n; ++i)
            {
                // i location in Y
                q = i;

                // compute columns
                if (  CsrUnSymRow(i,q,A,X,Y,X_index,ptr_idx,test_v_zero,prod_alpha) < 0 ) return -1;
            }
            return 0;
        }
};
template < typename PROD, typename SPARSEX, typename T, typename INDEXING >
class xPolicyGenericSparseMatrixGemvCsrUnSym < PROD,SPARSEX,T,INDEXING,xTraitRHSDistIndex >
{
    public:
        template < typename ID>
        static int  CsrUnSym(int n
                             , int m
                             , T *A
                             , const T *X
                             , T *Y
                             , const xPolicyGenericSparseMatrixGemvVindx < INDEXING,xTraitRHSDistIndex > &X_index
                             , const ID &ptr_idx
                             , const xDistIndex &idx_X
                             , const xDistIndex &idx_Y
                             , const SPARSEX &test_v_zero
                             , const PROD &prod_alpha)
        {

            // local
            int q;

            // loop on rows
            for (auto i : idx_Y)
            {
                // i in fortran indexing

                // i location in Y
                q = idx_Y.getPackedIndex(i);

                // compute columns
                if (  CsrUnSymRow(i-1,q,A,X,Y,X_index,ptr_idx,test_v_zero,prod_alpha) < 0 ) return -1;
            }
            return 0;
        }
};
}
#else
#error "you should not include xGenericSparseMatrixGemvPolicy.h directly"
#endif
#endif
