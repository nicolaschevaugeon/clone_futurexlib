/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */
#include <iostream>
#include <fstream>
#include <vector>
// xtool
#include "xMPIEnv.h"
#include "xDataType.h"
// xlinalg
#include "xGenericSparseMatrix.h"
#include "xGraphMatrix.h"
#include "xTraitsMatrix.h"

// setting: local vector/matrix dimension
#define SIZECL 13
#define SIZERL 7
#define BLOCKC 3
#define BLOCKR 2

using namespace std;
using namespace xtool;
using namespace xlinalg;

template < typename T >
void check1(int k,const char *s, ofstream &oc, vector < T > &R,vector < T > &S,int proc_id, MPI_Comm &world, bool verbose )
{
    if (verbose)
    {
        cout<<s<<"local test:"<<endl;
        for (int i = 0; i < k; ++i) cout<<" "<<R[i]<<","<<S[i];
        cout<<endl;
    }
    MPI_Allreduce(MPI_IN_PLACE, R.data(), k, xtool::xMPIDataType < T >(), MPI_SUM, world);
    if (!proc_id)
    {
        const T zero = xDataType < T >::zero();
        bool ko = false;
        for (int i = 0; i < k; ++i) ko += (( R[i]-S[i] ) != zero );
        oc<<s<<!ko<<endl;
    }
}
template < typename T >
void check2(int k, const char *s, ofstream &oc, vector < T > &RG, xDistVector < T > &R,vector < T > &S,int proc_id, MPI_Comm &world, bool verbose )
{
    if (verbose)
    {
        auto & ridx = R.getDistIndex();
        cout<<s<<"local test:"<<endl;
        for (auto i_f : ridx) cout<<" "<<R[ridx.getPackedIndex(i_f)]<<","<<S[i_f-1];
        cout<<endl;
    }
    R.gather(RG);
    if (!proc_id)
    {
        const T zero = xDataType < T >::zero();
        bool ko = false;
        for (int i = 0; i < k; ++i) ko += (( RG[i]-S[i] ) != zero );
        oc<<s<<!ko<<endl;
    }
}
void generateIdx(xDistIndex &idx,size_t size_per_proc,int proc_id, int nb_proc)
{
    // computation of local size
    // overlap 2
    assert(size_per_proc > 4);
    const size_t local_size = size_per_proc -2;
    const size_t offset = ( !proc_id ) ? 1 : local_size*proc_id+3;
    const size_t offset2 = offset-size_per_proc;

    size_t i;
    // add local index
    for (i = 0; i < local_size; ++i)
    {
        idx.insertIndex(i+offset);
    }
    // if proc 0 add last local index
    size_t start,stop,ls;
    if (!proc_id)
    {
        for (; i < size_per_proc; ++i)
        {
            idx.insertIndex(i+offset);
        }
        start = local_size;
        stop = size_per_proc;
        ls = size_per_proc;
    }
    // if proc > 0 add remote index (the last two from previous proc)
    else
    {
        for (; i < size_per_proc; ++i)
        {
            idx.insertIndex(i+offset2);
        }
        start = local_size-2;
        stop = local_size;
        ls = local_size;
    }
    // for all proc flag index remotely present in  proc+1 (if it exist)
    int proc_id_p_un = proc_id+1;
    if (proc_id_p_un < nb_proc)
    {
        for (i = start; i < stop; ++i)
        {
            idx.insertToFrom(i+offset,proc_id_p_un);
        }
    }
    int proc_id_m_un = proc_id-1;
    if (proc_id_m_un > -1)
    {
        for (i = local_size; i < size_per_proc; ++i)
        {
            idx.insertToFrom(i+offset2,proc_id_m_un);
        }
    }
    idx.finalize(ls,stop+offset-1,true);

}
template < typename T >
void test(const xDistIndex & idxc
          ,const xDistIndex & idxr
          ,MPI_Comm &worldt
          ,int nb_proct
          ,int proc_idt
          ,ofstream &oc
          ,bool verbose
          ,T C2
          ,T alpha
          ,T C34
          ,T C14
          )
{

    // constant
    T zero = xDataType < T >::zero();
    T one = xDataType < T >::one();
    T mone = -one;
    const int n = idxr.getGlobalIndexSize();
    const int m = idxc.getGlobalIndexSize();

    // sequential and distributed vector
    xDistVector < T > XD(idxc);
    xDistVector < T > XDS(idxc);
    xDistVector < T > XDSS(idxc);
    xDistVector < T > YDU(idxr);
    xDistVector < T > YDS(idxc);
    vector < T > X(m,one);
    vector < T > XS(m,zero);
    vector < T > YU(n,zero);
    vector < T > YS(m,zero);
    vector < T > SOLU(n,zero);
    vector < T > SOLUS(n,zero);
    vector < T > SOLS(m,zero);
    vector < T > SOLSS(m,zero);

    // graph construction
    xGraphMatrix gunsym(n,m,false);
    xGraphMatrix gsym(m,m,true);

    for (auto i : idxr)
    {
        --i;
        bool A = ( i/BLOCKR )%2;
        for (auto j : idxc)
        {
            bool B = ( --j/BLOCKC )%2;
            if (( A || B ) && !( A && B ))
                gunsym.add(i,j);
        }
    }
    gunsym.countNNZ();
    if (verbose) gunsym.print();

    for (auto i : idxc)
    {
        --i;
        bool A = !(( i/BLOCKC )%2 );
        for (auto j : idxc)
        {
            bool B = ( --j/BLOCKC )%2;
            if (( A || B ) && !( A && B ))
            {
                if (!( i < j )) gsym.add(i,j);
            }
        }
    }
    gsym.countNNZ();
    if (verbose) gsym.print();


    // matrix creation
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > Aus1(gunsym);
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixCindex > Aus2(gunsym);
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > Asy3(gsym);
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixFindex > Asy4(gsym);
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > Asy5(gsym);
    // matrix filling
    int ttc = idxc.getPackedLocalIndexSize();
    int ttr = idxr.getPackedLocalIndexSize();
    int tc = ( proc_idt+1 < nb_proct ) ? ttc-2 : idxc.getPackedIndexSize();
    int tr = ( proc_idt+1 < nb_proct ) ? ttr-2 : idxr.getPackedIndexSize();
    for (int j = 0; j < m; ++j)
    {
        int nbr = gunsym.getSizeCol(j);
        if (nbr)
        {
            int *col = gunsym.getCol(j);
            assert(col);
            T &valj = ( ( j/BLOCKC/2 )%2 ) ? mone : one;
            int j_f = j+1;
            int ct = idxc.getPackedIndex(j_f);
            bool A = ( ct < ttc );
            bool B = ( !( ct < tc ) && A );
            bool C = ( !A );
            if (C)
            {
                XD[ct] = C34;
                XDS[ct] = C34;
                XS[j] = one;
            }
            else if (B)
            {
                XD[ct] = C14;
                XDS[ct] = C14;
                XS[j] = one;
            }
            else
            {
                XD[ct] = one;
                XDS[ct] = zero;
            }
            for (int i = 0; i < nbr; ++i)
            {
                int i_f = col[i]+1;
                int rt = idxr.getPackedIndex(i_f);
                bool D = ( rt < ttr );
                T val = ( ( C && !D ) || ( B && ( !( rt < tr ) && D ))) ? valj/C2 : valj;
                Aus1.AddMatrix(i_f,j_f,val);
                Aus2.AddMatrix(i_f,j_f,val);
                SOLU[i_f-1] += val;
                if (C || B)
                    SOLUS[i_f-1] += val;
            }
        }
        nbr = gsym.getSizeCol(j);
        if (nbr)
        {
            int *col = gsym.getCol(j);
            assert(col);
            T &valj = ( ( j/BLOCKC/2 )%2 ) ? mone : one;
            int j_f = j+1;
            int ct = idxc.getPackedIndex(j_f);
            bool A = ( ct < ttc );
            bool B = ( !( ct < tc ) && A );
            bool C = ( !A );
            T sum = zero;
            T sums = zero;
            for (int i = 0; i < nbr; ++i)
            {
                int i_f = col[i]+1;
                int rt = idxc.getPackedIndex(i_f);
                bool D = ( rt < ttc );
                bool E = ( !( rt < tc ) && D );
                T val = ( ( C && !D ) || ( B && E)) ? valj/C2 : valj;
                Asy3.AddMatrix(i_f,j_f,val);
                Asy4.AddMatrix(i_f,j_f,val);
                Asy5.AddMatrix(i_f,j_f,val);
                SOLS[i_f-1] += val;
                if (i_f != j_f) sum += val;
                if (C || B)
                {
                    SOLSS[i_f-1] += val;
                    if (i_f != j_f && (!D || E) ) sums += val;
                }
            }
            SOLS[j_f-1] += sum;
            SOLSS[j_f-1] += sums;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, SOLU.data(), n, xtool::xMPIDataType < T >(), MPI_SUM, worldt);
    MPI_Allreduce(MPI_IN_PLACE, SOLUS.data(), n, xtool::xMPIDataType < T >(), MPI_SUM, worldt);
    MPI_Allreduce(MPI_IN_PLACE, SOLS.data(), m, xtool::xMPIDataType < T >(), MPI_SUM, worldt);
    MPI_Allreduce(MPI_IN_PLACE, SOLSS.data(), m, xtool::xMPIDataType < T >(), MPI_SUM, worldt);
    vector < T > SOLUA(SOLU);
    xCPPBlasDef < T >::scal(n, alpha, SOLUA.data());
    vector < T > SOLUSA(SOLUS);
    xCPPBlasDef < T >::scal(n, alpha, SOLUSA.data());
    vector < T > SOLSA(SOLS);
    xCPPBlasDef < T >::scal(m, alpha, SOLSA.data());

    // transposed
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > Aus1T(Aus1,true);
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixUnSym,xTraitMatrixSparceCSR,xTraitMatrixCindex > Aus2T(Aus2,true);
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCOO,xTraitMatrixFindex > &Asy3T = Asy3;
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixFindex > &Asy4T = Asy4;
    xGenericSparseMatrix < T,xTraitMatrixDefinitePositive,xTraitMatrixLowerSym,xTraitMatrixSparceCSC,xTraitMatrixCindex > &Asy5T = Asy5;

    if (verbose)
    {
        ofstream ofs("matrix_"+std::to_string(nb_proct)+"_"+std::to_string(proc_idt)+".txt",ios::out | ios::trunc);
        ofs<<"Aus1"<<endl;
        Aus1.printMatrixMarket(ofs);
        ofs<<"Aus2"<<endl;
        Aus2.printMatrixMarket(ofs);
        ofs<<"Aus1T"<<endl;
        Aus1T.printMatrixMarket(ofs);
        ofs<<"Aus2T"<<endl;
        Aus2T.printMatrixMarket(ofs);
        ofs<<"Asy3"<<endl;
        Asy3.printMatrixMarket(ofs);
        ofs<<"Asy3T"<<endl;
        Asy3T.printMatrixMarket(ofs);
        ofs<<"Asy4"<<endl;
        Asy4.printMatrixMarket(ofs);
        ofs<<"Asy4T"<<endl;
        Asy4T.printMatrixMarket(ofs);
        ofs<<"Asy5"<<endl;
        Asy5.printMatrixMarket(ofs);
        ofs<<"Asy5T"<<endl;
        Asy5T.printMatrixMarket(ofs);
        ofs.close();
    }
    // Gemv product with global vector  ==================================================================================================
    oc<<"Gemv product with global vector"<<endl;
    cout<<"Glob"<<endl;
    Aus1.gemv(m,n,one,zero,X.data(),YU.data(),false,false);
    check1(n,"GEMV UNSYM CSC C notranspose nosparseX alpha=1 beta=0: ",oc,YU,SOLU,proc_idt,worldt,verbose);

    Aus2.gemv(m,n,one,zero,X.data(),YU.data(),false,false);
    check1(n,"GEMV UNSYM CSR C notranspose nosparseX alpha=1 beta=0: ",oc,YU,SOLU,proc_idt,worldt,verbose);

    Asy3.gemv(m,m,one,zero,X.data(),YS.data(),false,false);
    check1(m,"GEMV SYM COO F notranspose nosparseX alpha=1 beta=0: ",oc,YS,SOLS,proc_idt,worldt,verbose);

    Asy4.gemv(m,m,one,zero,X.data(),YS.data(),false,false);
    check1(m,"GEMV SYM CSC F notranspose nosparseX alpha=1 beta=0: ",oc,YS,SOLS,proc_idt,worldt,verbose);

    Asy5.gemv(m,m,one,zero,X.data(),YS.data(),false,false);
    check1(m,"GEMV SYM CSC C notranspose nosparseX alpha=1 beta=0: ",oc,YS,SOLS,proc_idt,worldt,verbose);

    Aus1T.gemv(m,n,one,zero,X.data(),YU.data(),false,true);
    check1(n,"GEMV UNSYM CSC C transpose nosparseX alpha=1 beta=0: ",oc,YU,SOLU,proc_idt,worldt,verbose);

    Aus2T.gemv(m,n,one,zero,X.data(),YU.data(),false,true);
    check1(n,"GEMV UNSYM CSR C transpose nosparseX alpha=1 beta=0: ",oc,YU,SOLU,proc_idt,worldt,verbose);

    Asy3T.gemv(m,m,one,zero,X.data(),YS.data(),false,true);
    check1(m,"GEMV SYM COO F transpose nosparseX alpha=1 beta=0: ",oc,YS,SOLS,proc_idt,worldt,verbose);

    Asy4T.gemv(m,m,one,zero,X.data(),YS.data(),false,true);
    check1(m,"GEMV SYM CSC F transpose nosparseX alpha=1 beta=0: ",oc,YS,SOLS,proc_idt,worldt,verbose);

    Asy5T.gemv(m,m,one,zero,X.data(),YS.data(),false,true);
    check1(m,"GEMV SYM CSC C transpose nosparseX alpha=1 beta=0: ",oc,YS,SOLS,proc_idt,worldt,verbose);




    Aus1.gemv(m,n,alpha,zero,X.data(),YU.data(),false,false);
    check1(n,"GEMV UNSYM CSC C notranspose nosparseX alpha#1 beta=0: ",oc,YU,SOLUA,proc_idt,worldt,verbose);

    Aus2.gemv(m,n,alpha,zero,X.data(),YU.data(),false,false);
    check1(n,"GEMV UNSYM CSR C notranspose nosparseX alpha#1 beta=0: ",oc,YU,SOLUA,proc_idt,worldt,verbose);

    Asy3.gemv(m,m,alpha,zero,X.data(),YS.data(),false,false);
    check1(m,"GEMV SYM COO F notranspose nosparseX alpha#1 beta=0: ",oc,YS,SOLSA,proc_idt,worldt,verbose);

    Asy4.gemv(m,m,alpha,zero,X.data(),YS.data(),false,false);
    check1(m,"GEMV SYM CSC F notranspose nosparseX alpha#1 beta=0: ",oc,YS,SOLSA,proc_idt,worldt,verbose);

    Asy5.gemv(m,m,alpha,zero,X.data(),YS.data(),false,false);
    check1(m,"GEMV SYM CSC C notranspose nosparseX alpha#1 beta=0: ",oc,YS,SOLSA,proc_idt,worldt,verbose);

    Aus1T.gemv(m,n,alpha,zero,X.data(),YU.data(),false,true);
    check1(n,"GEMV UNSYM CSC C transpose nosparseX alpha#1 beta=0: ",oc,YU,SOLUA,proc_idt,worldt,verbose);

    Aus2T.gemv(m,n,alpha,zero,X.data(),YU.data(),false,true);
    check1(n,"GEMV UNSYM CSR C transpose nosparseX alpha#1 beta=0: ",oc,YU,SOLUA,proc_idt,worldt,verbose);


    Aus1.gemv(m,n,alpha,zero,X.data(),YU.data(),false,false);
    Aus1.gemv(m,n,alpha,mone,X.data(),YU.data(),false,false);
    Aus1.gemv(m,n,alpha,one,X.data(),YU.data(),false,false);
    check1(n,"GEMV UNSYM CSC C notranspose nosparseX alpha#1 beta=1: ",oc,YU,SOLUA,proc_idt,worldt,verbose);

    Aus1.gemv(m,n,one,zero,XS.data(),YU.data(),true,false);
    check1(n,"GEMV UNSYM CSC C notranspose sparseX alpha=1 beta=0: ",oc,YU,SOLUS,proc_idt,worldt,verbose);

    Aus1T.gemv(m,n,one,zero,XS.data(),YU.data(),true,true);
    check1(n,"GEMV UNSYM CSC C transpose sparseX alpha=1 beta=0: ",oc,YU,SOLUS,proc_idt,worldt,verbose);

    Aus1.gemv(m,n,alpha,zero,XS.data(),YU.data(),true,false);
    check1(n,"GEMV UNSYM CSC C notranspose sparseX alpha#1 beta=0: ",oc,YU,SOLUSA,proc_idt,worldt,verbose);

    Aus1T.gemv(m,n,alpha,zero,XS.data(),YU.data(),true,true);
    check1(n,"GEMV UNSYM CSC C transpose sparseX alpha#1 beta=0: ",oc,YU,SOLUSA,proc_idt,worldt,verbose);

    Asy3.gemv(m,m,one,zero,XS.data(),YS.data(),true,false);
    check1(m,"GEMV SYM COO F notranspose sparseX alpha=1 beta=0: ",oc,YS,SOLSS,proc_idt,worldt,verbose);

    Asy4.gemv(m,m,one,zero,XS.data(),YS.data(),true,false);
    check1(m,"GEMV SYM CSC F notranspose sparseX alpha=1 beta=0: ",oc,YS,SOLSS,proc_idt,worldt,verbose);

    Asy5.gemv(m,m,one,zero,XS.data(),YS.data(),true,false);
    check1(m,"GEMV SYM CSC C notranspose sparseX alpha=1 beta=0: ",oc,YS,SOLSS,proc_idt,worldt,verbose);

    try
    {
        Aus1.gemv(n,m,one,zero,X.data(),YU.data(),false,false);
    }
    catch (xGenericSparseMatrixException & e)
    {
        if (verbose) cout<<e.what()<<endl;
        oc<<"GEMV test dim: 1"<<endl;
    }
    try
    {
        Aus1.gemv(m,n,one,zero,X.data(),YU.data(),false,true);
    }
    catch (xGenericSparseMatrixException & e)
    {
        if (verbose) cout<<e.what()<<endl;
        oc<<"GEMV test dim transpose: 1"<<endl;
    }
    // Gemv product with dist  vector  ==================================================================================================
    oc<<"Gemv product with distributed vector"<<endl;
    cout<<"Dist"<<endl;
    Aus1.gemv(one,zero,XD,YDU,false,false);
    check2(n,"GEMV UNSYM CSC C notranspose nosparseX alpha=1 beta=0: ",oc, YU, YDU, SOLU,proc_idt,worldt,verbose);

    Aus2.gemv(one,zero,XD,YDU,false,false);
    check2(n,"GEMV UNSYM CSR C notranspose nosparseX alpha=1 beta=0: ",oc, YU, YDU, SOLU,proc_idt,worldt,verbose);

    Asy3.gemv(one,zero,XD,YDS,false,false);
    check2(m,"GEMV SYM COO F notranspose nosparseX alpha=1 beta=0: ",oc,YS,YDS,SOLS,proc_idt,worldt,verbose);

    Asy4.gemv(one,zero,XD,YDS,false,false);
    check2(m,"GEMV SYM CSC F notranspose nosparseX alpha=1 beta=0: ",oc,YS,YDS,SOLS,proc_idt,worldt,verbose);

    Asy5.gemv(one,zero,XD,YDS,false,false);
    check2(m,"GEMV SYM CSC C notranspose nosparseX alpha=1 beta=0: ",oc,YS,YDS,SOLS,proc_idt,worldt,verbose);

    Aus1T.gemv(one,zero,XD,YDU,false,true);
    check2(n,"GEMV UNSYM CSC C transpose nosparseX alpha=1 beta=0: ",oc, YU, YDU, SOLU,proc_idt,worldt,verbose);

    Aus2T.gemv(one,zero,XD,YDU,false,true);
    check2(n,"GEMV UNSYM CSR C transpose nosparseX alpha=1 beta=0: ",oc, YU, YDU, SOLU,proc_idt,worldt,verbose);

    Asy3T.gemv(one,zero,XD,YDS,false,true);
    check2(m,"GEMV SYM COO F transpose nosparseX alpha=1 beta=0: ",oc,YS,YDS,SOLS,proc_idt,worldt,verbose);

    Asy4T.gemv(one,zero,XD,YDS,false,true);
    check2(m,"GEMV SYM CSC F transpose nosparseX alpha=1 beta=0: ",oc,YS,YDS,SOLS,proc_idt,worldt,verbose);

    Asy5T.gemv(one,zero,XD,YDS,false,true);
    check2(m,"GEMV SYM CSC C transpose nosparseX alpha=1 beta=0: ",oc,YS,YDS,SOLS,proc_idt,worldt,verbose);



    Aus1.gemv(alpha,zero,XD,YDU,false,false);
    check2(n,"GEMV UNSYM CSC C notranspose nosparseX alpha#1 beta=0: ",oc, YU, YDU, SOLUA,proc_idt,worldt,verbose);

    Aus2.gemv(alpha,zero,XD,YDU,false,false);
    check2(n,"GEMV UNSYM CSR C notranspose nosparseX alpha#1 beta=0: ",oc, YU, YDU, SOLUA,proc_idt,worldt,verbose);

    Asy3.gemv(alpha,zero,XD,YDS,false,false);
    check2(m,"GEMV SYM COO F notranspose nosparseX alpha#1 beta=0: ",oc,YS,YDS,SOLSA,proc_idt,worldt,verbose);

    Asy4.gemv(alpha,zero,XD,YDS,false,false);
    check2(m,"GEMV SYM CSC F notranspose nosparseX alpha#1 beta=0: ",oc,YS,YDS,SOLSA,proc_idt,worldt,verbose);

    Asy5.gemv(alpha,zero,XD,YDS,false,false);
    check2(m,"GEMV SYM CSC C notranspose nosparseX alpha#1 beta=0: ",oc,YS,YDS,SOLSA,proc_idt,worldt,verbose);

    Aus1T.gemv(alpha,zero,XD,YDU,false,true);
    check2(n,"GEMV UNSYM CSC C transpose nosparseX alpha#1 beta=0: ",oc, YU, YDU, SOLUA,proc_idt,worldt,verbose);

    Aus2T.gemv(alpha,zero,XD,YDU,false,true);
    check2(n,"GEMV UNSYM CSR C transpose nosparseX alpha#1 beta=0: ",oc, YU, YDU, SOLUA,proc_idt,worldt,verbose);


    Aus1.gemv(alpha,zero,XD,YDU,false,false);
    Aus1.gemv(alpha,mone,XD,YDU,false,false);
    Aus1.gemv(alpha,one,XD,YDU,false,false);
    check2(n,"GEMV UNSYM CSC C notranspose nosparseX alpha#1 beta=1: ",oc, YU, YDU, SOLUA,proc_idt,worldt,verbose);

    Aus1.gemv(one,zero,XDS,YDU,true,false);
    check2(n,"GEMV UNSYM CSC C notranspose sparseX alpha=1 beta=0: ",oc, YU, YDU, SOLUS,proc_idt,worldt,verbose);

    Aus1T.gemv(one,zero,XDS,YDU,true,true);
    check2(n,"GEMV UNSYM CSC C transpose sparseX alpha=1 beta=0: ",oc, YU, YDU, SOLUS,proc_idt,worldt,verbose);

    Aus1.gemv(alpha,zero,XDS,YDU,true,false);
    check2(n,"GEMV UNSYM CSC C notranspose sparseX alpha#1 beta=0: ",oc, YU, YDU, SOLUSA,proc_idt,worldt,verbose);

    Aus1T.gemv(alpha,zero,XDS,YDU,true,true);
    check2(n,"GEMV UNSYM CSC C transpose sparseX alpha#1 beta=0: ",oc, YU, YDU, SOLUSA,proc_idt,worldt,verbose);

    Asy3.gemv(one,zero,XDS,YDS,true,false);
    check2(m,"GEMV SYM COO F notranspose sparseX alpha=1 beta=0: ",oc,YS,YDS,SOLSS,proc_idt,worldt,verbose);

    Asy4.gemv(one,zero,XDS,YDS,true,false);
    check2(m,"GEMV SYM CSC F notranspose sparseX alpha=1 beta=0: ",oc,YS,YDS,SOLSS,proc_idt,worldt,verbose);

    Asy5.gemv(one,zero,XDS,YDS,true,false);
    check2(m,"GEMV SYM CSC C notranspose sparseX alpha=1 beta=0: ",oc,YS,YDS,SOLSS,proc_idt,worldt,verbose);

    try
    {
        Aus1.gemv(alpha,zero,YDU,YDS,false,false);
    }
    catch (xGenericSparseMatrixException & e)
    {
        if (verbose) cout<<e.what()<<endl;
        oc<<"GEMV test dim: 1"<<endl;
    }
    try
    {
        Aus1.gemv(alpha,zero,YDS,YDU,false,true);
    }
    catch (xGenericSparseMatrixException & e)
    {
        if (verbose) cout<<e.what()<<endl;
        oc<<"GEMV test dim transpose: 1"<<endl;
    }
    return;
}
int main(int argc, char *argv[])
{
    // initialize mpi universe
    xtool::xMPIEnv::init(argc,argv);
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm world_g;
    MPI_Group general_group,selected_group;

    // general group
    MPI_Comm_group(world, &general_group );

    // init around id
    int proc_id,nb_proc;
    MPI_Comm_rank(world, &proc_id);
    MPI_Comm_size(world, &nb_proc);

    const bool verbose = true;

    //Redirection of all output to individual files
    string fo = "proc_"+std::to_string(proc_id)+"_output.txt";
    FILE * ok = freopen (fo.c_str(),"w",stdout);
    if (!ok) {std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << " compiled "<<__DATE__<<" at "<<__TIME__<< std::endl; throw; }


    // for group selection
    vector < int > select(nb_proc);
    for (int i = 0; i < nb_proc; ++i) select[i] = i;
    int nb[3] = {1,max(nb_proc/2,1),nb_proc};                 // work for one proc, but start to
                                                              // be 3 independent value
                                                              // with 4 proc


    // loop on increasing group size
    for (int pass = 0; pass < 3; ++pass)
    {

        // group selection and associated communicator
        MPI_Group_incl( general_group, nb[pass], &select[0], &selected_group );
        MPI_Comm_create( world, selected_group, &world_g );

        // Do the job
        if ( world_g != MPI_COMM_NULL)
        {
            int proc_id_g,nb_proc_g;
            MPI_Comm_rank(world_g, &proc_id_g);
            MPI_Comm_size(world_g, &nb_proc_g);

            // create index to be used with vector instances
            xlinalg::xDistIndex idxc(world_g);
            generateIdx(idxc,SIZECL,proc_id_g,nb_proc_g);
            xlinalg::xDistIndex idxr(world_g);
            generateIdx(idxr,SIZERL,proc_id_g,nb_proc_g);
            if (verbose)
            {
                cout <<"getGlobalIndexSize c : "<<idxc.getGlobalIndexSize()<<endl;
                cout <<"getPackedIndexSize c: "<<idxc.getPackedIndexSize()<<endl;
                cout <<"getPackedLocalIndexSize c: "<<idxc.getPackedLocalIndexSize()<<endl;
                cout <<"All index c: ";
                for (auto i : idxc) cout <<" "<<i;
                cout <<endl;
                cout <<"getGlobalIndexSize r : "<<idxr.getGlobalIndexSize()<<endl;
                cout <<"getPackedIndexSize r: "<<idxr.getPackedIndexSize()<<endl;
                cout <<"getPackedLocalIndexSize r: "<<idxr.getPackedLocalIndexSize()<<endl;
                cout <<"All index r: ";
                for (auto i : idxr) cout <<" "<<i;
                cout <<endl;
            }
            // check file
            ofstream oc;
            if (!proc_id_g)
            {
                oc.open("check_"+std::to_string(nb_proc_g)+".txt",ios::out | ios::trunc);
                oc<<"With "<<nb_proc_g<<" process"<<endl;
            }

            // do test for double
            if (!proc_id_g) oc<<"double"<<endl;
            test < double >(idxc,idxr,world_g,nb_proc_g,proc_id_g,oc,verbose,2.,25.,3./4.,1./4.);

            // do test for cplx
            if (!proc_id_g) oc<<"complex"<<endl;
            test < complex < double > >(idxc,idxr,world_g,nb_proc_g,proc_id_g,oc,verbose,complex < double >(2.,0.),complex < double >(6.,7.),complex < double >(3./4.,1./4.),complex < double >(1./4.,-1./4.));

            if (!proc_id_g) oc.close();
        }

        // syncro on all process: avoid call to group free by unused process
        MPI_Barrier(world);

        // cleaning xMPITag to avoid that if next created communicator in the loop
        // got the same address the xMPITag starts with non null counts
        if ( world_g != MPI_COMM_NULL)
        {
            xMPITag::resetTag(world_g);
            MPI_Comm_free( &world_g );
        }


        MPI_Group_free( &selected_group );
    }

    return xtool::xMPIEnv::finalize();
}
