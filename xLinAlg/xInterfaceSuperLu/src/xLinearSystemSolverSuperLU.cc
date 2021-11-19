/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include <fstream>

#include "xLinearSystemSolverSuperLUDataType.h"

//Set SuperLU Version. Remember to change also in xLinearSystemSolverSuperLU.h
// 3 is the only supported version
// 4 is experimental (some features do not work, like ConditionNumber)
// 5 is experimental for the moment

#include "SuperLU_Version.h"



// include for superlu 3.0 
#ifdef SSUPERLU
#include "ssp_defs.h"
#endif

#ifdef DSUPERLU
#if SUPERLU_VERSION == 3
    #include "dsp_defs.h"
#else
    #include "slu_ddefs.h"
#endif
#endif

#ifdef CSUPERLU
namespace superLuCmplx
{
#include "csp_defs.h"
} // end namespace
using namespace superLuCmplx;
#endif

#ifdef ZSUPERLU
#include "zsp_defs.h"
#endif

// include for superlu 4.0
// to be done : #include "slu_ddefs.h"
// all interface have to be revisited

#define _SOLVERSUPERLU_HEADER_SET 1

#include "xLinearSystemSolverSuperLU.h"


using std::cout;
using std::endl;

namespace xlinalg 
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverSuperLUBase class implementation ////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xLinearSystemSolverSuperLUBase::xLinearSystemSolverSuperLUBase() :
  verbose(0),keep(false),nrows(0),ncol(0),nzdata(0),A(nullptr),L(nullptr),U(nullptr),RHS(nullptr),SOL(nullptr),stat(nullptr),options(nullptr),Glu(nullptr), mem_usage(nullptr),
    perm_r(nullptr),perm_c(nullptr),ext_perm_c(nullptr),previous_ext_perm_c(false),diag_pivot_thresh(XLINEARSOLVERSUPERLU_EMPTY_OPTION*1.)
{
    // init SuperLU structure container
    try
    {
        A           = (void *)( new SuperMatrix );
        L           = (void *)( new SuperMatrix );
        U           = (void *)( new SuperMatrix );
        RHS         = (void *)( new SuperMatrix );
        SOL         = (void *)( new SuperMatrix );
        stat        = (void *)( new SuperLUStat_t );
        options     = (void *)( new superlu_options_t );
#if SUPERLU_VERSION == 5
	Glu         = (void *)(new GlobalLU_t);
#endif
        mem_usage   = (void *)( new mem_usage_t);

	
    }
    catch (std::bad_alloc & e)
    {
        cout << "Standard exception: " << e.what() << endl;
        std::ostringstream oss;
        oss << " SuperLU data structure where not allocated normaly ?!\n";
        throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }

    // set to null store part of SuperMatrix to avoid inconsistant free in cleanStore
    ((SuperMatrix *)A)->Store=(void *)nullptr;
    ((SuperMatrix *)L)->Store=(void *)nullptr;
    ((SuperMatrix *)U)->Store=(void *)nullptr;
    ((SuperMatrix *)RHS)->Store=(void *)nullptr;
    ((SuperMatrix *)SOL)->Store=(void *)nullptr;
    
    // set default option
    setDefaultParameter();

    // init saved default options
    for(int i=0;i<SIZE_ENUM; ++i) saved_option[i]=XLINEARSOLVERSUPERLU_EMPTY_OPTION;

}
xLinearSystemSolverSuperLUBase::~xLinearSystemSolverSuperLUBase()
{
    // delete SuperLU structure container
    try
    {
        if(A) 
        {
            cleanStore (A);
            delete (SuperMatrix *)( A );
        }
        if(L)           delete (SuperMatrix *)( L );
        if(U)           delete (SuperMatrix *)( U );
        if(RHS) 
        {
            cleanStore (RHS);
            delete (SuperMatrix *)( RHS );
        }
        if(SOL) 
        {
            cleanStore (SOL);
            delete (SuperMatrix *)( SOL );
        }
        if(stat)        delete (SuperLUStat_t *)( stat );
        if(options)     delete (superlu_options_t *)( options );
	#if SUPERLU_VERSION == 5
	if(Glu){
	  delete( GlobalLU_t *)( Glu);
        }
	#endif
        if(mem_usage)   delete (mem_usage_t *)( mem_usage );
    }
    catch (std::bad_alloc & e)
    {
        cout << "Standard exception: " << e.what() << endl;
        std::cerr << xLinearSystemSolverSuperLUException::createMessage(" SuperLU matrix structure where not deallocated normaly !\n",__FILE__,__LINE__,__DATE__,__TIME__) << std::endl;
        std::terminate();
    }
    if (verbose)
        cout << "Deleting SuperLU"<<endl;
    return;
} 

#define OPT(a,b,c) saved_option[a]=c;((superlu_options_t *)(options))->a=(b)(c);
template < >
void xLinearSystemSolverSuperLUBase :: setParam<int>( int param_tab, int param_index, int param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case INTER :
         {
             verbose = param_value;
             break;
         }
        case KEEP :
         {
             keep = (param_value=1);
             break;
         }
        case Fact:
         {
             OPT(Fact,fact_t,param_value)
             break;
         }
        case Equil:
         {
             OPT(Equil,yes_no_t,param_value)
             break;
         }
        case ColPerm:
         {
             OPT(ColPerm,colperm_t,param_value)
             break;
         }
        case Trans:
         {
             OPT(Trans,trans_t,param_value)
             break;
         }
        case IterRefine:
         {
             OPT(IterRefine,IterRefine_t,param_value)
             break;
         }
        case PrintStat:
         {
             OPT(PrintStat,yes_no_t,param_value)
             break;
         }
        case SymmetricMode:
         {
             OPT(SymmetricMode,yes_no_t,param_value)
             break;
         }
        case PivotGrowth:
         {
             OPT(PivotGrowth,yes_no_t,param_value)
             break;
         }
        case ConditionNumber:
         {
             OPT(ConditionNumber,yes_no_t,param_value)
             break;
         }
        case RowPerm:
         {
             OPT(RowPerm,rowperm_t,param_value)
             break;
         }
        case ReplaceTinyPivot:
         {
             OPT(ReplaceTinyPivot,yes_no_t,param_value)
             break;
         }
        case SolveInitialized:
         {
             OPT(SolveInitialized,yes_no_t,param_value)
             break;
         }
        case RefineInitialized:
         {
             OPT(RefineInitialized,yes_no_t,param_value)
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
#define OPTT(a,b) saved_option[a]=(int)b;((superlu_options_t *)(options))->a=(b);
template < >
void xLinearSystemSolverSuperLUBase :: setParam<yes_no_t>( int param_tab, int param_index, yes_no_t param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case Equil:
         {
             OPTT(Equil,param_value)
             break;
         }
        case PrintStat:
         {
             OPTT(PrintStat,param_value)
             break;
         }
        case SymmetricMode:
         {
             OPTT(SymmetricMode,param_value)
             break;
         }
        case PivotGrowth:
         {
             OPTT(PivotGrowth,param_value)
             break;
         }
        case ConditionNumber:
         {
             OPTT(ConditionNumber,param_value)
             break;
         }
        case ReplaceTinyPivot:
         {
             OPTT(ReplaceTinyPivot,param_value)
             break;
         }
        case SolveInitialized:
         {
             OPTT(SolveInitialized,param_value)
             break;
         }
        case RefineInitialized:
         {
             OPTT(RefineInitialized,param_value)
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
template < >
void xLinearSystemSolverSuperLUBase :: setParam<fact_t>( int param_tab, int param_index, fact_t param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case Fact:
         {
             OPTT(Fact,param_value)
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
template < >
void xLinearSystemSolverSuperLUBase :: setParam<colperm_t>( int param_tab, int param_index, colperm_t param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case ColPerm:
         {
             OPTT(ColPerm,param_value)
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
template < >
void xLinearSystemSolverSuperLUBase :: setParam<rowperm_t>( int param_tab, int param_index, rowperm_t param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case RowPerm:
         {
             OPTT(RowPerm,param_value)
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
template < >
void xLinearSystemSolverSuperLUBase :: setParam<trans_t>( int param_tab, int param_index, trans_t param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case Trans:
         {
             OPTT(Trans,param_value)
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
template < >
void xLinearSystemSolverSuperLUBase :: setParam<IterRefine_t>( int param_tab, int param_index, IterRefine_t param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case IterRefine:
         {
             OPTT(IterRefine,param_value)
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
template < >
void xLinearSystemSolverSuperLUBase :: setParam<double>( int param_tab, int param_index, double param_value)
{

    // depending on param_tab
    switch (param_tab)
    {
        case DiagPivotThresh:
         {
             diag_pivot_thresh = ((superlu_options_t *)(options))->DiagPivotThresh = param_value;
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of setParam use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_value is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}

template < >
void xLinearSystemSolverSuperLUBase :: getInternalData<int>( int param_tab,int param_index,int &param_ref)
{

    // depending on param_tab
    switch (param_tab)
    {
        case INTER :
         {
             param_ref = verbose;
             break;
         }
        case KEEP :
         {
             param_ref = keep?1:0;
             break;
         }
        case RefineSteps :
         {
             param_ref = ((SuperLUStat_t *)stat)->RefineSteps;
             break;
         }
        case Fact:
        case Equil:
        case ColPerm:
        case Trans:
        case IterRefine:
        case PrintStat:
        case SymmetricMode:
        case PivotGrowth:
        case ConditionNumber:
        case RowPerm:
        case ReplaceTinyPivot:
        case SolveInitialized:
        case RefineInitialized:
         {
             param_ref = saved_option[param_tab];
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Invocation of getInternal data use wrong param_tab ("<<param_tab<<") ; see xLinearSolverSuperLU.h enum_param_tab to use a good one\n Remeber that param_ref is typed\n";
             throw xLinearSystemSolverSuperLUException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }

}
void xLinearSystemSolverSuperLUBase::setDefaultParameter()
{
    set_default_options((superlu_options_t *)options);
}
void xLinearSystemSolverSuperLUBase::reInitParameter()
{
    for(int i=0;i<SIZE_ENUM; ++i) 
        if( saved_option[i]!=XLINEARSOLVERSUPERLU_EMPTY_OPTION)
        {
            setParam<int>(i,0,saved_option[i]);
        }
    if (diag_pivot_thresh>XLINEARSOLVERSUPERLU_EMPTY_OPTION*1.)
            setParam<double>(DiagPivotThresh,0,diag_pivot_thresh);
           
}

void xLinearSystemSolverSuperLUBase::initStat() 
{
    StatInit((SuperLUStat_t *)stat);
}
void xLinearSystemSolverSuperLUBase::clearStat()
{
    StatFree((SuperLUStat_t *)stat);
}
void xLinearSystemSolverSuperLUBase::printStat()
{
   if ( ((superlu_options_t *)(options))->PrintStat ) StatPrint((SuperLUStat_t *)stat);
}
void  xLinearSystemSolverSuperLUBase :: printPerf ()
{
    SCformat *Lstore;
    NCformat *Ustore;
    double   *utime;
    float    mega=1048576.;
    flops_t  *ops;
    SuperLUStat_t *statistic=(SuperLUStat_t *)stat;
    mem_usage_t *memory_usage=(mem_usage_t *)mem_usage;
    
    utime = statistic->utime;
    ops   = statistic->ops;
    
    if ( utime[FACT] != 0. )
        printf("\tFactor flops = %e\tMflops = %8.2f\n", ops[FACT],
               ops[FACT]*1e-6/utime[FACT]);
    printf("\tIdentify relaxed snodes     = %8.2f\n", utime[RELAX]);
    if ( utime[SOLVE] != 0. )
        printf("\tSolve flops = %.0f, Mflops = %8.2f\n", ops[SOLVE],
               ops[SOLVE]*1e-6/utime[SOLVE]);
    
    Lstore = (SCformat *) ((SuperMatrix *)L)->Store;
    Ustore = (NCformat *) ((SuperMatrix *)U)->Store;
    printf("\tNo of nonzeros in factor L = %d\n", Lstore->nnz);
    printf("\tNo of nonzeros in factor U = %d\n", Ustore->nnz);
    printf("\tNo of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);

#if SUPERLU_VERSION == 3
    printf("LU MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
           (memory_usage->for_lu)/mega, (memory_usage->total_needed)/mega,
           memory_usage->expansions);
#else
    printf("LU MB %.3f\ttotal MB needed %.3f\n",
           (memory_usage->for_lu)/mega, (memory_usage->total_needed)/mega);
#endif
        
    printf("\tFactor\tMflops\tSolve\tMflops\tEtree\tEquil\tRcond\tRefine\n");
    printf("PERF:%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
           utime[FACT], ops[FACT]*1e-6/utime[FACT],
           utime[SOLVE], ops[SOLVE]*1e-6/utime[SOLVE],
           utime[ETREE], utime[EQUIL], utime[RCOND], utime[REFINE]);

}
void  xLinearSystemSolverSuperLUBase :: cleanFactor ()
{
    if ( ((SuperMatrix *)L)->Store ) Destroy_SuperNode_Matrix((SuperMatrix *)L);
    if ( ((SuperMatrix *)U)->Store ) Destroy_CompCol_Matrix((SuperMatrix *)U);
}
void  xLinearSystemSolverSuperLUBase :: cleanStore (void * M)
{
    if (M)
    {
        SuperMatrix * pM=(SuperMatrix *)( M );
        if (pM->Store)
            Destroy_SuperMatrix_Store(pM);
    }
    return;
}
//  type dependant implementation
#ifdef SSUPERLU
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<float,xTraitMatrixSparceCSR>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, float *data)
{
    sCreate_CompRow_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NR, SLU_S, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<float,xTraitMatrixSparceCSC>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, float *data)
{
    sCreate_CompCol_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NC, SLU_S, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixDense<float>(void * D, const int n,const int m,const int ld, float *data)
{
    // little cleaning before allocating again 
    if (D) cleanStore(D);
    sCreate_Dense_Matrix((SuperMatrix *)D, n, m, data, ld, SLU_DN, SLU_S, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase :: gssvx < float >( void *options, void *A, int *perm_c,int *perm_r,int *etree,
            char *equed,void *R,void *C,void *L,void *U, void *work,int lwork,
            void *RHS,void *SOL, void *rpg,void *rcond,void *ferr,void *berr,
            void *mem_usage, void *stat, int *info)
{ 
    sgssvx((superlu_options_t *)options, (SuperMatrix *)A, perm_c,perm_r,etree,
            equed,(float *)R,(float *)C,(SuperMatrix *)L,(SuperMatrix *)U,
            work,lwork,(SuperMatrix *)RHS,(SuperMatrix *)SOL,
            (float *)rpg,(float *)rcond,(float *)ferr,(float *)berr,
            (mem_usage_t *)mem_usage, (SuperLUStat_t *)stat, info);
}
#endif
#ifdef DSUPERLU
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<double,xTraitMatrixSparceCSR>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, double *data)
{
    dCreate_CompRow_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NR, SLU_D, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<double,xTraitMatrixSparceCSC>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, double *data)
{
    dCreate_CompCol_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NC, SLU_D, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixDense<double>(void * D, const int n,const int m,const int ld, double *data)
{
    // little cleaning before allocating again 
    if (D) cleanStore(D);
    dCreate_Dense_Matrix((SuperMatrix *)D, n, m, data, ld, SLU_DN, SLU_D, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase :: gssvx < double >( void *options, void *A, int *perm_c,int *perm_r,int *etree,
            char *equed,void *R,void *C,void *L,void *U, void *work,int lwork,
            void *RHS,void *SOL, void *rpg,void *rcond,void *ferr,void *berr,
#if SUPERLU_VERSION == 5
	    void *Glu,						 
#endif 							 
            void *mem_usage, void *stat, int *info)
{ 
    dgssvx((superlu_options_t *)options, (SuperMatrix *)A, perm_c,perm_r,etree,
            equed,(double *)R,(double *)C,(SuperMatrix *)L,(SuperMatrix *)U,
            work,lwork,(SuperMatrix *)RHS,(SuperMatrix *)SOL,
            (double *)rpg,(double *)rcond,(double *)ferr,(double *)berr,
#if SUPERLU_VERSION == 5
	   (GlobalLU_t *)Glu,
#endif      
            (mem_usage_t *)mem_usage, (SuperLUStat_t *)stat, info);
}
#endif
#ifdef CSUPERLU
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<superLuCmplx::complex,xTraitMatrixSparceCSR>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, superLuCmplx::complex *data)
{
    cCreate_CompRow_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NR, SLU_C, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<superLuCmplx::complex,xTraitMatrixSparceCSC>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, superLuCmplx::complex *data)
{
    cCreate_CompCol_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NC, SLU_C, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixDense<superLuCmplx::complex>(void * D, const int n,const int m,const int ld, superLuCmplx::complex *data)
{
    // little cleaning before allocating again 
    if (D) cleanStore(D);
    cCreate_Dense_Matrix((SuperMatrix *)D, n, m, data, ld, SLU_DN, SLU_C, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase :: gssvx < superLuCmplx::complex >( void *options, void *A, int *perm_c,int *perm_r,int *etree,
            char *equed,void *R,void *C,void *L,void *U, void *work,int lwork,
            void *RHS,void *SOL, void *rpg,void *rcond,void *ferr,void *berr,
            void *mem_usage, void *stat, int *info)
{ 
    cgssvx((superlu_options_t *)options, (SuperMatrix *)A, perm_c,perm_r,etree,
            equed,(float *)R,(float *)C,(SuperMatrix *)L,(SuperMatrix *)U,
            work,lwork,(SuperMatrix *)RHS,(SuperMatrix *)SOL,
            (float *)rpg,(float *)rcond,(float *)ferr,(float *)berr,
            (mem_usage_t *)mem_usage, (SuperLUStat_t *)stat, info);
}
#endif
#ifdef ZSUPERLU
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<doublecomplex,xTraitMatrixSparceCSR>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, doublecomplex *data)
{
    zCreate_CompRow_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NR, SLU_Z, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<doublecomplex,xTraitMatrixSparceCSC>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, doublecomplex *data)
{
    zCreate_CompCol_Matrix((SuperMatrix *)M,n,m,nnz,data,index2,index1, SLU_NC, SLU_Z, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixDense<doublecomplex>(void * D, const int n,const int m,const int ld, doublecomplex *data)
{
    // little cleaning before allocating again 
    if (D) cleanStore(D);
    zCreate_Dense_Matrix((SuperMatrix *)D, n, m, data, ld, SLU_DN, SLU_Z, SLU_GE);
}
template < >
void xLinearSystemSolverSuperLUBase :: gssvx < doublecomplex > ( void *options, void *A, int *perm_c,int *perm_r,int *etree,
            char *equed,void *R,void *C,void *L,void *U, void *work,int lwork,
            void *RHS,void *SOL, void *rpg,void *rcond,void *ferr,void *berr,
            void *mem_usage, void *stat, int *info)
{ 
    zgssvx((superlu_options_t *)options, (SuperMatrix *)A, perm_c,perm_r,etree,
            equed,(double *)R,(double *)C,(SuperMatrix *)L,(SuperMatrix *)U,
            work,lwork,(SuperMatrix *)RHS,(SuperMatrix *)SOL,
            (double *)rpg,(double *)rcond,(double *)ferr,(double *)berr,
            (mem_usage_t *)mem_usage, (SuperLUStat_t *)stat, info);
}
#endif
/* for future if ever
 void xLinearSystemSolverSuperLUBase::setPointerSuperMatrixSparse<double,xTraitMatrixSparceDCSR>(void *M,const int n,const int m,const int nnz, int *index1, int *index2, double *data,.....)
 dCreate_CompRowLoc_Matrix_dist .......
 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverSuperLU class implementation ////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverSuperLUException class implementation ///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
std::string xLinearSystemSolverSuperLUException::createMessage(std::string info,std::string file,int Line,std::string date,std::string time){
    std::ostringstream oss;
    oss << "In file "<< file << " line " << Line << " compiled "<<date<<" at "<<time<<std::endl;
    oss << "xLinearSystemSolverSuperLUException : "<< info << std::endl;
    return oss.str();
}

xLinearSystemSolverSuperLUException::xLinearSystemSolverSuperLUException(std::string info,std::string file,int Line,std::string date,std::string time):msg(xLinearSystemSolverSuperLUException::createMessage(info, file, Line, date, time))
{}
/// general exception object : destructor
//xLinearSystemSolverSuperLUException :: ~xLinearSystemSolverSuperLUException() throw( ) = default;

/// mandatory what method
const char * xLinearSystemSolverSuperLUException::what() const throw( )
{
    return this->msg.c_str();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xLinearSystemSolverSuperLUException class implementation ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



} // end namespace





