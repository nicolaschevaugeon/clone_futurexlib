/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _SOLVERSUPERLU_H
#define _SOLVERSUPERLU_H

#include <iostream>
#include "SuperLU_Version.h"
//Set SuperLU Version. Remember to change also in xLinearSystemSolverSuperLU.cc
// 3 is the only supported version
// 4 is experimental (some features do not work, like ConditionNumber)
// 5 is experimental for the moment



//
#define XLINEARSOLVERSUPERLU_EMPTY_OPTION     -999999
//
/* code for status */
#define XLINEARSOLVERSUPERLU_INIT       1U
#define XLINEARSOLVERSUPERLU_ALLOC      2U
#define XLINEARSOLVERSUPERLU_SYMB       4U
#define XLINEARSOLVERSUPERLU_FACT       8U
#define XLINEARSOLVERSUPERLU_SOLV      16U

/// xLinearSystemSolverSuperLU is the interface to superLU.
//! It can be use as is, or with expert paramerter set.
//! <br/> by default following option setting is used :
//! <ul><li>    Fact = follow methode invocation 
//! </li><li>   Equil = YES;
//! </li><li>   ColPerm = COLAMD;
//! </li><li>   DiagPivotThresh = 1.0;
//! </li><li>   Trans = NOTRANS;
//! </li><li>   IterRefine = NOREFINE;
//! </li><li>   SymmetricMode = NO;
//! </li><li>   PivotGrowth = NO;
//! </li><li>   ConditionNumber = NO;
//! </li><li>   PrintStat = NO;
//! </li></ul>
//! <br/> This option members are accessed via setParam using enum enum_param_tab of xLinearSystemSolverSuperLUBase class with same names
//! <br/> All other API keyword in the documentation are available here as typedef are duplicate below.
//!       For example if you want to have reciprocal condition number you use setParam before computation
//!       like this :
//!             my_instance_of_xLinearSystemSolverSuperLU.setParam(xlinalg::xLinearSystemSolverSuperLUBase::ConditionNumber,0,YES);
//!       and after computation :
//!             my_sub_data_type rcond
//!             my_instance_of_xLinearSystemSolverSuperLU.getData(xlinalg::xLinearSystemSolverSuperLUBase::GET_RCOND,0,rcond);
//!
//! <br/> See documentation (and below) for complete list and usage of API keyord. Here are some of them :  
//! <ul><li>    yes_no_t : YES, NO
//! </li><li>   fact_t : DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED 
//! </li><li>   IterRefine_t : NOREFINE, SINGLE, DOUBLE, EXTRA
//! </li></ul>
//!
//!
//
// nota :
//  for now only one class is derived from xLinearSystemSolverSuperLUBase but
//  derived class may be created for // (SuperLU_MT and SuperLU_DIST) 
//  if done maybe more things might be put or supress from base class. For now after a short
//  analysys of SuperLu doc it look like data structure base (SuperMatrix) are comon 
//  betwen SuperLU SuperLU_MT SuperLU_DIST => data structure keeped in base class
//  There contain is filed in derived class methode as they looks different from doc reading ...
//  
// this include is here and not in imp file as xLinearSystemSolverSuperLUSubType is use here
// (see data_sub_type below)
#include "xLinearSystemSolverSuperLUTraitPolicy.h"

#ifndef  _SOLVERSUPERLU_HEADER_SET 
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=                                                                                                                                         ==
//= This part is a partial copy from util.h and [sdcz]sp_defs.h from SuperLu version 3.0                                                    ==
//= It's realy dangerous as if you link this interface with a other version bug will be very difficult to track                             ==
//= So please pay attention to what you are linking this interface with !!!!                                                                ==
//= In one hand it's dangerous in the other hand it give user the full usage of SuperLU with documentation Keywords transparently           ==
//= Same remark for complex type define here in xLinearSystemSolverSuperLUcplx                                                              ==
//=                                                                                                                                         ==
//=  Be careful also to interference with other library. For now those enum look to not interfer with other thing but as they               ==
//=  are in the public namspace they might interfear... to put in lagl namespace ?? or in a proper namespace ??? to thing about             ==
//=                                                                                                                                         ==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//
typedef enum {NO, YES}                                          yes_no_t;
typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;
typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD, MY_PERMC}colperm_t;
typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
typedef enum {NOREFINE, SINGLE=1, DOUBLE, EXTRA}                IterRefine_t;
typedef enum {LUSUP, UCOL, LSUB, USUB}                          MemType;
typedef enum {HEAD, TAIL}                                       stack_end_t;
typedef enum {SYSTEM, USER}                                     LU_space_t;


#include  "xLinearSystemSolverSuperLUcplx.h"


//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
//=                                                                                                                                         ==
//=   End dangerous part                                                                                                                    ==
//=                                                                                                                                         ==
//=WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING===WARNING==
#endif


namespace xlinalg
{
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverSuperLUException class //////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// interface derived class of standart exception
class xLinearSystemSolverSuperLUException : public std::exception
{
    public:
        xLinearSystemSolverSuperLUException(std::string,std::string,int,std::string,std::string);
//        ~xLinearSystemSolverSuperLUException() throw( ) override;
        const char * what() const throw( ) override;
    static std::string createMessage(std::string,std::string,int,std::string,std::string);

    private:
        std::string msg;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xLinearSystemSolverSuperLUException  class /////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
class xLinearSystemSolverSuperLUBase  
{
    public :
        ///Traits available for all derived class
        typedef xTraitMatrixCindex matrix_indexing;

        /// parameter access
        /// partially using naming convention from superlu_options_t in util.h SuperLU version 3.0
        //! see documentation to use
        enum enum_param_tab { INTER , KEEP,
          // option struture access
          Fact , Equil, ColPerm, Trans, IterRefine,
          PrintStat, SymmetricMode, DiagPivotThresh, PivotGrowth,
          ConditionNumber, RowPerm, ReplaceTinyPivot, SolveInitialized, RefineInitialized,
          // option to get internal data
          RefineSteps,
          // option to get data of interface
          GET_VECTC,GET_VECTR,GET_FERR,GET_BERR,GET_RCOND,GET_RPG,
          // to allocate table of correct size
          SIZE_ENUM };

        // constructor
        xLinearSystemSolverSuperLUBase ( );
        ~xLinearSystemSolverSuperLUBase ( );

        /// methode to controle SuperLU
        template< typename P>
        void setParam( int,int,P);

        /// methode to retreave information from SuperLu internal data
        // nota :
        //    for now not much is done. Must be completed to all stat and memusage data at least
        template< typename P>
        void getInternalData( int,int,P&);


    protected:

        /// setting superMatrix for sparse matrix 
        template < typename T, typename STORAGE_TYPE >
        void setPointerSuperMatrixSparse(void *M, const int n,const int m,const int nnz,int *index1, int *index2, T *data);
        
        /// setting superMatrix for dense matrix 
        template < typename T >
        void setPointerSuperMatrixDense(void *D, const int n,const int m,const int ld, T *data);

        /// initalise statistic container
        void initStat();
        /// print statistic container
        void printStat();
        /// clear statistic container
        void clearStat();

        /// clear Factor container (L and U)
        void cleanFactor();
        /// clear Store part of supermatrix container (A,...)
        void  cleanStore (void * );

        /// generic expert driver access
        template< typename T>
        void gssvx(void *options, void *A, int *perm_c,int *perm_r,int *etree,
            char *equed,void *R,void *C,void *L,void *U, void *work,int lwork,
            void *RHS,void *SOL, void *rpg,void *rcond,void *ferr,void *berr,
#if SUPERLU_VERSION == 5
	    void *Glu,
#endif		  
            void *mem_usage, void *stat, int *info);

        ///  performance printing function 
        void  printPerf();

        /// set default option from package
        void setDefaultParameter();
        /// re init options according to setParam setting
        void reInitParameter();

        /// verbosity of the interface itself
        int verbose;

        /// parameter which control  use of intermediate container to avoid
        /// any change of attached data during factorisation
        //! Warning this option imply copy of data values of the matrix into keep_xx container in derived class
        //! => more memory is used
        //! This option by default is false.
        //! It may be change with setParam with KEEP argument.
        //! if changed it have to be modified only befor attaching a matrix. It should not be modified after except
        //! if you attach a new container right after the modification
        bool keep;

        /// current number of rows an columns of the last matrix attached to the solver and number of non zero 
        int nrows,ncol,nzdata;

        /// container in SuperLU format
        void *A;        //real type is SuperMatrix
        void *L;        //real type is SuperMatrix
        void *U;        //real type is SuperMatrix
        void *RHS;      //real type is SuperMatrix
        void *SOL;      //real type is SuperMatrix
        void *stat;     //real type is SuperLUStat_t
        void *options;  //real type is superlu_options_t
	void *Glu;      //real type is GlobalLU_t
        void *mem_usage;//real type is mem_usage_t
	

        /// permutation vectors
        int *perm_r,*perm_c;
        int *ext_perm_c;
        bool previous_ext_perm_c;


        /// saved option vectors
        double diag_pivot_thresh;
        int saved_option[SIZE_ENUM];


};

template < typename T = double >
class xLinearSystemSolverSuperLU : public   xLinearSystemSolverSuperLUBase
{
    public :
        // Traits
        typedef T data_type;
        typedef typename xlinalg::xLinearSystemSolverSuperLUSubType<T>::sub_type data_sub_type;
        typedef xLinearSystemSolverSuperLUBase::matrix_indexing matrix_indexing;
        // for historical reason this interface is using  a CSR storage
        // SuperLu accept it but in reallity the solver transforme artificial the problem in it's transposed form
        // to deal with CSC storage in all case.
        // If you want to use CSC storage change coment below. You will need to have a compatible matrix container to attach to.
        // At time of this commit, for this initial version,  only xGenericSparceMatrix have this capability implemented
        typedef xTraitMatrixSparceCSR matrix_storage;
        //typedef xTraitMatrixSparceCSC matrix_storage;

        //constructor/destructor
        xLinearSystemSolverSuperLU () ;
        ~xLinearSystemSolverSuperLU () ;

        /// Connecting a Matrix to the solver
        template < typename M > void connectMatrix( M &matrix);
        /// Connecting a Matrix to the solver and give column permutation
        template < typename M > void connectMatrix( M &matrix , int *external_col_perm);

        /// do symbolique phase assuming that a previous connection exist with a matrix. connectMatrix  call this methode by default.
        //! This last point make this methode allmost anaivalable for now.  Lived in the public part but may become private
        //! with superLu no separation betwen sym and fact was implemented (to much work for now). symb is then empty in this implementation
        void symb();

        /// do numerical factorisation assuming that a previous symb methode was called
        void fact();
 
        /// do solving the linear system assuming that a previous symb methode was called at least
        //! if a previous call to fact was done the  numerical factorisation is not re-computed
        template < typename V  >
        void solve(const V &,  V & );

        /// methode to retreave information from SuperLu interface
        template< typename P>
        void getData( int,int,P&);

    private :

        
        // initialisation methodes
        template < typename M > void init ();

        // allocatting and dealocating SuperLU interface data space
        void allocateMemoryInterface(const int &n, const int &m, const int &nnz);
        void deallocateMemoryInterface();

        /// expert driver
        int expertDriver();

        /// advance status
        char status;

        /// Extra container for keep option
        T *keep_data;
        /// Extra pointer for keep option
        T *keep_orig;

        /// expert driver suplementary data
        int lwork;
        void *work;
        int     *etree; 
        data_sub_type  *R, *C;
        data_sub_type  *ferr, *berr;
        data_sub_type  rpg, rcond;
        char    equed[1];

};



} //end namespace

#include "xLinearSystemSolverSuperLU_imp.h"


#endif





