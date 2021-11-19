/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _SOLVERMUMPS_H
#define _SOLVERMUMPS_H


#include <iostream>
#include <cstdlib>
#include <cstring>
#include "mpi.h"
#include "xCSRVector.h"
#include "xCSRMatrix.h"
#include "xLinearSystemSolverMumpsTraitPolicy.h"

/* code for empty declaration of icntl control */
#define XLINEARSOLVERMUMPS_EMPTY_CONTROL     -999999
/* code for status */
#define XLINEARSOLVERMUMPS_INIT       1U
#define XLINEARSOLVERMUMPS_ALLOC      2U
#define XLINEARSOLVERMUMPS_SYMB       4U
#define XLINEARSOLVERMUMPS_FACT       8U
#define XLINEARSOLVERMUMPS_SOLV      16U


namespace lalg
{

class xLinearSystemSolverMumpsBase
{

    public:
        //Traits available for all derived class
        typedef xTraitMatrixFindex matrix_indexing;

        enum enum_param_tab { ICNTL,CNTL,PAR,INTER,WRITE_PROBLEM,INFOG,INFO,RINFOG,RINFO};

        xLinearSystemSolverMumpsBase (MPI_Comm communicator =  MPI_COMM_WORLD );

        /// methode to controle Mumps
        template< typename P>
        void setParam( int,int,P);
        
        /// methode to retreave information from Mumps internal data
        template< typename P>
        void getInternalData( int,int,P&);


    protected:


        // allocatting and dealocating Mumps structure
        template < typename T >
        void allocateMemoryStructMumps(void);
        template < typename S >
        void allocateMemoryStructMumpsPointer(void);
        template < typename T >
        void deallocateMemoryStructMumps(void);
        template < typename S >
        void deallocateMemoryStructMumpsPointer(void);

        template < typename T >
        void mumps(void);
         
        template < typename T >
        MPI_Datatype MPIDatatype(void);

        template < typename T >
        void setDefaultParameter(void);
        void setDefaultParameterFloat(void);
        void setDefaultParameterDouble(void);

        void reInitParSymCom(void);
        int * & getJobPtr(void);

        void setSym(xLinearSystemSolverMumpsBase* p,int sym_);

        // pointeur version of somme methode
        void reInitParSymCom(xLinearSystemSolverMumpsBase* p);
        int * & getJobPtr(xLinearSystemSolverMumpsBase* p);
        template < typename T > void mumps(xLinearSystemSolverMumpsBase* p);


        int icntl[40]; // control parameter to store setting permanantely (uggly as 40 may evolve form version :-( simplest for now)
        bool cntl[15]; // control parameter to store setting permanantely (uggly as 15 may evolve form version :-( simplest for now)
        double rcntl[15]; // Here a arbitrary choice is made. Whatever T is we feed rcntl and id_cntl with a double. This will
                          //  make no diference regaring real or complex arithmetic chosed as this array is alwayse of type REAL
                          //  but this REAL type differs from float to double depending on arimethic (float for s,c and double for d,z)
                          //  if using d or z no problem
                          //  if using s or c compiler may complaine about a possible loss of data as double value will be casted
                          //  to float. About To check/test
        char write_problem[256];

        void *mumps_struct; //real type is SMUMPS_STRUC_C/DMUMPS_STRUC_C/CMUMPS_STRUC_C/ZMUMPS_STRUC_C

        // pointer to simplify access
        void **id_a;
        void **id_a_loc;
        int *id_comm_fortran;
        int *id_icntl;
        void *id_cntl;
        int *id_infog;
        int *id_info;
        void *id_rinfog;
        void *id_rinfo;
        int **id_irn;
        int **id_irn_loc;
        int **id_jcn;
        int **id_jcn_loc;
        int *id_job;
        int *id_n;
        int *id_nz;
        int *id_nz_loc;
        int *id_par;
        int *id_nrhs;
        int *id_lrhs;
        void **id_rhs;
        int *id_sym;
        int *id_nz_rhs;
        int **id_irhs_sparse;
        int **id_irhs_ptr;
        void **id_rhs_sparse;
        char *id_write_problem;
        int *id_size_shur;
        int **id_listvar_schur;
        int *id_mblock;
        int *id_nblock;
        int *id_nprow;
        int *id_npcol;
        int *id_schur_mloc;
        int *id_schur_nloc;
        int *id_schur_lld;
        void **id_schur;
        void **id_redrhs;
        int *id_lredrhs;

        // verbosity of the interface itself
        int verbose;
        // default value
        int par,comm_fortran,sym;
        MPI_Comm univ;


};


/// Base classe for all xLinearSystemSolverMumpsMasterAndSlaves derived class ("all" here refere to template instanciations)
//! This class is in charge of managing acces to all instance of xLinearSystemSolverMumpsMasterAndSlaves for slavesActions
//! It shouldn't be used directely. Use xLinearSystemSolverMumpsMasterAndSlaves instead
class xLinearSystemSolverMumpsMasterAndSlavesBase : public xLinearSystemSolverMumpsBase
{
    public:
        // Traits
        // Intentionaly there is no traits here

        xLinearSystemSolverMumpsMasterAndSlavesBase(); 


    protected:

        /// interface instance id
        int num_interface_instance;

        /// send the right code to slaves to exit action loop => enter interface instance loop in slavesAction
        //! It is called only by the root
        void suspendSlavesActions(void);

        /// It is called only by the root to init action loop of the slave of that interface instance
        void initSlavesActions(void);
        
        /// methode to synchronise other process with main computation
        void checkInstance(void);

        /// give index of the first available interface pointeur in tab_interface_instance
        int getFirstAvailableInterfaceInstance(void);

        // protected static data : 
        // Implementing these here remove the burden of one static table per signature. This
        // lead to a unique table for all signature. This give a unique id for all instance 
        // independantly of it's signature 
        // In slavesActions it yould have been very penfull to establish a correspondance betwen
        // global instance number and index in approriate (with the good signature) static table.
        //
        /// number of active interface instance 
        static int nb_interface_instance;
        /// table of pointer to interface instance
        static std::vector < xLinearSystemSolverMumpsMasterAndSlavesBase * > tab_interface_instance;
        /// current interface instance in action in slavesActions
        static int current_interface_instance;
};

template < typename T = double  >
class xLinearSystemSolverMumpsMasterAndSlaves : public xLinearSystemSolverMumpsMasterAndSlavesBase
{

    public:
        // Traits
        typedef T data_type;
        typedef xLinearSystemSolverMumpsBase::matrix_indexing matrix_indexing;
        typedef xTraitMatrixSparceCOO matrix_storage;

        //constructor/destructor
        xLinearSystemSolverMumpsMasterAndSlaves ( void );
        ~xLinearSystemSolverMumpsMasterAndSlaves ( void );

        // Connecting a Matrix to the solver
        template < typename M > void connectMatrix( M &matrix);

        /// do symbolique phase assuming that a previous connection exist with a matrix. connectMatrix  call this methode by default.
        //! This last point make this methode allmost anaivalable for now.  Lived in the public part but may become private
        void symb(void);

        /// do numerical factorisation assuming that a previous symb methode was called
        void fact(void);
 
        /// do solving the linear system assuming that a previous symb methode was called at least
        //! if a previous call to fact was done the  numerical factorisation is not re-computed
        template < typename V  >
        void solve(const V &,  V & );

        /// this methode give (i,j) termes of the inverse of A
        void aijm1(int, int *, int *,T *);

        // methode to synchronise other slave process with main computation
        void slavesActions(void);


    protected:

        // initialisation methodes
        template < typename M > void init (void);
        void initDefaultParameter(void);

        char status;

};

template < typename T = double  >
class xLinearSystemSolverMumpsDistributed : public xLinearSystemSolverMumpsBase
{
    public:
        // Traits
        typedef T data_type;
        typedef xLinearSystemSolverMumpsBase::matrix_indexing matrix_indexing;
        typedef xTraitMatrixSparceCOO matrix_storage;
        // TODO : should be 
        // typedef xTraitMatrixSparceDCOO matrix_storage;
        // but work corresponding to DCOO not done wet.
        // for xGenericSparceMatix DCOO implamentation will be the exacte copy of COO one apriori
        // What might be intresting though is in xfem Assemble to use DCOO to reduce (MPI_reduce) right end side terms

        //constructor/destructor
        xLinearSystemSolverMumpsDistributed ( MPI_Comm communicator  = MPI_COMM_WORLD ); 
        ~xLinearSystemSolverMumpsDistributed ( void );

        /// Connecting a Matrix to the solver
        template < typename M > void connectMatrix( M &matrix);

        /// Connecting a Matrix with Shur complement to the solver
        template < typename M > void connectMatrixWithShur( M &matrix, int size_shur, int *list_var_shur, int &shur_mloc, int &shur_nloc, int &nprow=0, int &npcol=0, int &Mblock=0, int &Nblock=0);
        // Connecting a Shur complement matrix to the solver
        void connectShur( T *shur_matrix, int shur_lld);

        /// do symbolique phase assuming that a previous connection exist with a matrix. connectMatrix  call this methode by default.
        //! This last point make this methode allmost anaivalable for now.  Lived in the public part but may become private
        void symb(void);

        /// do numerical factorisation assuming that a previous symb methode was called
        void fact(void);
 
        /// do solving the linear system assuming that a previous symb methode was called at least
        //! if a previous call to fact was done the  numerical factorisation is not re-computed
        template < typename V  >
        void solve(const V &,  V & , bool brodcast=true);

        /// do solving the linear system assuming that a previous symb methode was called at least
        //! if a previous call to fact was done the  numerical factorisation is not re-computed
        //! This version supose that Shur complement matrix is/will be computed and you intend to do one of
        //! the tree resolution available in this context :
        //!       - solve internal problem
        //!       - solve to reduce rhs
        //!       - solve to extand
        template < typename V  >
        void solveWithShur(const V &,  V & , V & ,int icntl26, bool brodcast=false);

        /// Reset Mumps instance. This mainly remove factor allready computed if any.
        //! This diconnect solver from already connect Matrix if any
        void reset(void);

    protected:

        // initialisation methodes
        template < typename M > void init (void);
        void initDefaultParameter(void);

        int procid;

        bool schur_connected;

        char status;

};

// this the equivalent of master and slave (matrix is given by master proc) in distributed context (i.e. all proc call connect, solve, ...)
template < typename T = double  >
class xLinearSystemSolverMumpsCentralized : public xLinearSystemSolverMumpsBase
{
    public:
        // Traits
        typedef T data_type;
        typedef xLinearSystemSolverMumpsBase::matrix_indexing matrix_indexing;
        typedef xTraitMatrixSparceCOO matrix_storage;

        //constructor/destructor
         xLinearSystemSolverMumpsCentralized( void ); 
        ~xLinearSystemSolverMumpsCentralized ( void );

        // Connecting a Matrix to the solver
        template < typename M > void connectMatrix( M &matrix);

        /// do symbolique phase assuming that a previous connection exist with a matrix. connectMatrix  call this methode by default.
        //! This last point make this methode allmost anaivalable for now.  Lived in the public part but may become private
        void symb(void);

        /// do numerical factorisation assuming that a previous symb methode was called
        void fact(void);
 
        /// do solving the linear system assuming that a previous symb methode was called at least
        //! if a previous call to fact was done the  numerical factorisation is not re-computed
        template < typename V  >
        void solve(const V &,  V & , bool brodcast=true);

        /// Reset Mumps instance. This mainly remove factor allready computed if any.
        //! This diconnect solver from already connect Matrix if any
        void reset(void);

    protected:

        // initialisation methodes
        template < typename M > void init (void);
        void initDefaultParameter(void);

        int procid;

        char status;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverMumpsException class ////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// interface derived class of standart exception
class xLinearSystemSolverMumpsException : public std::exception
{
    public:
        xLinearSystemSolverMumpsException(std::string,std::string,int,std::string,std::string);
        ~xLinearSystemSolverMumpsException() throw( );
        virtual const char * what() const throw( );

    private:
        std::string msg;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xLinearSystemSolverMumpsException  class ///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







} //end namespace

#include "xLinearSystemSolverMumps_imp.h"


#endif
