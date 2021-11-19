/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _SOLVERMUMPS_H
#define _SOLVERMUMPS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mpi.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"
#include "xLinearSystemSolverMumpsTraitPolicy.h"

/* code for empty declaration of icntl control */
#define XLINEARSOLVERMUMPS_EMPTY_CONTROL -999999
/* code for status */
#define XLINEARSOLVERMUMPS_INIT 1U
#define XLINEARSOLVERMUMPS_ALLOC 2U
#define XLINEARSOLVERMUMPS_SYMB 4U
#define XLINEARSOLVERMUMPS_FACT 8U
#define XLINEARSOLVERMUMPS_SOLV 16U
/*
 * Maximum number of control parameter (ictnl size)
 *
 * Depend on version. The maximum among all version
 * is used.
 * Starting with  5.2.0 version 60 control parameters are available
 * */
#define XLINEARSOLVERMUMPS_ICNTL_MX_SIZE 60

namespace xlinalg
{
class xLinearSystemSolverMumpsBase
{
  public:
   // Traits available for all derived class
   typedef xTraitMatrixFindex matrix_indexing;

   enum enum_param_tab
   {
      ICNTL,
      CNTL,
      PAR,
      INTER,
      WRITE_PROBLEM,
      INFOG,
      INFO,
      RINFOG,
      RINFO
   };
   enum enum_param_inter
   {
      VERBOSE = 0,
      NBTHREAD,
      RESETMXNBTHREAD
   };

   xLinearSystemSolverMumpsBase(MPI_Comm communicator);

   /// Method to control Mumps
   //
   //! Since mumps version 5.1.2 open mp feature are available.
   //! Default behavior is one thread per interface instance (living in a communicator).
   //! If user do no modify explicitly this number of thread per interface instance a job
   //! launched on N MPI process will use N cores with mumps what ever OMP_NUM_THREAD environment
   //! variable setting is.
   //! If user set explicitly this number of thread per interface instance to M for example
   //! a N MPI process launched job will use N*M cores on all cluster (still independently of
   //! OMP_NUM_THREAD setting).
   //!
   //! Use setParam with INTER and NBTHREAD to set number of threads per interface instance.
   //! After any method that call mumps a reseting is done by the interface to fix the
   //! maximum number of possible threads to the value it have before the call.
   //! This value is obtained by the interface by calling omp_get_max_threads  at construction time or
   //! if user use setParam with INTER and RESETMXNBTHREAD. In this last case if user is using a
   //! zero value for last argument then omp_get_max_threads is used otherwise last argument
   //! is used (witch is not encouraged as it may hide a OMP_NUM_THREAD environment
   //! variable setting and perturb the rest of the application that may rely on that mechanism)
   //
   //! Note: Controlling number of thread like that permit to distinguish usage of mumps interface
   //! instance for a single instance that use many process (and maybe many thread per process) from
   //! a case where there is many interface that use one process or few process concurrently. In this
   //! last case maybe one thread is expected to be use per interface. Adding this control provide a way
   //! to customize those different instances. If you want to always follow OMP_NUM_THREAD environment
   //! variable setting,just after instantiating, use setParam with INTER and NBTHREAD with the last parameter
   //! set by omp_get_max_threads returned value.
   //
   template <typename P>
   void setParam(int, int, P);

   /// method to retrieve information from Mumps internal data
   template <typename P>
   void getInternalData(int, int, P &);

   /// specific method to show null pivot
   //! It is called automatically with a verbose level greater then 2 in derived class
   void printNullPivot();

  protected:
   // allocating and deallocating Mumps structure
   template <typename T>
   void allocateMemoryStructMumps();
   template <typename S>
   void allocateMemoryStructMumpsPointer();
   template <typename T>
   void deallocateMemoryStructMumps();
   template <typename S>
   void deallocateMemoryStructMumpsPointer();

   // method to call mumps with correct arithmetic
   template <typename T>
   void mumps();

   // general driver to call mumps. This intermediate method impose number of threads to be used
   // by mumps instance controlled by this instance.
   template <typename T>
   void mumpsOmp();

   template <typename T>
   void setDefaultParameter();
   void setDefaultParameterFloat();
   void setDefaultParameterDouble();

   void reInitParSymCom();
   int *&getJobPtr();

   void setSym(xLinearSystemSolverMumpsBase *p, int sym_);

   // pointeur version of somme methode
   void reInitParSymCom(xLinearSystemSolverMumpsBase *p);
   int *&getJobPtr(xLinearSystemSolverMumpsBase *p);
   template <typename T>
   void mumpsOmp(xLinearSystemSolverMumpsBase *p);

   int icntl[XLINEARSOLVERMUMPS_ICNTL_MX_SIZE];  // control parameter to store setting permanantely
   bool cntl[15];  // control parameter to store setting permanantely (uggly as 15 may evolve form version :-( simplest for now)
   double rcntl[15];  // Here a arbitrary choice is made. Whatever T is we feed rcntl and id_cntl with a double. This will
                      //  make no diference regaring real or complex arithmetic chosed as this array is alwayse of type REAL
                      //  but this REAL type differs from float to double depending on arimethic (float for s,c and double for
                      //  d,z) if using d or z no problem if using s or c compiler may complaine about a possible loss of data as
                      //  double value will be casted to float. About To check/test
   char write_problem[256];

   void *mumps_struct;  // real type is SMUMPS_STRUC_C/DMUMPS_STRUC_C/CMUMPS_STRUC_C/ZMUMPS_STRUC_C

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
   void **id_piv_null_list;

   // verbosity of the interface itself
   int verbose;
   // since V5.1.2 number of thread have to controlled: interface setting
   int nb_threads, mx_nb_threads;
   // default value
   int par, comm_fortran, sym;
   MPI_Comm univ;
   int procid;
   // ictnl size (version dependant)
   unsigned short icntl_size;
   // version
   unsigned short major_version, minor_version, subminor_version;
};

/// Base class for all xLinearSystemSolverMumpsMasterAndSlaves derived class ("all" here refere to template instanciations)
//! This class is in charge of managing acces to all instance of xLinearSystemSolverMumpsMasterAndSlaves for slavesActions
//! It shouldn't be used directely. Use xLinearSystemSolverMumpsMasterAndSlaves instead
class xLinearSystemSolverMumpsMasterAndSlavesBase : public xLinearSystemSolverMumpsBase
{
  public:
   // Traits
   // Intentionaly there is no traits here

   xLinearSystemSolverMumpsMasterAndSlavesBase(MPI_Comm communicator);

  protected:
   /// interface instance id
   int num_interface_instance;

   /// send the right code to slaves to exit action loop => enter interface instance loop in slavesAction
   //! It is called only by the root
   void suspendSlavesActions();

   /// It is called only by the root to init action loop of the slave of that interface instance
   void initSlavesActions();

   /// methode to synchronise other process with main computation
   void checkInstance();

   /// give index of the first available interface pointeur in tab_interface_instance
   int getFirstAvailableInterfaceInstance();

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
   static std::vector<xLinearSystemSolverMumpsMasterAndSlavesBase *> tab_interface_instance;
   /// current interface instance in action in slavesActions
   static int current_interface_instance;
};

template <typename T = double>
class xLinearSystemSolverMumpsMasterAndSlaves : public xLinearSystemSolverMumpsMasterAndSlavesBase
{
  public:
   // Traits
   typedef T data_type;
   typedef xLinearSystemSolverMumpsBase::matrix_indexing matrix_indexing;
   typedef xTraitMatrixSparceCOO matrix_storage;

   // constructor/destructor
   xLinearSystemSolverMumpsMasterAndSlaves(MPI_Comm communicator = MPI_COMM_WORLD);
   ~xLinearSystemSolverMumpsMasterAndSlaves() noexcept(false);

   // Connecting a Matrix to the solver
   template <typename M>
   void connectMatrix(M &matrix);

   /// do symbolique phase assuming that a previous connection exist with a matrix. connectMatrix  call this methode by default.
   //! This last point make this methode allmost anaivalable for now.  Lived in the public part but may become private
   void symb();

   /// do numerical factorisation assuming that a previous symb methode was called
   void fact();

   /// do solving the linear system assuming that a previous symb methode was called at least
   //! if a previous call to fact was done the  numerical factorisation is not re-computed
   template <typename V>
   void solve(const V &, V &);

   /// this methode give (i,j) termes of the inverse of A
   void aijm1(int, int *, int *, T *);

   // methode to synchronise other slave process with main computation
   void slavesActions();

  protected:
   // initialisation methodes
   template <typename M>
   void init();
   void initDefaultParameter();

   char status;

   // for compatibility with other interface
   std::vector<T> rhs_gather;
};

template <typename T = double>
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

   // constructor/destructor
   xLinearSystemSolverMumpsDistributed(MPI_Comm communicator = MPI_COMM_WORLD);
   ~xLinearSystemSolverMumpsDistributed() noexcept(false);

   /// Connecting a Matrix to the solver
   template <typename M>
   void connectMatrix(M &matrix, double ratio_reduce_comm_ = 0.);

   /// Connecting a Matrix with Shur complement to the solver
   template <typename M>
   void connectMatrixWithShur(M &matrix, int size_shur, int *list_var_shur, int &shur_mloc, int &shur_nloc, int &nprow,
                              int &npcol, int &Mblock, int &Nblock);
   // Connecting a Shur complement matrix to the solver
   void connectShur(T *shur_matrix, int shur_lld);

   /// do symbolique phase assuming that a previous connection exist with a matrix. connectMatrix  call this methode by default.
   //! This last point make this methode allmost anaivalable for now.  Lived in the public part but may become private
   void symb();

   /// do numerical factorisation assuming that a previous symb methode was called
   void fact();

   /// do solving the linear system assuming that a previous symb methode was called at least
   //! if a previous call to fact was done the  numerical factorisation is not re-computed
   template <typename V, typename U>
   void solve(const V &, U &, bool brodcast = true);

   /// do solving the linear system assuming that a previous symb methode was called at least
   //! if a previous call to fact was done the  numerical factorisation is not re-computed
   //! This version supose that Shur complement matrix is/will be computed and you intend to do one of
   //! the tree resolution available in this context :
   //!       - solve internal problem
   //!       - solve to reduce rhs
   //!       - solve to extand
   template <typename V, typename U, typename W>
   void solveWithShur(const V &, U &, W &, int icntl26, bool brodcast = false);

   /// Reset Mumps instance. This mainly remove factor allready computed if any.
   //! This diconnect solver from already connect Matrix if any
   //! It clean also reducing container
   void reset();

   /// Reduce connected Matrix
   void reduceMatrices();

  protected:
   // initialisation methodes
   template <typename M>
   void init(int n);
   void initDefaultParameter();

   bool schur_connected;

   // reducing members
   bool reduce_comm;
   int procid_reduced;
   double ratio_reduce_comm;
   MPI_Comm reduced_comm;
   MPI_Comm gathering_comm;
   T *unreduced_data;
   int nnz_loc;
   std::vector<int> nnz_loc_per_proc;
   std::vector<int> nnz_loc_disp;

   // instance status
   char status;

   // some rhs type need this extra location to gather on master proc
   std::vector<T> rhs_gather;
};

// this the equivalent of master and slave (matrix is given by master proc) in distributed context (i.e. all proc call connect,
// solve, ...)
template <typename T = double>
class xLinearSystemSolverMumpsCentralized : public xLinearSystemSolverMumpsBase
{
  public:
   // Traits
   typedef T data_type;
   typedef xLinearSystemSolverMumpsBase::matrix_indexing matrix_indexing;
   typedef xTraitMatrixSparceCOO matrix_storage;

   // constructor/destructor
   xLinearSystemSolverMumpsCentralized(MPI_Comm communicator = MPI_COMM_WORLD);
   ~xLinearSystemSolverMumpsCentralized();

   // Connecting a Matrix to the solver
   template <typename M>
   void connectMatrix(M &matrix);

   /// do symbolique phase assuming that a previous connection exist with a matrix. connectMatrix  call this methode by default.
   //! This last point make this methode allmost anaivalable for now.  Lived in the public part but may become private
   void symb();

   /// do numerical factorisation assuming that a previous symb methode was called
   void fact();

   /// do solving the linear system assuming that a previous symb methode was called at least
   //! if a previous call to fact was done the  numerical factorisation is not re-computed
   template <typename V, typename U>
   void solve(const V &, U &, bool brodcast = true);

   /// Reset Mumps instance. This mainly remove factor allready computed if any.
   //! This diconnect solver from already connect Matrix if any
   void reset();

  protected:
   // initialisation methodes
   template <typename M>
   void init();
   void initDefaultParameter();

   char status;

   // some rhs type need this extra location to gather on master proc
   std::vector<T> rhs_gather;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverMumpsException class //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// interface derived class of standart exception
class xLinearSystemSolverMumpsException : public std::exception
{
  public:
   xLinearSystemSolverMumpsException(std::string, std::string, int, std::string, std::string);
   ~xLinearSystemSolverMumpsException() throw() override;
   const char *what() const throw() override;

  private:
   std::string msg;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xLinearSystemSolverMumpsException  class /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace xlinalg

#include "xLinearSystemSolverMumps_imp.h"

#endif
