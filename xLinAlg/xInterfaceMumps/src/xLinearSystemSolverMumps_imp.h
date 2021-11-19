/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <complex>
#include <fstream>
#include <iostream>

#include "xMPIDataType.h"

namespace xlinalg
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverMumpsBase class implementation ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename S>
void xLinearSystemSolverMumpsBase ::allocateMemoryStructMumpsPointer()
{
   mumps_struct = (void *)(new S);
   S *id = (S *)(mumps_struct);
   id_a = (void **)(&(id->a));
   id_a_loc = (void **)(&(id->a_loc));
   id_comm_fortran = &(id->comm_fortran);
   id_icntl = &(id->icntl[0]);
   id_cntl = (void *)&(id->cntl[0]);
   id_infog = &(id->infog[0]);
   id_info = &(id->info[0]);
   id_rinfog = (void *)&(id->rinfog[0]);
   id_rinfo = (void *)&(id->rinfo[0]);
   id_irn = &(id->irn);
   id_irn_loc = &(id->irn_loc);
   id_jcn = &(id->jcn);
   id_jcn_loc = &(id->jcn_loc);
   id_job = &(id->job);
   id_n = &(id->n);
   id_nz = &(id->nz);
   id_nz_loc = &(id->nz_loc);
   id_par = &(id->par);
   id_nrhs = &(id->nrhs);
   id_lrhs = &(id->lrhs);
   id_rhs = (void **)(&(id->rhs));
   id_sym = &(id->sym);
   id_nz_rhs = &(id->nz_rhs);
   id_irhs_sparse = &(id->irhs_sparse);
   id_irhs_ptr = &(id->irhs_ptr);
   id_rhs_sparse = (void **)(&(id->rhs_sparse));
   id_write_problem = (char *)(id->write_problem);
   id_size_shur = &(id->size_schur);
   id_listvar_schur = &(id->listvar_schur);
   id_mblock = &(id->mblock);
   id_nblock = &(id->nblock);
   id_nprow = &(id->nprow);
   id_npcol = &(id->npcol);
   id_schur_mloc = &(id->schur_mloc);
   id_schur_nloc = &(id->schur_nloc);
   id_schur_lld = &(id->schur_lld);
   id_schur = (void **)(&(id->schur));
   id_redrhs = (void **)(&(id->redrhs));
   id_lredrhs = &(id->lredrhs);
   id_piv_null_list = (void **)&(id->pivnul_list);
   // ensure that pointer are real NULL; Mumps library do it but
   // this save the potential interogation
   id->irn = id->jcn = nullptr;
   id->a = nullptr;
   id->schur = nullptr;
   id->redrhs = nullptr;
}
template <typename S>
void xLinearSystemSolverMumpsBase ::deallocateMemoryStructMumpsPointer()
{
   delete ((S *)mumps_struct);
}

template <typename P>
void xLinearSystemSolverMumpsBase ::setParam(int param_tab, int param_index, P param_value)
{
   int n;

   // depending on param_tab
   switch (param_tab)
   {
      case CNTL:
      {
         // shift indexation as param_index is consistant with mumps doc (and fortran)
         n = param_index - 1;
         // WARNING WARNING WARNING :
         // here it's P which gives id_cntl type !
         // if type is double and P is float there will be certainely some problems
         // it's user responsability to give the argument type acording to interface instance type used
         ((P *)(id_cntl))[n] = param_value;
         rcntl[n] = (double)param_value;
         cntl[n] = true;
         break;
      }
      default:
      {
         std::ostringstream oss;
         oss << "Invocation of set_param use wrong param_tab ; see xLinearSolverMumps.h enum_param_tab to use a good one  \n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }
}
template <>
void xLinearSystemSolverMumpsBase ::setParam(int param_tab, int param_index, const char *param_value);
template <>
void xLinearSystemSolverMumpsBase ::setParam(int param_tab, int param_index, int param_value);

template <typename P>
void xLinearSystemSolverMumpsBase ::getInternalData(int tab_id, int tab_index, P &value)
{
   // shift indexation as param_index is consistant with mumps doc (and fortran)
   const int n = tab_index - 1;

   // depending on param_tab
   switch (tab_id)
   {
      case INFOG:
      {
         assert((major_version < 5 || (major_version == 5 && minor_version < 2)) ? n < 40 : n < 80);
         value = (P)(id_infog)[n];
         break;
      }
      case INFO:
      {
         assert((major_version < 5 || (major_version == 5 && minor_version < 2)) ? n < 40 : n < 80);
         value = (P)(id_info)[n];
         break;
      }
      case RINFOG:
      {
         assert(n < 40);
         value = (P)((P *)(id_rinfog))[n];
         break;
      }
      case RINFO:
      {
         assert(n < 40);
         value = (P)((P *)(id_rinfo))[n];
         break;
      }
      default:
      {
         std::ostringstream oss;
         oss << "Invocation of getInternalData use wrong tab_id ; see xLinearSolverMumps.h enum_param_tab to use a good one  \n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }
}

template <typename T>
void xLinearSystemSolverMumpsBase::mumpsOmp(xLinearSystemSolverMumpsBase *p)
{
   p->mumpsOmp<T>();
   return;
}

template <typename T>
void xLinearSystemSolverMumpsBase::mumpsOmp()
{
#ifdef _OPENMP
   omp_set_num_threads(nb_threads);
   mumps<T>();
   omp_set_num_threads(mx_nb_threads);
#else
   mumps<T>();
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverMumpsMasterAndSlaves class implementation /////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
xLinearSystemSolverMumpsMasterAndSlaves<T>::xLinearSystemSolverMumpsMasterAndSlaves(MPI_Comm communicator)
    : xLinearSystemSolverMumpsMasterAndSlavesBase(communicator), status(0)
{
   // init mumps struture
   try
   {
      allocateMemoryStructMumps<T>();
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " Mumps data structure was not allocated normaly ?!\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set default parameter
   // nota : a job -1 as not been done yet this call look yeard here but it permit to store icntl,cntl and write_problem
   // information in member of this interface. Information given in Mumps structure will be erased by
   // the job -1 call and re setted by the setDefaultParameter < T >();
   // This is a good way to fix default parameter for the wall life of this object. User may set is one
   // choice at any time. If it's before the first job -1 call it will be the same principal as here :
   //  stored in member interface and erased from mums struct and set again by setDefaultParameter
   initDefaultParameter();
}

template <typename T>
xLinearSystemSolverMumpsMasterAndSlaves<T>::~xLinearSystemSolverMumpsMasterAndSlaves() noexcept(false)
{
   int err;
   if (verbose) std::cout << "Deleting  Mumps instance (" << num_interface_instance << ")" << std::endl;

   // only master have to do this as slave receve this in there ifinite loop
   // doing that from slave may lead to race condition
   // this have to be checked has when a slave finish  it call a general cleaning
   // including a delete of xLinearSystemSolverMumpsMasterAndSlaves object ??? or MPI-reduce is done before ??
   // for security leave this test
   if (!procid)
   {
      // check that slave are using correct interface instance
      // infinite loop on slaves of that interface instance is assumed after this call
      checkInstance();

      // send the deleting  "MUMPS instance" code to slaves of that instance
      (*id_job) = -2;
      err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      // deleting  "MUMPS instance"
      mumpsOmp<T>();
   }

   // free structure of this instance
   try
   {
      deallocateMemoryStructMumps<T>();
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " Mumps matrix structure was not deallocated normaly !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // destroying this instance from interface instance table
   --nb_interface_instance;
   tab_interface_instance[num_interface_instance] =
       nullptr;  // the delete is from this call ! don't recursively delete the instance.

   if (!procid)
   {
      // this is the last instance to be destroyed so we have to finish the infinit loop on object instance
      if (nb_interface_instance < 1)
      {
         int task;
         // set and send the task to stop the instance loop
         task = -1;
         err = MPI_Bcast((void *)&(task), 1, MPI_INT, 0, univ);
         if (err != MPI_SUCCESS)
         {
            std::ostringstream oss;
            oss << " Prb with MPI_Bcast. Error returned : " << err << "\n";
            throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
         }
      }
      // otherwise we have to reset current instance status to unknown as we just remove the current one
      else
      {
         current_interface_instance = -1;
      }
   }
}

template <typename T>
template <typename M>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::init()
{
   // local
   int err, newsym;
   bool do_send_newsym = true;

   if (procid)
   {
      std::ostringstream oss;
      oss << "You shouldn't have call this methode from a process different from master. Here you are in process " << procid
          << "  \n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // synchronize slave
   checkInstance();

   // define sym from traits of matrix
   newsym = xTraitsLinearSystemSolverMumpsMatrixType<typename M::matrix_pattern, typename M::matrix_defined>::sym;

   // if Mumps instance have already been set we are in the case of a re - init
   if (status & XLINEARSOLVERMUMPS_INIT)
   {
      // see if sym change and prepare to send new version to slave if needed
      if (newsym == sym) do_send_newsym = false;

      // if user memory have already been  connected => connection is lost => status lose XLINEARSOLVERMUMPS_ALLOC
      // done naturaly below
      // memory have also to be cleaned in solver space : a -2 have to be done first
      //
      (*id_job) = -2;

      // send to slave
      int err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending job. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      // doing the work
      if (verbose) std::cout << "Begin reinit Mumps instance(" << num_interface_instance << ")" << std::endl;
      mumpsOmp<T>();

      // re-enter slave instance loop for this object
      initSlavesActions();
   }
   else
   {
      if (verbose) std::cout << "Begin init Mumps instance(" << num_interface_instance << ")" << std::endl;
   }

   if (do_send_newsym)
   {
      sym = newsym;
      reInitParSymCom();

      // send to slave sym and start initialisation of this instance
      (*id_job) = -4;

      err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending job. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      err = MPI_Bcast((void *)&(sym), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending sym. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // reset job to -1 for next call
      (*id_job) = -1;
   }
   else
   {
      // set job to -1 for next call
      (*id_job) = -1;
      // informe slaves of next action
      err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending job. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   // doing the work
   // initializing  "MUMPS instance"
   mumpsOmp<T>();

   // checking exit status
   if (id_infog[0] < 0)
   {
      std::ostringstream oss;
      oss << "Prb during (re)initialisation phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set default parameter
   setDefaultParameter<T>();

   // set status whatever it was
   status = XLINEARSOLVERMUMPS_INIT;
}

template <typename T>
template <typename M>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::connectMatrix(M &matrix)
{
   // temporary local
   int n, m, nnz;
   // temporary local pointeur
   int *index1;
   int *index2;
   T *data;
   matrix_storage storage;
   matrix_indexing indexing;

   // attache "matrix space" to the solver
   matrix.getMemoryAccess(n, m, nnz, &index1, &index2, &data, storage, indexing);

   // check dimension
   if (n != m)
   {
      std::ostringstream oss;
      oss << " n and m are note egale ? Mumps deal with square matrix only !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // initialise or re-initialise mumps instance
   init<M>();

   // verbose printing
   if (verbose)
   {
      std::cout << "Matrix size " << n << " x " << n << std::endl;
      std::cout << "Number of non zeros " << nnz << std::endl;
   }

   // set dimension and initial number of non zeros terms
   (*id_n) = n;
   (*id_nz) = nnz;

   // set pointer for matrix structure and values
   (*id_irn) = index1;
   (*id_jcn) = index2;
   (*id_a) = (void *)(data);

   // update status
   status |= XLINEARSOLVERMUMPS_ALLOC;

   // do by default the symbolic phase
   symb();

   // when verbose give extra basic information from sym
   if (verbose)
   {
      getInternalData(INFOG, 8, m);
      std::cout << "Matrix structural symmetry in % (100% = symmetric 0% = fully unsymetric) : " << m << std::endl;
      getInternalData(INFOG, 3, m);
      std::cout << "Estimed number of entries for factor : " << ((m > 0) ? m : ((-m) * 1000000)) << std::endl;
   }
}

template <typename T>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::symb()
{
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC))
   {
      // numeric and others may have already be done >> check
      if (status != (XLINEARSOLVERMUMPS_SYMB - 1))
      {
         std::ostringstream oss;
         oss << " Mumps symbolic phase is not intended for now to by used without re-initialisation of the mumps instance !\n";
         oss << " Reconnect a matrix will do that ! Otherwise you have to add your own implementation\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // check that slave are using correct instance
      checkInstance();

      (*id_job) = 1;

      // send to slave
      int err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending job. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // doing the work
      if (verbose) std::cout << "Begin symbolic factorisation with Mumps instance(" << num_interface_instance << ")" << std::endl;
      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << "Prb during symbolic phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      status |= XLINEARSOLVERMUMPS_SYMB;
   }
   else
   {
      std::ostringstream oss;
      oss << " Symbolic phase is not possible without a previous Mumps instance initialized and a matrix structure connected !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::fact()
{
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // numeric and others may have already be done >> clean
      if (status != (XLINEARSOLVERMUMPS_FACT - 1))
      {
         // to be done
      }

      // check that slave are using correct instance
      checkInstance();

      (*id_job) = 2;

      // send to slave
      int err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending job. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // doing the work
      if (verbose)
         std::cout << "Begin numerical factorisation with Mumps instance(" << num_interface_instance << ")" << std::endl;
      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << " Mumps numerical factorisation phase failed !\n\tid_infog[0]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      if (verbose > 2) printNullPivot();
      status |= XLINEARSOLVERMUMPS_FACT;
   }
   else
   {
      std::ostringstream oss;
      oss << "numerical factorisation phase is not possible without a previous Mumps matrix structure connect, and symbolicly "
             "factorized !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
template <typename V>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::solve(const V &B, V &X)
{
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // previous solve >> clean ??
      /*
         if (status!=(XLINEARSOLVERMUMPS_SOLV-1))
         {
         // to be done ??
         }
       */

      // check that slave are using correct instance
      checkInstance();

      if (status & XLINEARSOLVERMUMPS_FACT)
         (*id_job) = 3;
      else
         (*id_job) = 5;

      // send to slave
      int err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending job. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // duplicate B in X and attache to rhs in mumps structure
      (*id_rhs) = xPolicyLinearSystemSolverMumpsRHSTwoType<T, V, V>::copyBInX(B, X, rhs_gather);

      // set rhs depending on number of colonne
      (*id_nrhs) = xPolicyLinearSystemSolverMumpsRHSOneType<T, V>::nbColRhs(B);
      (*id_lrhs) = *(id_n);

      // doing the work
      if (verbose)
      {
         if (status & XLINEARSOLVERMUMPS_FACT)
            std::cout << "Begin solving with Mumps instance(" << num_interface_instance << ") and factor already calculated"
                      << std::endl;
         else
            std::cout << "Begin numerical factorisation and solving with Mumps instance(" << num_interface_instance << ")"
                      << std::endl;
      }

      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         if (status & XLINEARSOLVERMUMPS_FACT)
            oss << " Mumps solving phase failed !\n";
         else
            oss << "Mumps factorisation and solving phase failed infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      if (status & XLINEARSOLVERMUMPS_FACT)
         status |= XLINEARSOLVERMUMPS_SOLV;
      else
      {
         if (verbose > 2) printNullPivot();
         status |= XLINEARSOLVERMUMPS_SOLV | XLINEARSOLVERMUMPS_FACT;
      }
   }
   else
   {
      std::ostringstream oss;
      oss << " solving phase is not possible without a previous Mumps symbolic phase !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::aijm1(int nz_rhs, int *ptr_irhs_sparse, int *ptr_irhs_ptr, T *ptr_rhs_sparse)
{
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB) &&
       (status & XLINEARSOLVERMUMPS_FACT))
   {
      // previous solve >> clean ??
      if (status != (XLINEARSOLVERMUMPS_SOLV - 1))
      {
         // to be done ??
      }

      // ask for aij-1 termes
      setParam(ICNTL, 30, 1);

      // check that slave are using correct instance
      checkInstance();

      // solve request
      (*id_job) = 3;

      // send to slave
      int err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
      if (err != MPI_SUCCESS)
      {
         std::ostringstream oss;
         oss << " Prb with MPI_Bcast while sending job. Error returned : " << err << "\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // set rhs sparce matrix pointer
      (*id_nrhs) = *(id_n);
      (*id_nz_rhs) = nz_rhs;
      (*id_irhs_sparse) = ptr_irhs_sparse;
      (*id_irhs_ptr) = ptr_irhs_ptr;
      (*id_rhs_sparse) = (void *)ptr_rhs_sparse;

      // doing the work
      if (verbose)
         std::cout << "Begin solving to obtain aij-1 matrix terms with Mumps instance(" << num_interface_instance
                   << ") and factor already calculated" << std::endl;
      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << " Mumps aij-1 phase failed !\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      status |= XLINEARSOLVERMUMPS_SOLV;
   }
   else
   {
      std::ostringstream oss;
      oss << " for now aij-1 phase is not possible without a previous Mumps numerical factorisation !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

// slave treatment
// This is a special methode wich is used when mumps is the unique parallel part of the programe
// It is called only by the slaves. The solv and other metodes are only called by the root nodes.
//
template <typename T>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::slavesActions()
{
   bool do_loop = false;
   int task, err;
   const int n = tab_interface_instance.size();
   int *id_job_local;

   if (nb_interface_instance > 1)
   {
      // all object have to call this methode in the slave part but only one infinite loop have to be run
      // choose arbitrary the first created object
      if (getFirstAvailableInterfaceInstance() == num_interface_instance) do_loop = true;
   }
   else
      do_loop = true;

   if (do_loop)
   {
      // start "infinit loop" on solver object  : instance loop
      do
      {
         // reiceve task information from master
         err = MPI_Bcast((void *)&(task), 1, MPI_INT, 0, MPI_COMM_WORLD);
         if (err != MPI_SUCCESS)
         {
            std::ostringstream oss;
            oss << " Prb with MPI_Bcast while receiving tack  number . Error returned : " << err << "\n";
            throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
         }

         // check task number (index of tab_interface_instance); -1 treated below
         if (task < -1 || task >= n)
         {
            std::ostringstream oss;
            oss << " Task out of range !?  tab_interface_instance have a size=" << n
                << " with nb_interface_instance=" << nb_interface_instance << "\n";
            throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
         }

         // if task=-1 it finish the object loop properly
         // this append when the last instance is deleted
         if (task < 0) continue;

         // set the pointeur to the correct object
         xLinearSystemSolverMumpsMasterAndSlavesBase *ptr_solver = tab_interface_instance[task];
         id_job_local = getJobPtr(ptr_solver);

         // start "infinit loop" to perform mumps job asked by master : action loop for the curent instance
         do
         {
            // receive job action to performe
            err = MPI_Bcast((void *)(id_job_local), 1, MPI_INT, 0, univ);
            if (err != MPI_SUCCESS)
            {
               std::ostringstream oss;
               oss << " Prb with MPI_Bcast while receiving job number . Error returned : " << err << "\n";
               throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
            }

            // special case : jump out of the loop to perfom something else  :
            //       loopin again with a another instance
            // to be use with care as it doesn't clean up "MUMPS instance" properly (-2) and memory is still in place
            if ((*id_job_local) == -3) break;

            // special case : setting sym
            if ((*id_job_local) == -4)
            {
               // receive sym
               int local_sym;
               err = MPI_Bcast((void *)&(local_sym), 1, MPI_INT, 0, univ);
               if (err != MPI_SUCCESS)
               {
                  std::ostringstream oss;
                  oss << " Prb with MPI_Bcast while receiving sym value . Error returned : " << err << "\n";
                  throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
               }
               setSym(ptr_solver, local_sym);

               // sym is correctely fixed we can now initiate a Mumps instance
               (*id_job_local) = -1;
            }

            // in the case of a reinit set par,sym,comm_fortran
            if ((*id_job_local) == -1)
            {
               reInitParSymCom(ptr_solver);
            }

            // do the task
            mumpsOmp<T>(ptr_solver);

            // naturaly stop on that instance as a -2 destroy mumps internal instance
         } while ((*id_job_local) != -2);

         // stop and exit if no more instances are available
      } while (task > -1);
   }
}

template <typename T>
void xLinearSystemSolverMumpsMasterAndSlaves<T>::initDefaultParameter()
{
   // before setting any parameter initialise class member
   for (int i = 0; i < XLINEARSOLVERMUMPS_ICNTL_MX_SIZE; i++) icntl[i] = XLINEARSOLVERMUMPS_EMPTY_CONTROL;
   for (int i = 0; i < 15; i++) cntl[i] = false;
   strcpy(write_problem, "NAME_NOT_INITIALIZED");

   // even if it's their default, set control parameter
   // that folow to have them set for all instance
   // assembled centralized matrix as MumpsWithSlaves naturaly imply them
   setParam(ICNTL, 5, 0);
   setParam(ICNTL, 18, 0);

   // selec Metis as default ordering
   setParam(ICNTL, 7, 5);

   // by default no output from mumps
   setParam(ICNTL, 1, -1);
   setParam(ICNTL, 2, -1);
   setParam(ICNTL, 3, -1);
   setParam(ICNTL, 4, 0);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverMumpsDistributed     class implementation /////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
xLinearSystemSolverMumpsDistributed<T>::xLinearSystemSolverMumpsDistributed(MPI_Comm communicator)
    : xLinearSystemSolverMumpsBase(communicator),
      schur_connected(false),
      reduce_comm(false),
      procid_reduced(procid),
      ratio_reduce_comm(0.),
      reduced_comm(communicator),
      gathering_comm(MPI_COMM_NULL),
      unreduced_data(nullptr),
      nnz_loc(0),
      status(0)
{
   // init mumps struture
   try
   {
      allocateMemoryStructMumps<T>();
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " Mumps data structure was not allocated normaly ?!\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set default parameter
   // nota : a job -1 as not been done yet this call look yeard here but it permit to store icntl,cntl
   // information in member of this interface. Information given in Mumps structure will be erased by
   // the job -1 call and re setted by the setDefaultParameter < T >();
   // This is a good way to fix default parameter for the wall life of this object. User may set is one
   // choice at any time. If it's before the first job -1 call it will be the same principal as here :
   //  stored in member interface and erased from mums struct and set again by setDefaultParameter
   initDefaultParameter();
}

template <typename T>
xLinearSystemSolverMumpsDistributed<T>::~xLinearSystemSolverMumpsDistributed() noexcept(false)
{
   // Resetting Mumps instance to clear memory allocation
   reset();

   if (verbose) std::cout << "Deletting  mumps distributed instance" << std::endl;

   // free structure of this instance
   try
   {
      deallocateMemoryStructMumps<T>();
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " Mumps matrix structure was not deallocated normaly !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
template <typename M>
void xLinearSystemSolverMumpsDistributed<T>::init(int n)
{
   // define sym from traits of matrix
   sym = xTraitsLinearSystemSolverMumpsMatrixType<typename M::matrix_pattern, typename M::matrix_defined>::sym;

   // if Mumps instance have already been set we are in the case of a re - init
   if (status & XLINEARSOLVERMUMPS_INIT) reset();

   //
   int real_nb_proc, nb_proc_max;
   if (ratio_reduce_comm > 0.)
   {
      MPI_Comm_size(univ, &real_nb_proc);
      if (!procid) nb_proc_max = ratio_reduce_comm * n;
      MPI_Bcast(&nb_proc_max, 1, MPI_INT, 0, univ);
      reduce_comm = (nb_proc_max < real_nb_proc);
   }
   else
      reduce_comm = false;

   if (reduce_comm)
   {
      const int chunk = real_nb_proc / nb_proc_max;
      int shift = real_nb_proc % nb_proc_max;
      int over = 1;
      if (shift && procid < shift * (chunk + 1))
      {
         shift = 0;
      }
      else
         over = 0;
      int color = (procid - shift) / (chunk + over);
      // split comm with colors: 1 color per group of proc to be gathered
      MPI_Comm_split(univ, color, 0, &gathering_comm);

      // set now color to rank in  gathering_comm
      MPI_Comm_rank(gathering_comm, &color);

      // only rank 0 of each gathering group are retained
      if (color)
         color = MPI_UNDEFINED;
      else
         color = 1;
      MPI_Comm_split(univ, color, 0, &reduced_comm);
      if (reduced_comm != MPI_COMM_NULL)
         MPI_Comm_rank(reduced_comm, &procid_reduced);
      else
         procid_reduced = -1;

      if (verbose && !procid)
      {
         int red_size;
         MPI_Comm_size(reduced_comm, &red_size);
         std::cout << "Reducing comunicator for Mumps instance to " << red_size << " out of " << real_nb_proc << " processes"
                   << std::endl;
      }

      // now reduced_comm will be used for mumps
      comm_fortran = (int)MPI_Comm_c2f(reduced_comm);
   }

   if (verbose) std::cout << "Begin init Mumps instance" << std::endl;

   // reset/set sym,par,comm
   reInitParSymCom();

   // initializing  "MUMPS instance" : -1 call
   (*id_job) = -1;

   // doing the work
   if (reduced_comm != MPI_COMM_NULL) mumpsOmp<T>();

   // checking exit status
   if (id_infog[0] < 0)
   {
      std::ostringstream oss;
      oss << "Prb during (re)initialisation phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set default parameter
   setDefaultParameter<T>();

   // set status whatever it was
   status = XLINEARSOLVERMUMPS_INIT;
}
template <typename T>
void xLinearSystemSolverMumpsDistributed<T>::reset()
{
   // if Mumps instance have already been set we clear every thing
   if (status & XLINEARSOLVERMUMPS_INIT)
   {
      if (verbose) std::cout << "Resetting Mumps instance" << std::endl;

      if (reduce_comm)
      {
         if (reduced_comm != MPI_COMM_NULL)
         {
            // local memory allocation must be freed
            if ((*id_irn_loc)) delete ((*id_irn_loc));
            if ((*id_jcn_loc)) delete ((*id_jcn_loc));
            if ((*id_a_loc)) delete ((T *)(*id_a_loc));
         }
         // container have to be cleared
         unreduced_data = nullptr;
         nnz_loc_per_proc.clear();
         nnz_loc_disp.clear();
      }

      // cleanning connection
      (*id_irn_loc) = nullptr;
      (*id_jcn_loc) = nullptr;
      (*id_a_loc) = nullptr;

      // if user memory have already been  connected => connection is lost => status lose XLINEARSOLVERMUMPS_ALLOC
      // done naturaly below
      // memory have also to be cleaned in solver space : a -2 have to be done first
      //
      (*id_job) = -2;

      // doing the work
      if (reduced_comm != MPI_COMM_NULL) mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << "Prb during resetting phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // reset Shur computation
      setParam(ICNTL, 19, 0);
      (*id_size_shur) = 0;
      schur_connected = false;

      // reset comunicator if reduced one active
      if (reduce_comm)
      {
         if (reduced_comm != MPI_COMM_NULL)
         {
            // xtool::xMPITag::resetTag(reduced_comm);// uncoment if dataExchanger in use
            MPI_Comm_free(&reduced_comm);
         }

         assert(gathering_comm != MPI_COMM_NULL);
         // xtool::xMPITag::resetTag(gathering_comm);// uncoment if dataExchanger in use
         MPI_Comm_free(&gathering_comm);

         reduced_comm = univ;
         gathering_comm = MPI_COMM_NULL;
         procid_reduced = procid;
         comm_fortran = (int)MPI_Comm_c2f(univ);
         reduce_comm = false;
      }
   }

   status = 0;
   return;
}
template <typename T>
void xLinearSystemSolverMumpsDistributed<T>::reduceMatrices()
{
   // check
   if (reduce_comm)
   {
      if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC))
      {
         if (verbose) std::cout << "Reducing matrices on reduced communicator master" << std::endl;
         MPI_Gatherv(unreduced_data, nnz_loc, xtool::xMPIDataType<T>(), (*id_a_loc), &nnz_loc_per_proc[0], &nnz_loc_disp[0],
                     xtool::xMPIDataType<T>(), 0, gathering_comm);
      }
      else
      {
         std::ostringstream oss;
         oss << "Reducing  matrices is not possible without a previous Mumps matrix structure connection !\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }
}

template <typename T>
template <typename M>
void xLinearSystemSolverMumpsDistributed<T>::connectMatrix(M &matrix, double ratio_reduce_comm_)
{
   // temporary local
   int n_loc, m_loc;
   // temporary local pointeur
   int *index1;
   int *index2;
   T *data;
   matrix_storage storage;
   matrix_indexing indexing;

   // attache "matrix space" to the solver
   matrix.getMemoryAccess(n_loc, m_loc, nnz_loc, &index1, &index2, &data, storage, indexing);

   // no dimension check as in distributed paradigm local matrix may be of any shape
   // n_loc!=m_loc looks yeard but is possible

   // find dimention of matrix by inspecting max index of local irn and jcn
   int n_max = std::max(*std::max_element(index1, index1 + nnz_loc), *std::max_element(index2, index2 + nnz_loc));
   int n = n_max;

   // max on all proc
   int err = MPI_Reduce((void *)&n_max, (void *)&n, 1, MPI_INT, MPI_MAX, 0, univ);
   if (err != MPI_SUCCESS)
   {
      std::ostringstream oss;
      oss << " Prb with MPI_Reduce while finding matrix dimention. Error returned : " << err << "\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // change ratio
   ratio_reduce_comm = ratio_reduce_comm_;

   // initialise or re-initialise mumps instance
   init<M>(n);

   // set dimension
   (*id_n) = n;  // n is correct only for master but has it is only inpected by master no prb

   // gather data if reduce
   if (reduce_comm)
   {
      assert(gathering_comm != MPI_COMM_NULL);
      int nb_proc_grp;
      MPI_Comm_size(gathering_comm, &nb_proc_grp);
      int procid_grp;
      MPI_Comm_rank(gathering_comm, &procid_grp);

      // create nb per proc and disp
      nnz_loc_per_proc.resize((procid_grp) ? 1 : nb_proc_grp, 0);
      nnz_loc_disp.resize((procid_grp) ? 1 : nb_proc_grp + 1, 0);
      MPI_Gather(&nnz_loc, 1, MPI_INT, &nnz_loc_per_proc[0], 1, MPI_INT, 0, gathering_comm);
      if (!procid_grp)
         for (int p = 0; p < nb_proc_grp; ++p)
         {
            nnz_loc_disp[p + 1] = nnz_loc_disp[p] + nnz_loc_per_proc[p];
         }
      int nnz_loc_grp = (procid_grp) ? 1 : nnz_loc_disp[nb_proc_grp];

      // set for all (0 for non null proc of gathering_comm);
      (*id_nz_loc) = nnz_loc_grp;

      // only master of group (member of mumps comm) will allocated things
      if (!procid_grp)
      {
         // alloc buffer bluntly: contribution from all group process are stored
         // in the buffer without correct summation of common term. Mumps do it.
         // It would reduce final memory consumption to do the sum but add computation work to
         // remove redundant terms. Anyway memory at some point would be rather the same as
         // redundant information have to be exchanged. Here number of exchange is minimized
         // and thus receiving buffer have to be with a size including redundant terms.
         // May be changed in future if problematic.
         (*id_irn_loc) = new int[nnz_loc_grp];
         (*id_jcn_loc) = new int[nnz_loc_grp];
         (*id_a_loc) = (void *)(new T[nnz_loc_grp]);
      }

      // gather indexes into buffer
      MPI_Gatherv(index1, nnz_loc, MPI_INT, (*id_irn_loc), &nnz_loc_per_proc[0], &nnz_loc_disp[0], MPI_INT, 0, gathering_comm);
      MPI_Gatherv(index2, nnz_loc, MPI_INT, (*id_jcn_loc), &nnz_loc_per_proc[0], &nnz_loc_disp[0], MPI_INT, 0, gathering_comm);

      // save pointer
      unreduced_data = data;

      // verbose printing
      if (verbose)
      {
         if (!procid) std::cout << "Matrix size " << n << " x " << n << std::endl;
         std::cout << "Number of non zeros of the local matrix for proc " << procid << " : " << nnz_loc << std::endl;
         if (!procid_grp) std::cout << "Number of non zeros of the reduced local matrix: " << nnz_loc_grp << std::endl;
      }
   }
   else
   {
      // verbose printing
      if (verbose)
      {
         if (!procid) std::cout << "Matrix size " << n << " x " << n << std::endl;
         std::cout << "Number of non zeros of the local matrix for proc " << procid << " : " << nnz_loc << std::endl;
      }

      // set initial number of non zeros terms
      (*id_nz_loc) = nnz_loc;

      // set pointer for matrix structure and values
      (*id_irn_loc) = index1;
      (*id_jcn_loc) = index2;
      (*id_a_loc) = (void *)(data);
   }

   // update status
   status |= XLINEARSOLVERMUMPS_ALLOC;

   // do by default the reducing of current data
   reduceMatrices();

   // do by default the symbolic phase
   symb();

   // when verbose give extra basic information from symb
   if (verbose)
   {
      if (!procid)
      {
         int m;
         getInternalData(INFOG, 8, m);
         std::cout << "Matrix structural symmetry in % (100% = symmetric 0% = fully unsymetric)" << m << std::endl;
         getInternalData(INFOG, 3, m);
         std::cout << "Estimed number of entries for factor : " << ((m > 0) ? m : ((-m) * 1000000)) << std::endl;
      }
   }
}

template <typename T>
template <typename M>
void xLinearSystemSolverMumpsDistributed<T>::connectMatrixWithShur(M &matrix, int size_shur, int *list_var_shur, int &shur_mloc,
                                                                   int &shur_nloc, int &nprow, int &npcol, int &Mblock,
                                                                   int &Nblock)
{
   assert(!reduce_comm);
   // temporary local
   int n_loc, m_loc, nnz_loc;
   // temporary local pointeur
   int *index1;
   int *index2;
   T *data;
   matrix_storage storage;
   matrix_indexing indexing;

   // attache "matrix space" to the solver
   matrix.getMemoryAccess(n_loc, m_loc, nnz_loc, &index1, &index2, &data, storage, indexing);

   // no dimension check as in distributed paradigm local matrix may be of any shape
   // n_loc!=m_loc looks yeard but is possible

   // find dimention of matrix by inspecting max index of local irn and jcn
   int n_max = std::max(*std::max_element(index1, index1 + nnz_loc), *std::max_element(index2, index2 + nnz_loc));
   int n = n_max;

   // max on all proc
   int err = MPI_Reduce((void *)&n_max, (void *)&n, 1, MPI_INT, MPI_MAX, 0, univ);
   if (err != MPI_SUCCESS)
   {
      std::ostringstream oss;
      oss << " Prb with MPI_Reduce while finding matrix dimention. Error returned : " << err << "\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // initialise or re-initialise mumps instance
   init<M>(0);

   // verbose printing
   if (verbose)
   {
      if (!procid) std::cout << "Matrix size " << n << " x " << n << std::endl;
      std::cout << "Number of non zeros of the local matrix for proc " << procid << " : " << nnz_loc << std::endl;
   }

   // set dimension and initial number of non zeros terms
   (*id_n) = n;  // n is correct only for master but has it is only inpected by master no prb
   (*id_nz_loc) = nnz_loc;

   // set dimension and id's of Shur terms
   (*id_size_shur) = size_shur;  // same as n only needed on master
   (*id_listvar_schur) = list_var_shur;

   // set parameter for 2D block cyclic matrix
   (*id_nprow) = nprow;
   (*id_npcol) = npcol;
   (*id_nblock) = Nblock;
   (*id_mblock) = Mblock;

   // force config setting
   // Says to mumps that we want schur complement matrix
   /*
    * NOTE:  the real way to do that would have been to use this test:
      if (sym)
       setParam(ICNTL,19,2);
      else
      setParam(ICNTL,19,3);
    *
    * But this imply a lot of work if we want to merge common terms of
    * schur matrix in between domain. To simplify, as anyway memory
    * consumption is the same, we always ask for the full matrix (unsym).
    */
   setParam(ICNTL, 19, 3);

   // set pointer for matrix structure and values
   (*id_irn_loc) = index1;
   (*id_jcn_loc) = index2;
   (*id_a_loc) = (void *)(data);

   // update status
   status |= XLINEARSOLVERMUMPS_ALLOC;

   // do by default the symbolic phase
   symb();

   // set shur info
   shur_mloc = (*id_schur_mloc);
   shur_nloc = (*id_schur_nloc);
   nprow = (*id_nprow);
   npcol = (*id_npcol);
   Nblock = (*id_nblock);
   Mblock = (*id_mblock);

   // when verbose give extra basic information from symb
   if (verbose)
   {
      if (!procid)
      {
         int m;
         getInternalData(INFOG, 8, m);
         std::cout << "Matrix structural symmetry in % (100% = symmetric 0% = fully unsymetric)" << m << std::endl;
         getInternalData(INFOG, 3, m);
         std::cout << "Estimed number of entries for factor " << ((m > 0) ? m : ((-m) * 1000000)) << std::endl;
         std::cout << "Schur matrix global dimension " << size_shur << " x " << size_shur << std::endl;
         std::cout << "Schur matrix 2D cyclic parameter :" << std::endl;
         std::cout << "   nprow " << nprow << std::endl;
         std::cout << "   npcol " << npcol << std::endl;
         std::cout << "   Nblock " << Nblock << std::endl;
         std::cout << "   Mblock " << Mblock << std::endl;
      }
      std::cout << "Schur matrix local dimension " << shur_mloc << " x " << shur_nloc << std::endl;
   }
}

template <typename T>
void xLinearSystemSolverMumpsDistributed<T>::connectShur(T *shur_matrix, int shur_lld)
{
   assert(!reduce_comm);
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // When user force Schur matrix to be on root only the id_size_shur is null on other proc. So this test is only
      // done on root proc
      if (!procid && !(*id_size_shur))
      {
         std::ostringstream oss;
         oss << " You try to connect a Schur Matrix but you didn't previously define its structure during symbolic phase !\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      (*id_schur_lld) = shur_lld;
      (*id_schur) = (void *)(shur_matrix);
      schur_connected = true;
   }
   else
   {
      std::ostringstream oss;
      oss << "Connecting Shur matrix is not possible without a previous Mumps matrix structure connect, and symbolicly "
             "factorized !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverMumpsDistributed<T>::symb()
{
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC))
   {
      // numeric and others may have already be done >> check
      if (status != (XLINEARSOLVERMUMPS_SYMB - 1))
      {
         std::ostringstream oss;
         oss << " Mumps symbolic phase is not intended for now to by used without re-initialisation of the mumps instance !\n";
         oss << " Reconnect a matrix will do that ! Otherwise you have to add your own implementation\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      if (verbose) std::cout << "Begin symbolic factorisation " << std::endl;

      // doing the work
      (*id_job) = 1;
      if (reduced_comm != MPI_COMM_NULL) mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << "Prb during symbolic phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      status |= XLINEARSOLVERMUMPS_SYMB;
   }
   else
   {
      std::ostringstream oss;
      oss << " Symbolic phase is not possible without a previous Mumps instance initialized and a matrix structure connected !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverMumpsDistributed<T>::fact()
{
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // numeric and others may have already be done >> clean
      if (status != (XLINEARSOLVERMUMPS_FACT - 1))
      {
         // to be done
      }

      if ((*id_size_shur) != 0 && !schur_connected)
      {
         std::ostringstream oss;
         oss << " Factorisation phase is not possible without connecting Schur matrix as Schur complement asked !\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      if (verbose) std::cout << "Begin numerical factorisation" << std::endl;

      // doing the work
      (*id_job) = 2;
      if (reduced_comm != MPI_COMM_NULL) mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << " Mumps numerical factorisation phase failed !\n\tid_infog[0]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      if (verbose > 2) printNullPivot();
      status |= XLINEARSOLVERMUMPS_FACT;
   }
   else
   {
      std::ostringstream oss;
      oss << "numerical factorisation phase is not possible without a previous Mumps matrix structure connect, and symbolicly "
             "factorized !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
template <typename V, typename U>
void xLinearSystemSolverMumpsDistributed<T>::solve(const V &B, U &X, bool brodcast)
{
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // previous solve >> clean ??
      /*
         if (status!=(XLINEARSOLVERMUMPS_SOLV-1))
         {
         // to be done ??
         }
       */

      if (status & XLINEARSOLVERMUMPS_FACT)
         (*id_job) = 3;
      else
         (*id_job) = 5;

      // duplicate B in X and attache to rhs in mumps structure
      (*id_rhs) = xPolicyLinearSystemSolverMumpsRHSTwoType<T, V, U>::copyBInX(B, X, rhs_gather);

      // set rhs depending on number of colonne
      (*id_nrhs) = xPolicyLinearSystemSolverMumpsRHSOneType<T, V>::nbColRhs(B);
      (*id_lrhs) = *(id_n);

      // doing the work
      if (verbose)
      {
         if (status & XLINEARSOLVERMUMPS_FACT)
            std::cout << "Begin solving  with factor already calculated" << std::endl;
         else
            std::cout << "Begin numerical factorisation and solving " << std::endl;
      }

      if (reduced_comm != MPI_COMM_NULL) mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         if (status & XLINEARSOLVERMUMPS_FACT)
            oss << " Mumps solving phase failed !\n";
         else
            oss << "Mumps factorisation and solving phase failed infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      xPolicyLinearSystemSolverMumpsRHSTwoType<T, V, U>::brodcast(X, rhs_gather, *id_n, brodcast, xtool::xMPIDataType<T>(), univ);

      // release memory
      rhs_gather.clear();

      if (status & XLINEARSOLVERMUMPS_FACT)
         status |= XLINEARSOLVERMUMPS_SOLV;
      else
      {
         if (verbose > 2) printNullPivot();
         status |= XLINEARSOLVERMUMPS_SOLV | XLINEARSOLVERMUMPS_FACT;
      }
   }
   else
   {
      std::ostringstream oss;
      oss << " solving phase is not possible without a previous Mumps symbolic phase !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
template <typename V, typename U, typename W>
void xLinearSystemSolverMumpsDistributed<T>::solveWithShur(const V &B, U &X, W &R, int icntl26, bool brodcast)
{
   assert(!reduce_comm);
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // previous solve >> clean ??
      /*
         if (status!=(XLINEARSOLVERMUMPS_SOLV-1))
         {
         // to be done ??
         }
       */
      if (!procid && !(*id_size_shur))
      {
         std::ostringstream oss;
         oss << " Use solve and not this methode as Schur was not properly set !\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      else if (!schur_connected)
      {
         std::ostringstream oss;
         oss << " Solve phase is not possible without proper previous connection of a Schur matrix as Schur complement asked !\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      if (status & XLINEARSOLVERMUMPS_FACT)
         (*id_job) = 3;
      else
      {
         (*id_job) = 5;
      }

      // set type of resolution
      setParam(ICNTL, 26, icntl26);

      switch (icntl26)
      {
         case 1:
         case 2:
         {
            // set reduced/solution rhs array
            (*id_redrhs) = xPolicyLinearSystemSolverMumpsRHSOneType<T, W>::beginRhs(R);
            (*id_lredrhs) = xPolicyLinearSystemSolverMumpsRHSOneType<T, W>::nbRowRhs(R);
         }
         case 0:
         {
            // duplicate B in X and attache to rhs in mumps structure
            (*id_rhs) = xPolicyLinearSystemSolverMumpsRHSTwoType<T, V, U>::copyBInX(B, X, rhs_gather);

            // set rhs depending on number of column
            (*id_nrhs) = xPolicyLinearSystemSolverMumpsRHSOneType<T, V>::nbColRhs(B);
            (*id_lrhs) = *(id_n);
         }
      }

      // doing the work
      if (verbose)
      {
         if (status & XLINEARSOLVERMUMPS_FACT)
            std::cout << "Begin solving  with factor already calculated" << std::endl;
         else
            std::cout << "Begin numerical factorisation and solving " << std::endl;
      }

      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         if (status & XLINEARSOLVERMUMPS_FACT)
            oss << " Mumps solving phase failed !\n";
         else
            oss << "Mumps factorisation and solving phase failed infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      xPolicyLinearSystemSolverMumpsRHSTwoType<T, V, U>::brodcast(X, rhs_gather, *id_n, brodcast, xtool::xMPIDataType<T>(), univ);

      // release memory
      rhs_gather.clear();

      if (status & XLINEARSOLVERMUMPS_FACT)
         status |= XLINEARSOLVERMUMPS_SOLV;
      else
         status |= XLINEARSOLVERMUMPS_SOLV | XLINEARSOLVERMUMPS_FACT;

      // reset
      setParam(ICNTL, 26, 0);
   }
   else
   {
      std::ostringstream oss;
      oss << " solving phase is not possible without a previous Mumps symbolic phase !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverMumpsDistributed<T>::initDefaultParameter()
{
   // before setting any parameter initialise class member
   for (int i = 0; i < XLINEARSOLVERMUMPS_ICNTL_MX_SIZE; i++) icntl[i] = XLINEARSOLVERMUMPS_EMPTY_CONTROL;
   for (int i = 0; i < 15; i++) cntl[i] = false;
   strcpy(write_problem, "NAME_NOT_INITIALIZED");

   // Assembled form
   setParam(ICNTL, 5, 0);

   // Both structure and data distributed matrix as MumpsDistributed naturaly imply them
   setParam(ICNTL, 18, 3);

   // selec Metis as default ordering
   setParam(ICNTL, 7, 5);

   // by default no output from mumps
   setParam(ICNTL, 1, -1);
   setParam(ICNTL, 2, -1);
   setParam(ICNTL, 3, -1);
   setParam(ICNTL, 4, 0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  xLinearSystemSolverMumpsCentralized    class implementation /////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
xLinearSystemSolverMumpsCentralized<T>::xLinearSystemSolverMumpsCentralized(MPI_Comm communicator)
    : xLinearSystemSolverMumpsBase(communicator), status(0)
{
   // init mumps struture
   try
   {
      allocateMemoryStructMumps<T>();
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " Mumps data structure was not allocated normaly ?!\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set default parameter
   // nota : a job -1 as not been done yet this call look yeard here but it permit to store icntl,cntl
   // information in member of this interface. Information given in Mumps structure will be erased by
   // the job -1 call and re setted by the setDefaultParameter < T >();
   // This is a good way to fix default parameter for the wall life of this object. User may set is one
   // choice at any time. If it's before the first job -1 call it will be the same principal as here :
   //  stored in member interface and erased from mums struct and set again by setDefaultParameter
   initDefaultParameter();
}

template <typename T>
xLinearSystemSolverMumpsCentralized<T>::~xLinearSystemSolverMumpsCentralized()
{
   // Resseting Mumps instance to clear memory allocation
   reset();

   if (verbose) std::cout << "Deletting  mumps centralized instance" << std::endl;

   // free structure of this instance
   try
   {
      deallocateMemoryStructMumps<T>();
   }
   catch (std::bad_alloc &e)
   {
      std::cout << "Standard exception: " << e.what() << std::endl;
      std::ostringstream oss;
      oss << " Mumps matrix structure was not deallocated normaly !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
template <typename M>
void xLinearSystemSolverMumpsCentralized<T>::init()
{
   // define sym from traits of matrix
   sym = xTraitsLinearSystemSolverMumpsMatrixType<typename M::matrix_pattern, typename M::matrix_defined>::sym;

   // if Mumps instance have already been set we are in the case of a re - init
   if (status & XLINEARSOLVERMUMPS_INIT)
   {
      // if user memory have already been  connected => connection is lost => status lose XLINEARSOLVERMUMPS_ALLOC
      // done naturaly below
      // memory have also to be cleaned in solver space : a -2 have to be done first
      //
      (*id_job) = -2;

      if (verbose) std::cout << "Begin reinit Mumps instance" << std::endl;

      // doing the work
      mumpsOmp<T>();
   }
   else
   {
      if (verbose) std::cout << "Begin init Mumps instance" << std::endl;
   }

   // reset/set sym,par,comm
   reInitParSymCom();

   // initializing  "MUMPS instance" : -1 call
   (*id_job) = -1;

   // doing the work
   mumpsOmp<T>();

   // checking exit status
   if (id_infog[0] < 0)
   {
      std::ostringstream oss;
      oss << "Prb during (re)initialisation phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // set default parameter
   setDefaultParameter<T>();

   // set status whatever it was
   status = XLINEARSOLVERMUMPS_INIT;
}
template <typename T>
void xLinearSystemSolverMumpsCentralized<T>::reset()
{
   // if Mumps instance have already been set we clear every thing
   if (status & XLINEARSOLVERMUMPS_INIT)
   {
      if (verbose) std::cout << "Resetting Mumps instance" << std::endl;

      // if user memory have already been  connected => connection is lost => status lose XLINEARSOLVERMUMPS_ALLOC
      // done naturaly below
      // memory have also to be cleaned in solver space : a -2 have to be done first
      //
      (*id_job) = -2;

      // doing the work
      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << "Prb during resetting phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   status = 0;
   return;
}

template <typename T>
template <typename M>
void xLinearSystemSolverMumpsCentralized<T>::connectMatrix(M &matrix)
{
   // temporary local
   int n, m, nnz;
   // temporary local pointeur
   int *index1;
   int *index2;
   T *data;
   matrix_storage storage;
   matrix_indexing indexing;

   // only master really connect to given matrix
   if (!procid)
   {
      // attache "matrix space" to the solver
      matrix.getMemoryAccess(n, m, nnz, &index1, &index2, &data, storage, indexing);

      // check dimension
      if (n != m)
      {
         std::ostringstream oss;
         oss << " n and m are note egale ? Mumps deal with square matrix only !\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
   }

   // initialise or re-initialise mumps instance
   init<M>();

   if (!procid)
   {
      // verbose printing
      if (verbose)
      {
         std::cout << "Matrix size " << n << " x " << n << std::endl;
         std::cout << "Number of non zeros " << nnz << std::endl;
      }

      // set dimension and initial number of non zeros terms
      (*id_n) = n;
      (*id_nz) = nnz;

      // set pointer for matrix structure and values
      (*id_irn) = index1;
      (*id_jcn) = index2;
      (*id_a) = (void *)(data);
   }

   // update status
   status |= XLINEARSOLVERMUMPS_ALLOC;

   // do by default the symbolic phase
   symb();

   // when verbose give extra basic information from sym
   if (verbose)
   {
      getInternalData(INFOG, 8, m);
      std::cout << "Matrix structural symmetry in % (100% = symmetric 0% = fully unsymetric) : " << m << std::endl;
      getInternalData(INFOG, 3, m);
      std::cout << "Estimed number of entries for factor : " << ((m > 0) ? m : ((-m) * 1000000)) << std::endl;
   }
}

template <typename T>
void xLinearSystemSolverMumpsCentralized<T>::symb()
{
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC))
   {
      // numeric and others may have already be done >> check
      if (status != (XLINEARSOLVERMUMPS_SYMB - 1))
      {
         std::ostringstream oss;
         oss << " Mumps symbolic phase is not intended for now to by used without re-initialisation of the mumps instance !\n";
         oss << " Reconnect a matrix will do that ! Otherwise you have to add your own implementation\n";
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      if (verbose) std::cout << "Begin symbolic factorisation " << std::endl;

      // doing the work
      (*id_job) = 1;
      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << "Prb during symbolic phase infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      status |= XLINEARSOLVERMUMPS_SYMB;
   }
   else
   {
      std::ostringstream oss;
      oss << " Symbolic phase is not possible without a previous Mumps instance initialized and a matrix structure connected !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverMumpsCentralized<T>::fact()
{
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // numeric and others may have already be done >> clean
      if (status != (XLINEARSOLVERMUMPS_FACT - 1))
      {
         // to be done
      }

      if (verbose) std::cout << "Begin numerical factorisation" << std::endl;

      // doing the work
      (*id_job) = 2;
      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         oss << " Mumps numerical factorisation phase failed !\n\tid_infog[0]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }
      if (verbose > 2) printNullPivot();
      status |= XLINEARSOLVERMUMPS_FACT;
   }
   else
   {
      std::ostringstream oss;
      oss << "numerical factorisation phase is not possible without a previous Mumps matrix structure connect, and symbolicly "
             "factorized !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
template <typename V, typename U>
void xLinearSystemSolverMumpsCentralized<T>::solve(const V &B, U &X, bool brodcast)
{
   // check
   if ((status & XLINEARSOLVERMUMPS_INIT) && (status & XLINEARSOLVERMUMPS_ALLOC) && (status & XLINEARSOLVERMUMPS_SYMB))
   {
      // previous solve >> clean ??
      /*
         if (status!=(XLINEARSOLVERMUMPS_SOLV-1))
         {
         // to be done ??
         }
       */

      if (status & XLINEARSOLVERMUMPS_FACT)
         (*id_job) = 3;
      else
         (*id_job) = 5;

      // duplicate B in X and attache to rhs in mumps structure
      (*id_rhs) = xPolicyLinearSystemSolverMumpsRHSTwoType<T, V, U>::copyBInX(B, X, rhs_gather);

      // set rhs depending on number of colonne
      (*id_nrhs) = xPolicyLinearSystemSolverMumpsRHSOneType<T, V>::nbColRhs(B);
      (*id_lrhs) = *(id_n);

      // doing the work
      if (verbose)
      {
         if (status & XLINEARSOLVERMUMPS_FACT)
            std::cout << "Begin solving  with factor already calculated" << std::endl;
         else
            std::cout << "Begin numerical factorisation and solving " << std::endl;
      }

      mumpsOmp<T>();

      // checking exit status
      if (id_infog[0] < 0)
      {
         std::ostringstream oss;
         if (status & XLINEARSOLVERMUMPS_FACT)
            oss << " Mumps solving phase failed !\n";
         else
            oss << "Mumps factorisation and solving phase failed infog[1]=" << id_infog[0] << " infog[2]=" << id_infog[1];
         throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      xPolicyLinearSystemSolverMumpsRHSTwoType<T, V, U>::brodcast(X, rhs_gather, *id_n, brodcast, xtool::xMPIDataType<T>(), univ);

      // release memory
      rhs_gather.clear();

      if (status & XLINEARSOLVERMUMPS_FACT)
         status |= XLINEARSOLVERMUMPS_SOLV;
      else
      {
         if (verbose > 2) printNullPivot();
         status |= XLINEARSOLVERMUMPS_SOLV | XLINEARSOLVERMUMPS_FACT;
      }
   }
   else
   {
      std::ostringstream oss;
      oss << " solving phase is not possible without a previous Mumps symbolic phase !\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

template <typename T>
void xLinearSystemSolverMumpsCentralized<T>::initDefaultParameter()
{
   // before setting any parameter initialise class member
   for (int i = 0; i < XLINEARSOLVERMUMPS_ICNTL_MX_SIZE; i++) icntl[i] = XLINEARSOLVERMUMPS_EMPTY_CONTROL;
   for (int i = 0; i < 15; i++) cntl[i] = false;
   strcpy(write_problem, "NAME_NOT_INITIALIZED");

   // assembled centralized matrix as xLinearSystemSolverMumpsCentralized imply that
   setParam(ICNTL, 5, 0);
   setParam(ICNTL, 18, 0);

   // selec Metis as default ordering
   setParam(ICNTL, 7, 5);

   // by default no output from mumps
   setParam(ICNTL, 1, -1);
   setParam(ICNTL, 2, -1);
   setParam(ICNTL, 3, -1);
   setParam(ICNTL, 4, 0);
}

}  // namespace xlinalg
