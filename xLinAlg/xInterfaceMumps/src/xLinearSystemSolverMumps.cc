/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <iostream>
#include <sstream>
//#include <complex>

#include "xLinearSystemSolverMumps.h"

// as template lead to use all type of version
// all includes are used
#include "cmumps_c.h"
#include "dmumps_c.h"
#include "smumps_c.h"
#include "zmumps_c.h"

using namespace std;

namespace xtool
{
template <>
MPI_Datatype xMPIDataType<mumps_double_complex>()
{
   return MPI_DOUBLE_COMPLEX;
}

template <>
MPI_Datatype xMPIDataType<mumps_complex>()
{
   return MPI_COMPLEX;
}
}

namespace xlinalg
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// xLinearSystemSolverMumpsBase class implementation ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

xLinearSystemSolverMumpsBase::xLinearSystemSolverMumpsBase(MPI_Comm communicator)
    : mumps_struct(nullptr),
      verbose(1),
      nb_threads(1)
#ifdef _OPENMP
      ,
      mx_nb_threads(omp_get_max_threads())
#else
      ,
      mx_nb_threads(1)
#endif
      ,
      par(1),
      comm_fortran((int)MPI_Comm_c2f(communicator)),
      univ(communicator),
      procid(0),
      icntl_size(40),
      major_version(4),
      minor_version(7),
      subminor_version(0)
{
// extract version if possible, macro appears in 4.8.0 version
// If no version is given then default setting is 4.7.0
#ifdef MUMPS_VERSION
   {
      std::string stringvalues = MUMPS_VERSION;
      std::istringstream iss(stringvalues);
      char dot;
      iss >> major_version >> dot >> minor_version >> dot >> subminor_version;
   }
#endif
   // get procid
   int err = MPI_Comm_rank(univ, &procid);

   if (err != MPI_SUCCESS)
   {
      std::ostringstream oss;
      oss << " Prb with MPI_Comm_rank. Error returned : " << err << "\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   // tune icntl_size depending on version
   if (major_version > 5 || (major_version == 5 && minor_version > 1)) icntl_size = 60;
}

// [sdcz]mumps_c calling implementation
#ifdef USE_SMUMPS
template <>
void xLinearSystemSolverMumpsBase ::mumps<float>()
{
   smumps_c((SMUMPS_STRUC_C *)mumps_struct);
}
#endif
template <>
void xLinearSystemSolverMumpsBase ::mumps<double>()
{
   dmumps_c((DMUMPS_STRUC_C *)mumps_struct);
}
#ifdef USE_ZMUMPS
template <>
void xLinearSystemSolverMumpsBase ::mumps<mumps_double_complex>()
{
   zmumps_c((ZMUMPS_STRUC_C *)mumps_struct);
}

template <>
void xLinearSystemSolverMumpsBase ::mumps<std::complex<double>>()
{
   mumps<mumps_double_complex>();
}
#endif
#ifdef USE_CMUMPS
template <>
void xLinearSystemSolverMumpsBase ::mumps<mumps_complex>()
{
   cmumps_c((CMUMPS_STRUC_C *)mumps_struct);
}

template <>
void xLinearSystemSolverMumpsBase ::mumps<std::complex<float>>()
{
   mumps<mumps_complex>();
}
#endif

// allocateMemory implementation
template <>
void xLinearSystemSolverMumpsBase ::allocateMemoryStructMumps<float>()
{
   allocateMemoryStructMumpsPointer<SMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::allocateMemoryStructMumps<double>()
{
   allocateMemoryStructMumpsPointer<DMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::allocateMemoryStructMumps<mumps_complex>()
{
   allocateMemoryStructMumpsPointer<CMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::allocateMemoryStructMumps<mumps_double_complex>()
{
   allocateMemoryStructMumpsPointer<ZMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::allocateMemoryStructMumps<std::complex<double>>()
{
   allocateMemoryStructMumps<mumps_double_complex>();
}

// deallocateMemory implementation
template <>
void xLinearSystemSolverMumpsBase ::deallocateMemoryStructMumps<float>()
{
   deallocateMemoryStructMumpsPointer<SMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::deallocateMemoryStructMumps<double>()
{
   deallocateMemoryStructMumpsPointer<DMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::deallocateMemoryStructMumps<mumps_complex>()
{
   deallocateMemoryStructMumpsPointer<CMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::deallocateMemoryStructMumps<mumps_double_complex>()
{
   deallocateMemoryStructMumpsPointer<ZMUMPS_STRUC_C>();
}
template <>
void xLinearSystemSolverMumpsBase ::deallocateMemoryStructMumps<std::complex<double>>()
{
   deallocateMemoryStructMumps<mumps_double_complex>();
}

// setting default parameter
template <>
void xLinearSystemSolverMumpsBase ::setDefaultParameter<float>()
{
   setDefaultParameterFloat();
   return;
}

template <>
void xLinearSystemSolverMumpsBase ::setDefaultParameter<double>()
{
   setDefaultParameterDouble();
   return;
}

template <>
void xLinearSystemSolverMumpsBase ::setDefaultParameter<mumps_complex>()
{
   setDefaultParameterFloat();
   return;
}

template <>
void xLinearSystemSolverMumpsBase ::setDefaultParameter<mumps_double_complex>()
{
   setDefaultParameterDouble();
   return;
}

template <>
void xLinearSystemSolverMumpsBase ::setDefaultParameter<std::complex<double>>()
{
   setDefaultParameter<mumps_double_complex>();
   return;
}

void xLinearSystemSolverMumpsBase ::setDefaultParameterFloat()
{
   for (int i = 0; i < icntl_size; i++)
   {
      if (icntl[i] != XLINEARSOLVERMUMPS_EMPTY_CONTROL)
      {
         id_icntl[i] = icntl[i];
      }
   }
   for (int i = 0; i < 15; i++)
   {
      if (cntl[i])
      {
         ((float *)id_cntl)[i] = (float)rcntl[i];
      }
   }
   // if ( strcmp(write_problem,"NAME_NOT_INITIALIZED") ) ;
   strncpy(id_write_problem, write_problem, 256);
}
void xLinearSystemSolverMumpsBase ::setDefaultParameterDouble()
{
   for (int i = 0; i < icntl_size; i++)
   {
      if (icntl[i] != XLINEARSOLVERMUMPS_EMPTY_CONTROL)
      {
         id_icntl[i] = icntl[i];
      }
   }
   for (int i = 0; i < 15; i++)
   {
      if (cntl[i])
      {
         ((double *)id_cntl)[i] = rcntl[i];
      }
   }
   // if ( strcmp(write_problem,"NAME_NOT_INITIALIZED") ) ;
   strncpy(id_write_problem, write_problem, 256);
}

void xLinearSystemSolverMumpsBase ::reInitParSymCom()
{
   (*id_par) = par;
   (*id_sym) = sym;
   (*id_comm_fortran) = comm_fortran;
}

void xLinearSystemSolverMumpsBase ::reInitParSymCom(xLinearSystemSolverMumpsBase *p)
{
   p->reInitParSymCom();
   return;
}

int *&xLinearSystemSolverMumpsBase ::getJobPtr() { return id_job; }
int *&xLinearSystemSolverMumpsBase ::getJobPtr(xLinearSystemSolverMumpsBase *p) { return p->getJobPtr(); }

void xLinearSystemSolverMumpsBase ::setSym(xLinearSystemSolverMumpsBase *p, int sym_)
{
   p->sym = sym_;
   return;
}

template <>
void xLinearSystemSolverMumpsBase ::setParam(int param_tab, int param_index, const char *param_value)
{
   // depending on param_tab
   switch (param_tab)
   {
      case WRITE_PROBLEM:
      {
         // param_index is useless here
         // P must by definition be a char *
         strncpy(id_write_problem, param_value, 256);
         strncpy(write_problem, param_value, 256);
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
void xLinearSystemSolverMumpsBase ::setParam(int param_tab, int param_index, int param_value)
{
   int n;

   // depending on param_tab
   switch (param_tab)
   {
      case ICNTL:
      {
         // shift indexation as param_index is consistant with mumps doc (and fortran)
         n = param_index - 1;
         switch (n)
         {
            // rhs nature
            case 20:
            {
               if (param_value)
               {
                  std::ostringstream oss;
                  if (major_version > 4 && minor_version > 2)
                  {
                     if (param_value == 10 || param_value == 11)
                        oss << "ICNTL(20)=10 or 11 (distributed rhs) with this class will be available in future\n";
                     else
                        oss << "ICNTL(20) with this class can only be 0 (dense rhs) and in future 10,11 (distributed rhs)\n";
                  }
                  else
                     oss << "ICNTL(20) with this class can only be 0 (dense rhs)\n";
                  throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
               }
               break;
            }
            // sol nature
            case 21:
            {
               if (param_value)
               {
                  std::ostringstream oss;
                  oss << "ICNTL(21) with this class can only be 0 (centralized solution)\n";
                  throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
               }
               break;
            }
         }
         icntl[n] = id_icntl[n] = param_value;
         break;
      }
      case PAR:
      {
         (*id_par) = param_value;
         break;
      }
      case INTER:
      {
         switch (param_index)
         {
            case RESETMXNBTHREAD:
            {
#ifdef _OPENMP
               mx_nb_threads = (param_value > 0) ? param_value : omp_get_max_threads();
#else
               if (verbose) cout << "Using RESETMXNBTHREAD does no mean anything in non open MP context" << endl;
               mx_nb_threads = param_value;
#endif
               break;
            }
            case NBTHREAD:
            {
#ifndef _OPENMP
               if (verbose) cout << "Using NBTHREAD does no mean anything in non open MP context" << endl;
#endif
               nb_threads = param_value;
               break;
            }
            case VERBOSE:
            {
               verbose = param_value;
               break;
            }
            default:
            {
               std::ostringstream oss;
               oss << "Invocation of set_param use wrong param_inter ; see xLinearSolverMumps.h enum_param_inter to use a good "
                      "one  \n";
               throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
            }
         }
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
void xLinearSystemSolverMumpsBase ::printNullPivot()
{
   // null pivot
   const int ictnl24 = 24 - 1;
   if (!procid && (icntl[ictnl24] != XLINEARSOLVERMUMPS_EMPTY_CONTROL) && (icntl[ictnl24] != 0))
   {
      int nb_null_pivot;
      int *null_pivot_list = *reinterpret_cast<int **>(id_piv_null_list);
      getInternalData(INFOG, 28, nb_null_pivot);
      if (nb_null_pivot)
      {
         cout << "Null pivot detected:";
         for (int i = 0; i < nb_null_pivot; ++i)
         {
            cout << " " << null_pivot_list[i];
         }
         cout << endl;
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  xLinearSystemSolverMumpsMasterAndSlavesBase class implementation ////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// initilaisation of static members
int xLinearSystemSolverMumpsMasterAndSlavesBase ::nb_interface_instance = 0;
std::vector<xLinearSystemSolverMumpsMasterAndSlavesBase *> xLinearSystemSolverMumpsMasterAndSlavesBase::tab_interface_instance;
int xLinearSystemSolverMumpsMasterAndSlavesBase ::current_interface_instance = -1;

xLinearSystemSolverMumpsMasterAndSlavesBase ::xLinearSystemSolverMumpsMasterAndSlavesBase(MPI_Comm communicator)
    : xLinearSystemSolverMumpsBase(communicator), num_interface_instance(0)
{
   // register this object in table
   // To have a consistant structure betwen master and slave for all instances, all process
   // have to do it
   //
   // this is quite a sophistiquate implementation as it's almost imposible to delete a interface
   // instance and initiate a new one with the actual implementation of slaveAction.
   // maby in the future ....
   const int n = tab_interface_instance.size();
   while (num_interface_instance < n && tab_interface_instance[num_interface_instance] != nullptr) ++num_interface_instance;
   if (num_interface_instance < n)
      tab_interface_instance[num_interface_instance] = this;
   else
      tab_interface_instance.push_back(this);

   ++nb_interface_instance;
}

void xLinearSystemSolverMumpsMasterAndSlavesBase ::suspendSlavesActions()
{
   // local
   int err;

   // set the special id.job to force slave to get out of their instance loop
   (*id_job) = -3;

   // send special job action to performe
   err = MPI_Bcast((void *)(id_job), 1, MPI_INT, 0, univ);
   if (err != MPI_SUCCESS)
   {
      std::ostringstream oss;
      oss << " Prb with MPI_Bcast while sending special job number . Error returned : " << err << "\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}
void xLinearSystemSolverMumpsMasterAndSlavesBase ::initSlavesActions()
{
   // local
   int err;

   // send interface instance number to start the infinit action loop of this interface instance
   err = MPI_Bcast((void *)&(num_interface_instance), 1, MPI_INT, 0, univ);
   if (err != MPI_SUCCESS)
   {
      std::ostringstream oss;
      oss << " Prb with MPI_Bcast while sending slave number . Error returned : " << err << "\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }
}

void xLinearSystemSolverMumpsMasterAndSlavesBase ::checkInstance()
{
   // if interface instance is the current one, slave are in correct state and
   // there is nothing to do.
   // If not, change interface instance of slaves
   if (num_interface_instance != current_interface_instance)
   {
      // if allready used by a other interface instance
      if (current_interface_instance > -1) tab_interface_instance[current_interface_instance]->suspendSlavesActions();
      initSlavesActions();
      current_interface_instance = num_interface_instance;
   }
}

// give index of the first available interface pointeur in
int xLinearSystemSolverMumpsMasterAndSlavesBase ::getFirstAvailableInterfaceInstance()
{
   int i = 0;
   const int n = tab_interface_instance.size();
   while (i < n && tab_interface_instance[i] == nullptr) ++i;
   if (i == n)
   {
      std::ostringstream oss;
      oss << " Buuuuug !?!? enable to catch a available pointer in  tab_interface_instance of size=" << n
          << " with nb_interface_instance=" << nb_interface_instance << "\n";
      throw xLinearSystemSolverMumpsException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
   }

   return (i);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLinearSystemSolverMumpsException class implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xLinearSystemSolverMumpsException::xLinearSystemSolverMumpsException(std::string info, std::string file, int Line,
                                                                     std::string date, std::string time)
{
   std::ostringstream oss;
   oss << "In file " << file << " line " << Line << " compiled " << date << " at " << time << std::endl;
   oss << "xLinearSystemSolverMumpsException : " << info << std::endl;
   msg = oss.str();
}
/// general exception object : destructor
xLinearSystemSolverMumpsException ::~xLinearSystemSolverMumpsException() throw() = default;

/// mandatory what method
const char *xLinearSystemSolverMumpsException::what() const throw() { return this->msg.c_str(); }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xLinearSystemSolverMumpsException class implementation ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}  // end namespace
