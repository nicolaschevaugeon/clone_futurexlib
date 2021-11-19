/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef ___XMPIENV__H
#define ___XMPIENV__H
#include <exception>
#include <string>
#include "mpi.h"

// nota :
//   implementation is not thread safe as at least, finalize, if called by many thread may suffer of race
//   condition durring the finalized=true; operation

namespace xtool
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xMPIEnv class /////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class xMPIEnv
{
    public:
        // public members ////////////////
        enum enum_init_type { NO_THREAD,
                              THREAD_SINGLE,    // Only one thread will execute.
                              THREAD_FUNNELED,  // The process may be multi-threaded, but only the main thread will make MPI calls (all MPI calls are
                                                // funneled to the main thread).
                              THREAD_SERIALIZED, // The process may be multi-threaded, and multiple threads may make MPI calls, but only one at a time:
                                                 //      MPI calls are not made concurrently from two distinct threads (all MPI calls are serialized)
                              THREAD_MULTIPLE   // Multiple threads may call MPI, with no restrictions
        };

        // public destructor /////////////
        /// public destructor
        //! Note:
        //! If used it will finalize parallel environment and
        //! any attempt to reinitialize an environment will fail and give a null pointer
        //! Memory cleaning will be done after
        //!
        //! It is made public for this last point : cleaning at exit without
        //! asking explicitly to do it. For some case it might work. It's here for the lazy one who forget to clean things
        //! But user is highly encouraged to use finalize instead.
        //! Finalize do correctly close mpi parallel environment at the right moment at the right place.
        //! If not it's this destructor that will do it but we won't have any
        //! assurance that this will be done at the correct moment as it will
        //! depend on destructor calling sequence. For example  other object might be destructed after MPI_FINALIZE but
        //! need a communicator to do their own cleaning in between procs
        //!
        ~xMPIEnv();

        // public methodes ///////////////
        /// initiate mpi parallel environement and return 0 if OK or -1 if a finalize have already been done
        static int init(int &argc, char * * &argv,char ini_type = NO_THREAD);

        /// finalize mpi parallel environement
        static int finalize();

    private:
        // private constructor ///////////////
        xMPIEnv(int &argc, char * * &argv,char ini_type);
        // private members ////////////////
        static xMPIEnv* unique_instance;
        static bool finalized;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xMPIEnv class /////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xMPIEnvException class ////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// interface derived class of standart exception for xMPIEnv
class xMPIEnvException : public std::exception
{
    public:
        xMPIEnvException(std::string,std::string,int,std::string,std::string);
        ~xMPIEnvException() throw( ) override;
        const char * what() const throw( ) override;

    private:
        std::string msg;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xMPIEnvException  class ///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif
