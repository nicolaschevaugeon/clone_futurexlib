/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include <sstream>
#include "xMPIEnv.h"

namespace xtool
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xMPIEnv class implementation //////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////// Statics ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xMPIEnv * xMPIEnv::unique_instance = nullptr;
bool xMPIEnv::finalized = false;

/////////////////////////////////////// Constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xMPIEnv::xMPIEnv(int &argc, char * * &argv,char ini_type)
{
    bool do_init_thread = true;
    int required = 0,provided,ret;
    std::string name;

    // check MPI version
#ifndef FORCE_MPI2
    int version, subversion;
    MPI_Get_version( &version, &subversion);
    if (version < 3)
        throw xMPIEnvException("Exlibris library is conforming to MPI norm version 3.0 and above.\n You are not using a MPI implementation following this version. Choose an other MPI implementation.\n If really you have no other choice you can, at your own risk, force usage of version 2.X by recompiling every thing with FORCE_MPI2",__FILE__,__LINE__,__DATE__,__TIME__);

#endif

// supersede user choice
#ifdef WITH_DMUMPS
    ini_type = xMPIEnv::NO_THREAD;
#endif
// PTscotch requirement supersede previous
#ifdef WITH_SCOTCH6
    ini_type = xMPIEnv::THREAD_SINGLE;
#endif
// global directive requirement supersede previous
#ifdef WITH_MPI_THREAD_SINGLE
    ini_type = xMPIEnv::THREAD_SINGLE;
#endif
// global directive requirement supersede previous
#ifdef WITH_MPI_THREAD_FUNNELED
    ini_type = xMPIEnv::THREAD_FUNNELED;
#endif
// global directive requirement supersede previous
#ifdef WITH_MPI_THREAD_SERIALIZED
    ini_type = xMPIEnv::THREAD_SERIALIZED;
#endif
// pastix requirement supersede previous
#ifdef WITH_DPASTIX
    ini_type = xMPIEnv::THREAD_MULTIPLE;
#endif
// global directive requirement supersede previous
#ifdef WITH_MPI_THREAD_MULTIPLE
    ini_type = xMPIEnv::THREAD_MULTIPLE;
#endif

    // define what kind of environement is asked
    switch (ini_type)
    {
        case THREAD_SINGLE :
         {
             required = MPI_THREAD_SINGLE;
             name = "MPI_THREAD_SINGLE";
             break;
         }
        case THREAD_FUNNELED :
         {
             required = MPI_THREAD_FUNNELED;
             name = "MPI_THREAD_FUNNELED";
             break;
         }

        case THREAD_SERIALIZED :
         {
             required = MPI_THREAD_SERIALIZED;
             name = "MPI_THREAD_SERIALIZED";
             break;
         }
        case THREAD_MULTIPLE :
         {
             required = MPI_THREAD_MULTIPLE;
             name = "MPI_THREAD_MULTIPLE";
             break;
         }
        case NO_THREAD :
        default :
         {
             do_init_thread = false;
         }
    }

    // initialize paralle environement
    if (do_init_thread)
    {
        ret = MPI_Init_thread(&argc, &argv, required, &provided);
    }
    else
    {
        ret = MPI_Init(&argc,&argv);
    }

    // check returned value
    if (ret != MPI_SUCCESS)
    {
        std::ostringstream oss;
        oss << "Unsuccessful initialization of parallel environment !\n";
        if (do_init_thread)
            oss << "MPI_Init_thread with "<<name<<" return with error ="<<ret<<std::endl;
        else
            oss << "MPI_Init return with error ="<<ret<<std::endl;
        throw xMPIEnvException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);

    }

    // check returned value
    if (do_init_thread && required != provided)
    {
        std::ostringstream oss;
        oss << "ERROR : asked for "<<required<<" ("<<name<<") an only get "<<provided<<std::endl;
        oss<<"Note: Used MPI implementation may :"<<std::endl;
        oss<<"     -not provide such thread initiation"<<std::endl;
        oss<<"     -provide such thread initiation but was not compiled with appropriate option"<<std::endl;
        oss<<"     -provide such thread initiation but your application was not launched with appropriate option"<<std::endl;
        oss<<std::endl;
        oss<<"For example MVAPICH2-2.0 need MV2_ENABLE_AFFINITY=0 (env) to authorize THREAD_MULTIPLE otherwise it stays in MPI_THREAD_SERIALIZED"<<std::endl;
        throw xMPIEnvException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }


};
/////////////////////////////////////// End constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xMPIEnv::~xMPIEnv()
{
    // if not already done by the user
    if (!finalized)
    {
        std::cout<<"WARNING : You did not use finalize before the end of your application ! It may be a problem !"<<std::endl;
        finalize();
    }
}
/////////////////////////////////////// End Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////////
int xMPIEnv::init(int &argc, char ** &argv,char ini_type)
{
    if (finalized) return -1;
    if (!unique_instance)
    {
        unique_instance = new xMPIEnv(argc,argv,ini_type);
    }
    return 0;

}

int xMPIEnv::finalize()
{
    finalized = true;

    // finalize mpi parrallel environement
    int ret = MPI_Finalize();
    if ( ret != MPI_SUCCESS)
    {
        std::ostringstream oss;
        oss << "Problem with MPI_Finalize !\n error="<<ret<<std::endl;
        throw xMPIEnvException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }
    else return ret;

}
/////////////////////////////////////// End Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xMPIEnv class implementation //////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xMPIEnvException class implementation /////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// general exception used for all xMPIEnv throw
xMPIEnvException::xMPIEnvException(std::string info,std::string file,int Line,std::string date,std::string time)
{
    std::ostringstream oss;
    oss << "In file "<< file << " line " << Line << " compiled "<<date<<" at "<<time<<std::endl;
    oss << "xMPIEnvException : "<< info << std::endl;
    msg = oss.str();
}
/// general exception object : destructor
xMPIEnvException :: ~xMPIEnvException() throw( ) = default;;

/// mandatory what method
const char * xMPIEnvException::what() const throw( )
{
    return this->msg.c_str();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xMPIEnvException class implementation /////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



} // end of namespace
