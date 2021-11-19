/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */
#include <iostream>
#include <fstream>
#include <cstdio>
#include <typeinfo>
#include "parmetis.h"

/* NOTE: parmetis.h imply use of some mpi library (mpi.h).
 * But all what is needed by this program is only the 
 * header not the MPI runtime library. No parallelism
 * is needed. To simplify execution of this program by
 * cmake try_run on some plateform this program is
 * a pure sequential aplication (no MPI context introduced)
 */

int main(int argc, char *argv[])
{
    // the argument passed by cmake is the pass to source
    std::string path(argv[1]);
    std::string pathn(argv[1]);
    path += "/ParMetisInterface_types.h";
    pathn += "/.ParMetisInterface_types.h";
    std::ofstream ofsn (pathn.c_str(), std::ofstream::out);
    if (!( ofsn.is_open()))
    {
        std::cout << "Can't create temporary file " << pathn << std::endl;
        return -1;
    }

    ofsn<<"/*"<<std::endl<<" Don't change by hand ! Don't commit !"<<std::endl<<" This file is generated during cmake build"<<std::endl<<" */"<<std::endl;

#if ( PARMETIS_MAJOR_VERSION > 3 )
    if (typeid( idx_t ) == typeid( int ))
        ofsn<<"typedef int parmetis_indx_t;"<<std::endl;
    else if (typeid( idx_t ) == typeid( int32_t ))
        ofsn<<"typedef int32_t parmetis_indx_t;"<<std::endl;
    else if (typeid( idx_t ) == typeid( int64_t ))
        ofsn<<"typedef int64_t parmetis_indx_t;"<<std::endl;

    if (typeid( real_t ) == typeid( float ))
        ofsn<<"typedef float parmetis_real_t;"<<std::endl;
    else if (typeid( real_t ) == typeid( double ))
        ofsn<<"typedef double parmetis_real_t;"<<std::endl;
#else
    if (typeid( idxtype ) == typeid( int ))
        ofsn<<"typedef int parmetis_indx_t;"<<std::endl;
    if (typeid( idxtype ) == typeid( short ))
        ofsn<<"typedef short parmetis_indx_t;"<<std::endl;

    ofsn<<"typedef float parmetis_real_t;"<<std::endl;
#endif
    ofsn.close();
    std::cout << std::endl;
#if ( PARMETIS_MAJOR_VERSION > 3 )
    std::cout << "Identified ParMetis major version greater then 3 "<< std::endl;
#else
    std::cout << "Identified ParMetis major version less or equal to 3 "<< std::endl;
#endif
    std::cout << "temporary file created " << pathn << std::endl;
    // check now that it is new
    std::ifstream ifs (path.c_str(), std::ifstream::in);
    bool move = true;
    // first case : a file already present
    if (ifs.is_open())
    {
        std::cout << "old file " << path << " exist"<<std::endl;
        std::ifstream ifsn (pathn.c_str(), std::ifstream::in);
        if (!( ifsn.is_open()))
        {
            std::cout << "Tempory file " << pathn <<  "not available ???"<< std::endl;
            return -1;
        }
        std::string lo,ln;
        getline(ifs,lo);
        getline(ifsn,ln);
        bool same = true;
        do
        {
            same = ( lo == ln );
            getline(ifs,lo);
            getline(ifsn,ln);
        } while (!ifs.eof() && same);
        ifs.close();
        ifsn.close();
        if (same)
        {
            std::cout << "old file an temporary file are the same" <<std::endl;
            if ( std::remove( pathn.c_str() ) != 0 )
            {
                std::string msg("Tempory file "+pathn+"cannot be removed");
                std::perror(msg.c_str());
                return -1;
            }
            std::cout << "temporary file removed" <<std::endl;
            move = false;
        }
        else
            std::cout << "old file an temporary file are not the same" <<std::endl;
    }
    // second case: no file present => no doubt it is new, move set to true by default

    // if new file have to be moved
    if (move)
    {
        std::cout << "temporary file renamed into " << path << std::endl;
        if ( std::rename( pathn.c_str(), path.c_str() ) != 0 )
        {
            std::string msg("Tempory file "+pathn+"cannot be renamed in "+path);
            std::perror(msg.c_str());
        }
    }
    return 0;
}
