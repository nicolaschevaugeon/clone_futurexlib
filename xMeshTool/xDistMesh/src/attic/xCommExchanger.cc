/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include "xCommExchanger.h"

// std


using namespace std;


namespace distmesh
{
xCommExchangerException::xCommExchangerException(std::string info,std::string file,int Line,std::string date,std::string time)
{
    std::ostringstream oss;
    oss << "In file "<< file << " line " << Line << " compiled "<<date<<" at "<<time<<std::endl;
    oss << "xCommExchangerException : "<< info << std::endl;
    msg = oss.str();
}
/// general exception object : destructor
xCommExchangerException :: ~xCommExchangerException() throw( ) {}

/// mandatory what method
const char * xCommExchangerException::what() const throw( )
{
    return this->msg.c_str();
}

}         // end namspace
