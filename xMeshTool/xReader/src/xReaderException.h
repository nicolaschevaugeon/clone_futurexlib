/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */

#ifndef _READEREXP_XMESHTOOL
#define _READEREXP_XMESHTOOL
#include<sstream>

namespace xmeshtool
{
class xReaderException : public std::exception
{
    public:
        xReaderException(std::string info,std::string file,int Line,std::string date,std::string time)
        {
            std::ostringstream oss;
            oss << "In file "<< file << " line " << Line << " compiled "<<date<<" at "<<time<<std::endl;
            oss << "xReaderException : "<< info << std::endl;
            msg = oss.str();
        }
        ~xReaderException() throw( )override = default;
        const char * what() const throw( )override
        {
            return this->msg.c_str();
        }
    private:
        std::string msg;
};

} // end namespace

#endif

