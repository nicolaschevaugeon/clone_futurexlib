/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xGenericSparseMatrix.h"


namespace xlinalg
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xGenericSparseMatrixException class implementation /////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xGenericSparseMatrixException::xGenericSparseMatrixException(std::string info,std::string file,int Line,std::string date,std::string time){
        std::ostringstream oss;
        oss << "In file "<< file << " line " << Line << " compiled "<<date<<" at "<<time<<std::endl;
        oss << "xGenericSparseMatrixException : "<< info << std::endl;
        msg = oss.str();
}
/// general exception object : destructor
  xGenericSparseMatrixException :: ~xGenericSparseMatrixException() throw() = default;

/// mandatory what method
const char * xGenericSparseMatrixException::what() const throw()
{
   return this->msg.c_str();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xGenericSparseMatrixException class implementation /////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // end of namespace