/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef _XVECTORSCALARPROD_H_
#define _XVECTORSCALARPROD_H_


#include "xVector.h"
#include "xGenericOperations.h"


namespace xtensor{

using xVectorScalarProd =  typename xtool::xMult< xVector<double>, xVector<double>, double>;
using xVectorScalarProdf =  typename xtool::xMult< xVector<float>, xVectorf<float>, float>;

}


#endif
