/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
  
#include "xcValueJs.h"


namespace xcrack 
{


//we need to declare static members
int xcValueJs::mode_ = 0;
int xcValueJs::nb_modes = 0;
//we need to declare static members
int xcValueSifs::mode_ = 0;
int xcValueSifs::nb_modes = 0;
double xcValueSifs::young = 0.;
double xcValueSifs::poisson = 0.;
xcValueSifs::geom_t  xcValueSifs::geom = GEOM_3D;
//
xcValueJsAndSifs::field_t  xcValueJsAndSifs::field = JS;


}
