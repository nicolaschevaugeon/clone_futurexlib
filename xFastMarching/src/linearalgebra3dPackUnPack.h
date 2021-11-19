/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _linearalgebra3dPackUnPack_
#define _linearalgebra3dPackUnPack_

#include "linearalgebra3d.h"

// xtool
#include "xDataExchanger.h"
#include "xMPIDataType.h"




  template <typename T>
    inline void pack(const xfastmarching::vector2d<T>& vect,xtool::xMpiInputBuffer &buff)
    {
      buff.pack(&vect[0],2,xtool::xMPIDataType < T >());
    }
  template <typename T>
    inline void unPack(xfastmarching::vector2d<T>& vect,const xtool::xMpiOutputBuffer &buff)
    {
      buff.unPack(&vect[0],2,xtool::xMPIDataType < T >());
    }
  template <typename T>
    inline void pack(const xfastmarching::vector3d<T>& vect,xtool::xMpiInputBuffer &buff)
    {
      buff.pack(&vect[0],3,xtool::xMPIDataType < T >());
    }
  template <typename T>
    inline void unPack(xfastmarching::vector3d<T>& vect,const xtool::xMpiOutputBuffer &buff)
    {
      buff.unPack(&vect[0],3,xtool::xMPIDataType < T >());
    }


#endif
