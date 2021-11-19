/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/
#ifndef  _XASPECTRATIOAOMD_H
#define  _XASPECTRATIOAOMD_H

#include "xAspectRatio.h"

namespace AOMD
{
  class mTet;
  class mEntity;
}

namespace xtool
{
 
    //! xTool xAspectRatio class specialisation for mTet
    template < >
      class xAspectRatio < AOMD::mTet,enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >
      {
      public:
        static double aspectRatio (AOMD::mEntity * entity);
        static double aspectRatio (AOMD::mTet * entity);
        static double aspectRatio (AOMD::mEntity *v1, AOMD::mEntity *v2,AOMD::mEntity *v3, AOMD::mEntity *v4);
      };

  
} // end of namespace

#endif






