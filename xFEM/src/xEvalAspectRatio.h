/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef  _XEVALASPECTRATIO_H
#define  _XEVALASPECTRATIO_H

#include "xAspectRatioAOMD.h"
#include "xEval.h"

namespace AOMD
{
    class mEntity;
}

namespace xfem
{

//! evaluator to output aspect ratio
template <typename T, xtool::enum_xAspectRatio ART = xtool::enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >
class xEvalAspectRatio :public xEval < double >
{
    public:
        void operator() (const xGeomElem*  geo_appro, const xGeomElem* geo_integ, result_type & r)  const override
        {

            AOMD::mEntity* e = geo_appro->getEntity();
            assert(dynamic_cast < T * >( e ));
            r=xtool::xAspectRatio<T, ART >::aspectRatio(static_cast < T * >(e));
        }
};


} // end of namespace

#endif






