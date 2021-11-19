/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xStateOfValue.h"
#include "xValue.h"

namespace xfem
{

std::ostream& xStateOfValueNone::print(std::ostream& o) const { o << "xStateOfValueNone" << std::endl; return o; }

template <>
void xCloneState(const xStateOfValueFixed<double> * state, xValue<double> *v )
{
    v->setState(new xStateOfValueFixed<double>(v));
}

template <>
void xCloneState(const xStateOfValueDof * state, xValue<double> *v )
{
    v->setState(new xStateOfValueDof(state->Numdof));
}
template <>
void xCloneState(const xStateOfValueNone * state, xValue<double> *v )
{
    v->setState(new xStateOfValueNone());
}

} // end of namespace
