/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _XVALUE_IMP_H_
#define _XVALUE_IMP_H_

#ifndef _VALUE_H
#error Do NOT include xValue_imp.h alone
#endif

#include <typeinfo>

namespace xfem
{
//----------------------------------------------------------------------------------------------------------------
// xCloneValue specialization for xSingleValue in double aritmethic
template <>
xValue<double> *xCloneValue(const xSingleValue<double> *val, const std::map<xValue<double> *, xValue<double> *> &coresp);

}  // namespace xfem

#endif
