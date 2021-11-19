/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xValueLinearCombination.h"
#include "xValue.h"
#include "xStateOfValue.h"

namespace xfem {


template < >
xValue < double > * xCloneValue(const xValueLinearCombination<double> * val, const std::map < xValue < double > *,xValue < double > * > &coresp)
{
    // do copies
    xValueLinearCombination<double>::coeffs_t coef (val->coeffs);
    xValueLinearCombination<double>::xvalues_t vals (val->values);

    // remove last value which is always a fixed value created on the fly by the constructor
    vals.pop_back();

    // loop on value to get new ones via coresp
    xValueLinearCombination<double>::xvalues_t nvals;
    nvals.reserve(vals.size());
    auto ite = coresp.end();
    for (auto v : vals)
    {
        auto itf = coresp.find(v);
        if (itf != ite)
            nvals.push_back(itf->second);
        else
            break;
    }

    // if some values in coresp are missing we can't clone val
    if (nvals.size() != vals.size()) return nullptr;

    // get last coeff
    double coeff_last = coef.back();
    coef.pop_back();

    // create new value
    xValue < double >* new_val = new xValueLinearCombination<double>(coef,nvals,coeff_last);

    // state is create on demande by getState. Forced if original is set
    const xStateOfValue *state = val->getState();
    if (state)
        state = new_val->getState();
    return new_val;
}



}
