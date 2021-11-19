/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _XVALUE_LINEAR_COMBINATION_IMP_H_
#define _XVALUE_LINEAR_COMBINATION_IMP_H_

#ifndef _XVALUE_LINEAR_COMBINATION_H_
#error Do NOT include xValueLinearCombination_imp.h alone
#endif

namespace xfem {


//----------------------------------------------------------------------------------------------------------------
template <typename T>
xValueLinearCombinationBase<T>::xValueLinearCombinationBase(const coeffs_t& c, const xvalues_t& v, T coeff_last)
    : coeffs(c), values(v), v_last {new xSingleValue<T>}
{
    xSingleValue<T> *pv_last=v_last.get();
    v_last->setVal(xtool::xDataType<T>::one());
    v_last->setState(new xStateOfValueFixed<T>(pv_last));
    coeffs.push_back(coeff_last);
    values.push_back(pv_last);
}
template <typename T>
void xValueLinearCombinationBase<T>::setVal(T in) {}
template <typename T>
void xValueLinearCombinationBase<T>::setState(xStateOfValue* s) {}
template <typename T>
const xStateOfValue * xValueLinearCombinationBase<T>::getState()  const
{
    if (!state) { state = new xStateOfValueLinearCombination<T>(this); }
    return state;
}
template <typename T>
xStateOfValue * xValueLinearCombinationBase<T>::getState()
{
    if (!state) { state = new xStateOfValueLinearCombination<T>(this); }
    return state;
}
template <typename T>
size_t xValueLinearCombinationBase<T>::getSize()  const
{
    return  values.size();
}
template <typename T>
auto xValueLinearCombinationBase<T>::getValue(size_t i) const -> xvalue_t *
{
    assert(i<values.size());
    return values[i];
}
template <typename T>
std::ostream & xValueLinearCombinationBase<T>::printVal(std::ostream & o) const
{
    size_t n=values.size();
    o << "linear combination of the " << n << " values below " << std::endl;
    o << "the coefficients are "     << std::endl;
    for(size_t i=0;i<n;++i)
        o << this->getCoeff(i)<< std::endl;
    return o;
}
//----------------------------------------------------------------------------------------------------------------
template <typename T>
xValueLinearCombination<T>::xValueLinearCombination(const coeffs_t & c, const xvalues_t & v, T coeff_last)
    : xValueLinearCombinationBase<T>(c,v,coeff_last)
{
    const bool debug = xdebug_flag;
    if (debug) std::cout << "creating value average with " << values.size() << " values " << std::endl;
}
template <typename T>
xValueLinearCombination<T>::xValueLinearCombination(T c, xvalue_t * v, T coeff_last)
    : xValueLinearCombinationBase<T>(coeffs_t(1,c),xvalues_t(1,v),coeff_last)
{}

template <typename T>
T xValueLinearCombination<T>::getVal() const
{
    T v = xtool::xDataType<T>::zero();
    auto itc = coeffs.begin();
    for (auto val : values)
    {
        v += *itc * val->getVal();
        ++itc; 
    }
    return v;
}
template <typename T>
T xValueLinearCombination<T>::getCoeff(size_t i) const
{
    assert(i<values.size());
    return coeffs[i];
}
template <typename T>
std::ostream & xValueLinearCombination<T>::printVal(std::ostream & o) const
{
    xValueLinearCombinationBase<T>::printVal(o);
    o << "and the values are   " << std::endl;
    for (auto val : values)  val->print(o) << std::endl;
    return o;
}
//----------------------------------------------------------------------------------------------------------------
template <typename T>
xValueLinearCombinationParam<T>::xValueLinearCombinationParam(const coeffs_t & c, const xvalues_t& w, const xvalues_t & v, T coeff_last)
    : xValueLinearCombinationBase<T>(c,v,coeff_last),params(w)
{}
template <typename T>
T xValueLinearCombinationParam<T>::getVal() const
{
    T v = xtool::xDataType<double>::zero();
    auto itv = values.begin();
    auto itve = values.end();
    auto itc = coeffs.begin();
    for (auto param : params)
    { 
        v += *itc * ( *itv )->getVal()*param->getVal();
        ++itv; ++itc; 
    }
    for (; itv != itve; ++itv, ++itc) v += *itc * ( *itv )->getVal();
    return v;
}

template <typename T>
T xValueLinearCombinationParam<T>::getCoeff(size_t i) const
{
    assert(i<values.size());
    T v=coeffs[i];
    if (i<params.size())
        v*=params[i]->getVal();
    return v;
}

template <typename T>
std::ostream & xValueLinearCombinationParam<T>::printVal(std::ostream & o) const
{
    xValueLinearCombinationBase<T>::printVal(o);
    o << "the computed coefficients are "     << std::endl;
    auto itc = coeffs.begin();
    for (auto param : params)
    { 
        o<< "coef: "<<param->getVal()*(*itc)<<"= fixed coef: "<<*itc<<" * param value: ";
        param->print(o) << std::endl;
        ++itc; 
    }
    o << "and the values are   " << std::endl;
    for (auto val : values)  val->print(o) << std::endl;
    return o;
}
 
//----------------------------------------------------------------------------------------------------------------

}//End namespace


#endif
