/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _STATE_OF_VALUE_IMP_H
#define _STATE_OF_VALUE_IMP_H

#include <typeinfo>

namespace xfem
{

//-------------------------------------------------------------------
template <typename S, typename T>
void xCloneState(const S *state, xValue<T> * v)
{
    if (typeid(*state)==typeid(xStateOfValueFixed<T>))
        xCloneState(static_cast<const xStateOfValueFixed<T> *>(state),static_cast<xValue<double> *>(v));
    else if (typeid(*state)==typeid(xStateOfValueDof))
        xCloneState(static_cast<const xStateOfValueDof *>(state),static_cast<xValue<double> *>(v));
    else if (typeid(*state)==typeid(xStateOfValueNone))
    {
        if (typeid(T)==typeid(double))
            xCloneState(static_cast<const xStateOfValueNone *>(state),static_cast<xValue<double> *>(v));
        else
        {
            std::cout<<"state of type "<<typeid(*state).name()<<"for xValue of type "<<typeid(T).name()<<" is not yet clonable !"<<std::endl;
            throw -67845;
        }
    }
    else
    {
        std::cout<<"state of type "<<typeid(*state).name()<<" is not yet clonable !"<<std::endl;
        throw -12345;
    }

}


template <>
void xCloneState(const xStateOfValueFixed<double> * state, xValue<double> *v );
template <>
void xCloneState(const xStateOfValueDof   * state, xValue<double> *v );
template <>
void xCloneState(const xStateOfValueNone  * state, xValue<double> *v );


//-------------------------------------------------------------------

template<typename VT>
xStateOfValueFixed<VT>::xStateOfValueFixed(const xValue<VT>* v) :xStateOfValue(FIXED), val(v) {}

template<typename VT>
VT xStateOfValueFixed<VT>::getVal() const {return val->getVal();}

template<typename VT>
std::ostream& xStateOfValueFixed<VT>::print(std::ostream& o) const { o << "xStateOfValueFixed" << std::endl; return o; }


template<typename T>
xStateOfValueLinearCombination<T>::xStateOfValueLinearCombination(const xValueLinearCombinationBase<T>* l) :xStateOfValue(LINEARCOMBINATION), lin(l)
{
  const bool debug = xdebug_flag;
  if (debug) std::cout << "creating state average with "
                       << lin->getSize() << " values " << std::endl;
}

template<typename T>
size_t xStateOfValueLinearCombination<T>::size()                   const {return lin->getSize(); }

template<typename T>
T          xStateOfValueLinearCombination<T>::coeff(int i)         const {return lin->getCoeff(i); }

template<typename T>
xStateOfValue* xStateOfValueLinearCombination<T>::state(int i)     const {return lin->getValue(i)->getState(); }

template<typename T>
std::ostream& xStateOfValueLinearCombination<T>::print(std::ostream& o) const
        { o << "xStateOfValueLinearCombination "<< std::endl; return o; }

//-------------------------------------------------------------------
}

#endif
