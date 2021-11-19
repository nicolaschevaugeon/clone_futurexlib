/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#ifndef _XVALUE_LINEAR_COMBINATION_H_
#define _XVALUE_LINEAR_COMBINATION_H_

#include "xValue.h"
#include "xStateOfValue.h"



namespace xfem {


/// Base class concept of linear combination
template < typename T >
class xValueLinearCombinationBase : public xValue < T >
{
    public:
        typedef std::vector < T > coeffs_t;
        typedef xValue < T >  xvalue_t;
        typedef std::vector < xvalue_t * > xvalues_t;

        void   setVal(T in) override;
        void setState(xStateOfValue* s) override;
        const xStateOfValue * getState()  const override;
        xStateOfValue * getState() override;
        std::ostream & printVal(std::ostream & o) const override;

        size_t getSize() const;
        xvalue_t * getValue(size_t i) const;
        virtual T getCoeff(size_t i) const = 0;
    protected:

        xValueLinearCombinationBase(const coeffs_t & c, const xvalues_t & v, T coeff_last);
        coeffs_t coeffs;
        xvalues_t values;
        std::unique_ptr < xSingleValue < T > > v_last;
        using xValue < T >::state;
};

/*!
 * A value is tied to other values and a constant
 * v = sum_i c_i * v_i + b
 * c_i are the coefficients, v_i are the other value
 * b is a constant which by default is zero.
 *
 * in the implementation, b is an extra coefficient c_last
 * multiplying a fixed value v_last of value xtool::xDataType<T>::one()
 */

template < typename T >
class xValueLinearCombination : public xValueLinearCombinationBase < T >
{
    template < typename V, typename B >
    friend B * xCloneValue(const V * val, const std::map < B *,B * > &coresp);

    using xValueLinearCombinationBase < T >::values;
    using xValueLinearCombinationBase < T >::coeffs;

    public:
        typedef typename xValueLinearCombinationBase < T >::coeffs_t coeffs_t;
        typedef typename xValueLinearCombinationBase < T >::xvalues_t xvalues_t;
        typedef typename xValueLinearCombinationBase < T >::xvalue_t xvalue_t;

        xValueLinearCombination(const coeffs_t & c, const xvalues_t & v, T coeff_last = xtool::xDataType < T >::zero());
        xValueLinearCombination(T c, xvalue_t *v, T coeff_last = xtool::xDataType < T >::zero());
        T getVal() const override;
        T getCoeff(size_t i) const override;
        std::ostream & printVal(std::ostream & o) const override;
};

/*!
 * A value is tied to other values with help of associated values that
 * play the role of parametric coefficients. The parametric values and the
 * set of fixed coefficient give the set of coefficients that appears in the
 * final cinematic relation.
 * v = sum_i c_i v_i w_i sum_k c_k v_k + b
 * c_i,c_k are the fixed coefficients,
 * w_i are the parametric values considered at evaluation
 * v_i,v_k are the other value on which v depend linearly
 * b is a constant which by default is zero.
 *
 * Said differently v values is tied to other values via a parametric field
 *
 * in the implementation, b is an extra coefficient c_last
 * multiplying a fixed value v_last of value xtool::xDataType<T>::one()
 */

template < typename T >
class xValueLinearCombinationParam : public xValueLinearCombinationBase < T >
{
    using xValueLinearCombinationBase < T >::values;
    using xValueLinearCombinationBase < T >::coeffs;

    public:
        typedef typename xValueLinearCombinationBase < T >::coeffs_t coeffs_t;
        typedef typename xValueLinearCombinationBase < T >::xvalues_t xvalues_t;
        typedef typename xValueLinearCombinationBase < T >::xvalue_t xvalue_t;

        xValueLinearCombinationParam(const coeffs_t & c, const xvalues_t & w, const xvalues_t & v, T coeff_last = xtool::xDataType < T >::zero());
        T getVal() const override;
        T getCoeff(size_t i) const override;
        std::ostream & printVal(std::ostream & o) const override;

    private:
        xvalues_t params;
};

//Clone for double-valued specialization
template < >
xValue < double > * xCloneValue(const xValueLinearCombination < double > * val, const std::map < xValue < double > *,xValue < double > * > &coresp);

} //End namespace

#include "xValueLinearCombination_imp.h"

#endif
