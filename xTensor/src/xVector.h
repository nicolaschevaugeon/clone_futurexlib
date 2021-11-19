/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

//////////////////////////////////////////////////////////////////////

#ifndef _XVECTOR_H_
#define _XVECTOR_H_

#include <cmath>
#include <complex>
#include <functional>
#include <initializer_list>
#include <iostream>

#include "xDataType.h"
#include "xPoint.h"

namespace xtensor
{
template <typename VT>
class xTensor2;

/**
     Simple class for a Vector in 3d space.
  */
template <typename VT = double>
class xVector
{
   VT pos[3];

  public:
   using value_type = VT;
   inline xVector(VT x = xtool::xDataType<VT>::zero(), VT y = xtool::xDataType<VT>::zero(), VT z = xtool::xDataType<VT>::zero())
       : pos{x, y, z}
   {
   }
   inline xVector(const xVector &in) : pos{in(0), in(1), in(2)} {}
   inline xVector &operator=(const xVector &in)
   {
      std::copy(in.pos, in.pos + 3, pos);
      return *this;
   }
   inline xVector(const xPoint &from, const xPoint &to)
       : pos{static_cast<VT>(to(0) - from(0)), static_cast<VT>(to(1) - from(1)), static_cast<VT>(to(2) - from(2))}
   {
   }

   //! return a const pointer to the begining of an array  of 3 VT containing the data of the xVector:
   //! a0 a1 a2
   inline const VT *data() const { return pos; }
   //! return a pointer to the begining of an array  of 3 VT containing the data of the xVector:
   //! a0 a1 a2
   //! modifying the data inside the array modify the vector ...
   inline VT *data() { return pos; }

   inline const VT &operator()(int i) const;
   inline VT &operator()(int i);
   inline const VT &operator[](int i) const;
   inline VT &operator[](int i);

   inline xVector<VT> &operator+=(const xVector<VT> &other);
   inline xVector<VT> &operator-=(const xVector<VT> &other);
   inline xVector<VT> &operator*=(const VT &other);
   inline xVector<VT> &operator/=(const VT &other);
   xVector<VT> &operator*=(const xTensor2<VT> &other);
   /// Normalize the vector
   inline xVector<VT> &norm();
   /// Normalize and return the Norm as it was before normalization. If initial norm is 0., return 0. and leave the vector
   /// unchanged.
   inline VT normValue();
   /// Return the Magnitude (Euclidien Norm)
   inline VT mag() const;

   inline xVector<VT> operator-() const;
   inline xVector<VT> operator-(const xVector<VT> &other) const;
   inline xVector<VT> operator+(const xVector<VT> &other) const;
   /// return the scalar product of 2 xVector.
   inline VT operator*(const xVector<VT> &other) const
   {
      return pos[0] * other.pos[0] + pos[1] * other.pos[1] + pos[2] * other.pos[2];
   }

   // Cross-type case
   //    inline float operator * (const xVector<float> &other) const
   //    { return static_cast<float>(pos[0] * other.pos[0] + pos[1] * other.pos[1] + pos[2]*other.pos[2]);}

   /// LEFT mutiply of a xVector to a xTensor2.
   xVector<VT> operator*(const xTensor2<VT> &M) const;
   /// Compute the Cross Product of 2 vector.
   inline xVector<VT> operator%(const xVector<VT> &other) const;
   inline xVector<VT> operator*(const VT &other) const;

   inline xVector<VT> operator/(const VT &other) const;
   inline double angleRad(const xVector<VT> &v) const;
   inline double angleDeg(const xVector<VT> &v) const;

   inline ~xVector() = default;
   friend class xTensor2<VT>;
};

using xVectorDoubleComplex = xVector<std::complex<double>>;

#include "xVector_imp.h"

}  // namespace xtensor

namespace xtool
{
// Casts...Only explicit casts allowed...
struct xCastToDoubleVector : public std::unary_function<xtensor::xVector<float>, xtensor::xVector<double>>
{
  public:
   typedef std::unary_function<xtensor::xVector<float>, xtensor::xVector<double>> base_type;
   typedef typename base_type::argument_type argument_type;
   typedef typename base_type::result_type result_type;

   result_type operator()(argument_type &t) const
   {
      return xtensor::xVector<double>(static_cast<double>(t(0)), static_cast<double>(t(1)), static_cast<double>(t(2)));
   }
};
}  // namespace xtool

#endif
