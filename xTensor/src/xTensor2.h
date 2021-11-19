/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

//////////////////////////////////////////////////////////////////////
#ifndef _XTENSOR2_H_
#define _XTENSOR2_H_

#include <algorithm>
#include <cassert>
#include <cstring>
#include <functional>
#include <iostream>
#include <numeric>

#include "xVector.h"

namespace xtensor
{
// Forward declarations
template <typename VT>
class xTensor4;

class xTensor4Isotropic;
class xTensor4AnisoPlaneStrain;
class xTensor4AnisoPlaneStress;

/**
     Simple class for a 2nd order tensor in 3d which
     has a 3x3 matrix representation.
  */

template <typename VT = double>
class xTensor2
{
  public:
   /// Constructor : insert all term.
   inline xTensor2(const VT &a00, const VT &a01, const VT &a02, const VT &a10, const VT &a11, const VT &a12, const VT &a20,
                   const VT &a21, const VT &a22);
   /// Constructor :   xTensor2 row by row.
   inline xTensor2(const xVector<VT> &r0, const xVector<VT> &r1, const xVector<VT> &r2);
   /// Constructor : all terms equal to a. If not given a = 0.;
   inline xTensor2(const VT &a = xtool::xDataType<VT>::zero());
   /// A constructor that do not initialize the data in the xTensor2
   //     inline xTensor2(bool noinit){};
   /// Copy Constructor
   inline xTensor2(const xTensor2<VT> &in);
   /// Scaled Copy Constructor (return an xTensor2, which is a copy of in multiplyed by scal)
   inline xTensor2(const xTensor2<VT> &in, const VT &scal);
   /// Construct a xTensor2 from Trellis mTensor2
   // GILLES:    inline xTensor2(const Trellis_Util::mTensor2& in);
   /// Affectation operator. affect to this the data contained in other.
   inline xTensor2<VT> &operator=(const xTensor2<VT> &other);
   /// Return the value of the krneker symbol ( 0. if i!=j; 1 1. if i == j)
   inline VT delta(int i, int j) const;
   /// Return  the xTensor2 containing the contraction of the xTensor2 with the xTensor4 other.
   /*! res_ij  = this_kl*other_klij !*/
   xTensor2<VT> operator*(const xTensor4<VT> &other) const;
   xTensor2<VT> operator*(const xTensor4Isotropic &other) const;
   xTensor2<VT> operator*(const xTensor4AnisoPlaneStrain &other) const;
   xTensor2<VT> operator*(const xTensor4AnisoPlaneStress &other) const;
   /// return the xVector result of the simple contraction of this with the xVector other
   /*!  res_i = this_ij other_j */
   inline xVector<VT> operator*(const xVector<VT> &other) const;
   /// return the xTensor result of the simple contraction of this with the xTensor2 other
   /*!  res_ij = this_ik other_kj */
   inline xTensor2<VT> operator*(const xTensor2<VT> &other) const;
   /// Scaling this by multiplying all the term by scal.
   inline xTensor2<VT> &operator*=(const VT &scalar);
   /// Return a new xTensor2 that is the result of scaling this by scalar.
   inline xTensor2<VT> operator*(const VT &scalar) const;
   /// Scaling this by dividing  all the term by 1./scal
   inline xTensor2<VT> &operator/=(const VT &scalar);
   /// adding other to this.
   inline xTensor2<VT> &operator+=(const xTensor2<VT> &other);
   /// removing other to this.
   inline xTensor2<VT> &operator-=(const xTensor2<VT> &other);
   /// return a new tensor which is the result of adding other to this
   inline xTensor2<VT> operator+(const xTensor2<VT> &other) const;
   /// return a new xTensor2 which is the result of removing other to this.
   inline xTensor2<VT> operator-(const xTensor2<VT> &other) const;
   /// return a new xTensor2 which is the opposite of this
   inline xTensor2<VT> operator-() const;
   /// return the double contraction of this with other.
   /*! res = this_ij*other_ij */
   inline VT contract(const xTensor2<VT> &other) const;
   //! return a pointer to the first element of an array of 9 VT (3*3),
   //! containing the component of the xTensor2, stored in row major order.
   //! a00 a01 a02 a10 a11 a12 a20 a21 a22
   //! modifying data inside the array modify the tensor ...
   inline VT *data();
   //! return a const pointer to the first element of an array of 9 VT (3*3),
   //! containing the component of the xTensor2, stored in row major order.
   //! a00 a01 a02 a10 a11 a12 a20 a21 a22
   inline const VT *data() const;
   /// return reference to term ij of the xTensor2.
   inline VT &operator()(int i, int j);
   /// return const reference to term ij of the xTensor2.
   inline const VT &operator()(int i, int j) const;
   /// return an xVector res that is the solution of the linear system this *res = other
   xVector<VT> operator/(const xVector<VT> &other) const;
   /// Return the invert of this
   xTensor2<VT> invert() const;
   /// Return the trace of this (res = a00+a11+a22)
   VT trace() const;
   /// Return the second invariant of this
   VT trace2() const;
   /// Return the third invariant of this (determinant)
   VT det() const;
   long eigen(xVector<VT> e[3], VT v[3]) const;
   long eigen2d(xVector<VT> e[3], VT v[3]) const;
   /// This compute eigenvalues and theire first and  second derivative
   //! All those are calculated with analytical expression
   //! Compare to eigen methode it doesn't compute eigen vector but compute first and second derivative of the eigenvalues
   //! if dersec is a null pointeur only eigenvalues and theire first derivative are computed
   //! if dersec and der are null pointer only  eigenvalues are computed
   //! For now we can't obtain second derivative without computing first.
   //! if eignvalue equation lead to commplex solution a exception is thrown
   //! if eignvalue equation have a null polynomial discriminant computation of derivative throw a exception
   void getAnalyticalEigenvaluesAndDerivative(double lambda[3], double der[][3][3] = nullptr,
                                              double dersec[][3][3][3][3] = nullptr);
   inline xTensor2<VT> &symmetrize();
   xTensor2<VT> intersect2d(xTensor2<VT> &t) const;
   /// Transpose tensor this
   xTensor2<VT> &transpose();
   /// return the transpose of this.
   xTensor2<VT> operator!() const;
   /// return the vonMisesNorm of the tensor sqrt(3*a/2. - b*b/2  ) where a = this_ij*this_ij and b = this_ii
   VT vonMisesNorm() const;
   /// return a norm of the tensor : sqrt(a_ij*a_ij)
   VT normValue() const;
   /// set the export format used to push this to a stream.
   /*! if export format  == 1 the xTensor2  will be printed line by line with a line jump between each line?
     Otherwise , all the 9 values  are printed in one line.
   */
   static void setExportFormat(const bool export_format_) { export_format = export_format_; }

   static bool getExportFormat() { return export_format; }

  private:
   friend class xTensor4<VT>;
   friend class xTensor4Isotropic;
   template <typename UT>
   friend void axpy_tensor2(const UT &a, const xTensor2<UT> &x, xTensor2<UT> &y);

   template <typename UT>
   friend void copyvalues(const xTensor2<UT> &, VT *);

   template <typename UT>
   friend int sizedouble(const xTensor2<UT> &);

   VT pos[3][3];
   static bool export_format;
};

/// Compute the cubic root of the polynomial with coefficient given by coeff. array x contains the roots.
/*!  return the number of real root founded */
long FindCubicRoots(const double coeff[4], double x[3]);
/// Compute the null space of the 3x3 matrix pointed to by a
long NullSpace(const double *a, double *result, double eps, long n);
/// return a new xtensor 2 result of other*t.
template <typename VT>
inline xTensor2<VT> operator*(const VT &other, const xTensor2<VT> &t);
/// Output an xTensor2 in a stream
template <typename VT>
inline std::ostream &operator<<(std::ostream &s, const xTensor2<VT> &t);
/// Read the content of tensor t from a stream (row by row)
template <typename VT>
inline std::istream &operator>>(std::istream &s, xTensor2<VT> &t);
/// Return the xTensor2 res_ij = a_i*b_j
template <typename VT>
inline xTensor2<VT> tensor_product(const xVector<VT> &a, const xVector<VT> &b);
/// set xTensor2 y  = y + a*x ... Blas level 1 style interface
template <typename VT>
inline void axpy_tensor2(const VT &a, const xTensor2<VT> &x, xTensor2<VT> &y);
/// Copy the value of the xTensor2 in to the array pointed to by vals. need room for 9 VT ..
template <typename VT>
inline void copyvalues(const xTensor2<VT> &in, VT *vals);
/// Return the number of VT contain in an xTensor2 (9)
template <typename VT>
inline int sizedouble(const xTensor2<VT> &in);
/// Generate a rotation matrix in a xTensor2 container from a given vector.
//! Vector given, vn, is considered as a already normalized vector.
//! Transverse axes are arbitrary chosen to give the more accurate rotation matix.
//! Calling R the resulting Matrix, if V={x,y,z} represent a vector in a 3 dimentional cartesian coordinate system,
//! U = R.V = {x',y',z'} represente a vector in a 3 dimentional cartesian coordinate system having it's first axis colinear to vn.
//! This may be used for example to generate rotation matrix for xTensor4 methode rotate
template <typename VT = double>
inline xTensor2<VT> generateRotateFromVect(const xVector<VT> &vn);
/// Generate a rotation matrix in a xTensor2 container from a given vector.
//! Vector given, vn, is considered as a already normalized vector.
//! Calling R the resulting Matrix, if V={x,y,z} represents a vector in a 3 dimensional cartesian coordinate system,
//! U = R.V = {x',y',z'} represents a vector in a 3 dimensional cartesian coordinate system having it's first axis colinear to vn.
//! This may be used for example to generate matrix rotation for xTensor4 methode rotate
template <typename VT>
inline xTensor2<VT> generateRotateFromVect2D(const xVector<VT> &vn);
//! Cholesky Factorisation of a xTensor2
//! Only The lower triangle part of A is accessed (A_ij accessed if  j<= i)during the factorisation.
//!   Note that all the coefficients are accessed when the norm of A is computed.
//! Upon success it return true and the xTensor2 L such as A = LL^T and Lij =0 if j > i.
//! It is garantied to succeed if the A is numerically symmetric positive definite.
//! If A is only semi-definite (0 eigenvalues), Some pivot during the factorisation might become clause to zero.
//! if the |piv| < prec*norm(A), piv and the associated column of L are set to zero.
//! prec can be given as a optional parameter.
template <typename VT = double>  // Only coded for double
inline bool choleskyL(const xTensor2<VT> &A, xTensor2<VT> &L, const VT prec = 1.e-6);
//! Cholesky Factorisation of a xTensor2
//! Only The upper triangle part of A is accessed (A_ij accessed if  i<= j)during the factorisation.
//!   Note that all the coefficients are accessed when the norm of A is computed.
//! Upon success it return true and the xTensor2 U such as A = U^TU and Uij =0 if i > j.
//! It is garantied to succeed if the A is numerically symmetric positive definite.
//! If A is only semi-definite (0 eigenvalues), or close to semi-definite (small eigenvalues),
//! Some pivots during the factorisation might become close to zero.
//! if the |piv| < prec*norm(A), piv and the associated line of U are set to zero.
//! prec can be given as a optional parameter.
template <typename VT = double>  // Only coded for double
inline bool choleskyU(const xTensor2<VT> &A, xTensor2<VT> &U, const VT prec = 1.e-6);

using xTensor2DoubleComplex = xTensor2<std::complex<double>>;

#include "xTensor2_imp.h"

}  // namespace xtensor

#endif
