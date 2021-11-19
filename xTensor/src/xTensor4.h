/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

//////////////////////////////////////////////////////////////////////
#ifndef _XTensor4_H_
#define _XTensor4_H_

#include "xTensor2.h"

namespace xtensor {

/**
     Simple class for a 4th order tensor in 3d which
     has a 3x3x3x3 matrix representation.
  */

template<typename VT = double>
class xTensor4Base
{
public:
    virtual inline VT delta(int i, int j) const {return ((i == j)?1.0:0.0);}
    virtual xTensor2<VT> operator * (const xTensor2<VT> &other) const =0;
    virtual ~xTensor4Base()= default;
};


/// Class than specilize operation for xTensor4 that are "isotropic" : only 2 component needed.
class xTensor4Isotropic :public xTensor4Base<double>
{
public:
    double lam, mu;
    inline xTensor4Isotropic ();
    inline xTensor4Isotropic (const double& E, const double& nu);
    inline xTensor2<double> operator * (const xTensor2<double> &other) const override;
    inline xTensor4Isotropic operator * (const double &other) const;
    inline xTensor4Isotropic &operator *=(const double &other) ;
    inline double getLambda() const {return lam;}
    inline double getMu() const {return mu;}
    inline void setLambda(const double & _lam)  {lam = _lam;}
    inline void setMu(const double & _mu)  {mu = _mu;}
};



template<typename VT = double>
class xTensor4 : public xTensor4Base<VT>
{
public:
    /// Default Constructor. Note : data are not initialized.
    inline xTensor4 ();
    /// All the component of the constructed xTensor4 are set to init.
    inline xTensor4 (const VT& init);
    /// Give all the 81 component of the xTensor4 to onstruct it. Make it possible to use initializer lists.
    inline xTensor4 (const VT &a0000, const VT &a0001,  const VT &a0002,
                     const VT &a0010, const VT &a0011,  const VT &a0012,
                     const VT &a0020, const VT &a0021,  const VT &a0022,

                     const VT &a0100, const VT &a0101,  const VT &a0102,
                     const VT &a0110, const VT &a0111,  const VT &a0112,
                     const VT &a0120, const VT &a0121,  const VT &a0122,

                     const VT &a0200, const VT &a0201,  const VT &a0202,
                     const VT &a0210, const VT &a0211,  const VT &a0212,
                     const VT &a0220, const VT &a0221,  const VT &a0222,

                     const VT &a1000, const VT &a1001,  const VT &a1002,
                     const VT &a1010, const VT &a1011,  const VT &a1012,
                     const VT &a1020, const VT &a1021,  const VT &a1022,

                     const VT &a1100, const VT &a1101,  const VT &a1102,
                     const VT &a1110, const VT &a1111,  const VT &a1112,
                     const VT &a1120, const VT &a1121,  const VT &a1122,

                     const VT &a1200, const VT &a1201,  const VT &a1202,
                     const VT &a1210, const VT &a1211,  const VT &a1212,
                     const VT &a1220, const VT &a1221,  const VT &a1222,

                     const VT &a2000, const VT &a2001,  const VT &a2002,
                     const VT &a2010, const VT &a2011,  const VT &a2012,
                     const VT &a2020, const VT &a2021,  const VT &a2022,

                     const VT &a2100, const VT &a2101,  const VT &a2102,
                     const VT &a2110, const VT &a2111,  const VT &a2112,
                     const VT &a2120, const VT &a2121,  const VT &a2122,

                     const VT &a2200, const VT &a2201,  const VT &a2202,
                     const VT &a2210, const VT &a2211,  const VT &a2212,
                     const VT &a2220, const VT &a2221,  const VT &a2222
                     );
    inline xTensor4 (const VT c[3][3][3][3]);
    //! generate an isotropic xTensor4 usefull for mechanical application where E and nu are the young modulus and poisson ratio.
    inline xTensor4 (const VT E, const VT nu);
    //! generate some kind of cubic or isotropic transverse material by nullifying
    //! shear termes of the standard isotropic hook tensor
    inline xTensor4 (const VT E, const VT nu, const VT beta );
    //! Copy constructor
    inline xTensor4(const xTensor4& in);
    inline xTensor4 &operator = (const xTensor4 &other);
    //! convert Tensor4 from Trellis to xTensor4
    //GILLES:    inline xTensor4 (const Trellis_Util::mTensor4& in);
    //! return a pointer to the first element of an array of 81 VT (3*3*3*3),
    //! containing the component of the xTensor4, stored in row major order.
    //! a0000 a0001 a0002 a0010 a0011 a0012 a0020 a0021 a0022 .....
    //! modifying data inside the array modify the tensor ...
    inline VT *data();
    //! return a const pointer to the first element of an array of 81 VT (3*3*3*3),
    //! containing the component of the xTensor4, stored in row major order.
    //! a0000 a0001 a0002 a0010 a0011 a0012 a0020 a0021 a0022 .....
    inline const VT *data() const;
    inline const VT &operator () (int i, int j, int k, int l) const;
    inline VT & operator () (int i, int j, int k, int l);
    /// return a norm of the tensor : sqrt(a_ijkl*a_ijkl)
    inline VT normValue() const ;
    inline xTensor4 &operator += (const xTensor4 &other);
    inline xTensor4 &operator -= (const xTensor4 &other);
    inline xTensor4 operator -  (const xTensor4 &other) const;
    inline xTensor4 operator +  (const xTensor4 &other) const;
    inline xTensor4 operator *  (const VT &scalar) const;
    inline xTensor4 operator /  (const VT &scalar) const;
    inline xTensor2<VT> operator * (const xTensor2<VT> &other) const override;
    inline xTensor4 &operator /= (const VT &scalar);
    inline xTensor4 &operator *= (const VT &scalar);
    inline xTensor4 rotate(const xTensor2<VT>&);
    /// return the transpose of this (A^T_ijkl =  A_klij
    inline xTensor4 operator !() const;
    
private:
    VT pos[3][3][3][3];
};

/// Output an xTensor4 to a stream (usefull to print on screen for example)
template<typename VT>
inline std::ostream& operator<<(std::ostream& s, const xTensor4<VT>& t);

/// Product of 2 xTensor4 C =A *B where C_ijkl = A_ijmn * Bmnkl
template<typename VT>
inline xTensor4<VT> operator*(const xTensor4<VT> &A, const xTensor4<VT> &B);

//! Cholesky Factorisation of a xTensor4
//! Only The upper triangle part of A is accessed (A_ijkl accessed if  "3*i+j => 3*k+j) during the factorisation.
//!   Note that all the coefficients are accessed when the norm of A is computed.
//! Upon success it return true and the xTensor4 U such as A = U^TU and Uijkl =0 if 3*k+l < 3*i+j.
//!  Note that in the above equation multiplication is defined as A_ijkl = B_ijmn * C_mnkl and the
//!  Transposition refer to A^T_ijkl = A_klij
//! It is garantied to succeed if the A is numerically symmetric positive definite.
//! If A is only semi-definite (0 eigenvalues), Some pivot during the factorisation might become clause to zero.
//! if the |piv| < prec*norm(A), piv and the associated column of L are set to zero.
//! prec can be given as a optional parameter.
//Only coded for double
inline bool choleskyU(const xTensor4<double> &A, xTensor4<double> &U, const double prec = 1.e-6);

//Only coded for double
inline bool choleskyU(const xTensor4Isotropic &A4, xTensor4<double> &U, const double = 1.e-6);




//! Cholesky Factorisation of a xTensor4
//! Only The lower triangle part of A is accessed (A_ijkl accessed if  "3*i+j<= 3*k+j) during the factorisation.
//!   Note that all the coefficients are accessed when the norm of A is computed.
//! Upon success it return true and the xTensor4 L such as A = LL^T and Lijkl =0 if 3*k+l > 3*i+j.
//!  Note that in the above equation multiplication is defined as A_ijkl = B_ijmn * C_mnkl and the
//!  Transposition refer to A^T_ijkl = A_klij
//! It is garantied to succeed if the A is numerically symmetric positive definite.
//! If A is only semi-definite (0 eigenvalues), Some pivot during the factorisation might become clause to zero.
//! if the |piv| < prec*norm(A), piv and the associated column of L are set to zero.
//! prec can be given as a optional parameter.
//Only coded for double
inline bool choleskyL(const xTensor4<double> &A, xTensor4<double> &L, const double prec = 1.e-6);

//Only coded for double
inline bool choleskyL(const xTensor4Isotropic &A4, xTensor4<double> &L, const double prec = 1.e-6);


/// This class implements a tensor of order 4 for possibly anisotropic behaviour under the assumption of plane strain.
/// The coefficients of the tensor are stored as a xTensor2 T, which links the vector S=(sigma_11, sigma_22, sigma_12) to
/// the vector E=(epsilon_11, epsilon_22,   2 * epsilon_12 !!!! ) as S = T E
class xTensor4AnisoPlaneStrain :public xTensor4Base<double> {
public:
    xTensor2<double> coefs;
    double nu;
    double E;
    inline xTensor4AnisoPlaneStrain();
    inline xTensor4AnisoPlaneStrain(const xTensor2<double> & _coefs, double _nu);
    inline xTensor4AnisoPlaneStrain(const double _E, const double _nu);
    inline xTensor4AnisoPlaneStrain(const double _E, const double _nu, const double beta );
    inline xTensor2<double> operator * (const xTensor2<double> &other) const override;
    inline xTensor4AnisoPlaneStrain operator * (const double &other) const;
    inline xTensor4AnisoPlaneStrain rotate(const xTensor2<double>& rot) const;
    inline void updateCoeffs(const double beta = 1.);

};

class xTensor4AnisoPlaneStress :public xTensor4Base<double> {
public:
    xTensor2<double> coefs;
    double nu;
    double E;
    inline xTensor4AnisoPlaneStress();
    inline xTensor4AnisoPlaneStress(const xTensor2<double> & _coefs, double _nu);
    inline xTensor4AnisoPlaneStress(const double _E, const double _nu);
    inline xTensor2<double> operator * (const xTensor2<double> &other) const override;
    inline xTensor4AnisoPlaneStress operator * (const double &other) const;
    inline void updateCoeffs();
};



#include "xTensor4_imp.h"

} // end of namespace

#endif
