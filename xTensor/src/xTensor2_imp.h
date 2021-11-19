#ifndef _XTENSOR2_H_
#error Do NOT include xTensor2_imp.h alone
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// #########################################################################################
// IMPLEMENTATION of inlined function from class xTensor2
// #########################################################################################

template <typename VT>
inline xTensor2<VT>::xTensor2(const VT &a00, const VT &a01, const VT &a02, const VT &a10, const VT &a11, const VT &a12,
                              const VT &a20, const VT &a21, const VT &a22)
    : pos{{a00, a01, a02}, {a10, a11, a12}, {a20, a21, a22}}
{
}

template <typename VT>
inline xTensor2<VT>::xTensor2(const xVector<VT> &r0, const xVector<VT> &r1, const xVector<VT> &r2)
    : pos{{r0.pos[0], r0.pos[1], r0.pos[2]}, {r1.pos[0], r1.pos[1], r1.pos[2]}, {r2.pos[0], r2.pos[1], r2.pos[2]}}
{
}

template <typename VT>
inline xTensor2<VT>::xTensor2(const VT &a) : pos{{a, a, a}, {a, a, a}, {a, a, a}}
{
}

template <typename VT>
inline xTensor2<VT>::xTensor2(const xTensor2<VT> &in)
    : pos{{in.pos[0][0], in.pos[0][1], in.pos[0][2]},
          {in.pos[1][0], in.pos[1][1], in.pos[1][2]},
          {in.pos[2][0], in.pos[2][1], in.pos[2][2]}}
{
}

template <typename VT>
inline xTensor2<VT>::xTensor2(const xTensor2<VT> &in, const VT &scal)
    : pos{{in.pos[0][0] * scal, in.pos[0][1] * scal, in.pos[0][2] * scal},
          {in.pos[1][0] * scal, in.pos[1][1] * scal, in.pos[1][2] * scal},
          {in.pos[2][0] * scal, in.pos[2][1] * scal, in.pos[2][2] * scal}}
{
}

template <typename VT>
inline xTensor2<VT> &xTensor2<VT>::operator=(const xTensor2<VT> &other)
{
   VT *l = data();
   const VT *r = other.data();
   memcpy(l, r, 9 * sizeof(VT));
   return *this;
}

template <typename VT>
inline VT xTensor2<VT>::delta(int i, int j) const
{
   return ((i == j) ? 1.0 : 0.0);
}

template <typename VT>
inline xVector<VT> xTensor2<VT>::operator*(const xVector<VT> &other) const
{
   return xVector<VT>{pos[0][0] * other.pos[0] + pos[0][1] * other.pos[1] + pos[0][2] * other.pos[2],
                      pos[1][0] * other.pos[0] + pos[1][1] * other.pos[1] + pos[1][2] * other.pos[2],
                      pos[2][0] * other.pos[0] + pos[2][1] * other.pos[1] + pos[2][2] * other.pos[2]};
}

template <typename VT>
inline xTensor2<VT> xTensor2<VT>::operator*(const xTensor2<VT> &other) const
{
   return xTensor2<VT>{
       pos[0][0] * other.pos[0][0] + pos[0][1] * other.pos[1][0] + pos[0][2] * other.pos[2][0],
       pos[0][0] * other.pos[0][1] + pos[0][1] * other.pos[1][1] + pos[0][2] * other.pos[2][1],
       pos[0][0] * other.pos[0][2] + pos[0][1] * other.pos[1][2] + pos[0][2] * other.pos[2][2],

       pos[1][0] * other.pos[0][0] + pos[1][1] * other.pos[1][0] + pos[1][2] * other.pos[2][0],
       pos[1][0] * other.pos[0][1] + pos[1][1] * other.pos[1][1] + pos[1][2] * other.pos[2][1],
       pos[1][0] * other.pos[0][2] + pos[1][1] * other.pos[1][2] + pos[1][2] * other.pos[2][2],

       pos[2][0] * other.pos[0][0] + pos[2][1] * other.pos[1][0] + pos[2][2] * other.pos[2][0],
       pos[2][0] * other.pos[0][1] + pos[2][1] * other.pos[1][1] + pos[2][2] * other.pos[2][1],
       pos[2][0] * other.pos[0][2] + pos[2][1] * other.pos[1][2] + pos[2][2] * other.pos[2][2],
   };
}

template <typename VT>
inline xTensor2<VT> &xTensor2<VT>::operator*=(const VT &scalar)
{
   VT *l = data();
   for (unsigned int i = 0; i < 9; ++i)
   {
      (*l++) *= scalar;
   }
   return *this;
}

template <typename VT>
inline xTensor2<VT> xTensor2<VT>::operator*(const VT &scalar) const
{
   return xTensor2<VT>(*this, scalar);
}

template <typename VT>
inline xTensor2<VT> &xTensor2<VT>::operator/=(const VT &scalar)
{
   VT *l = data();
   for (unsigned int i = 0; i < 9; ++i)
   {
      *l++ /= scalar;
   }
   return *this;
}

template <typename VT>
inline xTensor2<VT> &xTensor2<VT>::operator+=(const xTensor2<VT> &other)
{
   VT *l = data();
   const VT *r = other.pos[0];
   for (unsigned int i = 0; i < 9; ++i) *l++ += *r++;
   return *this;
}

template <typename VT>
inline xTensor2<VT> &xTensor2<VT>::operator-=(const xTensor2<VT> &other)
{
   VT *l = data();
   const VT *r = other.data();
   for (unsigned int i = 0; i < 9; ++i)
   {
      *l -= *r;
      ++l;
      ++r;
   }
   return *this;
}

template <typename VT>
inline xTensor2<VT> xTensor2<VT>::operator+(const xTensor2<VT> &other) const
{
   xTensor2<VT> m(*this);
   m += other;
   return m;
}

template <typename VT>
inline xTensor2<VT> xTensor2<VT>::operator-(const xTensor2<VT> &other) const
{
   xTensor2<VT> m(*this);
   m -= other;
   return m;
}

template <typename VT>
inline xTensor2<VT> xTensor2<VT>::operator-() const
{
   xTensor2<VT> m;
   m -= *this;
   return m;
}

template <typename VT>
inline VT xTensor2<VT>::contract(const xTensor2<VT> &other) const
{
   const VT *l = data();
   const VT *r = other.data();
   VT x = *l++ * *r++;
   for (unsigned int i = 0; i < 8; ++i) x += *l++ * *r++;
   return x;
}

template <typename VT>
inline xTensor2<VT> &xTensor2<VT>::symmetrize()
{
   VT oneHalf = xtool::xDataType<VT>::one() / 2;
   pos[0][1] = pos[1][0] = oneHalf * (pos[0][1] + pos[1][0]);
   pos[0][2] = pos[2][0] = oneHalf * (pos[0][2] + pos[2][0]);
   pos[2][1] = pos[1][2] = oneHalf * (pos[2][1] + pos[1][2]);
   return *this;
}

template <typename VT>
inline VT *xTensor2<VT>::data()
{
   return pos[0];
}

template <typename VT>
inline const VT *xTensor2<VT>::data() const
{
   return pos[0];
}

template <typename VT>
inline VT &xTensor2<VT>::operator()(int i, int j)
{
   return pos[i][j];
}

template <typename VT>
inline const VT &xTensor2<VT>::operator()(int i, int j) const
{
   return pos[i][j];
}

// #########################################################################################
// IMPLEMENTATION of  inlined function
// #########################################################################################

template <typename VT>
inline xTensor2<VT> operator*(const VT &other, const xTensor2<VT> &t)
{
   xTensor2<VT> m(t * other);
   return m;
}

template <typename VT>
inline std::ostream &operator<<(std::ostream &s, const xTensor2<VT> &t)
{
   if (xTensor2<VT>::getExportFormat())
   {
      s << t(0, 0) << " " << t(0, 1) << " " << t(0, 2) << std::endl;
      s << t(1, 0) << " " << t(1, 1) << " " << t(1, 2) << std::endl;
      s << t(2, 0) << " " << t(2, 1) << " " << t(2, 2) << std::endl;
   }
   else
   {
      s << t(0, 0) << " " << t(0, 1) << " " << t(0, 2) << " " << t(1, 0) << " " << t(1, 1) << " " << t(1, 2) << " " << t(2, 0)
        << " " << t(2, 1) << " " << t(2, 2);
   }
   return s;
}

template <typename VT>
inline std::istream &operator>>(std::istream &s, xTensor2<VT> &t)
{
   return (s >> t(0, 0) >> t(0, 1) >> t(0, 2) >> t(1, 0) >> t(1, 1) >> t(1, 2) >> t(2, 0) >> t(2, 1) >> t(2, 2));
}

template <typename VT>
inline xTensor2<VT> tensor_product(const xVector<VT> &a, const xVector<VT> &b)
{
   return xTensor2<VT>{a(0) * b(0), a(0) * b(1), a(0) * b(2), a(1) * b(0), a(1) * b(1),
                       a(1) * b(2), a(2) * b(0), a(2) * b(1), a(2) * b(2)};
}

template <typename VT>
inline void axpy_tensor2(const VT &a, const xTensor2<VT> &x, xTensor2<VT> &y)
{
   VT *py = y.data();
   const VT *px = x.data();
   for (unsigned int i = 0; i < 9; ++i)
   {
      *py += a * (*px);
      ++py;
      ++px;
   }
}

template <typename VT>
inline void copyvalues(const xTensor2<VT> &in, VT *vals)
{
   const VT *r = in.data();
   for (unsigned int i = 0; i < 9; ++i)
   {
      *vals = *r;
      ++vals;
      ++r;
   }
}

template <typename VT>
inline int sizedouble(const xTensor2<VT> &in)
{
   return 9;
}

/// Generate a rotation matrix in a xTensor2 container from a given vector.
//! Vector given, vn, is considered as a already normalized vector.
//! Transverse axes are arbitrary chosen to give the more accurate rotation matix.
//! Calling R the resulting Matrix, if V={x,y,z} represent a vector in a 3 dimentional cartesian coordinate system,
//! U = R.V = {x',y',z'} represente a vector in a 3 dimentional cartesian coordinate system having it's first axis colinear to vn.
//! This may be used for example to generate rotation matrix for xTensor4 methode rotate
template <typename VT>
inline xTensor2<VT> generateRotateFromVect(const xVector<VT> &vn)
{
   const double v0 = vn(0);
   const double v1 = vn(1);
   const double v2 = vn(2);
   const double v02 = v0 * v0;
   const double v12 = v1 * v1;
   const double v22 = v2 * v2;

   double nt1 = v12 + v22;
   double nt2 = v02 + v22;
   double nt3 = v02 + v12;

   xVector<VT> t11;

   if (nt1 > nt3)
   {
      if (nt1 > nt2)
      {
         nt1 = sqrt(nt1);
         t11(1) = v2 / nt1;
         t11(2) = -v1 / nt1;
      }
      else
      {
         nt2 = sqrt(nt2);
         t11(0) = -v2 / nt2;
         t11(2) = v0 / nt2;
      }
   }
   else
   {
      if (nt3 > nt2)
      {
         nt3 = sqrt(nt3);
         t11(0) = v1 / nt3;
         t11(1) = -v0 / nt3;
      }
      else
      {
         nt2 = sqrt(nt2);
         t11(0) = -v2 / nt2;
         t11(2) = v0 / nt2;
      }
   }

   xVector<VT> t22 = vn % t11;
   t22.norm();
   return xTensor2<VT>(vn, t11, t22);
}

/// Generate a rotation matrix in a xTensor2 container from a given vector.
//! Vector given, vn, is considered as a already normalized vector.
//! Calling R the resulting Matrix, if V={x,y,z} represents a vector in a 3 dimensional cartesian coordinate system,
//! U = R.V = {x',y',z'} represents a vector in a 3 dimensional cartesian coordinate system having it's first axis colinear to vn.
//! This may be used for example to generate matrix rotation for xTensor4 methode rotate
template <typename VT>
inline xTensor2<VT> generateRotateFromVect2D(const xVector<VT> &vn)
{
   return xTensor2<VT>(vn, xVector<VT>(-vn(1), vn(0), xtool::xDataType<VT>::zero()),
                       xVector<VT>(xtool::xDataType<VT>::zero(), xtool::xDataType<VT>::zero(), xtool::xDataType<VT>::one()));
}

template <>
inline bool choleskyL(const xTensor2<double> &A, xTensor2<double> &L, const double prec)
{
   const double norm = A.normValue();
   const double eps = prec * norm;
   xVector<double> d = {0., 0., 0.};
   L = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
   for (unsigned int j = 0; j < 3; ++j)
   {
      d(j) = A(j, j);
      for (unsigned int k = 0; k < j; ++k)
      {
         d(j) -= L(j, k) * L(j, k) * d(k);
      }
      if (fabs(d(j)) > eps)
         for (unsigned int i = j + 1; i < 3; ++i)
         {
            L(i, j) = A(i, j);
            for (unsigned int k = 0; k < j; ++k)
            {
               L(i, j) -= L(i, k) * L(j, k) * d(k);
            }
            L(i, j) /= d(j);
         }
      else
      {
         d(j) = 0.;
         for (unsigned int i = j + 1; i < 3; ++i)
         {
            L(i, j) = 0.;
         }
      }
   }
   for (unsigned int j = 0; j < 3; ++j)
   {
      if (d(j) < 0.)
         return false;
      else
      {
         const double piv = sqrt(d(j));
         for (unsigned int i = j; i < 3; ++i)
         {
            L(i, j) *= piv;
         }
      }
   }
   return true;
}

template <>
inline bool choleskyU(const xTensor2<double> &A, xTensor2<double> &U, const double prec)
{
   const double norm = A.normValue();
   const double eps = prec * norm;
   xVector<double> d = {0., 0., 0.};
   U = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
   for (unsigned int i = 0; i < 3; ++i)
   {
      d(i) = A(i, i);
      for (unsigned int k = 0; k < i; ++k)
      {
         d(i) -= U(k, i) * U(k, i) * d(k);
      }
      if (fabs(d(i)) > eps)
         for (unsigned int j = i + 1; j < 3; ++j)
         {
            U(i, j) = A(i, j);
            for (unsigned int k = 0; k < i; ++k)
            {
               U(i, j) -= U(k, i) * U(k, j) * d(k);
            }
            U(i, j) /= d(i);
         }
      else
      {
         d(i) = 0.;
         for (unsigned int j = i + 1; j < 3; ++j)
         {
            U(i, j) = 0.;
         }
      }
   }
   for (unsigned int i = 0; i < 3; ++i)
   {
      if (d(i) < 0.)
         return false;
      else
      {
         const double piv = sqrt(d(i));
         for (unsigned int j = i; j < 3; ++j)
         {
            U(i, j) *= piv;
         }
      }
   }
   return true;
}

template <class charT, class traits, typename VT>
inline std::basic_istream<charT, traits> &read_tensor2(std::basic_istream<charT, traits> &strm, xTensor2<VT> &t)
{
   VT data[9];
   unsigned int nopen = 0;
   while ((strm.peek() == '{') || (strm.peek() == ' ') || (strm.peek() == '('))
   {
      if ((strm.peek() == '{') || (strm.peek() == '(')) nopen++;
      strm.ignore();
   }
   for (unsigned int i = 0; i < 9; ++i)
   {
      strm >> data[i];
      while ((strm.peek() == ',') || (strm.peek() == ' ') || (strm.peek() == ';'))
      {
         strm.ignore();
      }
   }
   while (nopen != 0 && strm.good())
   {
      if ((strm.peek() == '}') || (strm.peek() == ')')) nopen--;
      strm.ignore();
   }
   if (strm.fail())
   {
      std::cout << "a parentezes was open but never closed while reading a xTensor2" << std::endl;
      throw;
   }
   t = xTensor2<VT>(data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8]);
   return strm;
}

//-------------------------------------------------------------------------------------------------

template <typename VT>
VT xTensor2<VT>::trace() const
{
   return pos[0][0] + pos[1][1] + pos[2][2];
}

template <typename VT>
VT xTensor2<VT>::trace2() const
{
   double a00 = pos[0][0] * pos[0][0] + pos[1][0] * pos[0][1] + pos[2][0] * pos[0][2];
   double a11 = pos[1][0] * pos[0][1] + pos[1][1] * pos[1][1] + pos[1][2] * pos[2][1];
   double a22 = pos[2][0] * pos[0][2] + pos[2][1] * pos[1][2] + pos[2][2] * pos[2][2];

   return a00 + a11 + a22;
}

template <typename VT>
VT xTensor2<VT>::det() const
{
   return pos[0][0] * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]) -
          pos[0][1] * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]) +
          pos[0][2] * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
}

template <typename VT>
xTensor2<VT> &xTensor2<VT>::transpose()
{
   for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = i + 1; j < 3; j++)
      {
         double temp = pos[j][i];
         pos[j][i] = pos[i][j];
         pos[i][j] = temp;
      }
   return *this;
}

template <typename VT>
xTensor2<VT> xTensor2<VT>::operator!() const
{
   xTensor2 B(*this);
   B.transpose();
   return B;
}

template <typename VT>
VT xTensor2<VT>::vonMisesNorm() const
{
   VT vm;
   VT a = 0;
   VT b = trace();
   for (unsigned int i = 0; i < 3; ++i)
   {
      for (unsigned int j = 0; j < 3; ++j)
      {
         a += pos[i][j] * pos[i][j];
      }
   }
   vm = std::sqrt(1.5 * a - b * b * 0.5);
   return vm;
}

template <typename VT>
VT xTensor2<VT>::normValue() const
{
   xTensor2<VT> mat(*this);
   return std::sqrt(mat.contract(mat));
}

template <typename VT>
xTensor2<VT> xTensor2<VT>::intersect2d(xTensor2<VT> &t) const
{
   xVector<VT> e1[3], e2[3];
   VT v1[3], v2[3];

   eigen2d(e1, v1);
   t.eigen2d(e2, v2);

   if (fabs(v1[0]) > fabs(v2[0]))
   {
      xTensor2<VT> TR(e1[0], e1[1], e1[2]);
      xTensor2<VT> TRT = TR;
      TRT.transpose();
      xTensor2<VT> D(0.0);
      D(0, 0) = v1[0];
      D(2, 2) = 1.0;
      VT x = (t * e1[1]) * e1[1];
      if (fabs(x) > fabs(v1[1]))
         D(1, 1) = x;
      else
         D(1, 1) = v1[1];
      return (TR * D) * TRT;
   }
   else
   {
      xTensor2<VT> TR(e2[0], e2[1], e2[2]);
      xTensor2<VT> TRT = TR;
      TRT.transpose();
      xTensor2<VT> D(0.0);
      D(0, 0) = v2[0];
      D(2, 2) = 1.0;
      VT x = ((*this) * e2[1]) * e2[1];
      if (fabs(x) > fabs(v2[1]))
         D(1, 1) = x;
      else
         D(1, 1) = v2[1];
      return (TR * D) * TRT;
   }
}

template <typename VT>
xTensor2<VT> xTensor2<VT>::operator*(const xTensor4<VT> &other) const
{
   xTensor2 m(0.0);
   for (unsigned int i = 0; i < 3; i++)
   {
      for (unsigned int j = 0; j < 3; j++)
      {
         // m(i,j) = 0.0;
         for (unsigned int k = 0; k < 3; k++)
         {
            for (unsigned int l = 0; l < 3; l++)
            {
               m(i, j) += pos[k][l] * other(k, l, i, j);
            }
         }
      }
   }
   return m;
}

/* xTensor2 xTensor2<VT>::operator * (const xTensor4Isotropic &other) const
     {
     xTensor2 m(0.0);
     double t = trace();
     for(unsigned int i=0;i<3;i++){
     for(unsigned int j=0;j<3;j++){
     m(i,j) =  t * delta(i,j) * other.lam + other.mu * ( pos[i][j] + pos[j][i] );
     }
     }
     return m;
     }
  */
template <typename VT>
bool xTensor2<VT>::export_format = true;

//-------------------------------------------------------------------------------------------------
