/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xDenseMatrix.h"

#include <algorithm>
#include <exception>
#include <iomanip>
#include <sstream>
#include <string>

#include "xBlasDef.h"
#include "xCSRVector.h"

class xDenseMatrixException : public std::exception
{
  public:
   xDenseMatrixException(std::string type, int Line)
   {
      std::ostringstream oss;
      oss << "xDenseMatrix : " << type << " exception, line " << Line << std::endl;
      msg = oss.str();
   }
   ~xDenseMatrixException() throw() override = default;
   const char *what() const throw() override { return this->msg.c_str(); }

  private:
   std::string msg;
};

namespace xlinalg
{
xDenseMatrix::xDenseMatrix() : nl(0), nc(0), a(nullptr), col_ass(0) {}

xDenseMatrix::xDenseMatrix(const int _nl, const int _nc) : nl(_nl), nc(_nc), a(nullptr), col_ass(0)
{
   a = new double[nl * nc];
   std::fill(&a[0], &a[nl * nc], 0.);
}

xDenseMatrix::xDenseMatrix(const int _nl) : nl(_nl), nc(_nl), a(nullptr), col_ass(0)
{
   a = new double[nl * nc];
   std::fill(&a[0], &a[nl * nc], 0.);
}

xDenseMatrix::xDenseMatrix(const xDenseMatrix &in) : nl(in.nl), nc(in.nc), a(nullptr), col_ass(0)
{
   a = new double[nl * nc];
   std::copy(&in.a[0], &in.a[nl * nc], &a[0]);
}

xDenseMatrix::xDenseMatrix(std::istream &in, int option) : a(nullptr), col_ass(0)
{
   in >> nl >> nc;
   // std::cout << nl << " " << nc << std::endl;
   a = new double[nl * nc];
   for (int j = 0; j < nc; ++j)
   {
      for (int i = 0; i < nl; ++i) in >> (*this)(i, j);
   }
}

xDenseMatrix::~xDenseMatrix()
{
   if (a) delete[] a;
   a = nullptr;
}

void xDenseMatrix::resize(const int _nl, const int _nc, double v_init)
{
   nl = _nl;
   nc = _nc;
   if (a) delete[] a;
   a = new double[nl * nc];
   std::fill(&a[0], &a[nl * nc], v_init);
}

xDenseMatrix &xDenseMatrix::operator=(const xDenseMatrix &A)
{
   if (this == &A) return *this;
   if (!((nl * nc) == (A.nc * A.nl)))
   {
      if (a) delete[] a;
      a = new double[A.nl * A.nc];
   }
   nl = A.nl;
   nc = A.nc;
   std::copy(&A.a[0], &A.a[nl * nc], &a[0]);
   return *this;
}

// Other implementation
//   xDenseMatrix& xDenseMatrix::operator=(const xDenseMatrix &A){
//             xDenseMatrix temp(A);
//             nl = temp.nl;
//             nc = temp.nc;
//             std::swap(this->a,temp.a);
//             return *this;
//       }

void xDenseMatrix::operator+=(const xDenseMatrix &A)
{
   if ((nl == A.nl) && (nc == A.nc))
   {
      std::transform(&a[0], &a[nl * nc], &A.a[0], &a[0], std::plus<double>());
   }
   else
      throw xDenseMatrixException("Dimensions not compatible", __LINE__);
}

void xDenseMatrix::operator-=(const xDenseMatrix &A)
{
   if ((nl == A.nl) && (nc == A.nc))
   {
      std::transform(&a[0], &a[nl * nc], &A.a[0], &a[0], std::minus<double>());
   }
   else
      throw xDenseMatrixException("Dimensions not compatible", __LINE__);
}

void xDenseMatrix::operator*=(const double &alpha) { xCPPBlasDef<double>::scal(nl * nc, alpha, &a[0]); }
double xDenseMatrix::maxnorm() const
{
   int idmax = xCPPBlasDef<double>::iamax(nl * nc, a);
   return fabs(a[idmax - 1]);
}
const double &xDenseMatrix::operator()(const int i, const int j) const
{
#ifndef NDEBUG
   if ((i >= nl) || (j >= nc)) throw xDenseMatrixException("Out of bounds", __LINE__);
#endif
   return a[i + j * nl];
}

double &xDenseMatrix::operator()(const int i, const int j)
{
#ifndef NDEBUG
   if ((i >= nl) || (j >= nc)) throw xDenseMatrixException("Out of bounds", __LINE__);
#endif
   return a[i + j * nl];
}

void xDenseMatrix::AddMatrix(int i, int j, double val) { (*this)(i - 1, j - 1) += val; }

void xDenseMatrix::AddVal(const int &i, const double &val) { (*this)(i - 1, col_ass) += val; }

void xDenseMatrix::setAssCol(const int col_num)
{
#ifndef NDEBUG
   if ((col_num >= nc) || (col_num < 0)) throw xDenseMatrixException("Out of bounds", __LINE__);
#endif
   col_ass = col_num;
}

double xDenseMatrix::GetMatrix(int i, int j) const { return (*this)(i - 1, j - 1); }

xDenseMatrix operator*(const xDenseMatrix &A, const xDenseMatrix &B)
{
   if (A.nc != B.nl) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
   xDenseMatrix C(A.nl, B.nc);
   gemm(false, false, 1, A, B, 0., C);
   return C;
}

xCSRVector operator*(const xDenseMatrix &A, const xCSRVector &X)
{
   if (A.nc != X.size()) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
   xCSRVector Y(A.nl);
   gemv(false, 1., A, X, 0., Y);
   return Y;
}

void xDenseMatrix::print2screen() const
{
   for (int i = 0; i < nl; ++i)
   {
      for (int j = 0; j < nc; ++j)
      {
         std::cout << (*this)(i, j) << ' ';
      }
      std::cout << std::endl;
   }
}

void xDenseMatrix::export2Matlab(std::string fileName, std::string variableName) const
{
   std::ofstream fout(fileName.c_str());
   export2Matlab(fout, variableName);
   fout.close();
}

void xDenseMatrix::export2Matlab(std::ostream &out, std::string variableName) const
{
   out << variableName << "= [";
   for (int i = 0; i < nl; ++i)
   {
      for (int j = 0; j < nc; ++j)
      {
         out << (*this)(i, j) << ' ';
      }
      if (i != nl - 1)
         out << "; ";
      else
         out << "];" << std::endl;
   }
}

void xDenseMatrix::export2txt(std::string fileName, unsigned precision, int option) const
{
   std::ofstream out(fileName.c_str());
   if (option == 1)
   {
      out << nl << " " << nc << std::endl;
      ;
   }
   out.precision(precision);
   out << std::showpos << std::scientific;
   for (int i = 0; i < nl; ++i)
   {
      for (int j = 0; j < nc; ++j)
      {
         out << (*this)(i, j) << "  ";
      }
      out << std::endl;
   }
   out.close();
}

std::ostream &operator<<(std::ostream &output, const xDenseMatrix &A)
{
   const int nl = A.nl;
   const int nc = A.nc;

   output << std::setprecision(4) << "[";
   for (int i = 0; i < nl - 1; ++i)
   {
      for (int j = 0; j < nc - 1; ++j)
      {
         output << std::setw(8) << A(i, j) << " ";
      }
      output << std::setw(8) << A(i, nc - 1) << ";" << std::endl;
   }
   for (int j = 0; j < nc - 1; ++j)
   {
      output << std::setw(8) << A(nl - 1, j) << " ";
   }
   output << std::setw(8) << A(nl - 1, nc - 1) << "]" << std::endl;
   return output;  // for multiple << operators.
}

void gemm(bool TRANSA, bool TRANSB, double ALPHA, const xDenseMatrix &A, const xDenseMatrix &B, double BETA, xDenseMatrix &C)
{
   int M, N, K, LDA, LDB;
   char TA, TB;
   if (TRANSA)
   {
      M = A.ncolumn();
      K = A.nline();
      LDA = K;
      TA = 'T';
   }
   else
   {
      M = A.nline();
      K = A.ncolumn();
      LDA = M;
      TA = 'N';
   }
   if (TRANSB)
   {
      N = B.nline();
      LDB = N;
      TB = 'T';
      if (K != B.ncolumn()) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
   }
   else
   {
      N = B.ncolumn();
      LDB = K;
      TB = 'N';
      if (K != B.nline()) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
   }

   xlinalg::xCPPBlasDef<double>::gemm(TA, TB, M, N, K, ALPHA, A.ReturnValPointer(), LDA, B.ReturnValPointer(), LDB, BETA,
                                      C.ReturnValPointer(), M);
}

void syrk(const char &UPLO, const bool &TRANS, const double &ALPHA, const xDenseMatrix &A, const double &BETA, xDenseMatrix &C)
{
   int N, K, LDA, LDC;
   char T;
   N = C.nline();
   if (C.ncolumn() != N) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
   LDC = N;
   if (TRANS)
   {
      K = A.nline();
      if (A.ncolumn() != N) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
      LDA = K;
      T = 'T';
   }
   else
   {
      K = A.ncolumn();
      if (A.nline() != N) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
      LDA = N;
      T = 'N';
   }
   xlinalg::xCPPBlasDef<double>::syrk(UPLO, T, N, K, ALPHA, A.ReturnValPointer(), LDA, BETA, C.ReturnValPointer(), LDC);
}

void gemv(bool TRANSA, double ALPHA, const xDenseMatrix &A, const xCSRVector &x, double BETA, xCSRVector &y)
{
   int M = A.nline();
   int N = A.ncolumn();
   // T = Transpose // N = no Transpose
   char trans = 'N';         // No transpose
   if (TRANSA) trans = 'T';  // Transpose
   xlinalg::xCPPBlasDef<double>::gemv(trans, M, N, ALPHA, const_cast<double *>(A.ReturnValPointer()), M,
                                      const_cast<double *>(&x[0]), BETA, &y[0]);
}

void axpy(const double &alpha, const xDenseMatrix &A, xDenseMatrix &B)
{
   int N = A.nline() * A.ncolumn();
   int M = B.nline() * B.ncolumn();
   if (M != N) throw xDenseMatrixException("Dimensions not compatible", __LINE__);
   xlinalg::xCPPBlasDef<double>::axpy(N, alpha, const_cast<double *>(A.ReturnValPointer()), B.ReturnValPointer());
}

void transpose(const xDenseMatrix &A, xDenseMatrix &B)
{
   int nl = A.nline();
   int nc = A.ncolumn();
   B.resize(nc, nl);
   for (int i = 0; i < nl; ++i)
   {
      for (int j = 0; j < nc; ++j)
      {
         B(j, i) = A(i, j);
      }
   }
}

}  // namespace xlinalg
