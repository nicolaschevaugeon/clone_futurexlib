/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _BLASDEF_
#define _BLASDEF_
#include <complex>

const char TRANSPOSE = 'T';
const char NOTRANSPOSE = 'N';
const char UPPERTRIANGULAR = 'U';
const char LOWERTRIANGULAR = 'L';
const int INCONE = 1;

// BLAS L3 extern def
extern "C" void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, const double *alpha,
                       const double *A, const int *LDA, const double *B, const int *LDB, const double *beta, double *C,
                       const int *LDC);

extern "C" void sgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, const float *alpha,
                       const float *A, const int *LDA, const float *B, const int *LDB, const float *beta, float *C,
                       const int *LDC);

extern "C" void dsyrk_(const char *UPLO, const char *TRANS, const int *N, const int *K, const double *alpha, const double *A,
                       const int *LDA, const double *beta, double *C, const int *LDC);

extern "C" void dsyr2k_(const char *UPLO, const char *TRANS, const int *N, const int *K, const double *alpha, const double *A,
                        const int *LDA, const double *B, const int *LDB, const double *beta, double *C, const int *LDC);

extern "C" void dtrmm_(const char *SIDE, const char *UPLO, const char *TRANSA, const char *DIAG, const int *M, const int *N,
                       const double *alpha, const double *A, const int *LDA, double *B, const int *LDB);

// BLAS L2 extern def
extern "C" void dgemv_(const char *TRANS, const int *M, const int *N, const double *ALPHA, const double *A, const int *LDA,
                       const double *X, const int *INCX, const double *BETA, double *Y, const int *INCY);
extern "C" void dsymv_(const char *UPLO, const int *N, const double *ALPHA, const double *A, const int *LDA, const double *X,
                       const int *INCX, const double *BETA, const double *Y, const int *INCY);

extern "C" void dtrmv_(const char *UPLO, const char *TRANS, const char *DIAG, const int *N, const double *A, const int *LDA,
                       double *X, const int *INCX);

extern "C" void dtbsv_(const char *UPLO, const char *TRANS, const char *DIAG, int *N, int *K, double *A, int *LDA, double *X,
                       int *INCX);

// BLAS L1 extern def
extern "C" void daxpy_(const int *N, const double *a, const double *x, const int *incx, double *y, const int *incy);
extern "C" void zaxpy_(const int *N, const void *a, const void *x, const int *incx, void *y, const int *incy);
extern "C" double dnrm2_(const int *N, const double *x, const int *incx);
extern "C" float sdot_(const int *N, const float *x, const int *incx, const float *y, const int *incy);
extern "C" double ddot_(const int *N, const double *x, const int *incx, const double *y, const int *incy);
extern "C" double zdotc_(const int *N, const void *zx, const int *incx, const void *zy, const int *incy);
extern "C" double zdotu_(const int *N, const void *zx, const int *incx, const void *zy, const int *incy);
extern "C" void dscal_(const int *N, const double *alpha, double *x, const int *incx);
extern "C" void zscal_(const int *N, const void *alpha, void *x, const int *incx);
extern "C" int idamax_(const int *N, const double *x, const int *incx);
extern "C" int izamax_(const int *N, const void *x, const int *incx);

extern "C" void saxpy_(const int *N, const float *a, const float *x, const int *incx, float *y, const int *incy);
extern "C" void sscal_(const int *N, const float *alpha, float *x, const int *incx);
extern "C" int isamax_(const int *N, const float *x, const int *incx);

namespace xlinalg
{
template <typename T>
class xCPPBlasDef;

template <>
class xCPPBlasDef<float>
{
  public:
   static void scal(int N, float alpha, float *x, int incx = 1) { sscal_(&N, &alpha, x, &incx); }
   static void axpy(int N, float alpha, const float *x, float *y, int incx = 1, int incy = 1)
   {
      saxpy_(&N, &alpha, x, &incx, y, &incy);
   }
   static float dot(int N, const float *x, const float *y, int incx = 1, int incy = 1) { return sdot_(&N, x, &incx, y, &incy); }
   static int iamax(int N, const float *x, int incx = 1) { return isamax_(&N, x, &incx); }
   static void gemm(char TRANSA, char TRANSB, const int M, const int N, const int K, float alpha, const float *A, int LDA,
                    const float *B, int LDB, float beta, float *C, int LDC)
   {
      sgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
   }
};

template <>
class xCPPBlasDef<double>
{
  public:
   static void scal(int N, double alpha, double *x, int incx = 1) { dscal_(&N, &alpha, x, &incx); }
   static void axpy(int N, double alpha, const double *x, double *y, int incx = 1, int incy = 1)
   {
      daxpy_(&N, &alpha, x, &incx, y, &incy);
   }
   static double dot(int N, const double *x, const double *y, int incx = 1, int incy = 1)
   {
      return ddot_(&N, x, &incx, y, &incy);
   }
   static double nrm2(int N, const double *x, int incx = 1) { return dnrm2_(&N, x, &incx); }
   static int iamax(int N, double *x, int incx = 1) { return idamax_(&N, x, &incx); }

   static void gemv(char TRANS, int M, int N, double alpha, const double *A, int LDA, const double *X, double beta, double *Y,
                    int incx = 1, int incy = 1)
   {
      dgemv_(&TRANS, &M, &N, &alpha, A, &LDA, X, &incx, &beta, Y, &incy);
   }

   static void gemm(char TRANSA, char TRANSB, const int M, const int N, const int K, double alpha, const double *A, int LDA,
                    const double *B, int LDB, double beta, double *C, int LDC)
   {
      dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
   }
   static void syrk(char UPLO, char TRANS, int N, int K, double alpha, const double *A, int LDA, double beta, double *C, int LDC)
   {
      dsyrk_(&UPLO, &TRANS, &N, &K, &alpha, A, &LDA, &beta, C, &LDC);
   }
   static void syr2k(char UPLO, char TRANS, int N, int K, double alpha, const double *A, int LDA, const double *B, int LDB,
                     double beta, double *C, int LDC)
   {
      dsyr2k_(&UPLO, &TRANS, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
   }
   static void trmm(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double alpha, const double *A, int LDA, double *B,
                    int LDB)
   {
      dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &alpha, A, &LDA, B, &LDB);
   }

   static void symv(char UPLO, int N, double alpha, const double *A, int LDA, const double *X, double beta, const double *Y,
                    int incx = 1, int incy = 1)
   {
      dsymv_(&UPLO, &N, &alpha, A, &LDA, X, &incx, &beta, Y, &incy);
   }
};

template <>
class xCPPBlasDef<std::complex<double>>
{
  public:
   static void scal(int N, std::complex<double> alpha, std::complex<double> *x, int incx = 1) { zscal_(&N, &alpha, x, &incx); }
   static void axpy(int N, std::complex<double> alpha, const std::complex<double> *x, std::complex<double> *y, int incx = 1,
                    int incy = 1)
   {
      zaxpy_(&N, &alpha, x, &incx, y, &incy);
   }
   static const std::complex<double> dot(int N, const std::complex<double> *x, const std::complex<double> *y, int incx = 1,
                                         int incy = 1)
   {
      return dotc(N, x, y, incx, incy);
   }
   static const std::complex<double> dotc(int N, const std::complex<double> *x, const std::complex<double> *y, int incx = 1,
                                          int incy = 1)
   {
      return static_cast<std::complex<double>>(zdotc_(&N, x, &incx, y, &incy));
   }
   static const std::complex<double> dotu(int N, const std::complex<double> *x, const std::complex<double> *y, int incx = 1,
                                          int incy = 1)
   {
      return static_cast<std::complex<double>>(zdotu_(&N, x, &incx, y, &incy));
   }
   static int iamax(int N, std::complex<double> *x, int incx = 1) { return izamax_(&N, x, &incx); }
};

}  // namespace xlinalg

#endif
