/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xBasicIterativeSolvers_
#define _xBasicIterativeSolvers_
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>

#include "xDiagonalMatrix.h"

//! A very simple conjugate gradient implementation.
/*! Template of M which must be a model of Matrix
    and V which must be a model of Vector
    needed member function of type V:
    V::size() return the size of a Vector.
    V::V(int size) : vector constructor;
    V::V(const V & in) : vector copy constructor;

    needed function of V and M:
    void axpy(double alpha, const V &X, V& Y ) <- compute Y:= Y + alpha X. blas style Vector addition.
    void gemv(double alpha,  const M & A, const V & X , double beta,  V & Y )   <- Compute Y := beta Y +alpha AX. blas style
 Matrix vector product. double dot(const V & X, const V & Y)    <- return the scalar product of X and Y. void scal(double alpha,
 V& Y) <- compute Y :=  alpha *Y double nrm2(V& Y) <- compute the norm 2 of vector Y.
 !*/
template <class M, class V>
void conjgrad(const M& A, const V& b, V& x)
{
   double eps = 1.e-6;
   int N = x.size();
   V r(b);
   gemv(-1, A, x, 1, r);  // r = r -Ax;
   double drprp = dot(r, r);
   V p(r);
   V Ap(N);
   int k = 0;
   do
   {
      gemv(1, A, p, 0., Ap);  // Ap = A *p;
      double drr = drprp;
      double alpha = drr / dot(p, Ap);  //(r*r) / (p *Ap);
      axpy(alpha, p, x);                // x = x +a*p // x += alpha* p;
      axpy(-alpha, Ap, r);              // rp = r - alpha * Ap;
      if (nrm2(r) < eps) break;
      drprp = dot(r, r);
      double beta = drprp / drr;  // rp*rp / r*r
      scal(beta, p);
      axpy(1., r, p);  // p = rp +beta p
      ++k;
   } while (k < 10000);
   cout << "iter : " << k << endl;
}

//! A very simple conjugate gradient implementation, with a jacobi preconditioner
/*! Template of M which must be a model of Matrix
    and V which must be a model of Vector
    needed member function of type V:
    V::size() return the size of a Vector.
    V::V(int size) : vector constructor;
    V::V(const V & in) : vector copy constructor;
    M::operator()(int i, int j) -> acces the i,j term of the matrix. only use here to get thediagonal elements.
    needed functions of V and M:
    void axpy(double alpha, const V &X, V& Y ) <- compute Y:= Y + alpha X. blas style Vector addition.
    void gemv(double alpha,  const M & A, const V & X , double beta,  V & Y )   <- Compute Y := beta Y +alpha AX. blas style
 Matrix vector product. double dot(const V & X, const V & Y)    <- return the scalar product of X and Y. void scal(double alpha,
 V& Y) <- compute Y :=  alpha *Y double nrm2(V& Y) <- compute the norm 2 of vector Y.

 !*/
template <class M, class V>
void conjgradjacobi(const M& A, const V& b, V& x)
{
   xDiagonalMatrix<double> MP(A);
   double eps = 1.e-6;
   int N = x.size();
   V r(b);
   gemv(-1, A, x, 1, r);  // r = r -Ax;
   V z(N);
   solve(MP, r, z);
   // z=r;
   double drpzp = dot(r, z);
   V p(z);
   V Ap(N);
   int k = 0;
   do
   {
      gemv(1, A, p, 0., Ap);  // Ap = A *p;
      double drz = drpzp;
      double alpha = drz / dot(p, Ap);  //(r*r) / (p *Ap);
      axpy(alpha, p, x);                // x = x +a*p // x += alpha* p;
      axpy(-alpha, Ap, r);              // rp = r - alpha * Ap;
      if (nrm2(r) < eps) break;
      solve(MP, r, z);
      // z=r;
      drpzp = dot(r, z);
      double beta = drpzp / drz;  // rp*rp / r*r
      scal(beta, p);
      axpy(1., z, p);  // p = rp +beta p
      ++k;
   } while (k < 10000);
   cout << "iter : " << k << endl;
}

//! A very simple conjugate gradient implementation, with a ssor preconditioner
/*! Template of M which must be a model of Matrix
    and V which must be a model of Vector
    needed member function of type V:
    V::size() return the size of a Vector.
    V::V(int size) : vector constructor;
    V::V(const V & in) : vector copy constructor;
    M::operator()(int i, int j) -> acces the i,j term of the matrix. only use here to get thediagonal elements.
    needed functions of V and M:
    void axpy(double alpha, const V &X, V& Y ) <- compute Y:= Y + alpha X. blas style Vector addition.
    void gemv(double alpha,  const M & A, const V & X , double beta,  V & Y )   <- Compute Y := beta Y +alpha AX. blas style
 Matrix vector product. double dot(const V & X, const V & Y)    <- return the scalar product of X and Y. void scal(double alpha,
 V& Y) <- compute Y :=  alpha *Y double nrm2(V& Y) <- compute the norm 2 of vector Y. M and V must also conform to the
 prerequisite of the ssor template function provided below.
 !*/
template <class M, class V>
void conjgradssor(const M& A, const V& b, V& x)
{
   xDiagonalMatrix<double> D(A);
   double eps = 1.e-6;
   int N = x.size();
   V r(b);
   gemv(-1, A, x, 1, r);  // r = r -Ax;
   V z(N);
   ssor(A, D, r, z);
   // z=r;
   double drpzp = dot(r, z);
   V p(z);
   V Ap(N);
   int k = 0;
   double nr;
   do
   {
      gemv(1, A, p, 0., Ap);  // Ap = A *p;
      double drz = drpzp;
      double alpha = drz / dot(p, Ap);  //(r*r) / (p *Ap);
      axpy(alpha, p, x);                // x = x +a*p // x += alpha* p;
      axpy(-alpha, Ap, r);              // rp = r - alpha * Ap;
      nr = nrm2(r);
      if (nr < eps) break;
      ssor(A, D, r, z);
      // z=r;
      drpzp = dot(r, z);
      double beta = drpzp / drz;  // rp*rp / r*r
      scal(beta, p);
      axpy(1., z, p);  // p = rp +beta p
      ++k;
      // cout << "iter : "<< k  <<  " " << nr << endl;
   } while (k < 10000);
   cout << "iter : " << k << endl;
}

//! compute the "shift" of x by an ssor preconditioner
/*!
 A is supposed to be symetric
 A = d + L + L^T  where d is the diagonal of A  L is the lower triangular part of A.
 the function compute :
 x :=   (d+L)   d^(-1) (d+L^T) b
 Input : b, A, d. b is the vector to be shifted, A is the original matrix (supposed symmetric, but stored in full ), d is the
 diagonal of the original matrix. the function is template on M, the matrix type, D the type of the diagonal matrix, and V the
 type of the vectors. prerequist : void V::operator=(const V& in) <- copy. V V::V(int i) <- constructor. trsv(const char * uplow,
 const &M, & x) : compute the x := x + ( L+D) x  or x+ (L+D)^T x depending of the value of uplow ("L" or "U") (blas 2 type
 function)
*/
template <class M, class D, class V>
void ssor(const M& A, const D& d, const V& b, V& x)
{
   x = b;
   V y(A.ncolumns);
   trsv("L", A, x);
   gemv(1., d, x, 0., y);
   trsv("U", A, y);
   x = y;
}

#endif
