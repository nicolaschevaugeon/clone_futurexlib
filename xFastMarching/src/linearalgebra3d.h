/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _linearalgebra3d_
#define _linearalgebra3d_


#include <iostream>
#include <array>
#include <cmath>
#include <cassert>
#include <numeric>

namespace xfastmarching
{

  template <class T>
    class vector2d{
  public:
    vector2d()= default;;
  vector2d(T x, T y):data({x,y}){}
    inline T &operator()(const size_t &i){assert(i <2); return data[i];}
    inline const T &operator()(const size_t &i) const {assert(i <2); return data[i];}
    inline T &operator[](const size_t &i){assert(i <2); return data[i];}
    inline const T &operator[](const size_t &i) const {assert(i <2); return data[i];}
    inline T norm2() const{
            return sqrt(data[0]* data[0] + data[1]*data[1]) ;
          }
    template < class U >
      friend U dot(const vector2d<U> &a, const vector2d<U> &b );
  private:
    std::array<T, 2> data = {{0.,0.}};
  };

  template <class T>
    class _e3{
  public:
  _e3(T x=1.):l(x){}
    inline const T & operator()  (const size_t &i) const {
      if (i == 2) return l;
      throw;
    }
  private :
    T l;
  };

  static const _e3<double > e3 = {1.};

  template <class T>
    inline T externprod(const vector2d<T> &a, const vector2d<T> &b ){
    return a(0)*b(1) - a(1)*b(0);
  }

  template <class T>
    inline T dot (const T &a, const T &b){ return a*b;}


  template <class T>
    inline T dot(const vector2d<T> &a, const vector2d<T> &b ){
    return a.data[0]*b.data[0]+ a.data[1]*b.data[1];
    //return a(0)*b(0) + a(1)*b(1);
  }

  template <class T>
    inline vector2d<T> &  operator*=( vector2d<T> & a , T  alpha){
    a(0)*=alpha; a(1)*=alpha;
    return a;
  }

  template <class T>
    inline vector2d<T>   operator*( const vector2d<T> & a , T  alpha){
    vector2d<T> b(a);
    return b*=alpha;
  }

  template <class T>
    inline vector2d<T>   operator*(  T  alpha, const vector2d<T> & a){
    return a*alpha;
  }

  template <class T>
    inline vector2d<T>  scal( T  alpha,   vector2d<T> & a){
    return a*=alpha;
  }

  template <class T>
    inline T  nrm2( const vector2d<T> & a){
    return sqrt(a(0)*a(0) + a(1)*a(1));
  }

  template <class T>
    inline vector2d<T> & operator+=( vector2d<T> & a, const vector2d<T> &b){
    a(0) +=b(0); 
    a(1) +=b(1); 
    return a;
  }

  template <class T>
    inline vector2d<T> & operator-=( vector2d<T> & a, const vector2d<T> &b){
    a(0) -=b(0); 
    a(1) -=b(1); 
    return a;
  }

  template <class T>
    inline vector2d<T> operator+(const vector2d<T> & a, const vector2d<T> &b){
    vector2d<T > c(a);
    c+= b;
    return c;
  }

  template <class T>
    inline vector2d<T> operator-(const vector2d<T> & a, const vector2d<T> &b){
    vector2d<T > c(a);
    c-= b;
    return c;
  }

  template <class T>
    inline std::ostream& operator<<(std::ostream& os, const vector2d<T>& v){
    os << v(0)<< " " << v(1) << std::endl;
    return os;
  }

  ////
  template <class T>
    class vector3d{
  public:
    vector3d()= default;;
  vector3d(T x, T y, T z):data({x,y,z}){}
    inline T &operator()(const size_t &i){assert(i <3); return data[i];}
    inline const T &operator()(const size_t &i) const {assert(i <3); return data[i];}
    inline T &operator[](const size_t &i){assert(i <3); return data[i];}
    inline const T &operator[](const size_t &i) const {assert(i <3); return data[i];}
    typename::std::array<T, 3>::const_iterator begin() const{return data.begin();}
    typename::std::array<T, 3>::const_iterator end() const{return data.end();}
    inline T norm2() const{
            return sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] );
          }
  
  private:
    std::array<T, 3> data = {{0.,0.,0.}};
  };

  template <class T>
    inline vector3d<T> externprod(const vector3d<T> &a, const vector3d<T> &b );

  template <class T>
    inline T sum(const vector3d<T> &a ){
    return a(0)+a(1) +a(2);
  }

  template <class T>
    inline T dot(const vector3d<T> &a, const vector3d<T> &b );

  template <class T>
    inline vector3d<T> &  operator*=( vector3d<T> & a , T  alpha){
    a(0)*=alpha; a(1)*=alpha; a(2)*=alpha;
    return a;
  }

  template <class T>
    inline void axpy(T alpha, const vector2d<T> & a,  vector2d<T> &b){
    b(0) += alpha*a(0);
    b(1) += alpha*a(1);
  };

  template <class T>
    inline vector3d<T>   operator*( const vector3d<T> & a , T  alpha){
    vector3d<T> b(a);
    return b*=alpha;
  }

  template <class T>
    inline vector3d<T>   operator*(  T  alpha, const vector3d<T> & a){
    return a*alpha;
  }

  template <class T>
    inline vector3d<T>  scal( T  alpha,   vector3d<T> & a){
    return a*=alpha;
  }

  template <class T>
    inline T  nrm2( const vector3d<T> & a){
    return sqrt(a(0)*a(0) + a(1)*a(1) + a(2)*a(2));
  }

  template <class T>
    inline vector3d<T> & operator+=( vector3d<T> & a, const vector3d<T> &b);

  template <class T>
    inline vector3d<T> & operator-=( vector3d<T> & a, const vector3d<T> &b);

  template <class T>
    inline vector3d<T> operator+(const vector3d<T> & a, const vector3d<T> &b);

  template <class T>
    inline vector3d<T> operator-(const vector3d<T> & a, const vector3d<T> &b);

  template <class T>
    inline std::ostream& operator<<(std::ostream& os, const vector3d<T>& v);

  template <class T>
    class point3d{
  public:
    point3d()= default;;
  point3d(T x, T y, T  z):data({x,y,z}){}
    inline T &operator()(const size_t &i){assert(i <3); return data[i];}
    inline const T &operator()(const size_t &i) const {assert(i <3); return data[i];}
  private:
    std::array<T, 3> data = {0.,0.,0.};
  };

  template <class T>
    class point2d{
  public:
    point2d()= default;;
  point2d(T x, T y):data({x,y}){}
    inline T &operator()(const size_t &i){assert(i <2); return data[i];}
    inline const T &operator()(const size_t &i) const {assert(i <2); return data[i];}
  private:
    std::array<T, 3> data = {0.,0.};
  };

  template <class T>
    vector3d<T > fromto(const point3d<T> &from, const point3d<T> &to){
    return vector3d<T>{ to(0)-from(0), to(1)-from(1), to(2)-from(2)};
  }

  template <class T>
    vector2d<T > fromto(const point2d<T> &from, const point2d<T> &to){
    return vector2d<T>{ to(0)-from(0), to(1)-from(1)};
  }


  template <class T>
    class tensor2d{
  public:
    tensor2d()= default;;
    tensor2d (const T & x00, const T &x10, 
	      const T & x01, const T &x11):data({x00,x10,x01,x11}){}
  tensor2d(const tensor2d & in):data(in.data){}
    inline T &operator()( size_t  i,  size_t j){
      assert(i< 2 && j < 2);
      return data[i+j*2];
    }
    inline const T &operator()(size_t i, size_t j) const {
      assert(i< 2 && j < 2);
      return data[i+j*2];
    }
  private:
    std::array<T, 4> data = {{0.,0.,0., 0.}};
  };


  template <class T>
    class tensor3d2d{
  public:
    tensor3d2d()= default;;
    tensor3d2d (const T & x00, const T &x10, const T& x20,
		const T & x01, const T &x11, const T &x21):
    data{x00,x10, x20, x01, x11, x21 }{};
    tensor3d2d (const vector3d <T> & c0, const vector3d<T> & c1):data{ c0(0), c0(1), c0(2), c1(0), c1(1), c1(2)}{};
  
    tensor3d2d (const tensor3d2d & in)= default;;
    inline T &operator()( size_t  i,  size_t j){
      assert(i< 3 && j < 2);
      return data[i+j*3];
    }
    inline const T &operator()(size_t i, size_t j) const {
      assert(i< 3 && j < 2);
      return data[i+j*3];
    }
  private:
    std::array<T, 6> data = {{0.,0.,0., 0.,0.,0.} };
  };

  template <class T>
    inline vector3d<T> operator*( const tensor3d2d<T> &A, const vector2d<T> &X ){
    return vector3d<T> {  A(0,0)*X(0) + A(0,1)*X(1),
	A(1,0)*X(0) + A(1,1)*X(1),
	A(2,0)*X(0) + A(2,1)*X(1) };
  };

  template <class T>
    inline vector2d<T> operator*(const vector3d< T > & X,  const tensor3d2d<T> &A ){
    return vector2d<T> {  A(0,0)*X(0) + A(1,0)*X(1) + A(2,0)*X(2), 
	A(0,1)*X(0) + A(1,1)*X(1) + A(2,1)*X(2), 
	};
  };


  template <class T>
    std::ostream& operator<<(std::ostream& os, const tensor2d<T>& A);

  template <class T>
    inline T det(const tensor2d<T> & A);

  template <class T>
    tensor2d<T> invert(const tensor2d<T> & A );

  template <class T>
    inline void gemm( T alpha, const tensor2d<T> &A, const tensor2d<T> &B, T beta,  tensor2d<T> &C);

  template <class T>
    inline tensor2d<T> operator*( const tensor2d<T> &A, const tensor2d<T> &B);
  
  template <class T>
    inline void gemv( T alpha, const tensor2d<T> &A, const vector2d<T> &X, T beta,  vector2d<T> &Y);

  template <class T>
    inline vector2d<T> operator*( const tensor2d<T> &A, const vector2d<T> &X );

  template <class T>
    inline tensor2d<T> transpose( const tensor2d<T> &A ){
    return tensor2d<T>{A(0,0), A(0,1),  A(1,0), A(1,1)};

  }
  /// t3D
  template <class T>
    class tensor3d{
  public:
    tensor3d()= default;;
    tensor3d (const T & x00, const T &x10, const T &x20, 
	      const T & x01, const T &x11, const T &x21,
	      const T & x02, const T &x12, const T &x22):data({x00,x10,x20,x01,x11,x21,x02,x12,x22}){}

    tensor3d (const vector3d<T> & c0, const vector3d<T> & c1, const vector3d<T> &c2):data({c0(0), c0(1), c0(2), c1(0), c1(1), c1(2), c2(0), c2(1), c2(2)}){}

    tensor3d(const tensor3d<T > & in )= default;

    inline T &operator()( size_t  i,  size_t j){
      assert(i< 3 && j < 3);
      return data[i+j*3];
    }
    inline const T &operator()(size_t i, size_t j) const {
      assert(i< 3 && j < 3);
      return data[i+j*3];
    }
  public:
    std::array<T, 9> data = {{0.,0.,0., 0.,0.,0., 0.,0.,0.}};
  };


  template <class T>
    std::ostream& operator<<(std::ostream& os, const tensor3d<T>& A);

  template <class T>
    inline T det(const tensor3d<T> & A);

  template <class T>
    tensor3d<T> invert(const tensor3d<T> & A );

  template <class T>
    inline void gemm( T alpha, const tensor3d<T> &A, const tensor3d<T> &B, T beta,  tensor3d<T> &C);

  template <class T>
    inline tensor3d<T> operator*( const tensor3d<T> &A, const tensor3d<T> &B);
  
  template <class T>
    inline void gemv( T alpha, const tensor3d<T> &A, const vector3d<T> &X, T beta,  vector3d<T> &Y);

  template <class T>
    inline vector3d<T> operator*( const tensor3d<T> &A, const vector3d<T> &X );

  template <class T>
    inline tensor3d<T> transpose( const tensor3d<T> &A ){
    return tensor3d<T>{A(0,0), A(0,1), A(0,2), A(1,0), A(1,1), A(1,2), A(2,0), A(2,1), A(2,2)};

  }

  /// Start implementation here.

  template <class T>
    inline vector3d<T> externprod(const vector3d<T> &a, const vector3d<T> &b ){
    return vector3d<T> {
      a(1)*b(2) - a(2)*b(1),
	- a(0)*b(2) + a(2)*b(0),
	a(0)*b(1) - a(1)*b(0)
	};
  }

  template <class T>
    inline T dot(const vector3d<T> &a, const vector3d<T> &b ){
    return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
  }


  template <class T>
    inline vector3d<T> & operator+=( vector3d<T> & a, const vector3d<T> &b){
    a(0)+=b(0);  a(1)+=b(1);  a(2)+=b(2);
    return a;
  };

  template <class T>
    inline vector3d<T> & operator-=( vector3d<T> & a, const vector3d<T> &b){
    a(0)-=b(0);  a(1)-=b(1);  a(2)-=b(2);
    return a;
  };

  template <class T>
    inline vector3d<T> operator+(const vector3d<T> & a, const vector3d<T> &b){
    vector3d<T> c(a);
    return c+=b;
  };

  template <class T>
    inline vector3d<T> operator-(const vector3d<T> & a, const vector3d<T> &b){
    vector3d<T> c(a);
    return c-=b;
  };

  template <class T>
    inline void axpy(T alpha, const vector3d<T> & a,  vector3d<T> &b){
    b(0) += alpha*a(0);
    b(1) += alpha*a(1);
    b(2) += alpha*a(2);
  };

  template <class T>
    inline std::ostream& operator<<(std::ostream& os, const vector3d<T>& v)
    {
      os << v(0) << " " << v(1) << " " << v(2);
      return os;
    }

  // tensor2d implem start here
  template <class T>
    std::ostream& operator<<(std::ostream& os, const tensor2d<T>& A){
    os << A(0,0) << " " << A(0,1) << std::endl;
    os << A(1,0) << " " << A(1,1)  << std::endl;
    return os;
  }


  template <class T>
    T det(const tensor2d<T> & A)
    {
      return  A(0,0)*A(1,1) -A(0,1) *A(1,0);
    }

  template <class T>
    tensor2d<T> invert(const tensor2d<T> & A ){
    const T detA = det(A);
    if (detA == 0.0){
      std::cout << "tensor2d<T> A not invertible " << __FILE__ << __LINE__ << std::endl;
      throw -1;
    }
    const T idetA = 1./detA;
    return tensor2d<T>{
      idetA    *A(1,1), //invA(0,0)
	-idetA *A(1,0), //invA(1,0)
	-idetA *A(0,1), //invA(0,1)
	idetA  *A(0,0)};//invA(1,1)
  }

  template <class T>
    inline void gemm( T alpha, const tensor2d<T> &A, const tensor2d<T> &B, T  beta,  tensor2d<T> &C){
    C(0,0) = beta*C(0,0) + alpha*(A(0,0)*B(0,0) + A(0,1)*B(1,0));
    C(1,0) = beta*C(1,0) + alpha*(A(1,0)*B(0,0) + A(1,1)*B(1,0));
    C(0,1) = beta*C(0,1) + alpha*(A(0,0)*B(0,1) + A(0,1)*B(1,1));
    C(1,1) = beta*C(1,1) + alpha*(A(1,0)*B(0,1) + A(1,1)*B(1,1));
  }

  template <class T>
    tensor2d<T> operator*( const tensor2d<T> &A, const tensor2d<T> &B){
    tensor2d<T> C;
    gemm(1.,A, B, 0.,C);
    return C;
  }


  
  template <class T>
    inline void gemv( T alpha, const tensor2d<T> &A, const vector2d<T> &X, T beta,  vector2d<T> &Y){
    Y(0) = beta*Y(0) + alpha*(A(0,0)*X(0) + A(0,1)*X(1));
    Y(1) = beta*Y(1) + alpha*(A(1,0)*X(0) + A(1,1)*X(1));
  }

  template <class T>
    inline vector2d<T> operator*( const tensor2d<T> &A, const vector2d<T> &X ){
    vector2d<T> Y;
    gemv(1.,A, X, 0., Y);
    return Y;
  }

  template <class T>
    inline void solve( const tensor2d<T> &A, const vector2d<T> &b,  vector2d<T> &x){
    const double detA = det(A);
    if (detA == 0.0){
      throw -1;
    }
    const double idetA = 1. / (detA);
    x(0) = idetA*( A(1,1)*b(0) - A(0,1)*b(1));
    x(1) = idetA*(-A(1,0)*b(0) + A(0,0)*b(1));
  }


  /// tensor3d implementation start here :
  template <class T>
    std::ostream& operator<<(std::ostream& os, const tensor3d<T>& A){
    os << A(0,0) << " " << A(0,1) << " " << A(0,2)<< " " << std::endl;
    os << A(1,0) << " " << A(1,1) << " " << A(1,2)<< " " << std::endl;
    os << A(2,0) << " " << A(2,1) << " " << A(2,2)<< " " << std::endl;
    return os;
  }


  template <class T>
    T det(const tensor3d<T> & A)
    {
      return  A(0,0)*(A(1,1)*A(2,2) - A(1,2) *A(2,1))
	-A(0,1)*(A(1,0) * A(2,2) - A(1,2) * A(2,0))
	+A(0,2)*(A(1,0) * A(2,1) - A(1,1) * A(2,0));
    }

  template <class T>
    tensor3d<T> invert(const tensor3d<T> & A ){
    const T detA = det(A);
    if (detA == 0.0){
      std::cout << "tensor3d<T> A not invertible " << __FILE__ << __LINE__ << std::endl;
      throw -1;
    }
    /*
      const T idetA = 1./detA;
      return tensor3d<T>{
      idetA * (A(1,1) * A(2,2) - A(1,2) * A(2,1)),
      -idetA * (A(1,0) * A(2,2) - A(1,2) * A(2,0)),
      idetA * (A(1,0) * A(2,1) - A(1,1) * A(2,0)),
      -idetA * (A(0,1) * A(2,2) - A(0,2) * A(2,1)),
      idetA * (A(0,0) * A(2,2) - A(0,2) * A(2,0)),
      -idetA * (A(0,0) * A(2,1) - A(0,1) * A(2,0)),
      idetA * (A(0,1) * A(1,2) - A(0,2) * A(1,1)),
      -idetA * (A(0,0) * A(1,2) - A(0,2) * A(1,0)),
      idetA * (A(0,0) * A(1,1) - A(0,1) * A(1,0))};	
    */
    return tensor3d<T>{
      (A(1,1) * A(2,2) - A(1,2) * A(2,1))/detA,
	-(A(1,0) * A(2,2) - A(1,2) * A(2,0))/detA,
	(A(1,0) * A(2,1) - A(1,1) * A(2,0))/detA,
	-(A(0,1) * A(2,2) - A(0,2) * A(2,1))/detA,
	(A(0,0) * A(2,2) - A(0,2) * A(2,0))/detA,
	-(A(0,0) * A(2,1) - A(0,1) * A(2,0))/detA,
	(A(0,1) * A(1,2) - A(0,2) * A(1,1))/detA,
	-(A(0,0) * A(1,2) - A(0,2) * A(1,0))/detA,
	(A(0,0) * A(1,1) - A(0,1) * A(1,0))/detA};
  	
  }

  template <class T>
    tensor3d<T> dual(const tensor3d<T> & A ){
    const vector3d<T> a = {A(0,0), A(1,0), A(2,0)};
    const vector3d<T> b = {A(0,1), A(1,1), A(2,1)};
    const vector3d<T> c = {A(0,2), A(1,2), A(2,2)};
    const vector3d<T> tmpda = externprod(b, c);
    const vector3d<T> tmpdb = externprod(c, a);
    const vector3d<T> tmpdc = externprod(a, b);
    const T scala = dot(a,tmpda);
    const T scalb = dot(tmpdb, b);
    const T scalc = dot(tmpdc, c);
    return tensor3d<T>{
      tmpda(0)/scala, tmpda(1)/scala, tmpda(2)/scala,
	tmpdb(0)/scalb, tmpdb(1)/scalb, tmpdb(2)/scalb,
	tmpdc(0)/scalc, tmpdc(1)/scalc, tmpdc(2)/scalc
	};
  }


  template <class T>
    tensor3d<T> FTF(const tensor3d<T> & F ){
    const T C00 =  F(0,0)*F(0,0) +  F(1,0)*F(1,0) +  F(2,0)*F(2,0);
    const T C01 =  F(0,1)*F(0,0) +  F(1,1)*F(1,0) +  F(2,1)*F(2,0);
    const T C02 =  F(0,2)*F(0,0) +  F(1,2)*F(1,0) +  F(2,2)*F(2,0);
    const T C11 =  F(0,1)*F(0,1) +  F(1,1)*F(1,1) +  F(2,1)*F(2,1);
    const T C12 =  F(0,1)*F(0,2) +  F(1,1)*F(1,2) +  F(2,1)*F(2,2);
    const T C22 =  F(0,2)*F(0,2) +  F(1,2)*F(1,2) +  F(2,2)*F(2,2);
    return tensor3d<T>{
      C00, C01, C02,
	C01, C11, C12,
	C02, C12, C22};
  }

  template <class T>
    T sum(const tensor3d<T> & F ){
    return std::accumulate(F.data.begin(), F.data.end(), 0. );
  }


  template <class T>
    inline void gemm( T alpha, const tensor3d<T> &A, const tensor3d<T> &B, T  beta,  tensor3d<T> &C){
    C(0,0) = beta*C(0,0) + alpha*(A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0));
    C(1,0) = beta*C(1,0) + alpha*(A(1,0)*B(0,0) + A(1,1)*B(1,0) + A(1,2)*B(2,0));
    C(2,0) = beta*C(2,0) + alpha*(A(2,0)*B(0,0) + A(2,1)*B(1,0) + A(2,2)*B(2,0));
    C(0,1) = beta*C(0,1) + alpha*(A(0,0)*B(0,1) + A(0,1)*B(1,1) + A(0,2)*B(2,1));
    C(1,1) = beta*C(1,1) + alpha*(A(1,0)*B(0,1) + A(1,1)*B(1,1) + A(1,2)*B(2,1));
    C(2,1) = beta*C(2,1) + alpha*(A(2,0)*B(0,1) + A(2,1)*B(1,1) + A(2,2)*B(2,1));
    C(0,2) = beta*C(0,2) + alpha*(A(0,0)*B(0,2) + A(0,1)*B(1,2) + A(0,2)*B(2,2));
    C(1,2) = beta*C(1,2) + alpha*(A(1,0)*B(0,2) + A(1,1)*B(1,2) + A(1,2)*B(2,2));
    C(2,2) = beta*C(2,2) + alpha*(A(2,0)*B(0,2) + A(2,1)*B(1,2) + A(2,2)*B(2,2));
  }

  template <class T>
    tensor3d<T> operator*( const tensor3d<T> &A, const tensor3d<T> &B){
    tensor3d<T> C;
    gemm(1.,A, B, 0.,C);
    return C;
  }


  
  template <class T>
    inline void gemv( T alpha, const tensor3d<T> &A, const vector3d<T> &X, T beta,  vector3d<T> &Y){
    Y(0) = beta*Y(0) + alpha*(A(0,0)*X(0) + A(0,1)*X(1) + A(0,2)*X(2));
    Y(1) = beta*Y(1) + alpha*(A(1,0)*X(0) + A(1,1)*X(1) + A(1,2)*X(2));
    Y(2) = beta*Y(2) + alpha*(A(2,0)*X(0) + A(2,1)*X(1) + A(2,2)*X(2));
  }

  template <class T>
    inline vector3d<T> operator*( const tensor3d<T> &A, const vector3d<T> &X ){
    return vector3d<T>  {
      A(0,0)*X(0) + A(0,1)*X(1) + A(0,2)*X(2),
	A(1,0)*X(0) + A(1,1)*X(1) + A(1,2)*X(2),
	A(2,0)*X(0) + A(2,1)*X(1) + A(2,2)*X(2)
	}; 
  }

  template <class T>
    inline vector3d<T> operator*( const vector3d<T> &X , const tensor3d<T> &A){
    return vector3d<T>  {
      A(0,0)*X(0) + A(1,0)*X(1) + A(2,0)*X(2),
	A(0,1)*X(0) + A(1,1)*X(1) + A(2,1)*X(2),
	A(0,2)*X(0) + A(1,2)*X(1) + A(2,2)*X(2)
	};
  }


  template <class T>
    inline void solve( const tensor3d<T> &A, const vector3d<T> &b,  vector3d<T> &x){
    const double detA = det(A);
  
    if (detA == 0.0){
      throw -1;
    }
  
    //const double idetA = 1. / (detA);
    x(0) = (b(0) * (A(1,1) * A(2,2) - A(1,2) * A(2,1)) -
	    A(0,1) * (b(1) * A(2,2) - A(1,2) * b(2)) +
	    A(0,2) * (b(1) * A(2,1) - A(1,1) * b(2)))/detA;
  
    x(1) = (A(0,0) * (b(1) * A(2,2) - A(1,2) * b(2)) -
	    b(0) * (A(1,0) * A(2,2) - A(1,2) * A(2,0)) +
	    A(0,2) * (A(1,0) * b(2) - b(1) * A(2,0)))/detA;
  
    x(2) = (A(0,0) * (A(1,1) * b(2) - b(1) * A(2,1)) -
	    A(0,1) * (A(1,0) * b(2) - b(1) * A(2,0)) +
	    b(0) * (A(1,0) * A(2,1) - A(1,1) * A(2,0)))/detA;
  }



  template <class T >
    inline std::array<vector3d<T >, 3 > dual(const std::array<vector3d<T >, 3 > & in){
    std::array<vector3d<T >, 3 > out;
    for (int i=0; i < 3; ++i){
      const int& ip0 = i;
      const int ip1 = (i+1)%3;
      const int ip2 = (i+2)%3;
      vector3d<T> tmp = externprod(in[ip1], in[ip2]);
      const T scalp = dot (in[ip0], tmp); 
      assert (fabs(scalp) > 0.);
      const T iscalp = 1./scalp;
      out[ip0] = scal(iscalp,tmp);
    }
    return out;
  } 

  template <class T >
    inline std::array<vector3d<T>, 2 > dual(const std::array<vector3d<T>, 2 > & in ){
    std::array<vector3d<T>, 2 > out;
    vector3d<T> g3 = externprod(in[0], in[1]);
    const double nrm2g3 = nrm2(g3);
    assert (fabs(nrm2g3) > 0.);
    g3 = scal(1./nrm2g3, g3);
    vector3d<T> tmp = externprod(in[1], g3);
    T scalp = dot (in[0], tmp);
    assert (fabs(scalp) > 0.);
    out[0] = scal(((T)1.)/scalp,tmp);

    tmp = externprod(g3, in[0]);
    scalp = dot (in[1], tmp);
    assert (fabs(scalp) > 0.);
    out[1] = scal(((T)1.)/scalp,tmp);
    return out;
  }

  template <class T >
    inline tensor3d2d<T > dual(const tensor3d2d<T>  & in ){
    const vector3d<T> a = {in(0,0), in(1,0), in(2,0)};
    const vector3d<T> b = {in(0,1), in(1,1), in(2,1)};
    const vector3d<T> g3 = externprod(a,b);
    const vector3d<T> tmpda = externprod(b, g3);
    const vector3d<T> tmpdb = externprod(g3, a);
    const T scala = dot(a,tmpda);
    const T scalb = dot(tmpdb, b);
    return tensor3d2d<T>{
      tmpda(0)/scala, tmpda(1)/scala, tmpda(2)/scala,
	tmpdb(0)/scalb, tmpdb(1)/scalb, tmpdb(2)/scalb
	};
  };


  

  template <class T >
    std::array<vector2d<T >, 2 > dual(const std::array<vector2d<T >, 2 > & in ){
    std::array<vector2d<T >, 2 > out;
  
    const vector2d<T > & a= in[0];
    const vector2d<T > & b= in[1];
    const T D= (a(0)*b(1) - a(1)*b(0)); 
    assert (fabs(D)> 0.);
    const T invD = 1./D;
    out[0] = {b(1)*invD, - b(0)*invD};
    out[1] = {-a(1)*invD, a(0)*invD};
    return out;
  }

  template <class T >
    vector2d<T > dual(const vector2d<T > & in ){
    const T scalp = dot(in,in);
    assert (fabs(scalp)  > 0.);
    return vector2d<T >(in(0)/scalp, in(1)/scalp);
  }

  template <class T >
    std::array<vector2d<T >, 1 >  dual(const std::array<vector2d<T >, 1 > & in){
    std::array<vector2d<T >, 1 >  out = {dual(in[0])};
    return out;
  }

  template <class T >
    vector3d<T >   dual(const vector3d<T > & in ){
    const T scalp = dot(in,in); 
    assert (fabs(scalp)  > 0.);
    return vector3d<T >(in(0)/scalp, in(1)/scalp, in(2)/scalp);
  }

  template <class T >
    std::array<vector3d<T >, 1 > dual(const std::array<vector3d<T >, 1 > & in){
    std::array<vector3d<T >, 1 >  out = {dual(in[0])};
    return out;
  }

} //end of namespace
#endif
