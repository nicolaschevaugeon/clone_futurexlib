/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _GENERICOPERATION_
#define _GENERICOPERATION_
#include <cmath>
namespace xfem{
  template< typename T >
    struct xtool::xIdentity : public std::unary_function< T, T >
    {
    public:
      typedef std::unary_function< T, T >               base_type;
      typedef typename base_type::argument_type         argument_type;
      typedef typename base_type::result_type           result_type;
      
      T& operator()( T& t ) const { return t; }
      const T& operator()( const T& t ) const { return t; }
    };
  
  //T1 : Arg1type , T2: Arg2 type,   T3 : Result type
  template <class T1, class T2, class T3>
    class xtool::xMult : public std::binary_function<T1, T2, T3>
    {
    public:
      typedef typename std::binary_function<T1, T2, T3>::result_type result_type;
      typedef typename std::binary_function<T1, T2, T3>::first_argument_type first_argument_type;
      typedef typename std::binary_function<T1, T2, T3>::second_argument_type second_argument_type;
      
      result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f * s;  }
    };
  
 template <class T1, class T2, class T3>
    class xtool::xAdd : public std::binary_function<T1, T2, T3>
    {
    public:
      typedef typename std::binary_function<T1, T2, T3>::result_type result_type;
      typedef typename std::binary_function<T1, T2, T3>::first_argument_type first_argument_type;
      typedef typename std::binary_function<T1, T2, T3>::second_argument_type second_argument_type;
    xtool::xAdd(double _coeff=1.0):coeff(_coeff)
      {}
      result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f+s*coeff; }
      
    private :
      double coeff;
      
    };

 template <class T1, class T2, class T3>
    class xMod : public std::binary_function<T1, T2, T3>
    {
    public:
      typedef typename std::binary_function<T1, T2, T3>::result_type result_type;
      typedef typename std::binary_function<T1, T2, T3>::first_argument_type first_argument_type;
      typedef typename std::binary_function<T1, T2, T3>::second_argument_type second_argument_type;

      result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f % s; }
    };

  template <class T> 
class xtool::xScale : public std::unary_function< T, T>
  {
  public:
  xtool::xScale(double _a):a(_a){};
    T operator () ( T const &arg ) const
    {
      return arg*a;
    }
  private: double a;
  };


template < typename T > struct xAbs : std::unary_function<T, T>
{
T operator () ( T const &arg ) const
{
using std::abs;
return abs( arg );
}    
};

template < typename T > struct xChangeSign : std::unary_function<T, T>
{
T operator () ( T const &arg ) const
{
return -arg;
}    
};

template < typename T > struct xNegativePart : std::unary_function<T, T>
{
T operator () ( T const &arg ) const
{
using std::abs;
 return 0.5 * (arg - abs( arg ));
}    
};

template < typename T > struct xtool::xPositivePart : std::unary_function<T, T>
{
T operator () ( T const &arg ) const
{
using std::abs;
 return 0.5 * (arg + abs( arg ));
}    
};

struct xSqrt : std::unary_function<double, double>
{
double operator () ( double const &arg ) const
{
return sqrt( arg );
}    
};

struct xSquare : std::unary_function<double, double>
{
double operator () ( double const &arg ) const
{
return arg*arg;
}    
};

struct xInv : std::unary_function<double, double>
{
double operator () ( double const &arg ) const
{
return  1/arg;
}    
};
  
}


#endif
