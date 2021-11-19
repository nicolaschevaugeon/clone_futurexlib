/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _GENERICOPERATION_
#define _GENERICOPERATION_

#include <cmath>

namespace xtool{

  template< typename T >
    struct xIdentity : public std::unary_function< T, T >
    {
    public:
      typedef std::unary_function< T, T >               base_type;
      typedef typename base_type::argument_type         argument_type;
      typedef typename base_type::result_type           result_type;
      
      T& operator()( T& t ) const { return t; }
      const T& operator()( const T& t ) const { return t; }
    };


    template< typename T >
      struct xCastToDouble : public std::unary_function< T, double >
      {
      public:
        typedef std::unary_function< T, double >               base_type;
        typedef typename base_type::argument_type         argument_type;
        typedef typename base_type::result_type           result_type;

        result_type& operator()( argument_type& t ) const { return static_cast<result_type&>(t); }
        const result_type& operator()( const argument_type& t ) const { return static_cast<const result_type&>(t); }
      };
  


  //T1 : Arg1type , T2: Arg2 type,   T3 : Result type
  template <class T1, class T2, class T3>
    class xMult : public std::binary_function<T1, T2, T3>
    {
    public:
      typedef typename std::binary_function<T1, T2, T3>::result_type result_type;
      typedef typename std::binary_function<T1, T2, T3>::first_argument_type first_argument_type;
      typedef typename std::binary_function<T1, T2, T3>::second_argument_type second_argument_type;
      
      result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f * s;  }
    };
   

 


  template <class T1, class T2, class T3>
    class xAdd : public std::binary_function<T1, T2, T3>
    {
    public:
      typedef typename std::binary_function<T1, T2, T3>::result_type result_type;
      typedef typename std::binary_function<T1, T2, T3>::first_argument_type first_argument_type;
      typedef typename std::binary_function<T1, T2, T3>::second_argument_type second_argument_type;
    xAdd(double _coeff=1.0):coeff(_coeff)
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
    class xScale : public std::unary_function< T, T>
    {
    public:
    xScale(double _a):a(_a){};
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



  template < typename T > struct xPositivePart : std::unary_function<T, T>
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
