/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef  _xDispatcher_H
#define  _xDispatcher_H

#include <exception>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/begin.hpp>
#include <boost/mpl/next.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/back.hpp>

#include "xStateOfValue.h"


namespace xfem  {

struct bad_dispatch : std::exception {};

template <typename Executor,
          typename Sequence1,
          typename Sequence2 = Sequence1,
          typename BaseType1 = typename boost::mpl::back <Sequence1>::type,
          typename BaseType2 = typename boost::mpl::back <Sequence2>::type>
class xDispatcher
{
  public:
    typedef typename Executor::result_type   result_type;

    xDispatcher(const Executor& executor_ = Executor())
    : executor(executor_) {}
    
    result_type Dispatch(BaseType1& p1, BaseType2& p2)
    {
      FindType1(p1, p2, Begin1());
    }

  private:
    Executor executor;

    typedef typename Sequence1::type                 Types1;
    typedef typename boost::mpl::begin<Types1>::type Begin1;
    typedef typename boost::mpl::end  <Types1>::type End1;

    typedef typename Sequence2::type                 Types2;
    typedef typename boost::mpl::begin<Types2>::type Begin2;
    typedef typename boost::mpl::end  <Types2>::type End2;

    template<typename iterator>
    result_type FindType1(BaseType1& p1, BaseType2& p2, iterator)
    {
      typedef typename boost::mpl::deref<iterator>::type Type1;
      typedef typename boost::mpl::next<iterator>::type type;
      if (Type1* P1 = dynamic_cast<Type1*>(&p1)) {
          return FindType2<Type1>(*P1, p2, Begin2());
      } else {
          return FindType1(p1, p2, type());
      }
    }

    result_type FindType1(BaseType1& p1, BaseType2& p2, End1)
    {
      throw bad_dispatch();
    }

    template<typename Type1, typename iterator>
    result_type FindType2(Type1& P1, BaseType2& p2, iterator)
    {
      typedef typename boost::mpl::deref<iterator>::type Type2;
      typedef typename boost::mpl::next<iterator>::type type;
      if (Type2* P2 = dynamic_cast<Type2*>(&p2)) {
          return executor(P1, *P2);
      } else {
          return FindType2(P1, p2, type());
      }
    }

    template<typename Type1>
    result_type FindType2(Type1& p1, BaseType2& p2, End2)
    {
      throw bad_dispatch();
    }
};

////////////////////////////////////////////////////////: 

 template <typename Executor >
  class xDispatcherStateOfValue
{
  public:
    typedef typename Executor::result_type   result_type;

    xDispatcherStateOfValue(const Executor& executor_ = Executor())
    : executor(executor_) {}
    
    result_type Dispatch(const xStateOfValue& p1, const xStateOfValue& p2)
    {
      switch (p1.state){
      case (xStateOfValue::DOF) : {
	const xStateOfValueDof *P1  = static_cast< const xStateOfValueDof * > (&p1);
	return FindType2< xStateOfValueDof >(*P1, p2 );
      }
      case xStateOfValue::NONE :{
	const xStateOfValueNone *P1  = static_cast< const xStateOfValueNone * > (&p1);
	return FindType2< xStateOfValueNone >(*P1, p2 );
      }
      case xStateOfValue::FIXED :{
    const xStateOfValueFixed<typename Executor::value_t> *P1  = static_cast< const xStateOfValueFixed<typename Executor::value_t> * > (&p1);
    return FindType2< xStateOfValueFixed<typename Executor::value_t> >(*P1, p2 );
      }
      case xStateOfValue::LINEARCOMBINATION :{
    const xStateOfValueLinearCombination<typename Executor::value_t> *P1  = static_cast< const xStateOfValueLinearCombination<typename Executor::value_t> * > (&p1);
    return FindType2< xStateOfValueLinearCombination<typename Executor::value_t> >(*P1, p2 );
      }
      default :	throw bad_dispatch();
      }
    }

  private:
    Executor executor;
 
    template<typename Type1 >
      result_type FindType2(const Type1& P1, const xStateOfValue& p2) {
      switch (p2.state){
      case xStateOfValue::DOF :{
	const xStateOfValueDof *P2  = static_cast< const xStateOfValueDof * > (&p2);
	return executor(P1, *P2);
      }
      case xStateOfValue::NONE :{
	const xStateOfValueNone *P2  = static_cast< const xStateOfValueNone * > (&p2);
	return  executor(P1, *P2);
      }
      case xStateOfValue::FIXED :{
    const xStateOfValueFixed<typename Executor::value_t> *P2  = static_cast< const xStateOfValueFixed<typename Executor::value_t> * > (&p2);
	return  executor(P1, *P2);
      }
      case xStateOfValue::LINEARCOMBINATION :{
    const xStateOfValueLinearCombination<typename Executor::value_t> *P2  = static_cast< const xStateOfValueLinearCombination<typename Executor::value_t> * > (&p2);
	return  executor(P1, *P2);
      }
      default :	throw bad_dispatch();
      }
    }
};
} // namespace xfem

#endif
