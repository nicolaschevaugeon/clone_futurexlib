/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _STATE_OF_VALUEH
#define _STATE_OF_VALUEH
#include <iostream>
#include <string>

#include "xValKey.h"
#include "xValue.h"
#include "xValueManager.h"
#include "xDebug.h"

namespace xfem
{


class xStateOfValueDof : public xStateOfValue {
  public:
  xStateOfValueDof(int d) :xStateOfValue(DOF),  Numdof(d) {}
  std::ostream& print(std::ostream& o) const override  { o << " xStateOfValueDof NumDof is " << Numdof << " " << std::endl; return o; }
  int Numdof;
};


class xStateOfValueNone : public xStateOfValue 
{ 
 public :
 xStateOfValueNone() :xStateOfValue(NONE){};
 private:
  std::ostream& print(std::ostream& o) const override;
};


template<typename VT = double>
class xStateOfValueFixed : public xStateOfValue {
  public:
    xStateOfValueFixed(const xValue<VT>* v);
    VT getVal() const;
    std::ostream& print(std::ostream& o) const override;
  private:
    const xValue<VT>* val;
};

template<typename T>
class  xValueLinearCombinationBase;

template <typename T>
class xStateOfValueLinearCombination : public xStateOfValue {
public:
  typedef T coef_t;
  xStateOfValueLinearCombination(const xValueLinearCombinationBase<T>* l);
  size_t size() const;
  coef_t          coeff(int i) const;
  xStateOfValue* state(int i)  const;
  std::ostream& print(std::ostream& o) const override;
private:
  const xValueLinearCombinationBase<T>* lin;
};


class xTrue 
{
public:
 static  bool test(const xValKey& k)
    {
      return true; 
    } 
};
template <class Condition = xTrue, typename VT=double>
class xStateDofCreator : public Condition
{
public:
  xStateDofCreator(xValueManagerDist<VT>& val, const std::string& sub) : val_manager(&val), subset(sub) {}
  xValue<VT>* operator()(const xValKey& key, xValue<VT>* v)
    { 
      const bool debug = xdebug_flag;
      if(v->getState()) return v; 
      if (debug) std::cout << "test is " << this->test(key) << std::endl;
      if (!this->test(key)) return nullptr;
      int numdof = val_manager->size(subset) + 1;
      if (debug) 
	{
	  std::cout << "numdof in create " << numdof << " for " << std::endl;
	  key.getEnti()->print();
	}
      val_manager->add(v, subset);
      v->setState(new xStateOfValueDof(numdof)); 
      return v;
    }
private:
  xValueManagerDist<VT>* val_manager;
  std::string subset;
};

template<typename VT = double>
class xStateFixedCreator
{
public:
  xValue<VT>* operator()(const xValKey& k, xValue<VT>* v)
    { 
      if(v->getState()) return v; 
      v->setState(new xStateOfValueFixed<VT>(v));
      return v;
    }
private:
};

/*!
 * A function that clone given state and put it in given value v
 */
template <typename S, typename T>
void xCloneState(const S *state, xValue<T> * v);

} // end of namespace


#include "xStateOfValue_imp.h"
#endif
