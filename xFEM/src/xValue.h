/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _VALUE_H
#define _VALUE_H
// std
#include <algorithm>
#include <iosfwd>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
// mpi
#include "mpi.h"
// xtool
#include "xDataType.h"
// xfem
#include "xValKey.h"

namespace xfem
{
//! State of Value base class
//!
class xStateOfValue
{
  public:
   enum state_type
   {
      DOF,
      NONE,
      FIXED,
      LINEARCOMBINATION,
      UNDEFINED
   };
   xStateOfValue() : state(UNDEFINED) {}
   xStateOfValue(const state_type& _state) : state(_state) {}
   virtual ~xStateOfValue() = default;
   virtual std::ostream& print(std::ostream& o) const = 0;
   state_type state;
};

//! Value base class
//!
template <class T>
class xValue
{
  public:
   xValue();
   virtual ~xValue();
   virtual T getVal() const = 0;
   virtual void setVal(T in) = 0;
   virtual std::ostream& print(std::ostream& o) const;
   virtual std::ostream& printVal(std::ostream& o) const = 0;
   virtual std::ostream& printState(std::ostream& o) const;
   virtual const xStateOfValue* getState() const;
   virtual xStateOfValue* getState();
   virtual void setState(xStateOfValue* s);
   virtual void delState();
   virtual xValue<T>* getValPtr();

  protected:
   mutable xStateOfValue* state;

  private:
};

template <class T>
xValue<T>::xValue() : state(nullptr)
{
}
template <class T>
xValue<T>::~xValue()
{
   delete state;
}
template <class T>
std::ostream& xValue<T>::print(std::ostream& o) const
{
   printVal(o);
   o << std::endl;
   printState(o);
   return o;
}
template <class T>
std::ostream& xValue<T>::printState(std::ostream& o) const
{
   if (state) return state->print(o);
   return o;
}
template <class T>
const xStateOfValue* xValue<T>::getState() const
{
   return state;
}
template <class T>
xStateOfValue* xValue<T>::getState()
{
   return state;
}
template <class T>
void xValue<T>::setState(xStateOfValue* s)
{
   if (state) delState();
   state = s;
}
template <class T>
void xValue<T>::delState()
{
   if (state)
   {
      delete state;
      state = nullptr;
   }
}
template <class T>
xValue<T>* xValue<T>::getValPtr()
{
   return this;
}

template <typename VT>
class xSingleValue : public xValue<VT>
{
  public:
   xSingleValue() {}
   VT getVal() const override { return value; }
   void setVal(VT in) override { value = in; }
   std::ostream& printVal(std::ostream& o) const override
   {
      o << value;
      return o;
   }

  private:
   VT value{xtool::xDataType<VT>::zero()};
};

using xValueDouble = xSingleValue<double>;
using xValueFloat = xSingleValue<float>;
using xValueDoubleComplex = xSingleValue<std::complex<double>>;

template <typename T = double>
class xNValue : public xValue<T>
{
  public:
   xNValue(size_t N) : value(N, xtool::xDataType<T>::zero()) {}
   static void set_current(size_t i) { current = i; }
   T getVal() const override { return value[current]; }
   void setVal(T in) override { value[current] = in; }
   std::ostream& printVal(std::ostream& o) const override
   {
      o << value[current];
      return o;
   }

  private:
   std::vector<T> value;
   static size_t current;
};

template <typename T>
size_t xNValue<T>::current = 0;

class xValueError : public xValue<double>
{
  public:
   xValueError();
   double getVal() const override;
   void setVal(double in) override;
   std::ostream& printVal(std::ostream& o) const override { return o; }
   static void choice(const std::string& s);
   static double total(const std::string& s);
   static void finalize(MPI_Comm world);
   static void finalizeAgainABS(MPI_Comm world);
   static void clear();

  private:
   typedef enum
   {
      ENG,
      ABS2,
      ABS,
      REL,
      ENG_EXA,
      ABS2_EXA,
      ABS_EXA,
      REL_EXA,
      EFF_EXA,
      EFF,
      DNS,
      VOL
   } choice_t;
   static choice_t _choice;
   double abs2_loc;
   double abs2_exa_loc;
   double eng_loc;
   double eng_exa_loc;
   double vol_loc;
   static double abs2_tot;
   static double abs2_exa_tot;
   static double eng_tot;
   static double eng_exa_tot;
   static double vol_tot;
   static bool reduced;
};

/*!
 * A function that clone given value and return its copy
 * Cloning is the act of creating a equivalent value to the given
 * one that rely in its own memory space.
 */
template <typename V, typename B>
B* xCloneValue(const V* val, const std::map<B*, B*>& coresp);

}  // namespace xfem

#include "xValue_imp.h"
#endif
