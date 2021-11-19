/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __SOL_VISITOR__H
#define __SOL_VISITOR__H

#include "xValue.h"

namespace xlinalg
{
    class xCSRVector;
}

namespace xfem
{

// class xWriteSolutionVisitor  
// {
//  public:
//   xWriteSolutionVisitor(xlinalg::xCSRVector::const_iterator s) : current(s) {}
//   xValue<double>* operator()(xValue<double>* v) 
//   { v->setVal(*(current++)); return v; }
//  private:
//    xlinalg::xCSRVector::const_iterator current;
// };

  template <class T = double>
  class xWriteSolutionVisitorFromArray
  {
  public:
    xWriteSolutionVisitorFromArray(const T *s) : current(s) {}
    xValue<T>* operator()(xValue<T>* v)
    { v->setVal(*(current++)); return v; }
  private:
    const double *current;
  };

template <class VEC = xlinalg::xCSRVector>
class xWriteSolutionVisitor
{
 public:
  xWriteSolutionVisitor(typename VEC::const_iterator s) : current(s) {}
  xWriteSolutionVisitor(const VEC &V) : current(V.begin()) {}

  template <typename VT>
  xValue<VT>* operator()(xValue<VT>* v)
  { v->setVal(*(current++)); return v; }
 private:
   typename VEC::const_iterator  current;
};

template <class VEC = xlinalg::xCSRVector>
class xWriteScaledSolutionVisitor
{
 public:
  xWriteScaledSolutionVisitor(typename VEC::const_iterator s, double a) : current(s), fac(a) {}
  xValue<double>* operator()(xValue<double>* v)
  { v->setVal((*(current++))*fac); return v; }
 private:
   typename VEC::const_iterator  current;
  double fac;
};

template <class VEC = xlinalg::xCSRVector>
class xReadVectorVisitor 
{
 public:
  xReadVectorVisitor(typename VEC::iterator s) : current(s) {}
  xReadVectorVisitor(VEC &V) : current(V.begin()) {}
  xValue<double>*  operator()(xValue<double>* v) 
  { *(current++) = v->getVal(); return v; }
 private:
   typename VEC::iterator current;
};

template <class VEC = xlinalg::xCSRVector>
class xAddSolutionVisitor 
{
 public:
  xAddSolutionVisitor(typename VEC::const_iterator s) : current(s) {}
  xAddSolutionVisitor(const VEC &V) : current(V.begin()) {}
  xValue<double>* operator()(xValue<double>* v) 
  { v->setVal(v->getVal()+*(current++)); return v; }
 private:
   typename VEC::const_iterator current;
};

template <class VEC = xlinalg::xCSRVector>
class xAddScaledSolutionVisitor 
// dev. SLC for linear search in Newton-Raphson context,
// may be useful for other purpose.
{
 public:
  xAddScaledSolutionVisitor(typename VEC::const_iterator s, double a) : current(s), fac(a) {}
  xValue<double>* operator()(xValue<double>* v) 
  { v->setVal(v->getVal()+(*(current++))*fac); return v; }
 private:
   typename VEC::const_iterator current;
  double fac;
};

} // end of namespace

#endif
