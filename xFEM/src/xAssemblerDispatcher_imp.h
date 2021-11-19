/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _ASSEMBLER_DISPATCHER_IMP_H
#define _ASSEMBLER_DISPATCHER_IMP_H

#include <iostream>
#include <cassert>
#include "xDispatcher.h"
#include "xAssemblerDispatcher.h"
#include "xStateOfValue.h"
#include "xDebug.h"

namespace xfem {

template< typename M, typename V, typename S, typename T >
xMatVecScaAssembler<M,V,S,T>::xExecutor::xExecutor(T val_, S* sca_, V* vec_, M* mat_)
: val (val_), sca (sca_), vec (vec_), mat (mat_) {}

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueDof& dofi,const xStateOfValueDof& dofj) const
{
#ifndef NDEBUG
  const bool debug = xdebug_flag;
  if (debug) std::cout << "xStateOfValueDof& dofi,xStateOfValueDof& dofj detected" << std::endl;
#endif
  if (mat)
  {
     int i=dofi.Numdof;
     int j=dofj.Numdof;
     xPolicyAssemblerDispatcher<typename M::matrix_pattern>::reorder(i,j);
     mat->AddMatrix(i, j, val);
  }
}

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueDof& equ, const xStateOfValueFixed<S>& f) const
{
#ifndef NDEBUG
  const bool debug = xdebug_flag;
  if (debug) std::cout << "xStateOfValueDof&  equ, xStateOfValueFixed<>& f detected" << f.getVal() << std::endl;
#endif
  if (vec) vec->AddVal(equ.Numdof, -val * f.getVal() );
}

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueFixed<S>& f1, const xStateOfValueFixed<S>& f2) const
{
#ifndef NDEBUG
  const bool debug = xdebug_flag;
  if (debug) std::cout << "xStateOfValueFixed<>& f1, xStateOfValueFixed<>& f2 detected" << std::endl;
#endif
  AddValSca(sca,val * f1.getVal() * f2.getVal());
}

template< typename M, typename V, typename S, typename T >
void  xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueFixed<S>& f, const xStateOfValueDof& other) const
{
  if ( vec && xPolicyAssemblerDispatcher<typename M::matrix_pattern>::is_sym  )
  {
       vec->AddVal(other.Numdof, -val * f.getVal() );
  }
}

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValue& s, const xStateOfValueLinearCombination<S>& equ) const
{
  xExecutor other (*this);
  int n = equ.size();
  const double two = 2.0;
  for (int i = 0; i < n; ++i)
  {
    const xStateOfValue * p_state_valuei = equ.state(i);
    other.val = xPolicyAssemblerDispatcher<typename M::matrix_pattern>::valCoef(equ.coeff(i),this->val,&s,p_state_valuei,two);
    Dispatcher dispatcher (other);
    dispatcher.Dispatch(s, *p_state_valuei);
  }
}

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueLinearCombination<S>& equi,
                          const xStateOfValueLinearCombination<S>& equj) const
{
  xExecutor other (*this);
  int n = equi.size();
  int m = equj.size();
  if ( xPolicyAssemblerDispatcher<typename M::matrix_pattern>::is_sym  )
  {
     if ( (&equi) != (&equj) )
     {
        const double k=2.00E+0;
        for (int i = 0; i < n; ++i)
        {
          const xStateOfValue * p_state_valuei = equi.state(i);
          for (int j = 0; j < m; ++j)
          {
            const xStateOfValue * p_state_valuej = equj.state(j);
            if ( p_state_valuei == p_state_valuej )
                other.val = k*this->val * equi.coeff(i) * equj.coeff(j);
            else
                other.val = this->val * equi.coeff(i) * equj.coeff(j);
            Dispatcher dispatcher (other);
            dispatcher.Dispatch(*p_state_valuei, *p_state_valuej);
          }
        }
     }
     else
     {
        const double k=0.50E+0;
        for (int i = 0; i < n; ++i)
        {
          const xStateOfValue * p_state_valuei = equi.state(i);
          for (int j = 0; j < m; ++j)
          {
            const xStateOfValue * p_state_valuej = equj.state(j);
            if ( p_state_valuei == p_state_valuej )
                other.val = this->val * equi.coeff(i) * equj.coeff(j);
            else
                other.val = k*this->val * equi.coeff(i) * equj.coeff(j);
            Dispatcher dispatcher (other);
            dispatcher.Dispatch(*p_state_valuei, *p_state_valuej);
          }
        }
     }
  }
  else
  {
     for (int i = 0; i < n; ++i)
     {
       for (int j = 0; j < m; ++j)
       {
         other.val = this->val * equi.coeff(i) * equj.coeff(j);
         Dispatcher dispatcher (other);
         dispatcher.Dispatch(*equi.state(i), *equj.state(j));
       }
     }
  }
}

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueLinearCombination<S>& equ, const xStateOfValue& s) const
{
  xExecutor other (*this);
  int n = equ.size();
  const double two = 2.0;
  for (int i = 0; i < n; ++i)
  {
    const xStateOfValue * p_state_valuei = equ.state(i);
    other.val = xPolicyAssemblerDispatcher<typename M::matrix_pattern>::valCoef(equ.coeff(i),this->val,&s,p_state_valuei, two);
    Dispatcher dispatcher (other);
    dispatcher.Dispatch(*p_state_valuei,s);
  }
}

//for vector assembly

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueFixed<S>& equ, const xStateOfValueNone& no) const
{
}

template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueDof& equ, const xStateOfValueNone& no) const
{
  if (vec) vec->AddVal(equ.Numdof, val);
}
template< typename M, typename V, typename S, typename T >
void xMatVecScaAssembler<M,V,S,T>::xExecutor::operator()(const xStateOfValueNone& equ, const xStateOfValueNone& no) const
{
  //  *sca += val;
  AddValSca(sca,val);
}



template< typename M, typename V, typename S, typename T >
 void xMatVecScaAssembler<M,V,S,T>::assemble(const xStateOfValue* s1, const xStateOfValue* s2, S v,
                     M* mat, V* vec, S* sca)
 {
   xExecutor executor(v, sca, vec, mat);
   Dispatcher dispatcher (executor);
   dispatcher.Dispatch(*s1, *s2);
 }

template< typename M, typename V, typename S, typename T >
 void xMatVecScaAssembler<M,V,S,T>::assemble(const xStateOfValue* s, S v, V* vec, S* sca)
{
  xExecutor executor(v, sca, vec);
  Dispatcher dispatcher (executor);
  xStateOfValueNone none;
  dispatcher.Dispatch(*s, none);
}

template< typename M, typename V, typename S, typename T >
 void xMatVecScaAssembler<M,V,S,T>::assemble(S v, S* sca)
{
  xExecutor executor(v, sca);
  Dispatcher dispatcher (executor);
  xStateOfValueNone none;
  dispatcher.Dispatch(none, none);
}

template< typename M, typename V, typename S, typename T >
 void xMatVecScaAssembler<M,V,S,T>::AddMatrixMat(M* mat,int dofi, int dofj, T value)
{
  if (mat) mat->AddMatrix(dofi, dofj, value);
#ifndef NDEBUG
  const bool debug=false;
  if(debug)
    std::cout<<"Valeur ajoutee="<<value
      <<" coord matrice :"<<dofi<<" ,"<<dofj<<std::endl;
#endif
}

template< typename M, typename V, typename S, typename T >
 void xMatVecScaAssembler<M,V,S,T>::AddValVec(V* vec,int numdof, T value)
{
  if (vec) vec->AddVal(numdof, value);
}

template< typename M, typename V, typename S, typename T >
 void xMatVecScaAssembler<M,V,S,T>::AddValSca(S* sca, T value)
{
  if(sca) *sca += value;
}


} // namespace xfem

#endif
