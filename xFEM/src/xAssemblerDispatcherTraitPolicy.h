/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XTRAITPOLICYASSEMBLERDISPATCHER_H
#define _XTRAITPOLICYASSEMBLERDISPATCHER_H


#include "xTraitsMatrix.h"


namespace xfem
{

 class xStateOfValue;

// default unsymetrique : i,j pair don't have to be reorder and no special treatement for linear combination
// covers xTraitMatrixUnSym and xTraitMatrixUnSymPSym 
// covers also xTraitMatrixDiagonal as a side effect
template<typename PATTERN >
class  xPolicyAssemblerDispatcher
{
  public :
   const static bool is_sym=false;
   static inline void reorder( int& i, int& j ) 
   {
        return;
   }

   template<typename VT>
   static inline VT valCoef(VT coef,VT val,const xStateOfValue* dofi, const xStateOfValue* dofj, const double k)
   {
         return(val*coef);
   }
   
};

// Symetric matrix : case lower triangle filled
template<>
class  xPolicyAssemblerDispatcher<xTraitMatrixLowerSym>
{
  public :
   const static bool is_sym=true;
   static inline void reorder( int& i, int& j ) 
   {
        if (i<j)
        {
            int k=i; 
            i=j;
            j=k;
        }
        return;
   }

   template<typename VT>
   static inline VT valCoef(VT coef,VT val,const xStateOfValue* dofi, const xStateOfValue* dofj, const double k)
   {
        if (dofi == dofj ) return(k*val*coef);
        else return(val*coef);
   }
};

// Symetric matrix : case upper triangle filled
template<>
class  xPolicyAssemblerDispatcher<xTraitMatrixUpperSym>
{
  public :
   const static bool is_sym=true;
   static inline void reorder( int& i, int& j ) 
   {
        if (j<i)
        {
            int k=i; 
            i=j;
            j=k;
        }
        return;
   }

   template<typename VT>
   static inline VT valCoef(VT coef,VT val,const xStateOfValue* dofi,const  xStateOfValue* dofj, const double k)
   {
        if (dofi == dofj ) return(k*val*coef);
        else return(val*coef);
   }
};



} // end namespace

#endif
