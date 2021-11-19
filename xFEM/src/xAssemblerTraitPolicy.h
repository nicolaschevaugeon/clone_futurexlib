/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XTRAITPOLICYASSEMBLER_H
#define _XTRAITPOLICYASSEMBLER_H


#include "xTraitsMatrix.h"
#include "xValue.h"


namespace xfem
{

// Traits for xAssembler class  concerning the need to transpose
// terme when in this context (in xAssemblerBasicAndTranspose)
// default is no and the unique cases for now is when unsym matrix
template<typename PATTERN>
class  xTraitsAssembler
{
  public :
     const static bool do_transpose=false;
};

template<>
class  xTraitsAssembler<xTraitMatrixUnSym>
{
  public :
     const static bool do_transpose=true;
};
template<>
class  xTraitsAssembler<xTraitMatrixUnSymPSym>
{
  public :
     const static bool do_transpose=true;
};
//
/// Policies concerning assembly on zero or not for xAssembler class
template<typename ASSEMBLY_ON_ZERO>
class  xPolicyAssemblerOnZero;

template<>
class  xPolicyAssemblerOnZero<xTraitMatrixNoAssemblyOnZero>
{
  public :

    //Plus compact...efficacite ??

    template<typename VT>
    static inline  bool  isNull(VT & val)
    {
        return ( (val == xtool::xDataType<VT>::zero()) );
    }


//    //Moins compact...mais efficace ?
//  static inline  bool  isNull(double & val)
//  {
//       return ( (val == 0.0E+0) );
//  }

//  static inline  bool  isNull(float & val)
//  {
//      return ( (val == 0.f) );
//  }


};


template<>
class  xPolicyAssemblerOnZero<xTraitMatrixForcedAssemblyOnZero>
{
  public :
    template<typename VT>
   static inline  bool  isNull(VT & val)
   { 
        return ( false );
   }
};

// default symetrique : lower or upper are treated the same way as it's in the dispatcher where
// decision take place. Elementary matrice is considered to be upper.
template<typename PATTERN >
class  xPolicyAssemblerLoop
{
  public :
//   typedef  std::vector<xValue<double>*>::iterator  iterDof;

    template<typename ITERDOF>
   static inline  void  firstForLoop(ITERDOF & first, ITERDOF & firsti, ITERDOF &start, int i, int & j)
   {
        start = firsti;
        j = i;
        return;
   }
   template<typename ITERDOF>
   static inline  ITERDOF&  lastForLoop(ITERDOF & last, ITERDOF firsti)
   {
        return last;
   }
};

// Unsymetrique matrice : Elementary matrice is considered to be full. 
template<>
class  xPolicyAssemblerLoop<xTraitMatrixUnSym>
{
  public :
//   typedef  std::vector<xValue<double>*>::iterator  iterDof;

    template<typename ITERDOF>
   static inline  void  firstForLoop(ITERDOF & first, ITERDOF & firsti, ITERDOF &start, int i, int & j)
   {
        start = first;
        j = 0;
        return;
   }
   template<typename ITERDOF>
   static inline  ITERDOF&  lastForLoop(ITERDOF & last, ITERDOF firsti)
   {
        return last;
   }
};
template<>
class  xPolicyAssemblerLoop<xTraitMatrixUnSymPSym>
{
  public :
//   typedef  std::vector<xValue<double>*>::iterator  iterDof;
    template<typename ITERDOF>
   static inline  void  firstForLoop(ITERDOF & first, ITERDOF & firsti, ITERDOF &start, int i, int & j)
   {
        start = first;
        j = 0;
        return;
   }
   template<typename ITERDOF>
   static inline  ITERDOF&  lastForLoop(ITERDOF & last, ITERDOF firsti)
   {
        return last;
   }
};

// Diagonal matrice : Elementary matrice is watever it is. 
template<>
class  xPolicyAssemblerLoop<xTraitMatrixDiagonal>
{
  public :
//   typedef  std::vector<xValue<double>*>::iterator  iterDof;
    template<typename ITERDOF>
   static inline  void  firstForLoop(ITERDOF & first, ITERDOF & firsti, ITERDOF &start, int i, int & j)
   {
        start = first;
        j = 0;
        return;
   }
   template<typename ITERDOF>
   static inline  ITERDOF  lastForLoop(ITERDOF & last, ITERDOF firsti)
   {
        return ++firsti;
   }
};


} // end namespace

#endif
