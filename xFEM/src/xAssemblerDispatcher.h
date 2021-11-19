/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// -*- C++ -*-
#ifndef ASSEMBLERDISPATCHER_H
#define ASSEMBLERDISPATCHER_H

#include "xStateOfValue.h"
#include "xDispatcher.h"
#include "xAssemblerDispatcherTraitPolicy.h"

namespace xfem { 

template< typename M, typename V, typename S, typename T = S >
class xMatVecScaAssembler 
{    

public:

    class xExecutor
    {
    public:
        xExecutor (T val, S* sca, V* vec=nullptr, M* mat=nullptr);

        typedef void result_type;
        using value_t = T;

        //functions needed
        //for the matrix assembly
        void operator()(const xStateOfValueDof&,const xStateOfValueDof&) const;
        void operator()(const xStateOfValueDof&,const xStateOfValueFixed<S>&) const;
        void operator()(const xStateOfValueFixed<S>&,const xStateOfValueDof&) const;
        void operator()(const xStateOfValueFixed<S>&,const xStateOfValueFixed<S>&) const;
        void operator()(const xStateOfValue&,const xStateOfValueLinearCombination<S>&) const;
        void operator()(const xStateOfValueLinearCombination<S>&,const xStateOfValueLinearCombination<S>&) const;
        void operator()(const xStateOfValueLinearCombination<S>&,const xStateOfValue&) const;

        //for the vector assembly
        void operator()(const xStateOfValueDof&,     const xStateOfValueNone&) const;
        void operator()(const xStateOfValueFixed<S>&,   const xStateOfValueNone&) const;
        void operator()(const xStateOfValueNone& a, const xStateOfValueNone& b) const;

        void operator()(const xStateOfValue& a, const xStateOfValue& b) const {assert(0);}

    private:
        T       val;                  // current value to assemble
        S     * sca;
        V     * vec;
        M     * mat;
    };

    /*  typedef xDispatcher
    < xExecutor,
      boost::mpl::vector<const xStateOfValueDof, const xStateOfValueFixed,
      const xStateOfValueLinearCombination,
      const xStateOfValueNone, const xStateOfValue>
    >  Dispatcher;*/

    typedef xDispatcherStateOfValue < xExecutor  >  Dispatcher;


    static void assemble(const xStateOfValue* s1, const xStateOfValue* s2, S v, M* mat, V* vec, S* sca);
    static void assemble(const xStateOfValue* s, S v, V* vec, S* sca);
    static void assemble(S v, S* sca);
    static void AddMatrixMat(M* mat,int dofi, int dofj, T value);
    static void AddValVec(V* vec,int numdof, T value);
    static void AddValSca(S* sca, T value);

private:

};


} // end namespace

#include "xAssemblerDispatcher_imp.h"

#endif
