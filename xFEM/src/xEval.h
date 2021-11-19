/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _Eval_hh_
#define _Eval_hh_
#include <memory>

#include "xApproxFunction.h"
#include "xEntityToEntity.h"
#include "xGeomElem.h"
#include "xMapping.h"

namespace xfem
{
class xApproxFunction;

template <class T>
class xEval
{
  public:
   typedef T result_type;
   virtual ~xEval() = default;
   virtual void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, result_type& res) const = 0;

  protected:
  private:
};

template <class T>
class xEvalConstant : public xEval<T>
{
  public:
   typedef typename xEval<T>::result_type result_type;
   xEvalConstant(const result_type& v) : val(v) {}
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, result_type& res) const override { res = val; }

  private:
   result_type val;
};

template <class T>
class xEvalZero : public xEval<T>
{
  public:
   typedef typename xEval<T>::result_type result_type;
   xEvalZero() : val(0.) {}
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, result_type& res) const override { res = val; }

  private:
   result_type val;
};

template <class T, typename PointEval>
class xEvalFctAtPoint : public xEval<T>
{
  public:
   typedef typename xEval<T>::result_type result_type;
   xEvalFctAtPoint(PointEval& e) : eval(e) {}
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, result_type& res) const override
   {
      res = eval(geo_integ->getXYZ());
   }

  public:
   const PointEval& getPointEvalFunctor() const { return eval; };

  private:
   PointEval& eval;
};

/// xEvalUnary is a xEval buil on top of a given xEval and an unary operator applied to it.
template <class UnaryOperator>
class xEvalUnary : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

  private:
   const xEval<typename UnaryOperator::argument_type>& First;
   mutable typename UnaryOperator::argument_type first;
   UnaryOperator funct;

  public:
   xEvalUnary(const xEval<typename UnaryOperator::argument_type>& f) : First(f) {}

   xEvalUnary(UnaryOperator unOp, const xEval<typename UnaryOperator::argument_type>& f) : First(f), funct(unOp) {}

   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      First(appro, integ, first);
      result = funct(first);
   }
};

/// xEvalBinary is an xEval build on top of 2 xEval and a binary Operator.
template <class BinaryOperator>
class xEvalBinary : public xEval<typename BinaryOperator::result_type>
{
  public:
   typedef typename xEval<typename BinaryOperator::result_type>::result_type result_type;

  private:
   const xEval<typename BinaryOperator::first_argument_type>& First;
   const xEval<typename BinaryOperator::second_argument_type>& Second;
   mutable typename BinaryOperator::first_argument_type first;
   mutable typename BinaryOperator::second_argument_type second;
   BinaryOperator funct;

  public:
   xEvalBinary(const xEval<typename BinaryOperator::first_argument_type>& f,
               const xEval<typename BinaryOperator::second_argument_type>& s)
       : First(f), Second(s)
   {
   }

   xEvalBinary(BinaryOperator BinOp, const xEval<typename BinaryOperator::first_argument_type>& f,
               const xEval<typename BinaryOperator::second_argument_type>& s)
       : First(f), Second(s), funct(BinOp)
   {
   }

   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      First(appro, integ, first);
      Second(appro, integ, second);
      result = funct(first, second);
   }
};

template <class T>
class xEvalRedirectAppro : public xEval<T>
{
  public:
   xEvalRedirectAppro(const xEval<T>& e, xEntityToEntity i) : eval(e), integ2appro(i) {}

   void operator()(const xGeomElem* appro, const xGeomElem* integ, T& result) const
   {
      AOMD::mEntity* e_integ = integ->getEntity();
      AOMD::mEntity* e_appro = integ2appro(e_integ);
      xGeomElem geom_appro(e_appro);
      auto& mapping_appro{*appro->getMapping()};
      auto& mapping_integ{*integ->getMapping()};
      if ((geom_appro.getEntity() != integ->getEntity()) || (typeid(mapping_appro) != typeid(mapping_integ)))
         geom_appro.setUVWForXYZ(integ->getXYZ());
      else
         geom_appro.setUVW(integ->getUVW());
      eval(&geom_appro, integ, result);
   }

  private:
   const xEval<T>& eval;
   xEntityToEntity integ2appro;
};

template <class T>
class xEvalRedirectApproRedirectInteg : public xEval<T>
{
  public:
   xEvalRedirectApproRedirectInteg(const xEval<T>& e, xEntityToEntity _redirectappro, xEntityToEntity _redirectinteg)
       : eval(e), redirectappro(_redirectappro), redirectinteg(_redirectinteg)
   {
   }

  private:
   void operator()(const xGeomElem* appro, const xGeomElem* integ, T& result) const
   {
      AOMD::mEntity* e_integ = integ->getEntity();
      AOMD::mEntity* r_e_integ = redirectinteg(e_integ);
      AOMD::mEntity* e_appro = appro->getEntity();
      AOMD::mEntity* r_e_appro = redirectappro(e_appro);
      xGeomElem r_appro(r_e_appro);
      xGeomElem r_integ(r_e_integ);
      r_appro.setUVWForXYZ(appro->getXYZ());
      r_integ.setUVWForXYZ(integ->getXYZ());
      eval(&r_appro, &r_integ, result);
   }
   const xEval<T>& eval;
   xEntityToEntity redirectappro;
   xEntityToEntity redirectinteg;
};

template <class T>
class xEvalApproxFunction : public xEval<T>
{
  public:
   xEvalApproxFunction(const xApproxFunction& _approx) : approx(_approx) {}
   void operator()(const xGeomElem* appro, const xGeomElem* integ, T& result) const override
   {
      approx.getVal(appro, integ, result);
   }

  private:
   const xApproxFunction& approx;  // const
};

// Evaluateur gradient de la fonction d'enrichissement (K)
template <class T>
class xEvalGradApproxFunction : public xEval<T>
{
  public:
   xEvalGradApproxFunction(xApproxFunction& _appro) : appro(_appro){};
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& res) const override
   {
      appro.getGrad(geo_appro, geo_integ, res);
   };

  private:
   const xApproxFunction& appro;
};

// pour main test cas 1D
class xEvalBodyForce2D : public xEval<xtensor::xVector<>>
{
  public:
   xEvalBodyForce2D(double& _f) : f(_f){};
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, result_type& force) const override
   {
      const double x = geo_integ->getXYZ()(0);
      force = {f * pow(x / 2., 2.), 0., 0.};  //  f * (x/l)^2 , avec l = longueur de la plaque
   }

  private:
   const double f;
};

/// Eval from a given Lambda function
template <class T>
class xEvalFromLambda : public xEval<T>
{
  public:
   xEvalFromLambda(std::function<T(const xGeomElem*, const xGeomElem*)> lambda_) : lambda{lambda_} {}

   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& res) const
   {
      res = lambda(geo_appro, geo_integ);
      return;
   }

  private:
   std::function<T(const xGeomElem* geo_appro, const xGeomElem* geo_integ)> lambda;
};

/// Warper that provides same API as xEval but gives a way to update underlying evaluator.
//! Practical to store warper in complex context that stay the same while evaluator need to be changed.
//! Calling update reset underlying evaluator and then complex context may stay the same but it will use this
//! updated evaluator.
//! Be careful this warper takes ownership of the pointer given as argument at construction time. User MUST not gives
//! a pointer to stack variable or pointer deleted outside this class. This class will delete allocated heap memory
//! Usage is then: xEvalWarper<XXX> warper(new ....);
//! Warper instance if copied share the same evaluator with its copy till one of the two call update with a new evaluator.
//! The update method is not propagating to copies the changes. This makes possible duplication of complex context having
//! the same underlying evaluator and then choose afterward to modify it or not.
template <class T>
class xEvalWarper : public xEval<T>
{
  public:
   xEvalWarper(xEval<T>* eval_ = nullptr) : eval(eval_) {}
   inline void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, T& res) const override
   {
      assert(eval);
      (*eval)(geo_appro, geo_integ, res);
      return;
   }
   void update(xEval<T>* eval_)
   {
      eval.reset(eval_);
      return;
   }
   bool empty() const { return !(eval); }

  private:
   std::shared_ptr<xEval<T>> eval;
};

}  // namespace xfem

#endif
