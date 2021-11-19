/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _X_FIELD_H
#define _X_FIELD_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "xApproxFctPtr.h"
#include "xCommandOnGeomElem.h"
#include "xEval.h"
#include "xFieldStorage.h"
#include "xForm.h"
#include "xGeomElem.h"
#include "xSpace.h"
#include "xSpacePtr.h"
#include "xTensor2.h"
#include "xTensorsPtr.h"
#include "xValue.h"
#include "xValueManager.h"
#include "xVariabManager.h"

namespace xfem
{
//  class xFormZero;
/// base class for every field. Content yet to be determined !
class xFieldBase
{
  public:
};

/// basic class used in most algorithm. Know the relationship between values and spaces  a "field"
/*!
    Most important member :
    beginFcts,    endFcts,    beginValues,   endValues.
    !*/
template <typename VT = double>
class xField
{
  public:
   typedef VT value_t;
   typedef spacePtr_t spacePtr;
   typedef approxFctPtr_t shapeFctPtr;
   typedef std::vector<shapeFctPtr> approximations_type;
   typedef xValueManagerDist<VT> value_manager_t;
   typedef std::vector<typename value_manager_t::xvalue_t*> value_container_t;
   typedef std::vector<spacePtr> rep_type;
   typedef rep_type::const_iterator const_iterator;
   typedef rep_type::iterator iterator;

   ~xField();

   xField(xValueManagerDist<VT>* vm);

   /// This constructor build a new field the space pointer f is pushed in the field. No space constructed.
   xField(xValueManagerDist<VT>* vm, const spacePtr& f);

   /// This constructor build a new field containing  a  space of type T, created during construction of the xField
   template <class T>
   xField(xValueManagerDist<VT>* vm, const T& f);

   /// This constructor build a new field containing  a  space of type T1 and T2, each are created during construction of the
   /// xField
   template <class T1, class T2>
   xField(xValueManagerDist<VT>* vm, const T1& f1, const T2& f2);

   /// This constructor build a new field containing  a  space of type T1, T2 and T3 each are created during construction of the
   /// xField
   template <class T1, class T2, class T3>
   xField(xValueManagerDist<VT>* vm, const T1& f1, const T2& f2, const T3& f3);

   /// This function reset the fieldstorage .
   /*!
     Each time you look for fonctions or value inside the field, it goes throught an xFieldStorage member object, a 'policy', to
     check if this things are already known. This function destroy the previous xFieldStorage  and create a new one of type
     xFieldStorageElement, the default xFieldStoragePolicy.
     !*/
   void ResetStoragePolicy();
   /// This is a template version of ResetStoragePolicy. The template parameter must be a class derived from xFieldStorage.
   template <class S>
   void ResetStoragePolicy();
   /// this function return a pointer to the storage policy.
   /*! N.C. note :
     I think  this function is useless and ptentially dangerous. It should probably be removed
     !*/
   xFieldStorage<VT>* GetStoragePolicy();
   /// insert a new space of type T in the list of space of the field. Note that a new space is constructed
   template <class T>
   void insert(const T& f);
   /// insert a pointer to a space in the xField list of space. Note that no new space is created.
   void insert(const spacePtr& f);
   /// return the size of the spaces container.
   int size() const { return spaces.size(); }
   /// return an iterator to the begining of the container of spaces
   iterator begin() { return spaces.begin(); }
   /// return an iterator to the end of the container of spaces
   iterator end() { return spaces.end(); }
   /// return an iterator to the begining of the container of spaces
   const_iterator begin() const { return spaces.begin(); }
   /// return an iterator to the end of the container of spaces
   const_iterator end() const { return spaces.end(); }
   /// return the number of shape functions associated to entity e;
   int sizeFcts(AOMD::mEntity* e) const;
   /// return an iterator to the beginning of the container of  shape functions pointers associated to entity e;
   /*!  Note :
     depending of the storage policy, this can take more or less time.
     Also note that the returned iterator might be invalidated by a subsecant call to one of the beginFct, endFcts beginValues and
     endValues function.
     !*/
   approximations_type::iterator beginFcts(AOMD::mEntity* e) const;
   /// return an iterator to the end of the container of  shape functions pointers associated to entity e;
   approximations_type::iterator endFcts(AOMD::mEntity* e) const;
   /// return an iterator to the beginning of the container of values pointer associated to entity e;
   typename value_container_t::iterator beginValues(AOMD::mEntity* e) const;
   /// return an iterator to the end of the container of values pointer associated to entity e;
   typename value_container_t::iterator endValues(AOMD::mEntity* e) const;
   /// clear the container of space and reset the storage policy.
   void clear();

   /// Push the double associated to each xValue->getVal() associated to element e at the end of vector vals
   /*!
     N.C. note : this is probably too specific and should be I think removed. One can do the same with more flexibility using
     beginValues and end Values.
     !*/
   void getVals(AOMD::mEntity* e, std::vector<VT>& vals);
   /// set all the xValue->setVal() associated to eement e to value v.
   /*!
     N.C. note : this is probably too specific and should be I think removed. One can do the same with more flexibility using
     beginValues and endValues.
     !*/
   void setVal(AOMD::mEntity* e, const VT& v);

   /// return the double manager associated to the field.
   /*!
     N.C. note:
     This is currently in use in xAlgorithm. I think this is not very safe ... It would be better for it to desappear in the futur
     !*/
   value_manager_t* getValueManager() const { return value_manager; }

   /// This function evaluate all the shape functions at the integration point hidden in geo_appro/geo_integ; associated to the
   /// element pointed to geo_appro->getEntity(), then retrive the associated xValues and compute the approximation of the field
   /// at the integration point. the computed value is returned throught v and the return value ...
   template <class T>
   T& getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& v) const;

   /// This function evaluate all the shape functions gradient at the integration point hidden in geo_appro/geo_integ; associated
   /// to the element pointed to geo_appro->getEntity(), then retrive the associated xValues and compute the gradient of the
   /// approximation of the field at the integration point. the computed value is returned throught v and the return value ...
   template <class T>
   T& getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& v) const;

   /// This function evaluate all the shape functions gradient at the integration point hidden in geo_appro/geo_integ; associated
   /// to the element pointed to geo_appro->getEntity(), then retrive the associated xValues and compute the gradient of the
   /// approximation of the field at the integration point. the computed value is returned throught v and the return value ...
   /// This is the gradient with regard to the local coordinate of the element. (For this to work, all the approx function should
   /// have getgradlocal implemented, which is not necesserly the case ... )
   template <class T>
   T& getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& v) const;

   /// getFF, getGradFF, getGradFFEigen, getGradLocalFF and resetStorage are static member of xField.
   /*!
     N.C. note : I think this functions should be removed from xField: They make sense outside of the class.
    !*/
   template <class iterFct, class UnaryOperator>
   static void getFF(iterFct it, iterFct ite, std::vector<typename UnaryOperator::result_type>& ff, const xGeomElem* geo_appro,
                     const xGeomElem* geo_integ, UnaryOperator funct, bool resetStorage = false);

   template <class iterFct, class UnaryOperator>
   static void getGradFF(iterFct it, iterFct ite, std::vector<typename UnaryOperator::result_type>& grads,
                         const xGeomElem* geo_appro, const xGeomElem* geo_integ, UnaryOperator funct, bool resetStorage = false);

   template <class iterFct, class UnaryOperator, class CONTAINER>
   static void getGradFFEigen(iterFct it, iterFct ite, CONTAINER& grads, const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                              UnaryOperator funct, bool resetStorage = false);

   template <class iterFct, class UnaryOperator>
   static void getGradLocalFF(iterFct it, iterFct ite, std::vector<typename UnaryOperator::result_type>& grads,
                              const xGeomElem* geo_appro, const xGeomElem* geo_integ, UnaryOperator funct);

   /// resetStorage : This function reset the storage of the approxfunction iterator range. Note that this has noting to do with
   /// the xfield storage.
   template <class iterFct>
   static void resetStorage(iterFct it, iterFct ite, const xGeomElem* geo_appro, const xGeomElem* geo_integ);

  private:
   template <class T>
   class ptr_product : public std::binary_function<T, xValue<VT>*, T>
   {
     public:
      typedef typename std::binary_function<T, xValue<VT>*, T>::result_type result_type;
      result_type operator()(const T& d, xValue<VT>* v) const { return d * v->getVal(); }
   };

   template <class iterFct, class iterVal, class T>
   static void getVal(iterFct it, iterFct ite, iterVal first, iterVal last, int size_Fct_Val, const xGeomElem* geo_appro,
                      const xGeomElem* geo_integ, T& v);

   template <class iterFct, class iterVal>
   static void getVal(iterFct it, iterFct ite, iterVal first, iterVal last, int size_Fct_Val, const xGeomElem* geo_appro,
                      const xGeomElem* geo_integ, std::complex<double>& v);

   template <class iterFct, class iterVal>
   static void getVal(iterFct it, iterFct ite, iterVal first, iterVal last, int size_Fct_Val, const xGeomElem* geo_appro,
                      const xGeomElem* geo_integ, xtensor::xVectorDoubleComplex& v);

   template <class iterFct, class iterVal>
   static void getGrad(iterFct it, iterFct ite, iterVal first, iterVal last, int size_Fct_Val, const xGeomElem* geo_appro,
                       const xGeomElem* geo_integ, xtensor::xVector<>& v);

   template <class iterFct, class iterVal>
   static void getGrad(iterFct it, iterFct ite, iterVal first, iterVal last, int size_Fct_Val, const xGeomElem* geo_appro,
                       const xGeomElem* geo_integ, xtensor::xTensor2<>& v);

   template <class iterFct, class iterVal>
   static void getGrad(iterFct it, iterFct ite, iterVal first, iterVal last, int size_Fct_Val, const xGeomElem* geo_appro,
                       const xGeomElem* geo_integ, xtensor::xVectorDoubleComplex& v);

   template <class iterFct, class iterVal>
   static void getGrad(iterFct it, iterFct ite, iterVal first, iterVal last, int size_Fct_Val, const xGeomElem* geo_appro,
                       const xGeomElem* geo_integ, xtensor::xTensor2DoubleComplex& v);

   template <class iterFct, class iterVal>
   static void getGradLocal(iterFct it, iterFct ite, iterVal first, iterVal last, const xGeomElem* geo_appro,
                            const xGeomElem* geo_integ, xtensor::xTensor2<>& v);

   void newstorage();

   void clearstorage();

  private:
   rep_type spaces;
   int TensoDim;
   xValueManagerDist<VT>* value_manager;
   void _insert(spacePtr d);
   mutable std::vector<xValue<VT>*> vals_curr_ptr;
   typedef typename std::vector<xValue<VT>*>::iterator iter_ptr;
   xFieldStorage<VT>* storage;
};  // end xField declaration

/// allows to store pointwise values at each integration point. It requires a space to be given (xSpacePointwise),
/// as well as valid names for the variable to be stored and for the gradient if desired.It uses an xVariabManager to store the
/// couples (key,value).
class xFieldPointwise : public xFieldBase
{
  public:
   typedef std::shared_ptr<xSpacePointwise> spacePtr;
   typedef std::vector<spacePtr> rep_type;
   typedef rep_type::const_iterator const_iterator;
   typedef rep_type::iterator iterator;
   typedef xVariabManager value_manager_t;

   template <class T>
   xFieldPointwise(value_manager_t* vm, const T& f) : value_manager(vm), val("val"), grad("grad")
   {
      insert(f);
   }
   template <class T>
   xFieldPointwise(value_manager_t* vm, const T& f, const std::string& val_) : value_manager(vm), val(val_), grad("grad")
   {
      insert(f);
   }
   template <class T>
   xFieldPointwise(value_manager_t* vm, const T& f, const std::string& val_, const std::string& grad_)
       : value_manager(vm), val(val_), grad(grad_)
   {
      insert(f);
   }

   value_manager_t* getValueManager() const { return value_manager; }
   int size() const { return spaces.size(); }
   iterator begin() { return spaces.begin(); }
   iterator end() { return spaces.end(); }
   const_iterator begin() const { return spaces.begin(); }
   const_iterator end() const { return spaces.end(); }
   void clear() { spaces.clear(); }
   std::string& valstr() { return val; }
   std::string& gradstr() { return grad; }

   template <class T>
   T& getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& r) const
   {
      xValKey key;
      spaces[0]->getKey(geo_integ, key);
      xValue<tensorsPtr_t>* v = value_manager->find(key);
      tensorsPtr_t t = v->getVal();
      r = t->get(val, boost::mpl::identity<T>());
      return r;
   }

   template <class T>
   T& getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& r) const
   {
      xValKey key;
      spaces[0]->getKey(geo_integ, key);
      xValue<tensorsPtr_t>* v = value_manager->find(key);
      tensorsPtr_t t = v->getVal();
      r = t->get(grad, boost::mpl::identity<T>());
      return r;
   }

   template <class T>
   void setVal(xGeomElem* geo_appro, xGeomElem* geo_integ, T r)
   {
      xValKey key;
      spaces[0]->getKey(geo_integ, key);
      xValue<tensorsPtr_t>* v = value_manager->find(key);
      tensorsPtr_t t = v->getVal();
      t->get(val, boost::mpl::identity<T>()) = r;
   }

   template <class T>
   void setGrad(xGeomElem* geo_appro, xGeomElem* geo_integ, T r)
   {
      xValKey key;
      spaces[0]->getKey(geo_integ, key);
      xValue<tensorsPtr_t>* v = value_manager->find(key);
      tensorsPtr_t t = v->getVal();
      t->get(grad, boost::mpl::identity<T>()) = r;
   }

  private:
   template <class T>
   void insert(const T& f)
   {
      _insert(spacePtr(new T(f)));
   }

   void insert(const spacePtr& f) { _insert(f); }

   void _insert(spacePtr d) { spaces.push_back(d); }

   value_manager_t* value_manager;
   rep_type spaces;
   std::string val;
   std::string grad;
};

/// evaluates the field at any point (for an approximation)
/// or only at integration points (for xFieldPointwise)
template <class UnaryOperator, class Field = xField<>>
class xEvalField : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalField(const Field& f_) : f(f_) {}
   xEvalField(const Field& f_, const UnaryOperator& _funct) : f(f_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v;
      f.getVal(appro, integ, v);
      result = funct(v);
   }

  private:
   const Field& f;
   UnaryOperator funct;
};

/// evaluates the gradient of the field at any point (for an approximation)
/// or only at integration points if it is defined (for xFieldPointwise)
template <class UnaryOperator, class Field = xField<>>
class xEvalGradField : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalGradField(const Field& f_) : f(f_) {}
   xEvalGradField(const Field& f_, const UnaryOperator& _funct) : f(f_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v;
      f.getGrad(appro, integ, v);
      result = funct(v);
   }

  private:
   const Field& f;
   UnaryOperator funct;
};

/// evaluates the axisymmetric gradient of the field at any point (for an approximation)
/// or only at integration points if it is defined (for xFieldPointwise)
/// WARNING : this works only for vectorial fields
template <class UnaryOperator, class Field = xField<>>
class xEvalGradFieldAxisym : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalGradFieldAxisym(const Field& f_) : f(f_) {}
   xEvalGradFieldAxisym(const Field& f_, const UnaryOperator& _funct) : f(f_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const
   {
      // Not optimized at all ...
      typename UnaryOperator::argument_type v;
      xtensor::xVector<> val;
      f.getVal(appro, integ, val);

      auto uvw = integ->getXYZ();
      double r = uvw(0);

      f.getGrad(appro, integ, v);
      result = funct(v);
      if (r > 0.)
      {
         result(2, 2) = val(0) / r;
      }
      else
      {
         result(2, 2) = result(0, 0);
      }
   }

  private:
   const Field& f;
   UnaryOperator funct;
};

/// evaluates the field's mean value over each element
/// used primarily with xFieldPointwise
template <class UnaryOperator, class Field = xField<>>
class xEvalFieldMeanValue : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalFieldMeanValue(const Field& f_) : f(f_) {}
   xEvalFieldMeanValue(const Field& f_, const UnaryOperator& _funct) : f(f_), funct(_funct) {}
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v, vtot = typename UnaryOperator::argument_type();
      double wtot = 0.;
      xGeomElem appro_copy(*appro);
      xGeomElem integ_copy(*integ);
      /*		xGeomElem appro_copy(appro->getEntity());
      xGeomElem integ_copy(integ->getEntity());*/
      int Ng = integ_copy.GetNbIntegrationPoints();
      for (int gaussp = 0; gaussp < Ng; ++gaussp)
      {
         integ_copy.setUVW(gaussp);
         if (appro_copy.getEntity() != integ_copy.getEntity())
            appro_copy.setUVWForXYZ(integ_copy.getXYZ());
         else
            appro_copy.setUVW(integ_copy.getUVW());

         f.getVal(&appro_copy, &integ_copy, v);
         double w = integ_copy.GetWeight();
         vtot += v * w;
         wtot += w;
      }
      vtot /= wtot;
      result = funct(vtot);
   }

  private:
   const Field& f;
   UnaryOperator funct;
};

/// evaluates the gradient of the field's mean value over each element
/// used primarily with xFieldPointwise (when a gradient is available !)
template <class UnaryOperator, class Field = xField<>>
class xEvalGradFieldMeanValue : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalGradFieldMeanValue(const Field& f_) : f(f_) {}
   xEvalGradFieldMeanValue(const Field& f_, const UnaryOperator& _funct) : f(f_), funct(_funct) {}
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const
   {
      typename UnaryOperator::argument_type v, vtot = typename UnaryOperator::argument_type();
      double wtot = 0.;
      xGeomElem appro_copy(*appro);
      xGeomElem integ_copy(*integ);
      int Ng = integ->GetNbIntegrationPoints();
      for (int gaussp = 0; gaussp < Ng; ++gaussp)
      {
         integ_copy.setUVW(gaussp);
         if (appro_copy.getEntity() != integ->getEntity())
            appro_copy.setUVWForXYZ(integ->getXYZ());
         else
            appro_copy.setUVW(gaussp);

         f.getGrad(appro_copy, integ_copy, v);
         double w = integ_copy.GetWeight();
         vtot += v * w;
         wtot += w;
      }
      vtot /= wtot;
      result = funct(vtot);
   }

  private:
   const Field& f;
   UnaryOperator funct;
};

template <class UnaryOperator>
class xEvalGradFieldOptimized : public xEval<typename UnaryOperator::result_type>  //
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalGradFieldOptimized(const xField<>& f_) : f(f_) {}
   xEvalGradFieldOptimized(const xField<>& f_, const UnaryOperator& _funct) : f(f_), funct(_funct) {}
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const
   {
      const bool debug = xdebug_flag;
      typename UnaryOperator::argument_type v;
      f.getGradLocal(appro, integ, v);

      appro->PushBackRightTranspose(v);
      result = funct(v);
      if (debug) std::cout << result << std::endl;
   }

  private:
   const xField<>& f;
   UnaryOperator funct;
};

template <typename VT = double>
class xFillFieldFromZeroForm : public xIntegrateFormCommand
{
  public:
   xFillFieldFromZeroForm(xField<VT>& f, xFormZero<VT>& fo);
   void openApproxElem(xGeomElem* g_appro) override;
   void closeApproxElem(xGeomElem* g_appro) override;
   VT getTotal() { return total; }

  private:
   xField<VT>& field;
   xFormZero<VT>* form;
   VT total;
};

#include "xField_imp.h"

}  // namespace xfem

#endif
