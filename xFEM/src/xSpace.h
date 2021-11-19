/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __SPACE_H
#define __SPACE_H

#include <string>
#include <vector>

// aomd
#include "AOMDfwd.h"

// xfem
#include "xApproxFctPtr.h"
#include "xCommandOnGeomElem.h"
#include "xEval.h"
#include "xIntegrationRule.h"
#include "xSpacePtr.h"
#include "xTensors.h"
#include "xTensorsPtr.h"
#include "xValKey.h"
#include "xValue.h"
#include "xVariabManager.h"

namespace xfem
{
class xApproxFunction;
class xValKey;

class xSpaceBase
{
  public:
   typedef std::vector<xValKey> femKeys;
   virtual ~xSpaceBase() = default;
   virtual void getKeys(AOMD::mEntity*, femKeys*) = 0;
};

class xSpace : public xSpaceBase
{
  public:
   enum TensorialType_t
   {
      SCALAR,
      VECTOR_X,
      VECTOR_Y,
      VECTOR_Z
   };

  public:
   ~xSpace() override = default;
   //  typedef xSpaceBase::femKeys femKeys;
   typedef std::vector<xValKey> femKeys;
   typedef spacePtr_t spacePtr;
   typedef approxFctPtr_t shapeFctPtr;
   typedef femFcts_t femFcts;
   //  virtual void getKeys(AOMD::mEntity*, femKeys *) = 0;
   virtual void getKeysAndFcts(AOMD::mEntity*, femKeys*, femFcts*) = 0;
};

class xSpacePointwise : public xSpaceBase
{
  public:
   ~xSpacePointwise() override = default;
   xSpacePointwise(xIntegrationRule& integ_);
   virtual void getKey(const xGeomElem* geo_integ, xValKey& key) = 0;
   void set_GP_KEY_size(unsigned int i);  // Redefinit le nombre max de point d'integration
  protected:
   xIntegrationRule& integ;
   std::vector<xValKey::ids_size_t> GP_KEY;
};

class xSpacePwRegular : public xSpacePointwise
{
  public:
   xSpacePwRegular(const std::string& a, xIntegrationRule& integ_) : xSpacePointwise(integ_), Phys(xKeyInfo::getPhysId(a)) {}

   ~xSpacePwRegular() override = default;
   //  typedef std::vector<xValKey> femKeys;
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKey(const xGeomElem* geo_integ, xValKey& key) override;

  protected:
   xValKey::ids_size_t Phys;
};

class xSpaceRegular : public xSpace
{
  public:
   xSpaceRegular(const std::string& a, TensorialType_t c);

  protected:
   xValKey::ids_size_t Phys;
   TensorialType_t TensorialType;

  private:
};

class xSpaceFiltered : public xSpace
{
  public:
   typedef std::function<bool(AOMD::mEntity* e)> filter_t;

   ~xSpaceFiltered() override;

   template <class T, class F>
   xSpaceFiltered(const T& base, const F& f) : accept(f), base_space(new T(base))
   {
   }
   template <class F>
   xSpaceFiltered(spacePtr base, const F& f) : accept(f), base_space(base)
   {
   }
   // int  getTensoDim()  {return base_space->getTensoDim();}
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   filter_t accept;
   spacePtr base_space;
};

// the following is the same than xSpaceFiltered,
// except that the argument given to the filter is the key, not the key's AOMD::mEntity*
class xSpaceKeyFiltered : public xSpace
{
  public:
   typedef std::function<bool(xfem::xValKey&)> filter_t;

   ~xSpaceKeyFiltered() override;

   template <class T, class F>
   xSpaceKeyFiltered(const T& base, const F& f) : accept(f), base_space(new T(base))
   {
   }

   void getKeys(AOMD::mEntity* e, femKeys* keys) override;

   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   filter_t accept;
   spacePtr base_space;
};

// the following is the same than xSpaceFiltered,
// except that the filter is directly applied on the AOMD::mEntity* e,
// and not on "e->keys->getEnti"
// therefore, the keys are created if e pass the filter,
// and not always created and then possibly deleted as done in xSpaceFiltered...
class xDiscontinuousSpaceFiltered : public xSpace
{
  public:
   typedef std::function<bool(AOMD::mEntity* e)> filter_t;

   ~xDiscontinuousSpaceFiltered() override;

   template <class T, class F>
   xDiscontinuousSpaceFiltered(const T& base, const F& f) : accept(f), base_space(new T(base))
   {
   }

   void getKeys(AOMD::mEntity* e, femKeys* keys) override;

   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   filter_t accept;
   spacePtr base_space;
};

class xSpaceComposite : public xSpace
{
  public:
   xSpaceComposite() {}
   template <class T>
   xSpaceComposite(const T& t)
   {
      children.push_back(spacePtr(new T(t)));
   }
   template <class T1, class T2>
   xSpaceComposite(const T1& t1, const T2& t2)
   {
      children.push_back(spacePtr(new T1(t1)));
      children.push_back(spacePtr(new T2(t2)));
   }
   template <class T1, class T2, class T3>
   xSpaceComposite(const T1& t1, const T2& t2, const T3& t3)
   {
      children.push_back(spacePtr(new T1(t1)));
      children.push_back(spacePtr(new T2(t2)));
      children.push_back(spacePtr(new T3(t3)));
   }
   template <class T1, class T2, class T3, class T4>
   xSpaceComposite(const T1& t1, const T2& t2, const T3& t3, const T4& t4)
   {
      children.push_back(spacePtr(new T1(t1)));
      children.push_back(spacePtr(new T2(t2)));
      children.push_back(spacePtr(new T3(t3)));
      children.push_back(spacePtr(new T4(t4)));
   }
   template <class T>
   void insert(const T& t)
   {
      children.push_back(spacePtr(new T(t)));
   }
   void insert(const spacePtr& t) { children.push_back(t); }
   void clear()
   {
      // children.erase(children.begin(),children.end());// this doesn't call the destructors of children
      children.clear();
   }

   // int  getTensoDim()  {return (*children.begin())->getTensoDim();}

   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

   typedef std::vector<spacePtr>::const_iterator const_iterator;

  private:
   std::vector<spacePtr> children;
};

class xSpaceDifference : public xSpace
{
  public:
   template <class T1, class T2>
   xSpaceDifference(const T1& t1, const T2& t2)
   {
      first = spacePtr(new T1(t1));
      second = spacePtr(new T2(t2));
   }
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   spacePtr first, second;
};

class xSpaceConstant : public xSpaceRegular
{
  public:
   xSpaceConstant(const std::string& a) : xSpaceRegular(a, SCALAR), CONSTANT_ELEMENT(xKeyInfo::getGeomId("CONSTANT_ELEMENT")) {}
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   const xValKey::ids_size_t CONSTANT_ELEMENT;
};

class xSpaceConstantGlobal : public xSpaceRegular
{
  public:
   xSpaceConstantGlobal(const std::string& a, AOMD::mEntity* _global_entity)
       : xSpaceRegular(a, SCALAR), CONSTANT_ELEMENT(xKeyInfo::getGeomId("CONSTANT_ELEMENT")), global_entity(_global_entity)
   {
   }
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   const xValKey::ids_size_t CONSTANT_ELEMENT;
   AOMD::mEntity* global_entity;
};

class xSpaceLagrange : public xSpaceRegular
{
  public:
   enum lag_degree_t
   {
      DEGREE_ZERO,
      DEGREE_ONE,
      DEGREE_TWO,
      DEGREE_THREE
   };

  public:
   xSpaceLagrange(const std::string& a, TensorialType_t b, lag_degree_t d) : xSpaceRegular(a, b), Degree(d)
   {
      for (int i = 0; i < 25; ++i)
      {
         char str[10];
         std::sprintf(str, "%d", i + 1);
         std::string name("HIERARCHICAL_" + std::string(str));
         HIERARCHICAL_KEY[i] = xKeyInfo::getGeomId(name);
      }
   }
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  protected:
   lag_degree_t Degree;
   xValKey::ids_size_t HIERARCHICAL_KEY[25];
};

class xSpacePieceWiseLagrange : public xSpaceRegular
{
  public:
   enum lag_degree_t
   {
      DEGREE_ZERO
   };

  public:
   xSpacePieceWiseLagrange(const std::string& a, TensorialType_t b, lag_degree_t d) : xSpaceRegular(a, b), Degree(d)
   {
      for (int i = 0; i < 25; ++i)
      {
         char str[10];
         std::sprintf(str, "%d", i + 1);
         std::string name("HIERARCHICAL_PIECEWISE_" + std::string(str));
         HIERARCHICAL_KEY[i] = xKeyInfo::getGeomId(name);
      }
   }
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  protected:
   lag_degree_t Degree;
   xValKey::ids_size_t HIERARCHICAL_KEY[25];
};

class xSpaceXFEM : public xSpace
{
  public:
   template <class T, class E>
   xSpaceXFEM(const T& base, const E& enr, ValKeyModifier_t mk)
       : base_space(new T(base)), enrichment(new E(enr)), key_modifier(mk)
   {
   }

   template <class E>
   xSpaceXFEM(const spacePtr& base, const E& enr, ValKeyModifier_t mk)
       : base_space(base), enrichment(new E(enr)), key_modifier(mk)
   {
   }

   xSpaceXFEM() {}

   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  protected:
   spacePtr base_space;
   shapeFctPtr enrichment;
   ValKeyModifier_t key_modifier;
};

/// \brief{Shifted X-FEM space.}
class xSpaceXFEMShifted : public xSpace
{
  public:
   template <class T, class E>
   xSpaceXFEMShifted(const T& base, const E& enr, ValKeyModifier_t mk)
       : base_space(new T(base)), enrichment(new E(enr)), key_modifier(mk)
   {
   }

   xSpaceXFEMShifted() {}

   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  protected:
   spacePtr base_space;
   shapeFctPtr enrichment;
   ValKeyModifier_t key_modifier;
};

/// \brief{Corrected X-FEM space. (only works if the enrichment function is shifted)}
/// \param[in] ramp : Ramp function evaluator.
/// \param[in] gramp : Grad of ramp function evaluator.
class xSpaceXFEMShiftedCorrected : public xSpace
{
  public:
   template <class T, class E>
   xSpaceXFEMShiftedCorrected(const T& base, const E& enr, ValKeyModifier_t mk, xEval<double>* ramp_,
                              xEval<xtensor::xVector<>>* gramp_)
       : base_space(new T(base)), enrichment(new E(enr)), key_modifier(mk), ramp(ramp_), gramp(gramp_)
   {
   }

   xSpaceXFEMShiftedCorrected() {}

   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
   xtensor::xPoint getInterpolantUVW(xValKey& key, shapeFctPtr& func);  // ONLY CODED IN 2D

  protected:
   spacePtr base_space;
   shapeFctPtr enrichment;
   ValKeyModifier_t key_modifier;
   xEval<double>* ramp;
   xEval<xtensor::xVector<>>* gramp;
};

class xSpaceVectorXFEM : public xSpaceXFEM
{
  public:
   template <class T, class E>
   xSpaceVectorXFEM(const T& base, const E& enr, ValKeyModifier_t mk) : xSpaceXFEM(base, enr, mk)
   {
   }
   template <class E>
   xSpaceVectorXFEM(const spacePtr& base, const E& enr, ValKeyModifier_t mk) : xSpaceXFEM(base, enr, mk)
   {
   }
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
};

class xSpaceVector2XFEM : public xSpaceXFEM
{
  public:
   template <class T, class E>
   xSpaceVector2XFEM(const T& base, const E& enr, ValKeyModifier_t mk) : xSpaceXFEM(base, enr, mk)
   {
   }
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
};

class xGetKeysOnGeomElemCommand : public xCommandOnGeomElem
{
  public:
   xGetKeysOnGeomElemCommand(xSpaceBase::femKeys* keys_, const std::vector<xValKey::ids_size_t>& GP_KEY_,
                             const xValKey::ids_size_t& Phys_)
       : keys(keys_), GP_KEY(GP_KEY_), Phys(Phys_){};
   ~xGetKeysOnGeomElemCommand() override = default;

   void execute(xGeomElem* g_integ) override
   {
      const unsigned int nb = g_integ->GetNbIntegrationPoints();
      for (unsigned int k = 0; k < nb; k++)
      {
         keys->push_back(xValKey(Phys, GP_KEY[k], g_integ->getEntity()));
      }
   }

  private:
   xSpaceBase::femKeys* keys;
   const std::vector<xValKey::ids_size_t>& GP_KEY;
   const xValKey::ids_size_t& Phys;
};

class xPointwiseVisitor
{
  public:
   virtual void setInfos(xGeomElem* integ, xGeomElem* appro)
   {
      geo_integ = integ;
      geo_appro = appro;
   }
   virtual void Visit(xValue<tensorsPtr_t>& v) = 0;
   virtual ~xPointwiseVisitor() = default;

  protected:
   xGeomElem* geo_integ;
   xGeomElem* geo_appro;
};

class xPointwiseCommand : public xCommandOnGeomElem
{
  public:
   xPointwiseCommand(xPointwiseVisitor& vi, xSpacePointwise& s, xVariabManager& v) : visitor(vi), space(s), var_manager(v) {}
   void execute(xGeomElem* geo_integ) override
   {
      // const bool debug=false;
      const int nb = geo_integ->GetNbIntegrationPoints();
      for (int k = 0; k < nb; k++)
      {
         geo_integ->setUVW(k);
         if (geom_appro != geo_integ)
            geom_appro->setUVWForXYZ(geo_integ->getXYZ());
         else
            geom_appro->setUVW(k);
         xValKey key;
         space.getKey(geo_integ, key);
         xValue<tensorsPtr_t>* v = var_manager.find(key);
         visitor.setInfos(geo_integ, geom_appro);
         visitor.Visit(*v);
      }
   }

  private:
   xPointwiseVisitor& visitor;
   xSpacePointwise& space;
   xVariabManager& var_manager;
};

template <typename T>
class xSetFieldPointwiseVisitor : public xPointwiseVisitor
{
  public:
   xSetFieldPointwiseVisitor(std::string& variable_, const xEval<T>& given_) : variable(variable_), given(given_) {}

   void Visit(xValue<tensorsPtr_t>& v) override
   {
      tensorsPtr_t tensors = v.getVal();
      T& r = tensors->get(variable, boost::mpl::identity<T>());
      given(geo_appro, geo_integ, r);
   }

  private:
   std::string& variable;
   const xEval<T>& given;
};

template <typename T>
class xSetFieldPointwiseCheckChangeVisitor : public xPointwiseVisitor
{
  public:
   xSetFieldPointwiseCheckChangeVisitor(std::string& variable_, const xEval<T>& given_)
       : variable(variable_), given(given_), number_of_changes(0), total_number(0)
   {
   }

   void Visit(xValue<tensorsPtr_t>& v) override
   {
      tensorsPtr_t tensors = v.getVal();
      T& r_previous = tensors->get(variable, boost::mpl::identity<T>());
      T temp = r_previous;
      given(geo_appro, geo_integ, r_previous);
      if (r_previous != temp) number_of_changes++;
      total_number++;
   }
   int getNumberOfChanges() { return number_of_changes; }
   double getChangeRatio()
   {
      double ratio = static_cast<double>(number_of_changes);
      ratio /= static_cast<double>(total_number);
      return ratio;
   }

  private:
   std::string& variable;
   const xEval<T>& given;
   int number_of_changes, total_number;
};

template <typename T>
class xAddToFieldPointwiseVisitor : public xPointwiseVisitor
{
  public:
   xAddToFieldPointwiseVisitor(std::string& variable_, const xEval<T>& given_) : variable(variable_), given(given_) {}

   void Visit(xValue<tensorsPtr_t>& v) override
   {
      tensorsPtr_t tensors = v.getVal();
      T& r = tensors->get(variable, boost::mpl::identity<T>());
      T res;
      given(geo_appro, geo_integ, res);
      r += res;
   }

  private:
   std::string& variable;
   const xEval<T>& given;
};

}  // namespace xfem

#endif
