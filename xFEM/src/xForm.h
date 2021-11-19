/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _FORMS_BASE_H
#define _FORMS_BASE_H

#include <iostream>
#include <vector>

#include "xApproxFctPtr.h"
#include "xDebug.h"
#include "xFemMatrix.h"
#include "xGeomElem.h"
#include "xTensor2.h"
#include "xTensorOperations.h"
#include "xVector.h"

namespace xfem
{
class xValKey;
class xApproxFunction;

namespace xFormProducts
{
/**
   product (vector vector)
 */
template <typename VT>
inline static VT product(const xtensor::xVector<VT>& a, const xtensor::xVector<VT>& b)
{
   return a * b;
}

//  //Cross-type product: double-float
inline static float product(const xtensor::xVector<double>& a, const xtensor::xVector<float>& b)
{
   return static_cast<float>(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/**
   product (double double)
 */
inline static double product(double a, double b) { return a * b; }

inline static std::complex<double> product(std::complex<double> a, std::complex<double> b) { return a * b; }

// Cross product double-complex
inline static std::complex<double> product(double a, std::complex<double> b) { return a * b; }
/**
   product (double double double)
 */
template <typename VT>
inline static VT product(VT a, VT b, VT c)
{
   return a * b * c;
}
/**
   product (vector double vector)
 */
template <typename VT>
inline static VT product(const xtensor::xVector<VT>& a, VT b, const xtensor::xVector<VT>& c)
{
   return (a * b) * c;
}
/**
   product (double double double)
 */
template <typename VT>
inline static VT product(const xtensor::xVector<VT>& a, const xtensor::xTensor2<VT>& b, const xtensor::xVector<VT>& c)
{
   const bool debug = false;
   if (debug) std::cout << "a " << a << std::endl;
   if (debug) std::cout << "b " << b << std::endl;
   if (debug) std::cout << "c " << c << std::endl;
   return a * (b * c);
}
/**
   product (vector vector double)
 */
template <typename VT>
inline static double product(const xtensor::xVector<VT>& a, const xtensor::xVector<VT>& b, double c)
{
   return (a * b) * c;
}
/**
   product (double vector vector)
 */
template <typename VT>
inline static double product(double a, const xtensor::xVector<VT>& b, const xtensor::xVector<VT>& c)
{
   return a * (b * c);
}

template <typename VT>
inline static double product(const xtensor::xTensor2<VT>& a, const xtensor::xTensor2<VT>& c)
{
   return a.contract(c);
   // return (a(0,0)*c(0,0) + a(1,1)*c(1,1) + a(2,2)*c(2,2) + a(0,1)*c(0,1) + a(1,0)*c(1,0)
   //          +	a(0,2)*c(0,2) + a(2,0)*c(2,0) + a(1,2)*c(1,2) + a(2,1)*c(2,1));
}

/**
   product (double xTensor2 xTensor2)
 */
template <typename VT>
inline static double product(const xtensor::xTensor2<VT>& a, double b, const xtensor::xTensor2<VT>& c)
{
   return b * product(a, c);
}

template <typename VT>
inline static VT product(const xtensor::xTensor2Sym<VT>& a, const xtensor::xTensor2Sym<VT>& c)
{
   return a.contract(c);
}

template <typename VT>
inline static double product(const xtensor::xTensor2<VT>& a, const xtensor::xTensor4<VT>& b, const xtensor::xTensor2<VT>& c)
{
   //    return product((b*a),1.0,c);

   return product(a, 1.0, (b * c));
   // std::cout << " bla " << std::endl;
   // return product(a,b*c);
}

template <typename VT>
inline static double product(const xtensor::xTensor2<VT>& a, const xtensor::xTensor4Isotropic& b, const xtensor::xTensor2<VT>& c)
{
   //    return product((b*a),1.0,c);

   return product(a, 1.0, (b * c));
   // std::cout << " bla " << std::endl;
   // return product(a,b*c);
}

template <typename VT>
inline static double product(const xtensor::xTensor2<VT>& a, const xtensor::xTensor2<VT>& b, const xtensor::xTensor2<VT>& c)
{
   return product((a * b), 1.0, c);
}
}  // namespace xFormProducts

bool integ_is_appro(xGeomElem* geo_integ, xGeomElem* geo_appro);

inline void set_appro_uvw_integ_is_same(xGeomElem* geo_integ, xGeomElem* geo_appro) { geo_appro->setUVW(geo_integ->getUVW()); };

inline void set_appro_uvw_integ_is_not_same(xGeomElem* geo_integ, xGeomElem* geo_appro)
{
   geo_appro->setUVWForXYZ(geo_integ->getXYZ());
};

//! Form base class
class xForm
{
  public:
   typedef std::vector<xValKey> femKeys;
   typedef approxFctPtr_t shapeFctPtr;
   typedef femFcts_t femFcts;

  public:
   virtual ~xForm();
   virtual void accumulate(xGeomElem* geo_integ) = 0;

  protected:
};

//! FormBilinear base class
///  \class xFormBilinear xForm.h "xFormBilinear base class"
///  \fn init : initialize geo_appro_test, geo_appro_trial, test et trial, la taille de Matrix.
///  \fn accumulate : loop over gauss points and call accumulate_pnt to compute the matrix and accumulate it. (called by the
///  integrator, via a commandOnGeomElem)

template <typename VT = double>
class xFormBilinear : public xForm
{
  public:
   xFormBilinear() {}
   ~xFormBilinear() override = 0;
   virtual void init(xGeomElem* geo_appro_test, femFcts* test, xGeomElem* geo_appro_trial, femFcts* trial);
   virtual void init(xGeomElem* geo_appro, femFcts* test);
   void accumulate(xGeomElem* geo_integ) override;
   virtual void finalize() = 0;
   const xFemMatrix<VT>& getMatrix() const { return Matrix; }
   xFemMatrix<VT>& getMatrix() { return Matrix; }
   // virtual bool isSym(femFcts* test_, femFcts* trial_){return false;}
  protected:
   virtual void accumulate_pnt(xGeomElem* geo_integ) = 0;
   femFcts* test;
   femFcts* trial;
   xFemMatrix<VT> Matrix;
   xGeomElem *geo_appro_test, *geo_appro_trial;

  private:
};

template <typename OperatorLeft, typename MaterialLaw, typename OperatorRight, typename VT = double>
class xFormBilinearWithLaw : public xFormBilinear<VT>
{
  public:
   xFormBilinearWithLaw(MaterialLaw& law) : Left(OperatorLeft()), Right(OperatorRight()), Law(law) {}
   xFormBilinearWithLaw(OperatorLeft& l, MaterialLaw& law) : Left(l), Law(law), Right(OperatorRight()) {}
   xFormBilinearWithLaw(OperatorLeft& l, MaterialLaw& law, OperatorRight& r) : Left(l), Law(law), Right(r) {}
   ~xFormBilinearWithLaw() override = default;
   void finalize() override {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      Law(this->geo_appro_test, geo_integ, phys);
      std::vector<typename OperatorLeft::result_type> values_left;
      std::vector<typename OperatorRight::result_type> values_right;
      const int test_size = this->test->size();
      const int trial_size = this->trial->size();
      values_left.reserve(test_size);
      values_right.reserve(trial_size);
      Left.eval(this->test, this->geo_appro_test, geo_integ, values_left);
      Right.eval(this->trial, this->geo_appro_trial, geo_integ, values_right);
      const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();

      for (size_t i = 0; i < test_size; i++)
      {
         const typename OperatorRight::result_type tmp = values_left[i] * phys;
         for (size_t j = 0; j < trial_size; ++j)
         {
            VT x = xFormProducts::product(tmp, values_right[j]) * wdet;
            this->Matrix.add(i, j, x);
         }
      }
   }

  private:
   const OperatorLeft Left;
   MaterialLaw& Law;
   const OperatorRight Right;
   typename MaterialLaw::result_type phys;
};

template <typename OperatorLeft, typename MaterialLaw, typename VT>
class xFormBilinearWithLaw<OperatorLeft, MaterialLaw, OperatorLeft, VT> : public xFormBilinear<VT>
{
  public:
   xFormBilinearWithLaw(MaterialLaw& law, bool sym = true) : Left(OperatorLeft()), Right(OperatorLeft()), Law(law), sym_ope(sym)
   {
   }
   xFormBilinearWithLaw(OperatorLeft& l, MaterialLaw& law, bool sym = true) : Left(l), Right(l), Law(law), sym_ope(sym) {}
   xFormBilinearWithLaw(OperatorLeft& l, MaterialLaw& law, OperatorLeft& r) : Left(l), Right(r), Law(law), sym_ope(false) {}
   ~xFormBilinearWithLaw() override = default;
   void finalize() override
   {
      if (this->test == this->trial && sym_ope) this->Matrix.symmetrize();
   }

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Law(this->geo_appro_test, geo_integ, phys);
      if (this->test == this->trial && sym_ope)
      {
         const size_t test_size = this->test->size();
         std::vector<typename OperatorLeft::result_type> values_left;
         values_left.reserve(test_size);
         Left.eval(this->test, this->geo_appro_test, geo_integ, values_left);
         const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();

         for (size_t i = 0; i < test_size; ++i)
         {
            const typename OperatorLeft::result_type tmp = values_left[i] * phys;
            if (debug) std::cout << "LEFT\n";
            if (debug) std::cout << values_left[i] << std::endl;
            for (size_t j = i; j < test_size; ++j)
            {
               VT x = xFormProducts::product(tmp, values_left[j]) * wdet;
               if (debug) std::cout << "On ajoute " << x << " En (" << i << "," << j << ")" << std::endl;
               this->Matrix.addSym(i, j, x);
               if (debug) std::cout << "Valeur=" << x << std::endl;
            }
         }
      }
      else
      {
         std::vector<typename OperatorLeft::result_type> values_left;
         std::vector<typename OperatorLeft::result_type> values_right;
         const size_t test_size = this->test->size();
         const size_t trial_size = this->trial->size();
         values_left.reserve(test_size);
         values_right.reserve(trial_size);
         Left.eval(this->test, this->geo_appro_test, geo_integ, values_left);
         Right.eval(this->trial, this->geo_appro_trial, geo_integ, values_right);

         const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
         for (size_t i = 0; i < test_size; ++i)
         {
            const typename OperatorLeft::result_type tmp = values_left[i] * phys;
            for (size_t j = 0; j < trial_size; ++j)
            {
               VT x = xFormProducts::product(tmp, values_right[j]) * wdet;
               this->Matrix.add(i, j, x);
            }
         }
      }
   }

  private:
   const OperatorLeft Left, Right;
   MaterialLaw& Law;
   typename MaterialLaw::result_type phys;
   bool sym_ope;
};

template <typename OperatorLeft, typename OperatorRight, typename VT = double>
class xFormBilinearWithoutLaw : public xFormBilinear<VT>
{
  public:
   xFormBilinearWithoutLaw() : Left(OperatorLeft()), Right(OperatorRight()) {}
   xFormBilinearWithoutLaw(OperatorLeft& l, OperatorRight& r) : Left(l), Right(r) {}
   xFormBilinearWithoutLaw(OperatorLeft& l) : Left(l) {}
   xFormBilinearWithoutLaw(OperatorRight& r) : Right(r) {}
   ~xFormBilinearWithoutLaw() override = default;
   void finalize() override {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      std::vector<typename OperatorLeft::result_type> values_left;
      std::vector<typename OperatorRight::result_type> values_right;
      Left.eval(this->test, this->geo_appro_test, geo_integ, values_left);
      Right.eval(this->trial, this->geo_appro_trial, geo_integ, values_right);
      const size_t test_size = this->test->size();
      const size_t trial_size = this->trial->size();
      const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();

      for (size_t i = 0; i < test_size; ++i)
      {
         const typename OperatorRight::result_type tmp = values_left[i];
         for (size_t j = 0; j < trial_size; ++j)
         {
            VT x = xFormProducts::product(tmp, values_right[j]) * wdet;
            this->Matrix.add(i, j, x);
         }
      }
   }

  private:
   const OperatorLeft Left;
   const OperatorRight Right;
};

template <typename OperatorLeft, typename VT>
class xFormBilinearWithoutLaw<OperatorLeft, OperatorLeft, VT> : public xFormBilinear<VT>
{
  public:
   xFormBilinearWithoutLaw() : Left(OperatorLeft()), Right(OperatorLeft()), OpleftequalOpRight(true) {}
   xFormBilinearWithoutLaw(OperatorLeft& l) : Left(l), Right(l), OpleftequalOpRight(true) {}
   xFormBilinearWithoutLaw(OperatorLeft& l, OperatorLeft& r) : Left(l), Right(r), OpleftequalOpRight(false) {}
   ~xFormBilinearWithoutLaw() override = default;
   void finalize() override
   {
      if (this->test == this->trial && OpleftequalOpRight) this->Matrix.symmetrize();
   }

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      if ((this->test == this->trial) && OpleftequalOpRight)
      {
         const size_t test_size = this->test->size();
         std::vector<typename OperatorLeft::result_type> values_left;
         values_left.reserve(test_size);
         Left.eval(this->test, this->geo_appro_test, geo_integ, values_left);
         const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
         for (size_t i = 0; i < test_size; ++i)
         {
            const typename OperatorLeft::result_type tmp = values_left[i];
            for (size_t j = i; j < test_size; ++j)
            {
               VT x = xFormProducts::product(tmp, values_left[j]) * wdet;
               this->Matrix.addSym(i, j, x);
            }
         }
      }
      else
      {
         const size_t test_size = this->test->size();
         const size_t trial_size = this->trial->size();
         std::vector<typename OperatorLeft::result_type> values_left;
         std::vector<typename OperatorLeft::result_type> values_right;
         values_left.reserve(test_size);
         values_right.reserve(trial_size);
         Left.eval(this->test, this->geo_appro_test, geo_integ, values_left);
         Right.eval(this->trial, this->geo_appro_trial, geo_integ, values_right);
         const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
         for (size_t i = 0; i < test_size; ++i)
         {
            const typename OperatorLeft::result_type tmp = values_left[i];
            for (size_t j = 0; j < trial_size; ++j)
            {
               VT x = xFormProducts::product(tmp, values_right[j]) * wdet;
               this->Matrix.add(i, j, x);
            }
         }
      }
   }

  private:
   const OperatorLeft Left, Right;
   bool OpleftequalOpRight;
};

template <typename VT = double>
class xFormLinear : public xForm
{
  public:
   ~xFormLinear() override = default;
   void init(xGeomElem* geo_appro, femFcts* te);
   const xFemVector<VT>& getVector() const { return Vector; }
   xFemVector<VT>& getVector() { return Vector; }
   void accumulate(xGeomElem* geo_integ) override;
   void finalize() {}

  protected:
   virtual void accumulate_pnt(xGeomElem* geo_integ) = 0;
   femFcts* test;
   xGeomElem* geo_appro;
   xFemVector<VT> Vector;

  private:
};

template <typename OperatorTest, typename PhysicalLoad, typename VT = double>
class xFormLinearWithLoad : public xFormLinear<VT>
{
  public:
   xFormLinearWithLoad(const PhysicalLoad& l) : Test(OperatorTest()), Load(l) {}
   xFormLinearWithLoad(OperatorTest& t, const PhysicalLoad& l) : Test(t), Load(l) {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Load(this->geo_appro, geo_integ, load);
      if (debug) std::cout << "load is " << load << std::endl;
      std::vector<typename OperatorTest::result_type> values_test;
      Test.eval(this->test, this->geo_appro, geo_integ, values_test);
      const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
      const size_t test_size = this->test->size();
      for (size_t i = 0; i < test_size; ++i)
      {
         VT x = xFormProducts::product(values_test[i], load) * wdet;  // Implicit cast ??
         if (debug) std::cout << "x in linear form with load is " << x << std::endl;
         if (debug) std::cout << "value test[" << i << "]=" << std::endl << values_test[i] << std::endl;
         this->Vector.add(i, x);
      }
   }

  private:
   const OperatorTest Test;
   const PhysicalLoad& Load;
   typename PhysicalLoad::result_type load;
};

template <typename VT = double>
class xFormZero : public xForm
{
  public:
   ~xFormZero() override = default;
   const xFemScalar<VT>& getScalar() const { return Scalar; }
   xFemScalar<VT>& getScalar() { return Scalar; }
   void init(xGeomElem* geo_appro);
   void accumulate(xGeomElem* geo_integ) override;
   void finalize() {}

  protected:
   virtual void accumulate_pnt(xGeomElem* geo_integ) = 0;
   xFemScalar<VT> Scalar;
   xGeomElem* geo_appro;

  private:
};

template <typename VT = double>
class xFormZeroUnit : public xFormZero<VT>
{
  public:
   xFormZeroUnit() {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const VT val = xtool::xDataType<VT>::one();
      double x = geo_integ->GetWeight() * val * geo_integ->GetDetJac();
      this->Scalar.add(x);
   }

  private:
};

template <typename Eval, typename VT = double>
class xFormZeroEvalLinearForm : public xFormZero<VT>
{
  private:
   Eval& ev;
   typename Eval::result_type res;

  public:
   xFormZeroEvalLinearForm(Eval& ev_) : ev(ev_) {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      ev(this->geo_appro, geo_integ, res);
      this->Scalar.add(geo_integ->GetWeight() * res * geo_integ->GetDetJac());
   }
};

template <typename EvalLeft, typename MaterialLaw, typename OperatorTest>
class xFormLinearEvalBilinearFormWithLaw : public xFormLinear<double>
{
  private:
   const OperatorTest Test;
   EvalLeft& Left;
   typename EvalLeft::result_type left;
   MaterialLaw& Law;
   typename MaterialLaw::result_type phys;

  public:
   xFormLinearEvalBilinearFormWithLaw(EvalLeft& left, MaterialLaw& law, OperatorTest& t) : Left(left), Law(law), Test(t) {}

   xFormLinearEvalBilinearFormWithLaw(EvalLeft& left, MaterialLaw& law) : Left(left), Law(law), Test(OperatorTest()) {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Law(geo_appro, geo_integ, phys);
      Left(geo_appro, geo_integ, left);
      const typename OperatorTest::result_type tmp = left * phys;
      std::vector<typename OperatorTest::result_type> values_test;
      Test.eval(test, geo_appro, geo_integ, values_test);
      const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
      const int test_size = test->size();
      for (size_t i = 0; i < test_size; ++i)
      {
         double x = xFormProducts::product(tmp, values_test[i]) * wdet;
         if (debug) std::cout << "x in linear form with load is " << x << std::endl;
         if (debug) std::cout << "value test[" << i << "]=" << std::endl << values_test[i] << std::endl;
         Vector.add(i, x);
      }
   }
};

template <typename EvalLeft, typename OperatorTest>
class xFormLinearEvalBilinearFormWithoutLaw : public xFormLinear<double>
{
  private:
   const OperatorTest Test;
   EvalLeft& Left;
   typename EvalLeft::result_type left;

  public:
   xFormLinearEvalBilinearFormWithoutLaw(EvalLeft& left, OperatorTest& t) : Left(left), Test(t) {}

   xFormLinearEvalBilinearFormWithoutLaw(EvalLeft& left) : Left(left), Test(OperatorTest()) {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Left(geo_appro, geo_integ, left);
      const typename OperatorTest::result_type tmp = left;
      std::vector<typename OperatorTest::result_type> values_test;
      Test.eval(test, geo_appro, geo_integ, values_test);
      const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
      const int test_size = test->size();
      for (size_t i = 0; i < test_size; ++i)
      {
         double x = xFormProducts::product(tmp, values_test[i]) * wdet;
         if (debug) std::cout << "x in linear form with load is " << x << std::endl;
         if (debug) std::cout << "value test[" << i << "]=" << std::endl << values_test[i] << std::endl;
         Vector.add(i, x);
      }
   }
};

template <typename EvalLeft, typename MaterialLaw, typename EvalRight, typename VT = double>
class xFormZeroEvalBilinearFormWithLaw : public xFormZero<VT>
{
  private:
   EvalLeft& Left;
   typename EvalLeft::result_type left;
   MaterialLaw& Law;
   typename MaterialLaw::result_type phys;
   EvalRight& Right;
   typename EvalRight::result_type right;

  public:
   xFormZeroEvalBilinearFormWithLaw(EvalLeft& left, MaterialLaw& law, EvalRight& right) : Left(left), Law(law), Right(right) {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Law(this->geo_appro, geo_integ, phys);
      Left(this->geo_appro, geo_integ, left);
      Right(this->geo_appro, geo_integ, right);
      VT val = xFormProducts::product(left, phys, right);
      if (debug) std::cout << "left " << left << std::endl;
      if (debug) std::cout << "phys " << phys << std::endl;
      if (debug) std::cout << "right " << right << std::endl;
      VT x = geo_integ->GetWeight() * val * geo_integ->GetDetJac();
      this->Scalar.add(x);
   }
};

template <typename EvalLeft, typename MaterialLaw, typename VT>
class xFormZeroEvalBilinearFormWithLaw<EvalLeft, MaterialLaw, EvalLeft, VT> : public xFormZero<VT>
{
  private:
   EvalLeft &Left, &Right;
   typename EvalLeft::result_type left, right;
   MaterialLaw& Law;
   typename MaterialLaw::result_type phys;
   bool symm;

  public:
   xFormZeroEvalBilinearFormWithLaw(EvalLeft& left, MaterialLaw& law) : Left(left), Right(left), Law(law), symm(true) {}

   xFormZeroEvalBilinearFormWithLaw(EvalLeft& left, MaterialLaw& law, EvalLeft& right)
       : Left(left), Right(right), Law(law), symm(false)
   {
   }

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Law(this->geo_appro, geo_integ, phys);
      Left(this->geo_appro, geo_integ, left);
      VT val;
      if (symm)
      {
         val = xFormProducts::product(left, phys, left);
         if (debug) std::cout << "left " << left << std::endl;
      }
      else
      {
         Right(this->geo_appro, geo_integ, right);
         val = xFormProducts::product(left, phys, right);
         if (debug) std::cout << "left " << left << std::endl;
         if (debug) std::cout << "right " << right << std::endl;
      }
      this->Scalar.add(geo_integ->GetWeight() * val * geo_integ->GetDetJac());
   }
};

template <typename EvalLeft, typename EvalRight, typename VT = double>
class xFormZeroEvalBilinearFormWithoutLaw : public xFormZero<VT>
{
  private:
   EvalLeft& Left;
   typename EvalLeft::result_type left;
   EvalRight& Right;
   typename EvalRight::result_type right;

  public:
   xFormZeroEvalBilinearFormWithoutLaw(EvalLeft& left, EvalRight& right) : Left(left), Right(right) {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Left(this->geo_appro, geo_integ, left);
      Right(this->geo_appro, geo_integ, right);
      VT val = xFormProducts::product(left, right);
      if (debug) std::cout << "left " << left << std::endl;
      if (debug) std::cout << "right " << right << std::endl;
      this->Scalar.add(geo_integ->GetWeight() * val * geo_integ->GetDetJac());
   }
};

template <typename EvalLeft, typename VT>
class xFormZeroEvalBilinearFormWithoutLaw<EvalLeft, EvalLeft, VT> : public xFormZero<VT>
{
  private:
   EvalLeft &Left, &Right;
   typename EvalLeft::result_type left, right;
   bool symm;

  public:
   xFormZeroEvalBilinearFormWithoutLaw(EvalLeft& left) : Left(left), Right(left), symm(true) {}

   xFormZeroEvalBilinearFormWithoutLaw(EvalLeft& left, EvalLeft& right) : Left(left), Right(right), symm(false) {}

  protected:
   void accumulate_pnt(xGeomElem* geo_integ) override
   {
      const bool debug = false;
      Left(this->geo_appro, geo_integ, left);
      VT val;
      if (symm)
      {
         val = xFormProducts::product(left, left);
         if (debug) std::cout << "left " << left << std::endl;
      }
      else
      {
         Right(this->geo_appro, geo_integ, right);
         val = xFormProducts::product(left, right);
         if (debug) std::cout << "left " << left << std::endl;
         if (debug) std::cout << "right " << right << std::endl;
      }
      this->Scalar.add(geo_integ->GetWeight() * val * geo_integ->GetDetJac());
   }
};

//----------------------------------------------------------------------

//   template<typename VT>
//   void xFormBilinear::accumulate(xGeomElem*  geo_integ) {
//     //const bool debug = false;
//     typedef std::function < void (xGeomElem * , xGeomElem *  )>   set_appro_t;
//     set_appro_t set_appro_test_uvw =
//       integ_is_appro(geo_integ, geo_appro_test)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
//     set_appro_t set_appro_trial_uvw =
//       integ_is_appro(geo_integ, geo_appro_trial)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;

//     const int nb = geo_integ->GetNbIntegrationPoints();
//     for(int k = 0; k < nb; k++){
//       geo_integ->setUVW(k);
//       set_appro_test_uvw(geo_integ, geo_appro_test);
//       set_appro_trial_uvw(geo_integ, geo_appro_trial);
//       accumulate_pnt(geo_integ);
//     }
//     return;
//   }

template <typename VT>
void xFormLinear<VT>::init(xGeomElem* g, femFcts* te)
{
   this->geo_appro = g;
   this->test = te;
   this->Vector.resize(this->test->size());
}

template <typename VT>
void xFormLinear<VT>::accumulate(xGeomElem* geo_integ)
{
   const bool debug = false;
   const unsigned int nb = geo_integ->GetNbIntegrationPoints();
   typedef std::function<void(xGeomElem*, xGeomElem*)> set_appro_t;
   set_appro_t set_appro_uvw =
       integ_is_appro(geo_integ, geo_appro) ? set_appro_uvw_integ_is_same : set_appro_uvw_integ_is_not_same;
   for (unsigned int k = 0; k < nb; k++)
   {
      geo_integ->setUVW(k);
      set_appro_uvw(geo_integ, geo_appro);
      accumulate_pnt(geo_integ);
   }
   if (debug) std::cout << "Vector of the Linear form" << std::endl;
   if (debug) std::cout << Vector << std::endl;
}

template <typename VT>
void xFormBilinear<VT>::init(xGeomElem* g_te, femFcts* te, xGeomElem* g_tr, femFcts* tr)
{
   this->geo_appro_test = g_te;
   this->test = te;
   this->geo_appro_trial = g_tr;
   this->trial = tr;
   this->Matrix.resize(this->test->size(), this->trial->size());
}

template <typename VT>
void xFormBilinear<VT>::init(xGeomElem* g, femFcts* t)
{
   this->geo_appro_test = g;
   this->geo_appro_trial = g;
   this->test = t;
   this->trial = t;
   this->Matrix.resize(this->test->size(), this->trial->size());
}

template <typename VT>
xFormBilinear<VT>::~xFormBilinear() = default;

template <typename VT>
void xFormBilinear<VT>::accumulate(xGeomElem* geo_integ)
{
   // const bool debug = false;
   typedef std::function<void(xGeomElem*, xGeomElem*)> set_appro_t;
   set_appro_t set_appro_test_uvw =
       integ_is_appro(geo_integ, geo_appro_test) ? set_appro_uvw_integ_is_same : set_appro_uvw_integ_is_not_same;
   set_appro_t set_appro_trial_uvw =
       integ_is_appro(geo_integ, geo_appro_trial) ? set_appro_uvw_integ_is_same : set_appro_uvw_integ_is_not_same;

   const int nb = geo_integ->GetNbIntegrationPoints();
   for (int k = 0; k < nb; k++)
   {
      geo_integ->setUVW(k);
      set_appro_test_uvw(geo_integ, geo_appro_test);
      set_appro_trial_uvw(geo_integ, geo_appro_trial);
      accumulate_pnt(geo_integ);
   }
   return;
}

template <typename VT>
void xFormZero<VT>::init(xGeomElem* g)
{
   geo_appro = g;
   Scalar.setVal(0.0);
}

template <typename VT>
void xFormZero<VT>::accumulate(xGeomElem* geo_integ)
{
   const int nb = geo_integ->GetNbIntegrationPoints();
   typedef std::function<void(xGeomElem*, xGeomElem*)> set_appro_t;
   set_appro_t set_appro_uvw =
       integ_is_appro(geo_integ, geo_appro) ? set_appro_uvw_integ_is_same : set_appro_uvw_integ_is_not_same;
   for (int k = 0; k < nb; k++)
   {
      geo_integ->setUVW(k);
      set_appro_uvw(geo_integ, geo_appro);
      accumulate_pnt(geo_integ);
   }
}

//----------------------------------------------------------------------

}  // namespace xfem

#endif
