/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef  _FORMS_OPTIMIZED_H
#define  _FORMS_OPTIMIZED_H

#include "xForm.h"
#include "xEnv.h"
#include "xSpace.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xValue.h"
#include "xField.h"
#include "xTensor4.h"
#include "xDenseMatrix.h"

namespace xfem{
  //! This template function is meant to give a uniform way to obtain a pointer to the data contain in a double, a xtensor::xVector, 
  //!    a xtensor::xTensor2 or a xtensor::xTensor4. more type of data could be added by adding new template specialization.
  template <class T > double * getData(T& in);
  
  //! This template function is meant to give a uniform way to obtain a const pointer to the data contain in 
  //! a const  double, a  const xtensor::xVector,    a const  xtensor::xTensor2 or a xtensor::xTensor4.
  //!  more type of data could be added by adding new template specialization.
  template <class T > const double * getData(const T& in);
  
  //! This template function is meant to give a uniform way to obtain a the size of the data contained in a
  //! a double, a   xtensor::xVector,    a   xtensor::xTensor2 or  xtensor::xTensor4.
  //!  more type of data could be added by adding new template specialization.
  template <class T > size_t getSize();
  
  //! This function is meant to fill up the upper triangle of a dense matrix with 
  //! the lower triangle, to obtain a symmetric matrix
  inline void matrixcopylowtriangletouptriangle(xlinalg::xDenseMatrix & Matrix){
    for (int i =0; i < Matrix.nline();++i){
      for (int j =0; j <  i; ++j){
	Matrix(j,i) = Matrix(i,j);
      }
    }
  };

  //! This class is meant to be a base class for xFormBilinearOptimizedWithoutLaw and xFormBilinearOptimizedWithLaw
  //! It implementent the common function and data structures used by the two previous classes.
  //! It should not be use directly that's why the constructor is protected.
  template < typename OPERATORTEST, typename OPERATORTRIAL >
    class xFormBilinearOptimizedBase : public xForm{
  public:
    typedef xlinalg::xDenseMatrix matrix_t;
    //! This function initialize the form in cases where test and trial are the same.
    void init(xGeomElem* _geo_appro_test, femFcts* _test);  
    //! This function initialize the form in cases where test and trial are differents.
    void init(xGeomElem* _geo_appro_test, femFcts* _test, xGeomElem* _geo_appro_trial , femFcts* _trial);
    //! This function return the internal matrix. meant to be called after a call to finalise.
    const matrix_t& getMatrix() const;
    //! This function return the internal matrix. meant to be called after a call to finalise.
    matrix_t& getMatrix();
  protected:
    xFormBilinearOptimizedBase(OPERATORTEST& o);
    xFormBilinearOptimizedBase(OPERATORTEST& otest,  OPERATORTRIAL &otrial);
    void accumulate_test_is_trial_no_law(xGeomElem*  geo_integ);  
    void accumulate_test_is_not_trial_no_law(xGeomElem*  geo_integ);
    template < typename MATERIALLAW > 
      void accumulate_test_is_trial_law(xGeomElem*  geo_integ, const MATERIALLAW &Law);
    template < typename MATERIALLAW >
      void accumulate_test_is_not_trial_law(xGeomElem*  geo_integ, const MATERIALLAW &Law); 
    //! the vector of current test functions
    femFcts* test;
    //! the vector of current trial functions
    femFcts* trial;
    //! The elementary matrix.
    matrix_t Matrix;
    //! contain the element on which the test space is defined.
    xGeomElem* geo_appro_test; 
    //! contain the element on which the trial space  is defined.
    xGeomElem* geo_appro_trial; 
    //! Operator applied to test function.
    const OPERATORTEST   OperatorTest;
    //! Operator Applied to trial function.
    const OPERATORTRIAL  OperatorTrial;
    //! A bool set at construction which is true if both operator are the same.
    bool optestisoptrial;
  };

  //! This class represent a bilinear form meant to be used by the Assemble algorithm.
  //! It does the same as xFormBilinearWithoutLaw but faster in cases
  //! where the space and the number of integration point is large enought,
  //! by using a BLAS implementation of the integration loop called by accumulate.
  //! Note that this version is not tested yet.
  template < typename OPERATORTEST,  typename OPERATORTRIAL >
    class xFormBilinearOptimizedWithoutLaw : public xFormBilinearOptimizedBase< OPERATORTEST, OPERATORTRIAL >  {
  public :
    xFormBilinearOptimizedWithoutLaw(OPERATORTEST& o);
    xFormBilinearOptimizedWithoutLaw(OPERATORTEST& otest, OPERATORTRIAL &otrial);
    void accumulate(xGeomElem*  geo_integ);
    void finalize();
  };


  //! This class represent a bilinear form meant to be used by the Assemble algorithm.
  //! It does the same as xFormBilinearLaw but faster in cases
  //! where the space and the number of integration point is large enought,
  //! by using a BLAS implementation of the integration loop called by accumulate.
  //! This version is tested in Xtest/xfem-seq-test/OptimizedForm.
  template < typename OPERATORTEST, typename MATERIALLAW, typename OPERATORTRIAL >
    class xFormBilinearOptimizedWithLaw : public xFormBilinearOptimizedBase< OPERATORTEST, OPERATORTRIAL >  {
  public:
    xFormBilinearOptimizedWithLaw(OPERATORTEST& otest,  MATERIALLAW& l, OPERATORTRIAL &otrial);
    xFormBilinearOptimizedWithLaw(OPERATORTEST& otest,  MATERIALLAW& l);
    void accumulate(xGeomElem*  geo_integ) override;
    void finalize();
  private :
    MATERIALLAW &Law;
  };

  
  //! This class is a specialized version of xFormBilinearWithLaw : it expect same space and same operator
  //! to right and left, and the Law must be symmetric positive.
  //! Op (u)T D Op (u')  (for example Grad uT : D : Grad u' ) is replaced by (Op(u):U)T : (U:Op(u))
  //! Which is coded using a rank k update (syrk) which demand less operation Than
  //! Op( u) *(D * Op( u) ), provided that we have access to U where U^TU = D. (A cholesky factorisation)
  //! This version is tested in Xtest/xfem-seq-test/OptimizedForm.
  template < typename  OPERATOR, typename MATERIALLAW >
    class xFormBilinearOptimizedSymmetricWithLawFactored : public xForm{
  public : 
    typedef  xlinalg::xDenseMatrix matrix_t;
    //! This is the constructor. Note that only 1 operator is accepted (This form is meant to deal ONLY with symmetric bilinear forms.)
    //! The law must be an xEval that return the factor U of the tensor that define the law. (Upper version : D= UTU) 
    xFormBilinearOptimizedSymmetricWithLawFactored(OPERATOR& o, MATERIALLAW& l);
    //! Calling this init will throw an exception : Again this bilinearForm is meant to treat only symmetric problems.
    //! It is there only so that the Assemble Algorithm Compile with this bilinear form.
    void init(xGeomElem* _geo_appro_test, femFcts* _test, xGeomElem* _geo_appro_trial, femFcts* _trial);  
    //! intitialize the Bilinear form when test and trial are the same.
    void init(xGeomElem* _geo_appro_test, femFcts* _test);   
    void accumulate(xGeomElem*  geo_integ) override;
    //! Finalize the assembly (In particular it fill up the lower part of the local finite element matrix Which is supposed to be symmetric, but only the lower part is computed when calling accumulate)
    void finalize();  
    const matrix_t& getMatrix() const;
    matrix_t& getMatrix();
  private:
    femFcts* test;
    matrix_t  Matrix;
    xGeomElem* geo_appro_test;  
    const OPERATOR  Operator;
    MATERIALLAW &  Law;
    //! setUGMatk_XX are private helper function that are meant to set up all the data before the call to syrk.
    //! So far there are 3 differents versions, With differents tread off for optimization... 
    //! by default it's currently _v1 that is used.
    void setUGMatk_v0 (const size_t k, const size_t &nb, const size_t &size_tensor, const  size_t &test_size, 
		       const std::vector<typename OPERATOR::result_type> &values,
		       const double &sqrtwdet,   const typename MATERIALLAW::result_type &uphys,  matrix_t &UGMat);
    
    void setUGMatk_v1 (const size_t k, const size_t &nb, const size_t &size_tensor, const  size_t &test_size, 
		       const std::vector<typename OPERATOR::result_type> &values,
		       const double &sqrtwdet,   const typename MATERIALLAW::result_type &uphys,  matrix_t &UGMat);
    
    void setUGMatk_v2 (const size_t k, const size_t &nb, const size_t &size_tensor, const  size_t &test_size, 
		       const std::vector<typename OPERATOR::result_type> &values,
		       const double &sqrtwdet,   const typename MATERIALLAW::result_type &uphys,  matrix_t &UGMat);
    
  };
  
  
  // ##########################################################
  // implementation of the template classes start here.
  // ##########################################################
  template <>
    double * getData<> (double  &in){
    return  &in;
  }; 
  
  template <>
    double * getData<> (xtensor::xVector<> &in){
    return in.data();
  };
  
  template <>
  double * getData<> (xtensor::xTensor2<> &in){
    return  in.data();
  };
  
  template <>
  double *getData<>(xtensor::xTensor4<> &in){
    return in.data();
  };
  
  template <>
    const double * getData<> (const double  &in){
    return  &in;
  }; 
  
  template <>
    const double * getData<> (const xtensor::xVector<> &in){
    return in.data();
  };
  
  template <>
  const double * getData<> (const xtensor::xTensor2<> &in){
    return  in.data();
  };
  
  template <>
  const double *getData<>(const xtensor::xTensor4<> &in){
    return in.data();
  };
  
  template<>
    size_t getSize< double >() {return 1;};
  
  template<>
    size_t getSize< xtensor::xVector<> >() {return 3; };
  
  template<>
  size_t getSize< xtensor::xTensor2<> >() {return 9; };
  
  template<>
  size_t getSize< xtensor::xTensor4<> >() {return 81; };
  
  //##################################################################################################################
  // implementation of the template classe xFormBilinearOptimizedBase  start here.
  //##################################################################################################################
  template < typename OPERATORTEST, typename OPERATORTRIAL >
    xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::xFormBilinearOptimizedBase(OPERATORTEST& o):  OperatorTest(o),  OperatorTrial(o), optestisoptrial(true) {}
  
  template < typename OPERATORTEST, typename OPERATORTRIAL >
    xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >:: xFormBilinearOptimizedBase(OPERATORTEST& otest,  OPERATORTRIAL &otrial):  OperatorTest(otest), OperatorTrial(otrial),  optestisoptrial(&otest== &otrial) {}
  
  template < typename OPERATORTEST, typename OPERATORTRIAL >
    void  xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    init(xGeomElem* _geo_appro_test, femFcts* _test){ 
    if (optestisoptrial){
      geo_appro_test = _geo_appro_test ;
      test = _test;
      const int n = test->size();
      Matrix.resize(n, n);
      std::fill(&Matrix(0,0), (&Matrix(0,0))+n*n, 0. );
    }
    else {
      init(_geo_appro_test, _test, _geo_appro_test, _test );
    }
  }
  
  template < typename OPERATORTEST, typename OPERATORTRIAL > 
    void xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    init(xGeomElem* _geo_appro_test, femFcts* _test, xGeomElem* _geo_appro_trial , femFcts* _trial)
    {
      geo_appro_test = _geo_appro_test ;
      geo_appro_trial = _geo_appro_trial ;
      test = _test;
      trial = _trial;
      const int m = test->size();
      const int n = trial->size();
      Matrix.resize(m, n);
      std::fill(&Matrix(0,0), (&Matrix(0,0))+m*n, 0. );
    }

  template < typename OPERATORTEST, typename OPERATORTRIAL > 
    const typename xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >
    ::matrix_t& xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    getMatrix() const {return Matrix; }
  
  template < typename OPERATORTEST, typename OPERATORTRIAL > 
    typename xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    matrix_t& xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    getMatrix() {return Matrix; } 
  
  template < typename OPERATORTEST, typename OPERATORTRIAL > 
    void  xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    accumulate_test_is_trial_no_law(xGeomElem*  geo_integ)
    {
      std::function < void (xGeomElem * , xGeomElem *  )> 
	set_appro_uvw = integ_is_appro(geo_integ, geo_appro_test)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
      const int nb = geo_integ->GetNbIntegrationPoints(); 
      const size_t test_size = test->size();
      const size_t size_tensor = getSize<typename OPERATORTEST::result_type>() ;
      matrix_t ValuesMat0(size_tensor*nb, test_size ) ;
      for(int k = 0; k < nb; k++){ 
	geo_integ->setUVW(k);
	set_appro_uvw(geo_integ, geo_appro_test);
	std::vector<typename OPERATORTEST::result_type> values;
	values.reserve(test_size);
	OperatorTest.eval(test, geo_appro_test, geo_integ, values);
	const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
	const double sqrtwdet = sqrt(wdet);
	for(size_t i=0;i<test_size;++i){
	  typename OPERATORTEST::result_type v = values[i]*sqrtwdet;
	  const double *datav = getData(v);
	  double *pValuesMat0 =  &ValuesMat0(k*size_tensor, i);
	  memcpy(pValuesMat0, datav, size_tensor*sizeof(double) );
	}
      } 
      syrk('L', true, 1.,  ValuesMat0, 1., Matrix); 
      return;
    }

  template < typename OPERATORTEST, typename OPERATORTRIAL > 
    void  xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    accumulate_test_is_not_trial_no_law(xGeomElem*  geo_integ){
    std::function < void (xGeomElem * , xGeomElem *  )> 
      set_appro_test_uvw = integ_is_appro(geo_integ, geo_appro_test)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
    std::function < void (xGeomElem * , xGeomElem *  )> 
      set_appro_trial_uvw = integ_is_appro(geo_integ, geo_appro_trial)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
    const int nb = geo_integ->GetNbIntegrationPoints(); 
    const size_t test_size  = test->size();
    const size_t trial_size = trial->size();
    const size_t size_tensor_test  = getSize<typename OPERATORTEST::result_type>();
    const size_t size_tensor_trial = getSize<typename OPERATORTRIAL::result_type>();
    matrix_t ValuesMatTest(size_tensor_test*nb, test_size ) ;
    matrix_t ValuesMatTrial(size_tensor_trial*nb, trial_size ) ;
    for(int k = 0; k < nb; k++){ 
      geo_integ->setUVW(k);
      set_appro_test_uvw(geo_integ, geo_appro_test);
      set_appro_trial_uvw(geo_integ, geo_appro_trial);
      std::vector<typename OPERATORTEST::result_type> values_test;
      values_test.reserve(test_size);
      std::vector<typename OPERATORTRIAL::result_type> values_trial;
      values_trial.reserve(trial_size);
      OperatorTest.eval(test,   geo_appro_test, geo_integ, values_test);
      OperatorTrial.eval(trial, geo_appro_trial, geo_integ, values_trial);
      const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
      for(size_t i=0;i<test_size;++i){
	const double *data_vtest = getData(values_test[i]);
	double *pValuesMatTest =  &ValuesMatTest(k*size_tensor_test, i);
	memcpy(pValuesMatTest, data_vtest, size_tensor_test*sizeof(double) );
      }
      for(size_t i=0;i<trial_size;++i){
	typename OPERATORTEST::result_type Triali = values_trial[i]*wdet;
	const double *data_vtrial = getData(Triali);
	double *pValuesMatTrial =  &ValuesMatTrial(k*size_tensor_trial, i);
	memcpy(pValuesMatTrial, data_vtrial, size_tensor_trial*sizeof(double) );
      } 
    }
    gemm(true, false, 1., ValuesMatTest, ValuesMatTrial, 1., Matrix );
    return;
  }


  template < typename OPERATORTEST, typename OPERATORTRIAL > 
    template < typename MATERIALLAW > 
    void  xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    accumulate_test_is_trial_law(xGeomElem*  geo_integ, const MATERIALLAW &Law)
    {
      std::function < void (xGeomElem * , xGeomElem *  )> 
	set_appro_uvw = integ_is_appro(geo_integ, geo_appro_test)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
      const int nb = geo_integ->GetNbIntegrationPoints(); 
      const size_t test_size = test->size();
      const size_t size_tensor = getSize<typename OPERATORTEST::result_type>() ;
      matrix_t ValuesMat0(size_tensor*nb, test_size ) ;
      matrix_t ValuesMat1(size_tensor*nb, test_size ) ;
      for(int k = 0; k < nb; k++){ 
	geo_integ->setUVW(k);
	set_appro_uvw(geo_integ, geo_appro_test);
	typename MATERIALLAW::result_type phys;
	Law(geo_appro_test, geo_integ, phys);
	std::vector<typename OPERATORTEST::result_type> values;
	values.reserve(test_size);
	OperatorTest.eval(test, geo_appro_test, geo_integ, values);
	const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
	for(size_t i=0;i<test_size;++i){
	  typename OPERATORTEST::result_type &v = values[i];
	  typename OPERATORTEST::result_type pv =  (phys*v)*wdet;
	  const double *datav = getData(v);
	  const double *datapv = getData(pv);
	  double *pValuesMat0 =  &ValuesMat0(k*size_tensor, i);
	  double *pValuesMat1 =  &ValuesMat1(k*size_tensor, i);
	  memcpy(pValuesMat0, datav, size_tensor*sizeof(double) );
	  memcpy(pValuesMat1, datapv, size_tensor*sizeof(double) );
	}
      }
      gemm(true, false, 1., ValuesMat0, ValuesMat1, 1., Matrix );
      return;
    }

  template < typename OPERATORTEST, typename OPERATORTRIAL > 
    template < typename MATERIALLAW > 
    void  xFormBilinearOptimizedBase <OPERATORTEST,  OPERATORTRIAL >::
    accumulate_test_is_not_trial_law(xGeomElem*  geo_integ, const MATERIALLAW &Law){
    std::function < void (xGeomElem * , xGeomElem *  )> 
      set_appro_test_uvw = integ_is_appro(geo_integ, geo_appro_test)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
    std::function < void (xGeomElem * , xGeomElem *  )> 
      set_appro_trial_uvw = integ_is_appro(geo_integ, geo_appro_trial)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
    const int nb = geo_integ->GetNbIntegrationPoints(); 
    const size_t test_size  = test->size();
    const size_t trial_size = trial->size();
    const size_t size_tensor_test  = getSize<typename OPERATORTEST::result_type>();
    const size_t size_tensor_trial = getSize<typename OPERATORTRIAL::result_type>();
    matrix_t ValuesMatTest(size_tensor_test*nb, test_size ) ;
    matrix_t ValuesMatTrial(size_tensor_trial*nb, trial_size ) ;
    for(int k = 0; k < nb; k++){ 
      geo_integ->setUVW(k);
      set_appro_test_uvw(geo_integ, geo_appro_test);
      set_appro_trial_uvw(geo_integ, geo_appro_trial);
      typename MATERIALLAW::result_type phys;
      Law(geo_appro_test, geo_integ, phys);
      std::vector<typename OPERATORTEST::result_type> values_test;
      values_test.reserve(test_size);
      std::vector<typename OPERATORTRIAL::result_type> values_trial;
      values_trial.reserve(trial_size);
      OperatorTest.eval(test,   geo_appro_test, geo_integ, values_test);
      OperatorTrial.eval(trial, geo_appro_trial, geo_integ, values_trial);
      const double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();
      for(size_t i=0;i<test_size;++i){
	const double *data_vtest = getData(values_test[i]);
	double *pValuesMatTest =  &ValuesMatTest(k*size_tensor_test, i);
	memcpy(pValuesMatTest, data_vtest, size_tensor_test*sizeof(double) );
      }
      for(size_t i=0;i<trial_size;++i){
	typename OPERATORTEST::result_type LawTriali = (phys*values_trial[i])*wdet;
	const double *data_vlawtrial = getData(LawTriali);
	double *pValuesMatTrial =  &ValuesMatTrial(k*size_tensor_trial, i);
	memcpy(pValuesMatTrial, data_vlawtrial, size_tensor_trial*sizeof(double) );
      } 
    }
    gemm(true, false, 1., ValuesMatTest, ValuesMatTrial, 1., Matrix );
    return;
  }
  
  
  //##################################################################################################################
  // implementation of the template classe xFormBilinearOptimizedWithoutLaw  start here.
  //##################################################################################################################
  template < typename OPERATORTEST,  typename OPERATORTRIAL >
    xFormBilinearOptimizedWithoutLaw<OPERATORTEST, OPERATORTRIAL >::
    xFormBilinearOptimizedWithoutLaw(OPERATORTEST& o)
    :xFormBilinearOptimizedBase< OPERATORTEST, OPERATORTRIAL >(o) {}
 
  template < typename OPERATORTEST,  typename OPERATORTRIAL >
    xFormBilinearOptimizedWithoutLaw<OPERATORTEST, OPERATORTRIAL >::
    xFormBilinearOptimizedWithoutLaw(OPERATORTEST& otest, OPERATORTRIAL &otrial)
    :xFormBilinearOptimizedBase< OPERATORTEST, OPERATORTRIAL > (otest, otrial){}

  template < typename OPERATORTEST,  typename OPERATORTRIAL >
    void xFormBilinearOptimizedWithoutLaw<OPERATORTEST, OPERATORTRIAL >::
    accumulate(xGeomElem*  geo_integ)
    {
      if (this->optestisoptrial) this->accumulate_test_is_trial_no_law(geo_integ);
      else this->accumulate_test_is_not_trial_nolaw(geo_integ);
    }
 
  template < typename OPERATORTEST,  typename OPERATORTRIAL >
    void xFormBilinearOptimizedWithoutLaw<OPERATORTEST, OPERATORTRIAL >::
    finalize()
    {
      if (this->optestisoptrial) matrixcopylowtriangletouptriangle(this->Matrix);
    }
  
  //##################################################################################################################
  // implementation of the template classe xFormBilinearOptimizedWithLaw  start here.
  //##################################################################################################################
  template < typename OPERATORTEST, typename MATERIALLAW, typename OPERATORTRIAL >
    xFormBilinearOptimizedWithLaw< OPERATORTEST, MATERIALLAW, OPERATORTRIAL > ::
    xFormBilinearOptimizedWithLaw(OPERATORTEST& otest,  MATERIALLAW& l, OPERATORTRIAL &otrial)
    :xFormBilinearOptimizedBase< OPERATORTEST, OPERATORTRIAL > (otest, otrial), Law(l){}
  
  template < typename OPERATORTEST, typename MATERIALLAW, typename OPERATORTRIAL >
    xFormBilinearOptimizedWithLaw< OPERATORTEST, MATERIALLAW, OPERATORTRIAL > ::
    xFormBilinearOptimizedWithLaw(OPERATORTEST& otest,  MATERIALLAW& l)
    :xFormBilinearOptimizedBase< OPERATORTEST, OPERATORTRIAL > (otest), Law(l) {}
  
  template < typename OPERATORTEST, typename MATERIALLAW, typename OPERATORTRIAL >
    void xFormBilinearOptimizedWithLaw< OPERATORTEST, MATERIALLAW, OPERATORTRIAL > ::
    accumulate(xGeomElem*  geo_integ)
    { 
      if (this->optestisoptrial) this->accumulate_test_is_trial_law(geo_integ, Law);
      else this->accumulate_test_is_not_trial_law(geo_integ, Law);
    }

  template < typename OPERATORTEST, typename MATERIALLAW, typename OPERATORTRIAL >
    void xFormBilinearOptimizedWithLaw< OPERATORTEST, MATERIALLAW, OPERATORTRIAL > ::
    finalize(){ }
  
  //##################################################################################################################
  // implementation of the template classe xFormBilinearOptimizedSymmetricWithLawFactored  start here.
  //##################################################################################################################
  template < typename  OPERATOR, typename MATERIALLAW >
    xFormBilinearOptimizedSymmetricWithLawFactored< OPERATOR, MATERIALLAW >
    ::xFormBilinearOptimizedSymmetricWithLawFactored(OPERATOR& o, MATERIALLAW& l):  Operator(o), Law(l){}
  
  template < typename  OPERATOR, typename MATERIALLAW >
    void xFormBilinearOptimizedSymmetricWithLawFactored< OPERATOR, MATERIALLAW >::
    init(xGeomElem* _geo_appro_test, femFcts* _test, xGeomElem* _geo_appro_trial, femFcts* _trial)
    {
      std::cout  << "The Bilinear Form " << " xFormBilinearOptimizedSymmetricWithLawFactored " << " expect same space on right and left "<< __FILE__ << __LINE__ << std::endl;
      throw;
    }

  template < typename  OPERATOR, typename MATERIALLAW >
    void xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    init(xGeomElem* _geo_appro_test, femFcts* _test)
    {
      geo_appro_test = _geo_appro_test ;
      test = _test;
      const int n = test->size();
      Matrix.resize(n, n);
      std::fill(&Matrix(0,0), (&Matrix(0,0))+n*n, 0. );
    }
  
  template < typename  OPERATOR, typename MATERIALLAW >  
    void xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    accumulate(xGeomElem*  geo_integ)
    {
      std::function < void (xGeomElem * , xGeomElem *  )> 
	set_appro_uvw = integ_is_appro(geo_integ, geo_appro_test)? set_appro_uvw_integ_is_same:set_appro_uvw_integ_is_not_same;
      const size_t nb = geo_integ->GetNbIntegrationPoints(); 
      const size_t test_size = test->size();
      const size_t size_tensor = getSize<typename OPERATOR::result_type>() ;
      matrix_t UGMat(size_tensor*nb, test_size ) ;
      std::vector<typename OPERATOR::result_type> values;
      values.reserve(test_size);
      
      typename MATERIALLAW::result_type uphys;
      for(size_t k = 0; k < nb; k++){ 
	geo_integ->setUVW(k);
	set_appro_uvw(geo_integ, geo_appro_test);
	Law(geo_appro_test, geo_integ, uphys);
	Operator.eval(test, geo_appro_test, geo_integ, values);
	const double sqrtwdet = sqrt(geo_integ->GetWeight() * geo_integ->GetDetJac());
	setUGMatk_v1 (  k, nb, size_tensor, test_size, values, sqrtwdet, uphys, UGMat);
	values.clear();
      }
      syrk('L', true, 1.,  UGMat, 1., Matrix);
      return;
    }
 
  template < typename  OPERATOR, typename MATERIALLAW >  
    void xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    finalize()
    {
      matrixcopylowtriangletouptriangle(Matrix);
    }
  
  template < typename  OPERATOR, typename MATERIALLAW >  
    const typename xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    matrix_t& xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    getMatrix() const {return Matrix; }

  template < typename  OPERATOR, typename MATERIALLAW >  
    typename xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    matrix_t& xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    getMatrix() {return Matrix; } 
  
  template < typename  OPERATOR, typename MATERIALLAW >  
    void xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    setUGMatk_v0 (const size_t k, const size_t &nb, const size_t &size_tensor, const  size_t &test_size, 
		  const std::vector<typename OPERATOR::result_type> &values,
		  const double &sqrtwdet,   const typename MATERIALLAW::result_type &uphys,  matrix_t &UGMat)
    {
      typename OPERATOR::result_type        UGki;
      const double  *pUGki =  getData(UGki);
      for(size_t i=0;i<test_size;++i){
	UGki = (uphys*values[i])*sqrtwdet;;
	double *pUGMatki =  &UGMat(k*size_tensor, i);
	memcpy(pUGMatki, pUGki, size_tensor*sizeof(double) );
      }
    }
   
  template < typename  OPERATOR, typename MATERIALLAW >  
    void xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>
    ::setUGMatk_v1 (const size_t k, const size_t &nb, const size_t &size_tensor, const  size_t &test_size, 
		    const std::vector<typename OPERATOR::result_type> &values,
		    const double &sqrtwdet,   const typename MATERIALLAW::result_type &uphys,  matrix_t &UGMat)
    {
      matrix_t GMat(size_tensor, test_size ) ;
      for(size_t i=0;i<test_size;++i){
	const double * pGki = getData(values[i]);
	double *pGMatki =  &GMat(0, i);
	memcpy(pGMatki, pGki, size_tensor*sizeof(double) );	
      }
    
      const char TRANSA='T';
      const char TRANSB='N';
      const int M = size_tensor;
      const int N = test_size;
      const int K = size_tensor;
      const double alpha = sqrtwdet;
      const double * A = getData(uphys);
      const int LDA = size_tensor;
      const double * B = &GMat(0, 0);
      const int LDB = size_tensor;
      const double beta = 0.;
      double * C = &UGMat(k*size_tensor, 0); 
      const int LDC = size_tensor*nb;
      xlinalg::xCPPBlasDef<double>::gemm(TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    }
  
  template < typename  OPERATOR, typename MATERIALLAW >  
    void xFormBilinearOptimizedSymmetricWithLawFactored<OPERATOR, MATERIALLAW>::
    setUGMatk_v2 (const size_t k, const size_t &nb, const size_t &size_tensor, const  size_t &test_size, 
		  const std::vector<typename OPERATOR::result_type> &values,
		  const double &sqrtwdet,   const typename MATERIALLAW::result_type &uphys,  matrix_t &UGMat)
    {
      for(size_t i=0;i<test_size;++i){
	const double *pUGki = getData(values[i]);
	double *pUGMatki =  &UGMat(k*size_tensor, i);
	memcpy(pUGMatki, pUGki, size_tensor*sizeof(double) );
      }
      
      const char SIDE = 'L';
      const char UPLO = 'L';
      const char TRANSA='T';
      const char DIAG  ='N';
      const int M = size_tensor;
      const int N = test_size;
      const double alpha = sqrtwdet;
      const double * A = getData(uphys);
      const int LDA = size_tensor;
      double * B = &UGMat(k*size_tensor, 0);
      const int LDB = size_tensor*nb;
      xlinalg::xCPPBlasDef<double>::trmm(SIDE, UPLO, TRANSA, DIAG, M, N, alpha,  A, LDA, B, LDB);
    }
 
}// end of namespace xfem.
#endif
