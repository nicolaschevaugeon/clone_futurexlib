/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __MATERIAL_SENSITIVITY_H
#define __MATERIAL_SENSITIVITY_H
#include <iostream>
#include <cassert>
#include <string>
#include "AOMDfwd.h"
#include "xTensor2.h"
#include "xTensors.h"
#include "xTensor4.h"
#include "xMaterial.h"
#include "xMaterialManager.h"
#include "xEval.h"


namespace xfem
{


class  xEvalMaterialValueTensor2 : public xEval<xtensor::xTensor2<> > {
  public :
    typedef xEval<xtensor::xTensor2<> >::result_type result_type;
    
  xEvalMaterialValueTensor2(const std::string& _value_name, const xEval<xTensors> & _tensors_variables_eval):value_name(_value_name), tensors_variables_eval(_tensors_variables_eval)
    {
    };
    
    void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, result_type& value) const override{ 
      xMaterial &mat = *xMaterialManagerSingleton::instance().getMaterial(geo_integ); 
      xTensors variables;
      variables.setSignature (mat.getVariablesSignature()) ;
      tensors_variables_eval(geo_appro, geo_integ, variables);
      mat.setVariables(variables );
      value= mat.getValues().tensor2(value_name);
    }
  private:
    std::string value_name;
    const xEval<xTensors> & tensors_variables_eval ;
  
    
  };

class  xEvalMaterialValues : public xEval< xTensors  > {
  public :
    typedef xEval<xTensors>::result_type result_type;
    
  xEvalMaterialValues( const xEval<xTensors> & _tensors_variables_eval): tensors_variables_eval(_tensors_variables_eval)
    {
    };
    
    void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, result_type& value) const override{ 
      xMaterial &mat = *xMaterialManagerSingleton::instance().getMaterial(geo_integ); 
      xTensors variables;
      variables.setSignature (mat.getVariablesSignature()) ;
      tensors_variables_eval(geo_appro, geo_integ, variables);
      mat.setVariables(variables ); 
      value.setSignature (mat.getValuesSignature()) ;
      value = mat.getValues();
    }
  private:
    const xEval<xTensors> & tensors_variables_eval ;
      
  };

template <typename T>
class  xUniformMaterialSensitivity : public xEval<T> {

public:
typedef typename xEval<T>::result_type result_type;
~xUniformMaterialSensitivity() override = default;
xUniformMaterialSensitivity(const std::string& variable_)
  : mat_manager(xMaterialManagerSingleton::instance()), variable(variable_)   {}
void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, result_type& curr_law) const override
{
    const bool debug = xdebug_flag;
    if(debug)
      {
        std::cout<<"Appro=\n";
	geo_appro->getEntity()->print();
        std::cout<<"Integ=\n";
	geo_integ->getEntity()->print();
      }
    if (debug) std::cout << "in setcurrentelem zone_id is " << geo_integ->getZone() << std::endl;
    xMaterial * mat = mat_manager.getMaterial(geo_integ);
    if (debug) std::cout << "after mat_manager.getMaterial" << std::endl;
    if (debug) assert(mat != nullptr);
    mat->sensitivityTo(variable, curr_law);
}  
private:
 xMaterialManager& mat_manager;
 std::string variable;
};

template <typename T>
class  xUniformMaterialParameter : public xEval<T> {

public:
typedef typename xEval<T>::result_type result_type;
~xUniformMaterialParameter() override = default;
xUniformMaterialParameter(const std::string& parametername_)
  //: mat_manager(xMaterialManagerSingleton::instance()), 
  :parametername(parametername_)   {}
void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, result_type& curr_para) const override
{
    xMaterial * mat = xMaterialManagerSingleton::instance().getMaterial(geo_integ);
    curr_para =  mat->getProperties()->scalar(parametername);
}  
private:
//xMaterialManager& mat_manager;
 std::string parametername;
};

  template <typename T >
class  xEvalDependantMaterialSensitivity : public xEval<T> {

public:
typedef typename xEval<T>::result_type result_type;
virtual ~xEvalDependantMaterialSensitivity() = default;
  xEvalDependantMaterialSensitivity(const std::string& variable_, const xEval< double >& _temperature)
    : mat_manager(xMaterialManagerSingleton::instance()), variable(variable_), temperature(_temperature)   {}
void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, result_type& curr_law) const
{
    const bool debug = xdebug_flag;
    xMaterial * mat = mat_manager.getMaterial(geo_integ);
    if (debug) std::cout << "after mat_manager.getMaterial" << std::endl;
    if (debug) assert(mat != nullptr);
    double Temp;
    temperature( geo_appro, geo_integ, Temp);
    mat->setTemperature(Temp);
    mat->sensitivityTo(variable, curr_law);
}  
private:
 xMaterialManager& mat_manager;
 std::string variable;
 const xEval< double > &temperature; 
 
};




} // end of namespace

#endif





