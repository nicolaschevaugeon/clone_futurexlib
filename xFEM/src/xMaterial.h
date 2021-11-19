/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __MATERIAL_H
#define __MATERIAL_H
#include <cassert>
#include <iostream>
#include <map>
#include <string>

#include "AOMDfwd.h"
#include "xEval.h"
#include "xTensor2.h"
#include "xTensor3.h"
#include "xTensor4.h"
#include "xTensors.h"
#include "xTensorsPtr.h"

namespace xfem
{
class xMaterial
{
  public:
   virtual ~xMaterial() = default;

   virtual void checkProperties() = 0;

   virtual void sensitivityTo(const std::string& phys_token, xtensor::xTensor2<>& sensitivity)
   {
      ErrorMessageInFunction("sensitivityToVector");
      return;
   }

   virtual void sensitivityTo(const std::string& phys_token, xtensor::xTensor3<>& sensitivity)
   {
      ErrorMessageInFunction("sensitivityToVector2");
      return;
   }

   virtual void sensitivityTo(const std::string& phys_token, xtensor::xTensor4<>& sensitivity)
   {
      ErrorMessageInFunction("sensitivityToMatrix");
      return;
   }

   virtual void sensitivityTo(const std::string& phys_token, xtensor::xTensor4Isotropic& sensitivity)
   {
      ErrorMessageInFunction("sensitivityToMatrix");
      return;
   }

   virtual void sensitivityTo(const std::string& phys_token, xtensor::xTensor4AnisoPlaneStrain& sensitivity)
   {
      ErrorMessageInFunction("sensitivityToMatrix");
      return;
   }

   virtual void sensitivityTo(const std::string& phys_token, xtensor::xTensor4AnisoPlaneStress& sensitivity)
   {
      ErrorMessageInFunction("sensitivityToMatrix");
      return;
   }

   virtual void sensitivityTo(const std::string& phys_token, double& val)
   {
      ErrorMessageInFunction("sensitivityToScalar");
      return;
   }

   // Added SLC : for non linear purpose
   virtual void sensitivityTo(const std::string& phys_token, xtensor::xVector<>& sensitivity)
   {
      ErrorMessageInFunction("sensitivityToScalar2");
      return;
   }

   virtual void computeCurrentState()
   {
      ErrorMessageInFunction("ComputeCurrentState");
      return;
   }

   virtual void computeCurrentState(std::string phys_to_update)
   {
      (*this).computeCurrentState();
      return;
   }

   const xTensors* getProperties() const { return &properties; }
   xTensors* getProperties() { return &properties; }

   void readProperties(const std::string& filename)
   {
      properties.read(filename);
      fname = filename;
   }
   xTensorsSignature* getVariablesSignature() { return &variables_signature; }
   xTensorsSignature* getValuesSignature() { return &values_signature; }

   const xTensors& getVariables() const { return variables; }
   const xTensors& getValues()
   {
      computeCurrentState();
      return values;
   }

   void setVariables(const xTensors& _variables)
   {
      variables = _variables;
      uptodate = false;
   };

   void setOldVariables(tensorsPtr_t o) { old = o; }
   void setCurrentVariables(tensorsPtr_t c) { curr = c; }
   // this is a very ugly hack, just for now
   virtual void setTemperature(double T)
   {
      std::cout << "not implemented" << std::endl;
      throw;
   }

  protected:
   // variables
   xTensorsSignature variables_signature;
   xTensorsSignature values_signature;
   xTensorsSignature properties_signature;
   tensorsPtr_t curr;
   tensorsPtr_t old;
   // properties

   xTensors properties;
   xTensors variables;
   xTensors values;

   bool uptodate;
   std::string fname;

  private:
   void ErrorMessageInFunction(const char* fct) const
   {
      std::cerr << fct << " not coded " << std::endl;
      assert(0);
      return;
   }
};

class xConductive : virtual public xMaterial
{
  public:
   xConductive();
   // Definition of the pure virtual functions
   void checkProperties() override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor2<>& sensitivity) override;
   inline static double delta(int i, int j) { return ((i == j) ? 1.0 : 0.0); }
   static void SetThermicConductivityIsotropic(double k, xtensor::xTensor2<>& conductivity);

  protected:
   // Useful tensors derived from the above
   // set up in the constructor
   xtensor::xTensor2<> thermic_conductivity;
};

class xNewtonianFluid : virtual public xMaterial
{
  public:
   xNewtonianFluid();
   void checkProperties() override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4<>& sensitivity) override;

  private:
   double delta(int i, int j) { return ((i == j) ? 1.0 : 0.0); }

  protected:
   xtensor::xTensor4<> viscosity_tensor;
};

class xEvalElasticVariables : public xEval<xTensors>
{
  public:
   xEvalElasticVariables(const xEval<xtensor::xTensor2<>>& _eval_strain);
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xTensors& tensors) const override;

  private:
   const xEval<xtensor::xTensor2<>>& eval_strain;
};

class xElastic : virtual public xMaterial
{
  public:
   xElastic();
   void checkProperties() override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4<>& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4Isotropic& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4AnisoPlaneStrain& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4AnisoPlaneStress& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, double& sensitivity) override;
   void computeCurrentState() override;
   /////
   static inline double delta(int i, int j) { return ((i == j) ? 1.0 : 0.0); }
   static void SetElasticStiffnessIsotropic(double E, double nu, xtensor::xTensor4<>& tensor);
   static void SetElasticStiffnessIsotropicPlate(double E, double nu, xtensor::xTensor4<>& tensor);
   static void SetElasticComplianceIsotropicPlate(double E, double nu, xtensor::xTensor4<>& tensor);
   static void SetElasticStiffnessIsotropicFromLame(double lam, double mu, xtensor::xTensor4<>& tensor);
   static void SetElasticComplianceIsotropic(double E, double nu, xtensor::xTensor4<>& tensor);
   static void SetElasticStiffnessTransverselyIsotropic(double E1, double E3, double G3, double nu, double nu3,
                                                        xtensor::xTensor4<>& tensor);
   static void SetElasticComplianceTransverselyIsotropic(double E1, double E3, double G3, double nu, double nu3,
                                                         xtensor::xTensor4<>& tensor);
   static void SetElasticStiffnessOrthotropic(double E1, double E2, double E3, double nu12, double nu23, double nu31, double mu12,
                                              double mu23, double mu31, xtensor::xTensor4<>& tensor);
   static void SetShearStiffnessIsotropic(double E, double nu, xtensor::xTensor4<>& tensor);
   static void SetInvLameOneIsotropic(double E, double nu, double& val);
   static void SetPlateBendingStiffnessIsotropic(double E, double nu, double t, xtensor::xTensor4<>& tensor);
   static void SetPlateShearStiffnessIsotropic(double E, double nu, double t, xtensor::xTensor2<>& tensor);

  protected:
   double lam, mu;
   double E, nu;
   xtensor::xTensor4<> elastic_stiffness;
   xtensor::xTensor4<> elastic_compliance;
};

class xEvalThermoElasticVariables : public xEval<xTensors>
{
  public:
   xEvalThermoElasticVariables(const xEval<double>& _eval_temperature, const xEval<xtensor::xTensor2<>>& _eval_strain);
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xTensors& tensors) const override;

  private:
   const xEval<double>& eval_temperature;
   const xEval<xtensor::xTensor2<>>& eval_strain;
};

class xThermoElastic : public xMaterial
{
  public:
   xThermoElastic();
   void checkProperties() override;
   void computeCurrentState() override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4Isotropic& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4<>& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, double& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor2<>& sensitivity) override;
   inline double delta(int i, int j) const { return ((i == j) ? 1.0 : 0.0); };

  protected:
   double E, nu, lam, mu, k, alpha, T0;
   xtensor::xTensor4<> elastic_stiffness;
   xtensor::xTensor2<> thermic_conductivity;
   void setThermicConductivityIsotropic(double k);
   void setElasticStiffnessIsotropicFromLame(xtensor::xTensor4<>& stiffness);
   void SetElasticComplianceIsotropic(double E, double nu, xtensor::xTensor4<>& compliance) const;
};

class xThermoElasticTemperatureDependant : public xMaterial
{
  public:
   xThermoElasticTemperatureDependant();
   void checkProperties() override;
   void setTemperature(double T) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4Isotropic& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor4<>& sensitivity) override;

   void sensitivityTo(const std::string& phys_token, double& sensitivity) override;
   void sensitivityTo(const std::string& phys_token, xtensor::xTensor2<>& sensitivity) override;
   static inline double delta(int i, int j) { return ((i == j) ? 1.0 : 0.0); }

  protected:
   double E, nu, lam, mu, k, alpha;
   xtensor::xTensor4<> elastic_stiffness;
   xtensor::xTensor2<> thermic_conductivity;
   void setThermicConductivityIsotropic(double k);
   void setElasticStiffnessIsotropicFromLame(xtensor::xTensor4<>& stiffness);

  private:
   xPieceWiseLinear Et, nut, alphat;
   double temperature;
};

}  // namespace xfem

#endif
