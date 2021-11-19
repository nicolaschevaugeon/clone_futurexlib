#include <unordered_map>

#include "xApproxFunction.h"
#include "xAttachedDataManagerAOMD.h"
#include "xSpaceMultiXFEM.h"
#include "xSupportComponent.h"

namespace xfem
{
class xApproxFunctionBySupportComponentGenerator : public xGeneratorApproxFunctionEnrichedXFEMWith
{
  public:
   xApproxFunctionBySupportComponentGenerator(
       std::function<const xfem::xSupportComponent*(AOMD::mEntity* e)> _getSupportComponent);
   void generateKeysAndFcts(AOMD::mEntity* e, xfem::xValKey& key, xfem::ValKeyModifierPtr key_modifier, approxFctPtr_t f,
                            femKeys* keys, femFcts_t* fcts) override;
   void generateKeys(xfem::xValKey& key, xfem::ValKeyModifierPtr key_modifier, femKeys* keys) override;
   void clear();

  protected:
   virtual void getEnrichmentFunctions(AOMD::mEntity* enrichedentity, femFcts_t& functions) = 0;
   xPartition getPartition(const AOMD::mEntity& e) const;
   std::function<const xfem::xSupportComponent*(AOMD::mEntity* e)> getSupportComponent;
   xinterface::aomd::xAttachedDataManagerAOMD<femFcts_t> generated_functions;
};

class xApproxFunctionConstantBySupportComponentGenerator : public xApproxFunctionBySupportComponentGenerator
{
  public:
   using xApproxFunctionBySupportComponentGenerator::xApproxFunctionBySupportComponentGenerator;

  private:
   void getEnrichmentFunctions(AOMD::mEntity* enrichedentity, femFcts_t& functions) override;
};

class xApproxFunctionConstantAndRampBySupportComponentGenerator : public xApproxFunctionBySupportComponentGenerator
{
  public:
   using xApproxFunctionBySupportComponentGenerator::xApproxFunctionBySupportComponentGenerator;

  private:
   void getEnrichmentFunctions(AOMD::mEntity* enrichedentity, femFcts_t& functions) override;
};

class xHeavisideBySub : public xfem::xApproxFunction
{
  public:
   xHeavisideBySub(std::unordered_set<const AOMD::mEntity*> _data);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   const std::unordered_set<const AOMD::mEntity*> data;
};

class xLinearBySub : public xfem::xApproxFunction
{
  public:
   xLinearBySub(std::unordered_map<const AOMD::mVertex*, double> _data);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   std::unordered_map<const AOMD::mVertex*, double> data;
   std::vector<double> vals(const AOMD::mEntity& e_integ) const;
};
}  // namespace xfem
