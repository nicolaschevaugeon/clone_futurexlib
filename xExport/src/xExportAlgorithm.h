#ifndef XEXPORTALGORITHM_H
#define XEXPORTALGORITHM_H
#include "xApplyCommandOnIntegrationRule.h"
#include "xExport.h"
#include "xLevelSet.h"
#include "xVectorLevelSet.h"
//#include "xAlgorithm.h"

namespace xexport
{
template <typename T, class RANGE>
void Export(const xfem::xEval<T>& eval, xexport::xExport& pexport, const std::string& fieldName,
            const xfem::xIntegrationRule& integration_rule, RANGE range, const bool simplex = true,
            xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   Export(eval, pexport, fieldName, integration_rule, range.begin(), range.end(), simplex, integ2appro);
}

template <typename T, class ITER>
void Export(const xfem::xEval<T>& eval, xexport::xExport& pexport, const std::string& fieldName,
            const xfem::xIntegrationRule& integration_rule, ITER it, ITER end, const bool simplex = true,
            xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   if (!pexport.processStarted())
   {
      std::string fs = fieldName;
      pexport.openFile(fs.c_str());
      pexport.startView(fieldName.c_str());
      xexport::xPlotCommand<T> plot(eval, pexport, simplex);
      xfem::ApplyCommandOnIntegrationRule(plot, integration_rule, it, end, integ2appro);
      pexport.endView();
      pexport.closeFile();
   }
   else
   {
      pexport.startView(fieldName.c_str());
      xexport::xPlotCommand<T> plot(eval, pexport, simplex);
      xfem::ApplyCommandOnIntegrationRule(plot, integration_rule, it, end, integ2appro);
      pexport.endView();
   }
   return;
}

template <typename T, class RANGE>
void Export(const xfem::xEval<T>& eval, xexport::xExport& pexport, const std::string& fieldName,
            const xfem::xIntegrationRule& integration_rule, RANGE range, xfem::xEntityToEntity integ2appro)
{
   Export(eval, pexport, fieldName, integration_rule, range.begin(), range.end(), integ2appro);
}

template <typename T, class ITER>
void Export(const xfem::xEval<T>& eval, xexport::xExport& pexport, const std::string& fieldName,
            const xfem::xIntegrationRule& integration_rule, ITER it, ITER end, xfem::xEntityToEntity integ2appro)
{
   Export(eval, pexport, fieldName, integration_rule, it, end, true, integ2appro);
}

template <typename ITER>
void Export(const xfem::xVectorLevelSet& vls, xexport::xExport& pexport, ITER it, ITER end, const std::string& fieldName,
            const bool simplex = true)
{
   xfem::xEvalVectorLevelSetVector<xtool::xIdentity<xtensor::xVector<>>> eval_vls(vls);
   Export(eval_vls, pexport, fieldName + "_vector", xfem::xIntegrationRuleBasic(0), it, end, simplex);
   xfem::xEvalVectorLevelSetInOut<xtool::xIdentity<double>> eval_vlsInOut(vls);
   Export(eval_vlsInOut, pexport, fieldName + "_inout", xfem::xIntegrationRuleBasic(0), it, end, simplex);
}

void Export(const xfem::xLevelSet& ls, xexport::xExport& pexport, const std::string& fieldName, const bool simplex = true);

void Export(const xfem::xMesh*, xexport::xExport& pexport, const std::string& name);
void Export(const xfem::xMesh&, xexport::xExport& pexport, const std::string& name);

/// To export fields defined at integration points
template <class T, class ITER>
void Export(const xfem::xEvalField<T, xfem::xFieldPointwise>& eval, xexport::xExport& pexport, const std::string& fieldName,
            const xfem::xIntegrationRule& integration_rule, ITER it, ITER end, const bool simplex = true,
            xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   if (!pexport.processStarted())
   {
      std::string fs = fieldName;
      pexport.openFile(fs.c_str());
      pexport.startView(fieldName.c_str());
      // xexport::xPlotPWCommand<T> plot(field_PW, pexport, integration_rule);
      xexport::xPlotPWCommand<T> plot(eval, pexport, integration_rule);
      ApplyCommandOnIntegrationRule(plot, integration_rule, it, end, integ2appro);
      pexport.endView();
      pexport.closeFile();
   }
   else
   {
      pexport.startView(fieldName.c_str());
      // xexport::xPlotPWCommand<T> plot(field_PW, pexport,integration_rule);
      xexport::xPlotPWCommand<T> plot(eval, pexport, integration_rule);
      ApplyCommandOnIntegrationRule(plot, integration_rule, it, end, integ2appro);
      pexport.endView();
   }
   return;
}

}  // namespace xexport
#endif  // XEXPORTALGORITHM_H
