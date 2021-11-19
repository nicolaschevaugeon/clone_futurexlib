/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xcInteractionIntegrals_
#define _xcInteractionIntegrals_
#include <string>

#include "mTensor2.h"
#include "xEnv.h"
#include "xEval.h"
#include "xField.h"
#include "xForm.h"
#include "xLevelSet.h"
#include "xLevelSetOperators.h"
#include "xOperators.h"
#include "xParseData.h"
#include "xRegion.h"
#include "xTensor2.h"
#include "xValueManager.h"
#include "xVector.h"
#include "xcFront.h"

namespace xfem
{
class xIntegrationRule;
class xSpaceRegular;
}  // namespace xfem
using namespace xfem;

namespace xcrack
{
template <template <class> class DATAMANAGER>
class xAcceptLocalPart
{
  public:
   xAcceptLocalPart(const DATAMANAGER<AOMD::mEntity *> &_was_created_by) : was_created_by(_was_created_by) {}
   bool operator()(AOMD::mEntity *e)
   {
      if (was_created_by.getData(*e)) return true;
      return false;
   }

  private:
   const DATAMANAGER<AOMD::mEntity *> &was_created_by;
};

class xcCrackBase;
class xcEvalJintTerms;

class xcRhsForInteractionsIntegrals
{
  public:
   xcRhsForInteractionsIntegrals(const std::string &_name, const xEval<xtensor::xTensor2<>> &_P,
                                 const xEval<xtensor::xVector<>> &_s)
       : name(_name), P(_P), s(_s)
   {
   }
   const std::string name;
   const xEval<xtensor::xTensor2<>> &P;
   const xEval<xtensor::xVector<>> &s;
};

/// This class take a lCrack at construction, separate the crack in simply connected parts and compute sifs along the crack
/*!
  parameters should contain :
  int sifs_nb_modes
  double sifs_rho_geo
  int sifs_nb_layers_core
  int sifs_nb_layers_cylinder
  int sifs_verbosity
  !*/

class xcInteractionIntegralsOnCrack
{
  public:
   struct resultsfrontpartpoint
   {
      double s, x, y, z;
      std::map<std::string, double> vals;
   };
   xcInteractionIntegralsOnCrack(const xcCrackBase &crack, const xParseData &parameters,
                                 const std::vector<xcRhsForInteractionsIntegrals *> &rhs);
   ~xcInteractionIntegralsOnCrack();
   void computeSIFS(const xField<> &disp_l, const xIntegrationRule &integrator_vol, const xIntegrationRule &integrator_bnd);
   void postResults(std::map<xcFrontPart *, std::vector<resultsfrontpartpoint>> &results) const;
   void setVerbosity(bool v) { verbose = v; };
   void exportSIFS();

  private:
   // void computeSIFSOnFrontPart(const xField<> &disp_l, const   xIntegrationRule& integrator_vol, const
   //				xIntegrationRule& integrator_bnd, const xcFrontPart& front_part) const;
   void exportSIFSOnFrontPart(const xcFrontPart &front_part) const;

   template <class ASSEMBLER>
   void assembleRhsforFrontPart(const xField<> &J_field, const xIntegrationRule &integrator_vol,
                                const xIntegrationRule &integrator_bnd, const xcFrontPart &front_part, ASSEMBLER &assembler_rhs,
                                const xEval<xtensor::xTensor2<>> &P, const xEval<xtensor::xVector<>> &s) const;
   const xcCrackBase &crack;
   xcFront front;
   bool verbose;
   const xParseData &parameters;
   double young, poisson, Estar;
   int dim_mesh;
   xValueManagerDist<double> double_manager;
   std::map<string, xField<> *> fields_I;
   const std::vector<xcRhsForInteractionsIntegrals *> &rhss;
};

// class xcSetRadialFunction : public xLevelSetModifier {

// public:
//   xcSetRadialFunction() {}
//   void visit(xLevelSet& f, xRegion target);
// private:
// };

// Rhs = int_omega (div P :  q) ++ int_domega  ( q P n)
template <class ASSEMBLER>
void xcInteractionIntegralsOnCrack::assembleRhsforFrontPart(const xField<> &J_field, const xIntegrationRule &integrator_vol,
                                                            const xIntegrationRule &integrator_bnd, const xcFrontPart &front_part,
                                                            ASSEMBLER &assembler_rhs, const xEval<xtensor::xTensor2<>> &P,
                                                            const xEval<xtensor::xVector<>> &s) const
{
   const xEvalUnary<xtool::xChangeSign<xtensor::xTensor2<>>> mP = (P);
   xFormLinearWithLoad<xGradOperator<xtool::xIdentity<xtensor::xTensor2<>>>, xEvalUnary<xtool::xChangeSign<xtensor::xTensor2<>>>>
       mPgradq(mP);

   xEvalNormal eval_normal;

   xEvalBinary<xtool::xMult<xtensor::xTensor2<>, xtensor::xVector<>, xtensor::xVector<>>> Pn(P, eval_normal);
   xFormLinearWithLoad<xValOperator<xtool::xIdentity<xtensor::xVector<>>>, xEval<xtensor::xVector<>>> qPn(Pn);
   xRegion region_for_int(&front_part.mesh, front_part.subset_for_int_label);

   Assemble(mPgradq, assembler_rhs, integrator_vol, J_field, region_for_int.begin(), region_for_int.end());
   int dim_region = region_for_int.dim();
   xFilteredRegion<xIter, xAcceptOnClassifiedBoundary> domain_for_integral_ext_bnd(
       region_for_int.begin(dim_region - 1), region_for_int.end(dim_region - 1), xAcceptOnClassifiedBoundary());
   Assemble(qPn, assembler_rhs, integrator_bnd, J_field, domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end(),
            xUpperAdjacency());

   if (!(dynamic_cast<const xEvalZero<xtensor::xVector<>> *>(&s)))
   {
      xFormLinearWithLoad<xValOperator<xtool::xIdentity<xtensor::xVector<>>>, xEval<xtensor::xVector<>>> sq(s);
      Assemble(sq, assembler_rhs, integrator_vol, J_field, region_for_int.begin(), region_for_int.end());
   }

   /*
     xexport::xExportGmshAscii pexport;
     xIntegrationRuleBasic iexp;
     Export(eval_normal, pexport, "normal", iexp, domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end());
     Export(eval_normal, pexport, "normal2", iexp,domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end(),
     true, xUpperAdjacency());
     Export(Pn, pexport, "boundaryterm", iexp,domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end(),
     true, xUpperAdjacency());
   */
}

}  // namespace xcrack
#endif
