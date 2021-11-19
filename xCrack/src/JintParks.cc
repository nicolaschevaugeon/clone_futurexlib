/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#include "JintParks.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "xcValueJs.h"

// xfem include
#include "lCrack.h"
#include "xAlgorithm.h"
#include "xApproxFunction.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"
#include "xCommandOnGeomElem.h"
#include "xElement.h"
#include "xEnv.h"
#include "xExportGmsh.h"
#include "xFemMatrix.h"
#include "xField.h"
#include "xGeomElem.h"
#include "xLinearSystemSolverSuperLU.h"
#include "xMaterial.h"
#include "xMaterialManager.h"
#include "xMaterialSensitivity.h"
#include "xMesh.h"
#include "xMeshCut.h"
#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
#include "xSolVisitor.h"
#include "xSpace.h"
#include "xStateOfValue.h"
#include "xTensors.h"
#include "xValue.h"
#include "xValueCreators.h"
#include "xValueManager.h"
#include "xVariabManager.h"
#include "xVectorField.h"
#include "xZone.h"
#include "xcFormLinearEnergyRelease.h"

// Trellis include

#include "mEdge.h"
#include "mEntity.h"
#include "mVertex.h"

#define NEW_LEVEL_SET

using namespace AOMD;
using namespace xfem;
// using namespace xlinalg;
using namespace std;
namespace xcrack
{
std::string int2string(int i)
{
   char str[30];
   sprintf(str, "%d", i);
   return std::string(str);
}

void lCrack::printSifs(std::string filename, const xValueManagerDist<double>& double_manager,
                       const std::function<double(const xtensor::xPoint&)>& p_to_s) const
{
   cout << "entering printSifs" << endl;
   xIsPhys test_1d("J_front");
   set<double> s_done;
   FILE* OUTPUT2 = fopen(filename.c_str(), "w");
   fprintf(OUTPUT2, "x y z s K1 K2 K3\n");

   for (const auto& e_val : double_manager)
   {
      if (test_1d(e_val.first.getPhys()))
      {
         mVertex* v = static_cast<mVertex*>(e_val.first.getEnti());
         xtensor::xPoint p = v->point();
         double s = p_to_s(p);
         if (s_done.find(s) == s_done.end())
         {
            xcValueJsAndSifs::sifs(1);
            double K1 = e_val.second->getVal();
            xcValueJsAndSifs::sifs(2);
            double K2 = e_val.second->getVal();
            xcValueJsAndSifs::sifs(3);
            double K3 = e_val.second->getVal();
            fprintf(OUTPUT2, "%e %e %e %e %e %e %e \n", p(0), p(1), p(2), s, K1, K2, K3);
            s_done.insert(s);
         }
      }
   }
   // sortie de la map

   fclose(OUTPUT2);
}

void lCrack::getJint3DParks(xField<>& disp_l, xIntegrationRule& integrator_vol, xIntegrationRule& integrator_bnd,
                            const double& rho, const int& nb_internodes, std::string filename,
                            const std::function<double(const xtensor::xPoint&)>& p_to_s)
{
   assert(dim() == 3);
   const bool debug = true;
   const bool K_zero = false;
   if (debug) cout << "using Parks method " << endl;

   // mesh_crack_front->modifyAllState();
   AOMD::mMesh& mmesh = mesh_crack_front->getMesh();
   mmesh.modifyState(3, 2, true);
   mmesh.modifyState(3, 1, true);
   mmesh.modifyState(3, 0, true);
   mmesh.modifyState(2, 1, true);
   mmesh.modifyState(2, 0, true);
   mmesh.modifyState(1, 0, true);
   mmesh.modifyState(0, 1, true);
   mmesh.modifyState(0, 2, true);
   mmesh.modifyState(0, 3, true);
   mmesh.modifyState(1, 2, true);
   mmesh.modifyState(1, 3, true);
   mmesh.modifyState(2, 3, true);

   // création maillage 1D pour les dofs
   xMesh mesh_crack_front_dofs_owner;
   xMesh* mesh_crack_front_dofs = &mesh_crack_front_dofs_owner;

   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> e0d_s_loc;
   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> e1d_dof_length;
   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey> e0d_dof_e0d;
   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey> e0d_e1d_dof;

   extractMeshCrackFrontDofsForJint3DParks(mesh_crack_front_dofs, e0d_s_loc, e1d_dof_length, e0d_dof_e0d, e0d_e1d_dof,
                                           nb_internodes);

   cout << "nb of dofs on the front " << e0d_dof_e0d.size() << endl;
   AOMD_Util::Instance()->ex_port("mesh_crack_front.msh", &mesh_crack_front->getMesh());
   AOMD_Util::Instance()->ex_port("mesh_crack_front_dofs.msh", &mesh_crack_front_dofs->getMesh());

   // declaration du systeme
   xValueManagerDist<double> double_manager;

   xSpaceLagrange J_front("J_front", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
   xField<> J_front_l(&double_manager, J_front);

   mEntity* e = *mesh->begin(3);
   xfem::xGeomElem geo(e);
   xMaterial* mat = xMaterialManagerSingleton::instance().getMaterial(&geo);
   const xTensors* properties = mat->getProperties();
   double young = properties->scalar("YOUNG_MODULUS");
   double poisson = properties->scalar("POISSON_RATIO");
   double Estar = young / (1. - poisson * poisson);

   xcValueJs::setNbModes(1);
   xcValueSifs::setNbModes(3);
   xcValueSifs::setGeom(xcValueSifs::GEOM_3D);
   xcValueSifs::setYoungAndPoisson(young, poisson);

   xValueCreator<xcValueJsAndSifs> creator_double;
   DeclareInterpolation(J_front_l, creator_double, mesh_crack_front_dofs->begin(1), mesh_crack_front_dofs->end(1));
   xStateDofCreator<> snh(double_manager, "dofs");
   DeclareState(J_front_l, snh, mesh_crack_front_dofs->begin(1), mesh_crack_front_dofs->end(1));
   ValueCreatorFront creator_linear_combination(&double_manager, e0d_s_loc, e1d_dof_length, e0d_e1d_dof);
   DeclareInterpolation(J_front_l, creator_linear_combination, mesh_crack_front->begin(1), mesh_crack_front->end(1));

   //
   double_manager.PrintForDebug("dcl_J_front.dbg");

   // creation de l'ensemble des elements utiles pour l'integrale J : "domain_for_integral"
   const xLevelSet& lsn = *getFieldn();
   const xLevelSet& lst = *getFieldt();

   xSubMesh& domain_for_integral = mesh->createSubMesh("domain_for_integral");
   for (mEntity* v : mesh->range(0))
   {
      double a = lsn(v);
      double b = lst(v);
      double r = sqrt(a * a + b * b);
      if (r <= rho)
      {
         for (int i = 0; i < v->size(3); ++i)
         {
            mEntity* e = v->get(3, i);
            domain_for_integral.add(e);
            for (int ii = 0; ii < 3; ++ii)
            {
               for (int j = 0; j < e->size(ii); ++j)
               {
                  mEntity* vv = e->get(ii, j);
                  domain_for_integral.add(vv);
               }
            }
         }
      }
   }
   if (domain_for_integral.size(3) == 0)
   {
      cerr << " you must increase the size of rho for the domain integral because the domain is empty " << endl;
      assert(0);
      abort();
   }
   // xRegion domain_for_integral(mesh, "domain_for_integral");

   xFilteredRegion<xIter, xAcceptOnBoundary> domain_for_integral_ext_bnd(domain_for_integral.begin(2), domain_for_integral.end(2),
                                                                         xAcceptOnBoundary());

   //
   xSubMesh& domain_for_integral_outer_layer = mesh->createSubMesh("domain_for_integral_outer_layer");
   for (mEntity* v : domain_for_integral.range(0))
   {
      double a = lsn(v);
      double b = lst(v);
      double r = sqrt(a * a + b * b);
      if (r >= rho)
      {
         for (int i = 0; i < v->size(3); ++i)
         {
            mEntity* e = v->get(3, i);
            if (domain_for_integral.find(e)) domain_for_integral_outer_layer.add(e);
         }
      }
   }

   // xRegion domain_for_integral_outer_layer(mesh, "domain_for_integral_outer_layer");

   // if (debug) mesh->printSubsetEntities();

   // for all nodes in "domain_for_integral", we set the corresponding s_loc and e1d

   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> e0d3d_s_loc;
   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey> e0d3d_e1d_dof;

   set3DElementsToFrontRelation(xRegion(&domain_for_integral), e0d_s_loc, e1d_dof_length, e0d_e1d_dof, e0d3d_s_loc,
                                e0d3d_e1d_dof);
   double_manager.PrintForDebug("dcl_inside.dbg");

   // level set creation.
#ifdef OLD_LEVEL_SET
   // old version
   mesh->allocateSubsetEntities("domain_for_sloc_level_set");
   xRegion domain_for_sloc_level_set(mesh, "domain_for_sloc_level_set");
   xLevelSet sloc_level_set(domain_for_sloc_level_set);
   bool sloc_ok = true;
   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey> e0d3d_e1d_dof_other;

   setSlocLevelSet(domain_for_integral, sloc_level_set, e0d3d_s_loc, e1d_dof_length, e0d3d_e1d_dof, e0d3d_e1d_dof_other, sloc_ok);

   cout << "sloc_ok is :" << sloc_ok << endl;
   if (!sloc_ok) return;
   sloc_level_set.exportGmsh("sloc_level_set");
   if (debug) cout << "before line xcut::xPhysSurf transitions " << endl;
   cutElementsOnTransitions(sloc_level_set);
   xcut::xPhysSurf transitions(sloc_level_set, 0);
   if (debug) cout << "adter  line xcut::xPhysSurf transitions " << endl;

   // J_3d
   xSpaceLagrange lagrange_3d("J_3d", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
   xField<> J_3d_l(&double_manager, lagrange_3d);

   xScalarFunctionDerivDiscXFEM enrichment_function(transitions);
   xSpaceFiltered::filter_t filter(bind1st(mem_fun(&xcut::xPhysSurf::boundary_strict), &transitions));

   xValKeyExtend key_modifier_center("_CENTER");
   xSpaceXFEM space_full_center(lagrange_3d, enrichment_function, key_modifier_center);
   xSpaceFiltered enriched_center(space_full_center, filter);
   J_3d_l.insert(enriched_center);

   xValKeyExtend key_modifier_left("_LEFT");
   xSpaceXFEM space_full_left(lagrange_3d, enrichment_function, key_modifier_left);
   xSpaceFiltered enriched_left(space_full_left, filter);
   J_3d_l.insert(enriched_left);

   xValKeyExtend key_modifier_right("_RIGHT");
   xSpaceXFEM space_full_right(lagrange_3d, enrichment_function, key_modifier_right);
   xSpaceFiltered enriched_right(space_full_right, filter);
   J_3d_l.insert(enriched_right);

   ValueCreatorJ3dOld creator_linear_combination_3d(&double_manager, e0d3d_s_loc, e1d_dof_length, e0d3d_e1d_dof,
                                                    e0d3d_e1d_dof_other);

   DeclareInterpolation(J_3d_l, creator_linear_combination_3d, domain_for_integral.begin(), domain_for_integral.end());
   double_manager.PrintForDebug("dcl_J_3d.dbg");

#endif

// a level set is defined for each e0d_dof
#ifdef NEW_LEVEL_SET

   std::unordered_map<mEntity*, xLevelSet*, EntityHashKey, EntityEqualKey> e0d_dof_sloc_level_set;

   setAllSlocLevelSets(xRegion(&domain_for_integral), e0d_dof_sloc_level_set, e0d3d_s_loc, e1d_dof_length, e0d3d_e1d_dof);

   if (debug)
   {
      for (const auto& pe_pls : e0d_dof_sloc_level_set)
      {
         xexport::xExportGmshAscii pexport;
         Export(*pe_pls.second, pexport, "sloc_level_set_" + int2string(pe_pls.first->getId()));
      }
   }

   xSpaceLagrange lagrange_3d("J_3d", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
   xField<> J_3d_l(&double_manager, lagrange_3d);
   std::unordered_map<mEntity*, xcut::xPhysSurf*, EntityHashKey, EntityEqualKey> e0d_dof_phys_surf;

   for (const auto& pe_pls : e0d_dof_sloc_level_set)
   {
      int e0d_dof_id = pe_pls.first->getId();
      if (debug) cout << "enriching for node " << e0d_dof_id << endl;
      xLevelSet& sloc_level_set = *pe_pls.second;

      if (debug)
      {
         xRegion r = sloc_level_set.getSupport();
         cout << " region name is " << r.getSubName() << endl;
      }

      xFitToVertices fit(1.e-2);
      sloc_level_set.accept(fit);

      cutElementsOnTransitions(sloc_level_set);

      // Problem here : this xcut::xPhysSurf Constructor does not exist anymore ...
      //	xcut::xPhysSurf* transitions = new xcut::xPhysSurf(sloc_level_set, 0);
      xcut::xPhysSurf* transitions = nullptr;

      e0d_dof_phys_surf.insert(make_pair(pe_pls.first, transitions));

      // xcut::xPhysSurf transitions(sloc_level_set, 0);

      // xScalarFunctionDerivDiscXFEM enrichment_function(transitions);
      xScalarFunctionDerivDiscXFEM enrichment_function(sloc_level_set);
      xSpaceFiltered::filter_t filter(bind1st(mem_fun(&xcut::xPhysSurf::boundary_strict), transitions));

      xcValKeyExtendAndSetRefe key_modifier_center("_CENTER", e0d_dof_id);
      xSpaceXFEM space_full_center(lagrange_3d, enrichment_function, key_modifier_center);
      xSpaceFiltered enriched_center(space_full_center, filter);
      J_3d_l.insert(enriched_center);

      xcValKeyExtendAndSetRefe key_modifier_left("_LEFT", e0d_dof_id);
      xSpaceXFEM space_full_left(lagrange_3d, enrichment_function, key_modifier_left);
      xSpaceFiltered enriched_left(space_full_left, filter);
      J_3d_l.insert(enriched_left);

      xcValKeyExtendAndSetRefe key_modifier_right("_RIGHT", e0d_dof_id);
      xSpaceXFEM space_full_right(lagrange_3d, enrichment_function, key_modifier_right);
      xSpaceFiltered enriched_right(space_full_right, filter);
      J_3d_l.insert(enriched_right);
   }
   if (debug) cout << "done with the enrichment" << endl;

   ValueCreatorJ3d creator_linear_combination_3d(&double_manager, e0d3d_s_loc, e1d_dof_length, e0d3d_e1d_dof,
                                                 mesh_crack_front_dofs);

   DeclareInterpolation(J_3d_l, creator_linear_combination_3d, domain_for_integral.begin(), domain_for_integral.end());
   double_manager.PrintForDebug("dcl_J_3d.dbg");
#endif

   // fonction fr pour Parks
   xcEvalFrParks eval_fr(*this, rho, xRegion(&domain_for_integral), xRegion(&domain_for_integral_outer_layer));
   xcEvalGradFrParks eval_grad_fr(*this, rho, xRegion(&domain_for_integral), xRegion(&domain_for_integral_outer_layer));

   // q_vector est la direction de grad_lst partout
   // sauf sur les bords ou on ne garde que la partie
   // tangente Ã  la surface
   xVectorField q_dir(&domain_for_integral);
   if (debug) cout << "before setQdir" << endl;
   setQDirectionParks(q_dir);
   if (debug) cout << "after setQdir" << endl;
   if (debug) q_dir.exportGmsh("q_dir");
   if (debug) cout << "after export Qdir" << endl;

   xexport::xExportGmshAscii pexport;

   // q global
   //
   xcEvalQGlobalVector eval_q_global(q_dir, eval_fr);
   xcEvalGradQGlobalVector eval_grad_q_global(q_dir, eval_fr, eval_grad_fr);

   Export(eval_q_global, pexport, "q_global", integrator_vol, domain_for_integral.begin(), domain_for_integral.end());
   Export(eval_q_global, pexport, "q_global_bnd_ext", integrator_bnd, domain_for_integral_ext_bnd.begin(),
          domain_for_integral_ext_bnd.end(), xUpperAdjacency());

   Export(eval_q_global, pexport, "q_global_front", integrator_vol, mesh_crack_front->begin(1), mesh_crack_front->end(1),
          xUpperCreatorRecursive());

   if (debug)
   {
      xcValueJsAndSifs::js(1);
      for (mEntity* v : mesh_crack_front_dofs->range(0))
      {
         string domain_loc = "domain_for_sloc_level_set_" + int2string(v->getId());

         xValKey key_debug(xKeyInfo::getPhysId("J_front"), xKeyInfo::getGeomId("HIERARCHICAL_1"), v);
         xValue<double>* val = double_manager.find(key_debug);
         val->setVal(1.0);

         xEvalField<xtool::xIdentity<double>> eval_approx(J_3d_l);
         xEvalBinary<xtool::xMult<xtensor::xVector<>, double, xtensor::xVector<>>> eval_approx_vec(eval_q_global, eval_approx);
         Export(eval_approx_vec, pexport, "approx_vec_" + int2string(v->getId()), integrator_vol, domain_for_integral.begin(),
                domain_for_integral.end());

         Export(eval_approx, pexport, "approx_sca_" + int2string(v->getId()), integrator_vol, domain_for_integral.begin(),
                domain_for_integral.end());

         xEvalGradField<xtool::xIdentity<xtensor::xVector<>>> eval_grad_approx(J_3d_l);
         Export(eval_grad_approx, pexport, "approx_grad_sca_" + int2string(v->getId()), integrator_vol,
                domain_for_integral.begin(), domain_for_integral.end());

         Export(eval_approx_vec, pexport, "approx_vec_on_front_trace_" + int2string(v->getId()), integrator_vol,
                mesh_crack_front->begin(1), mesh_crack_front->end(1), xUpperCreatorRecursive());

         Export(eval_approx, pexport, "approx_sca_on_front_trace_" + int2string(v->getId()), integrator_vol,
                mesh_crack_front->begin(1), mesh_crack_front->end(1), xUpperCreatorRecursive());

         xEvalField<xtool::xIdentity<double>> eval_approx_front(J_front_l);
         Export(eval_approx_front, pexport, "approx_sca_on_front_" + int2string(v->getId()), integrator_vol,
                mesh_crack_front->begin(1), mesh_crack_front->end(1));

         val->setVal(0.0);
      }
      for (const auto& pe_pls : e0d_dof_sloc_level_set)
      {
         mEntity* v = pe_pls.first;
         string domain_loc = "domain_for_sloc_level_set_" + int2string(v->getId());

         xValKey key_debug(xKeyInfo::getPhysId("J_front"), xKeyInfo::getGeomId("HIERARCHICAL_1"), v);
         xValue<double>* val = double_manager.find(key_debug);
         val->setVal(1.0);

         xEvalField<xtool::xIdentity<double>> eval_approx(J_3d_l);

         // Export(eval_approx, pexport,"approx_sca_loc_" + int2string(v->getId()), integrator_vol,
         //	   mesh->begin_sub(3, domain_loc), mesh->end_sub(3, domain_loc) );

         xEvalGradField<xtool::xIdentity<xtensor::xVector<>>> eval_grad_approx(J_3d_l);

         // Export(eval_grad_approx, pexport,"approx_grad_sca_loc_" + int2string(v->getId()), integrator_vol,
         //	   mesh->begin_sub(3, domain_loc), mesh->end_sub(3, domain_loc) );

         val->setVal(0.0);
      }
   }

   // codage de la matrice
   xlinalg::xCSRVector RHS(double_manager.size("dofs"));
   xlinalg::xCSRVector sol(double_manager.size("dofs"));
   xlinalg::xCSRMatrix M(double_manager.size("dofs"));

   xlinalg::xLinearSystemSolverSuperLU<> solver;
   // xLinearSystem system(&M, &solver);

   // xAssemblerBasic<> assembler_mat(M);
   xAssemblerLumped<> assembler_mat(M);

   xIntegrationRuleBasic integration_rule_MAT(2);
   xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<double>>, xValOperator<xtool::xIdentity<double>>> L2_bilinear;
   if (!K_zero)
      Assemble(L2_bilinear, assembler_mat, integration_rule_MAT, J_front_l, J_front_l, mesh_crack_front->begin(1),
               mesh_crack_front->end(1));
   //    if (!K_zero) Assemble(L2_bilinear, assembler_mat, integration_rule_MAT, J_3d_l, J_3d_l,
   //			  mesh_crack_front->begin(1),
   //			  mesh_crack_front->end(1), xUpperCreatorRecursive() );

   xAssemblerBasic<> assembler_rhs(RHS);

   xEvalGradField<xtool::xIdentity<xtensor::xTensor2<>>> grad_disp(disp_l);
   xUniformMaterialSensitivity<xtensor::xTensor4<>> hooke("strain");
   xEvalField<xtool::xIdentity<double>> eval_sol(J_3d_l);
   xEvalBinary<xtool::xMult<xtensor::xVector<>, double, xtensor::xVector<>>> eval_sol_vec(eval_q_global, eval_sol);

#if 1
   // J computation
   cout << " solving for J energy release " << endl;
   xcValueJsAndSifs::js(1);
   xcEvalEshelbyHPP eval_eshelby(grad_disp, hooke);
   xEvalConstant<xtensor::xVector<>> eval_source_term(xtensor::xVector<>(0.));

   // test purpose solution must be one
   // xEvalConstant<double> load(1.0);
   // xFormLinearWithLoad<xValOperator<xtool::xIdentity<double> >,  xEval<double> > unit_load_form(load);;
   // Assemble(unit_load_form, assembler_rhs, integrator_vol, J_front_l, mesh_crack_front->begin(1), mesh_crack_front->end(1));

   xcFormLinearEnergyRelease form_linear_J(eval_q_global, eval_grad_q_global, eval_eshelby, eval_source_term);
   if (!K_zero)
      Assemble(form_linear_J, assembler_rhs, integrator_vol, J_3d_l, domain_for_integral.begin(), domain_for_integral.end());

   if (!K_zero)
   {
      solver.connectMatrix(M);
      solver.solve(RHS, sol);
   }
   Visit(xWriteSolutionVisitor<>(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));

   Export(eval_sol_vec, pexport, "J_vec", integrator_vol, mesh_crack_front->begin(1), mesh_crack_front->end(1),
          xUpperCreatorRecursive());
   Export(eval_sol, pexport, "J_sca", integrator_vol, mesh_crack_front->begin(1), mesh_crack_front->end(1),
          xUpperCreatorRecursive());
#endif

   // Sifs computation

   for (int mode = 1; mode <= 3; ++mode)
   {
      cout << " solving sifs for mode " << mode << endl;
      xcValueJsAndSifs::sifs(mode);
      xcEvalEshelbyHPPInteraction eval_eshelby_interaction(grad_disp, hooke, *this, mode);

      RHS.ZeroArray();
      xcFormLinearEnergyRelease form_linear_interaction(eval_q_global, eval_grad_q_global, eval_eshelby_interaction,
                                                        eval_source_term);
      if (!K_zero)
         Assemble(form_linear_interaction, assembler_rhs, integrator_vol, J_3d_l, domain_for_integral.begin(),
                  domain_for_integral.end());

      // xcFormLinearEnergyRelease form_linear_interaction(eval_q_global, eval_grad_q_global, eval_eshelby_interaction);
      // if (!K_zero) Assemble(form_linear_interaction, assembler_rhs, integrator_vol, J_3d_l,
      // domain_for_integral_outer_layer.begin(), domain_for_integral_outer_layer.end());

      xcFormLinearEnergyReleaseBoundary form_linear_interaction_boundary(eval_q_global, eval_eshelby_interaction);
      if (!K_zero)
         Assemble(form_linear_interaction_boundary, assembler_rhs, integrator_bnd, J_3d_l, domain_for_integral_ext_bnd.begin(),
                  domain_for_integral_ext_bnd.end(), xUpperAdjacency());

      sol.ZeroArray();
      if (!K_zero)
      {
         solver.connectMatrix(M);
         solver.solve(RHS, sol);
      }
      transform(sol.begin(), sol.end(), sol.begin(), bind2nd(multiplies<double>(), 0.5 * Estar));

      Visit(xWriteSolutionVisitor<>(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));

      Export(eval_sol_vec, pexport, "sif_parks_mode_" + int2string(mode) + "_vec", integrator_vol, mesh_crack_front->begin(1),
             mesh_crack_front->end(1), xUpperCreatorRecursive());
      Export(eval_sol, pexport, "sif_parks_mode_" + int2string(mode) + "_sca", integrator_vol, mesh_crack_front->begin(1),
             mesh_crack_front->end(1), xUpperCreatorRecursive());
   }

   // Err export
   xcValueJsAndSifs::err();
   Export(eval_sol_vec, pexport, "K_err_vec", integrator_vol, mesh_crack_front->begin(1), mesh_crack_front->end(1),
          xUpperCreatorRecursive());

   //   std::ofstream out(filename.c_str());
   printSifs(filename, double_manager, p_to_s);
   //   out.close();

#ifdef NEW_LEVEL_SET
   for (std::unordered_map<mEntity*, xLevelSet*, EntityHashKey, EntityEqualKey>::const_iterator it =
            e0d_dof_sloc_level_set.begin();
        it != e0d_dof_sloc_level_set.end(); ++it)
   {
      delete it->second;
   }
   for (std::unordered_map<mEntity*, xcut::xPhysSurf*, EntityHashKey, EntityEqualKey>::const_iterator it =
            e0d_dof_phys_surf.begin();
        it != e0d_dof_phys_surf.end(); ++it)
   {
      delete it->second;
   }
#endif

   return;
}

xValue<double>* ValueCreatorFront::operator()(const xValKey& key) const  // ou xValueDouble
{
   const bool debug = false;
   mEntity* v = key.getEnti();
   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator it = e0d_e1d_dof.find(v);
   if (it == e0d_e1d_dof.end())
   {
      cerr << " node e0d " << v->getId() << " has no associated e1d_dof " << endl;
      assert(it != e0d_e1d_dof.end());
      abort();
   }
   mEntity* e1d_dof = it->second;
   mEntity* v_prev_dof = e1d_dof->get(0, 0);
   mEntity* v_next_dof = e1d_dof->get(0, 1);

   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator it_s = e0d_s_loc.find(v);
   assert(it_s != e0d_s_loc.end());
   double sloc = it_s->second;
   if (debug) cout << " s_loc of e0d is " << sloc << " for the node " << endl;
   if (debug) v->print();

   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator it_length = e1d_dof_length.find(e1d_dof);
   assert(it_length != e1d_dof_length.end());
   double length_e1d = it_length->second;
   if (debug) cout << " length of e1D_dof is " << length_e1d << " for the edge " << endl;
   if (debug) e1d_dof->print();

   xValKey key_prev = key;
   key_prev.setEnti(v_prev_dof);

   xValKey key_next = key;
   key_next.setEnti(v_next_dof);

   std::vector<xValue<double>*> values_t;
   xValue<double>* val_prev = double_manager->find(key_prev);
   xValue<double>* val_next = double_manager->find(key_next);
   values_t.push_back(val_prev);
   values_t.push_back(val_next);

   std::vector<double> coeffs_t;
   coeffs_t.push_back(1. - sloc / length_e1d);
   coeffs_t.push_back(sloc / length_e1d);

   return new xValueLinearCombination<double>(coeffs_t, values_t);  // a voir
}

xValue<double>* ValueCreatorJ3dOld::operator()(const xValKey& key) const
{
   const bool debug = false;
   if (debug) cout << " key in ValueCreatorJ3d" << key << endl;
   mEntity* v = key.getEnti();
   xIsGeom test("HIERARCHICAL_1");
   xIsGeom test_left("HIERARCHICAL_1_LEFT");
   xIsGeom test_center("HIERARCHICAL_1_CENTER");
   xIsGeom test_right("HIERARCHICAL_1_RIGHT");
   assert(test(key) || test_left(key) || test_center(key) || test_right(key));

   if (test(key))
   {
      xValKey key_tmp = key;
      key_tmp.setPhys(xKeyInfo::getPhysId("J_front"));
      return value_creator_front(key_tmp);
   }

   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator it = e0d3d_e1d_dof.find(v);
   if (it == e0d3d_e1d_dof.end())
   {
      cerr << " node e0d3d " << v->getId() << " has no associated e1d_dof " << endl;
      assert(it != e0d3d_e1d_dof.end());
      abort();
   }
   mEntity* e1d_dof = it->second;

   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator ito = e0d3d_e1d_dof_other.find(v);
   if (ito == e0d3d_e1d_dof_other.end())
   {
      cerr << " node e0d3d " << v->getId() << " has no associated e1d_dof_other " << endl;
      assert(ito != e0d3d_e1d_dof_other.end());
      abort();
   }
   mEntity* e1d_dof_other = ito->second;

   mEdge* e1d_dof_left = (mEdge*)e1d_dof;
   mEdge* e1d_dof_right = (mEdge*)e1d_dof_other;
   lCrack::orderEdges(e1d_dof_left, e1d_dof_right);

   mEntity* e0d_dof_left = e1d_dof_left->get(0, 0);
   mEntity* e0d_dof_center = e1d_dof_left->get(0, 1);
   mEntity* e0d_dof_right = e1d_dof_right->get(0, 1);

   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator itl = e1d_dof_length.find(e1d_dof_left);
   double l_left = itl->second;
   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator itr = e1d_dof_length.find(e1d_dof_right);
   double l_right = itr->second;

   double coeff_left = -1. / (2. * l_left);
   double coeff_center = 1. / (2. * l_left) + 1. / (2. * l_right);
   double coeff_right = -1. / (2. * l_right);

   xValue<double>* ret = nullptr;
   if (test_left(key))
   {
      xValKey key_left = key;
      key_left.setEnti(e0d_dof_left);
      key_left.setPhys(phys_id);
      key_left.setGeom(geom_id);
      xValue<double>* val_left = double_manager->find(key_left);
      ret = new xValueLinearCombination<double>(coeff_left, val_left);
   }
   if (test_center(key))
   {
      xValKey key_center = key;
      key_center.setEnti(e0d_dof_center);
      key_center.setPhys(phys_id);
      key_center.setGeom(geom_id);
      xValue<double>* val_center = double_manager->find(key_center);
      ret = new xValueLinearCombination<double>(coeff_center, val_center);
   }
   if (test_right(key))
   {
      xValKey key_right = key;
      key_right.setEnti(e0d_dof_right);
      key_right.setPhys(phys_id);
      key_right.setGeom(geom_id);
      xValue<double>* val_right = double_manager->find(key_right);
      ret = new xValueLinearCombination<double>(coeff_right, val_right);
   }

   return ret;
}

xValue<double>* ValueCreatorJ3d::operator()(const xValKey& key) const
{
   const bool debug = false;
   if (debug) cout << " key in ValueCreatorJ3d" << key << endl;
   mEntity* v = key.getEnti();
   if (debug) cout << " location is " << ((mVertex*)v)->point() << endl;
   if (debug) cout << " node is     " << ((mVertex*)v)->getId() << endl;

   xIsGeom test("HIERARCHICAL_1");
   xIsGeom test_left("HIERARCHICAL_1_LEFT");
   xIsGeom test_center("HIERARCHICAL_1_CENTER");
   xIsGeom test_right("HIERARCHICAL_1_RIGHT");
   assert(test(key) || test_left(key) || test_center(key) || test_right(key));

   if (test(key))
   {
      xValKey key_tmp = key;
      key_tmp.setPhys(xKeyInfo::getPhysId("J_front"));
      return value_creator_front(key_tmp);
   }

   int e0d_dof_id = key.getRefe();
   if (debug) cout << " in  ValueCreatorJ3d considering node " << e0d_dof_id << " on the front " << endl;
   mVertex* e0d_dof = mesh_front_dofs->getMesh().getVertex(e0d_dof_id);

   std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator it = e0d3d_e1d_dof.find(v);
   if (it == e0d3d_e1d_dof.end())
   {
      cerr << " node e0d3d " << v->getId() << " has no associated e1d_dof " << endl;
      assert(it != e0d3d_e1d_dof.end());
      abort();
   }
   mEntity* e1d_dof = it->second;
   if (debug)
   {
      cout << " e1d_dof after declaration is " << endl;
      e1d_dof->print();
   }

   // std::unordered_map<mEntity*,mEntity*,EntityHashKey,EntityEqualKey>::const_iterator ito = e0d3d_e1d_dof_other.find(v);
   // if (ito == e0d3d_e1d_dof_other.end())
   //  {
   //	cerr << " node e0d3d " << v->getId() << " has no associated e1d_dof_other " << endl;
   //	assert(ito != e0d3d_e1d_dof_other.end());
   //    abort();
   //
   //}
   mEntity* e1d_dof_other;
   if (debug) cout << " e0d_dof->size(1) is " << e0d_dof->size(1) << endl;
   assert(e0d_dof->size(1) == 2);

   mEntity* eda = e0d_dof->get(1, 0);
   mEntity* edb = e0d_dof->get(1, 1);
   if (debug)
   {
      cout << " eda is " << endl;
      eda->print();
      cout << " edb is " << endl;
      edb->print();
   }

   if (eda == e1d_dof)
      e1d_dof_other = edb;
   else
      e1d_dof_other = eda;

   if (debug)
   {
      cout << " e1d_dof is " << endl;
      e1d_dof->print();
      cout << " e1d_dof_other is " << endl;
      e1d_dof_other->print();
   }

   mEdge* e1d_dof_left = (mEdge*)e1d_dof;
   mEdge* e1d_dof_right = (mEdge*)e1d_dof_other;
   if (debug) cout << " before order edges " << endl;
   lCrack::orderEdges(e1d_dof_left, e1d_dof_right);

   mEntity* e0d_dof_left = e1d_dof_left->get(0, 0);
   mEntity* e0d_dof_center = e1d_dof_left->get(0, 1);
   mEntity* e0d_dof_right = e1d_dof_right->get(0, 1);

   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator itl = e1d_dof_length.find(e1d_dof_left);
   double l_left = itl->second;
   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator itr = e1d_dof_length.find(e1d_dof_right);
   double l_right = itr->second;

   double coeff_left = -1. / (2. * l_left);
   double coeff_center = 1. / (2. * l_left) + 1. / (2. * l_right);
   double coeff_right = -1. / (2. * l_right);

   xValue<double>* ret = nullptr;
   if (test_left(key))
   {
      xValKey key_left(phys_id, geom_id, e0d_dof_left);
      // key_left.setEnti(e0d_dof_left);
      // key_left.setPhys(phys_id);
      // key_left.setGeom(geom_id);
      xValue<double>* val_left = double_manager->find(key_left);
      ret = new xValueLinearCombination<double>(coeff_left, val_left);
   }
   if (test_center(key))
   {
      xValKey key_center(phys_id, geom_id, e0d_dof_center);
      // key_center.setEnti(e0d_dof_center);
      // key_center.setPhys(phys_id);
      // key_center.setGeom(geom_id);
      xValue<double>* val_center = double_manager->find(key_center);
      ret = new xValueLinearCombination<double>(coeff_center, val_center);
   }
   if (test_right(key))
   {
      xValKey key_right(phys_id, geom_id, e0d_dof_right);
      // key_right.setEnti(e0d_dof_right);
      // key_right.setPhys(phys_id);
      // key_right.setGeom(geom_id);
      xValue<double>* val_right = double_manager->find(key_right);
      ret = new xValueLinearCombination<double>(coeff_right, val_right);
   }

   return ret;
}

double lCrack::getClosestSegmentFromTo(mVertex* v, xMesh* mesh1d, mEntity*& edge_closest, double& percent_position) const
{
   const bool debug = false;
   xtensor::xPoint p = v->point();
   bool first = true;
   double mem1 = -1., alpha, mem2 = 0.;
   assert(mesh1d->size(1) != 0);
   for (mEntity* edge : mesh1d->range(1))
   {
      //  ... on va chercher les deux noeuds du segment
      mVertex* v1 = (mVertex*)edge->get(0, 0);
      mVertex* v2 = (mVertex*)edge->get(0, 1);
      //  ... on a besoin des coordonnees des deux noeuds
      xtensor::xPoint p1 = v1->point();
      xtensor::xPoint p2 = v2->point();

      xtensor::xVector<> P(p1, p);
      xtensor::xVector<> AB(p1, p2);
      double mag = AB.mag();
      alpha = (AB * P) / (mag * mag);
      if (alpha < 0.) alpha = 0.;
      if (alpha > 1.) alpha = 1.;
      xtensor::xVector<> PP = P - AB * alpha;
      double dist = PP.mag();
      if (first || (dist < mem1))
      {
         mem1 = dist;
         edge_closest = edge;
         mem2 = alpha;
         first = false;
      }
   }
   percent_position = mem2;
   if (debug)
   {
      cout << " in getClosestSegmentFromTo before return " << endl;
      cout << " results for node " << v->getId() << endl;
      cout << " at position " << v->point() << endl;
      cout << " closest edge is " << endl;
      edge_closest->print();
      mEntity* v1 = edge_closest->get(0, 0);
      mEntity* v2 = edge_closest->get(0, 1);
      cout << " with nodes " << endl;
      v1->print();
      v2->print();
      cout << "percent_position " << percent_position << endl;
   }

   return mem1;
}

void lCrack::cutElementsOnTransitions(const xLevelSet& sloc_level_set) const
{
   /*const bool debug = false;
   xMesh *interface_transition =new xMesh;
   interface_transition->createInterfaceFromLevelSet(sloc_level_set, xMesh::get_was_created_by(), xMesh::get_r_on_edge());
   if (debug) AOMD_Util::Instance()->ex_port("interface_transition.msh", interface_transition);
   mesh->cutAlongInterfaceRecursive(interface_transition, sloc_level_set);
   // warning no delete
   */
   const bool debug = false;
   xMesh interface_transition;
   xcut::createIsoZeroMeshFromLevelSet(sloc_level_set, interface_transition, xMesh::get_was_created_by(), xMesh::get_r_on_edge());
   if (debug) AOMD_Util::Instance()->ex_port("interface_transition.msh", &interface_transition.getMesh());
   xcut::cutAlongInterfaceRecursive(interface_transition, sloc_level_set, xMesh::get_const_was_created_by(),
                                    xMesh::get_partition(), xMesh::get_is_duplicated_in(), xMesh::get_was_duplicated_from(),
                                    xMesh::get_is_in_partition_of());
}

void lCrack::setSlocLevelSet(const xRegion& domain_for_integral, xLevelSet& sloc_level_set,
                             const std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e0d3d_s_loc,
                             const std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e1d_dof_length,
                             const std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>& e0d3d_e1d_dof,
                             std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>& e0d3d_e1d_dof_other,
                             bool& sloc_ok) const
{
   /*
   const bool debug = false;
   std::unordered_set<mEntity*,EntityHashKey,EntityEqualKey> edges;
   for (xIter it = domain_for_integral.begin(); it != domain_for_integral.end(); ++it)
     {
       mEntity* e = *it;
       if (debug)
         {
           cout << "treating element in setSlocLevelSet " << endl;
           e->print();
         }
       edges.clear();
       for (int i = 0; i < e->size(0); ++i)
         {
           mEntity* v = e->get(0,i);
           std::unordered_map<mEntity*,mEntity*,EntityHashKey,EntityEqualKey>::const_iterator ite = e0d3d_e1d_dof.find(v);
           assert(ite != e0d3d_e1d_dof.end());
           edges.insert(ite->second);
         }
       assert(edges.size() > 0);
       if (edges.size() == 1)
         {
         }
       else if (edges.size() == 2)
         {
           mesh->add_sub(e,"domain_for_sloc_level_set");
           for (int i=0; i < 3; ++i)
             {
               for (int j=0; j < e->size(i); ++j)
                 {
                   mEntity* vv = e->get(i,j);
                   mesh->add_sub(vv, "domain_for_sloc_level_set");
                 }
             }

           mEdge * ed1 = (mEdge*) *edges.begin();
           mEdge * ed2 = (mEdge*) *(++edges.begin());
           orderEdges(ed1, ed2);
           for (int i = 0; i < e->size(0); ++i)
             {
               mEntity* v = e->get(0,i);
               std::unordered_map<mEntity*,mEntity*,EntityHashKey,EntityEqualKey>::const_iterator ite = e0d3d_e1d_dof.find(v);
               mEntity* ed = ite->second;
               std::unordered_map<mEntity*,double,EntityHashKey,EntityEqualKey>::const_iterator its = e0d3d_s_loc.find(v);
               double s_loc = its->second;
               double val;
               if (ed == ed1)
                 {
                   std::unordered_map<mEntity*,double,EntityHashKey,EntityEqualKey>::const_iterator itl = e1d_dof_length.find(ed);
                   double length = itl->second;
                   val = s_loc - length;
                   e0d3d_e1d_dof_other[v] = ed2;
                 }
               else
                 {
                   val = s_loc;
                   e0d3d_e1d_dof_other[v] = ed1;
                 }
               if (debug)
                 {
                   cout << " for node  " << v->getId() << endl;
                   mVertex* vv = (mVertex*) v;
                   cout << " at position " << vv->point() << endl;
                   cout << " sloc  is " << s_loc << endl;
                   cout << " level set  is " << val << endl;
                 }
               if ( sloc_level_set.isDefinedAt(v) && sloc_level_set(v) != val)
                 {
                   cerr << " A node is being given two different sloc_level_set values " << endl;
                   cerr << " the values are " << sloc_level_set(v) << " and " << val << endl;
                   cerr << " please increase nb_internodes value " << endl;
                   sloc_ok = false;
                   cout<<"sloc_ok 1 is :"<<sloc_ok<<endl;
                   return;
 //		  assert(0);
 //		  abort();
                 }
               sloc_level_set(v) = val;
             }
         }
       else if (edges.size() > 2)
         {
           cerr << " mesh too coarse for SIF computation, please raise nb_internodes" << endl;
           sloc_ok = false;
           cout<<"sloc_ok 2 is :"<<sloc_ok<<endl;
           return;
 //	  assert(0);
 //	  abort();
         }
  //     cout<<"sloc_ok 3 is :"<<sloc_ok<<endl;

     }


   //post-treatment for robustness
   xFitToVertices fit(1.e-2);
   sloc_level_set.accept(fit);

   //to avoid elements in between the slices with a change of sign
   const xRegion& domain_for_sloc_level_set = sloc_level_set.getSupport();
   for (xIter it = domain_for_integral.begin(); it != domain_for_integral.end(); ++it)
     {
       mEntity* e = *it;
       if (!domain_for_sloc_level_set.IsInRegion(e) && sloc_level_set.isDefinedOnElement(e))
         {
           //check is there is no change of sign
           std::vector<double> sloc_ls = sloc_level_set.getVals(e);
           if (*std::min_element(sloc_ls.begin(), sloc_ls.end()) * *std::max_element(sloc_ls.begin(), sloc_ls.end()) < 0.)
             {
               cerr << " problem is sloc_level_set definition, increase nb_internodes " << endl;
               assert(0);
               abort();
             }
         }
     }
   */
}

void lCrack::setAllSlocLevelSets(const xRegion& domain_for_integral,
                                 std::unordered_map<mEntity*, xLevelSet*, EntityHashKey, EntityEqualKey>& e0d_dof_sloc_level_set,
                                 const std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e0d3d_s_loc,
                                 const std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e1d_dof_length,
                                 const std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>& e0d3d_e1d_dof) const
{
   const bool debug = false;
   std::unordered_set<mEntity*, EntityHashKey, EntityEqualKey> edges;
   for (mEntity* e : domain_for_integral)
   {
      if (debug)
      {
         cout << "treating element in setAllSlocLevelSet " << endl;
         e->print();
      }
      edges.clear();
      for (int i = 0; i < e->size(0); ++i)
      {
         mEntity* v = e->get(0, i);
         std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator ite = e0d3d_e1d_dof.find(v);
         assert(ite != e0d3d_e1d_dof.end());
         edges.insert(ite->second);
      }
      assert(edges.size() > 0);
      if (debug)
      {
         cout << " nb edges is  " << edges.size() << endl;
         cout << "edges are " << endl;
         for (mEntity* pe : edges) pe->print();
      }

      if (edges.size() == 1)
      {
      }
      else if (edges.size() == 2)
      {
         mEdge* ed1 = (mEdge*)*edges.begin();
         mEdge* ed2 = (mEdge*)*(++edges.begin());
         orderEdges(ed1, ed2);
         // retrieve the correct level set.
         mVertex* e0d_dof = ed1->commonVertex(ed2);
         xSubMesh& domain_for_sloc_level_set_e0d_dof =
             mesh->createSubMesh("domain_for_sloc_level_set_" + int2string(e0d_dof->getId()));
         // int id_dof = e0d_dof->getId();
         std::unordered_map<mEntity*, xLevelSet*, EntityHashKey, EntityEqualKey>::const_iterator ite =
             e0d_dof_sloc_level_set.find(e0d_dof);
         xLevelSet* sloc_level_set;
         if (ite == e0d_dof_sloc_level_set.end())
         {
            if (debug) cout << " adding a new level set for dof " << e0d_dof->getId() << endl;
            // to make things easy all level set are defined in teh "boudin" and are
            // set to zero intially

            xRegion reg(mesh, "domain_for_sloc_level_set_" + int2string(e0d_dof->getId()));
            sloc_level_set = new xLevelSet(reg);
            e0d_dof_sloc_level_set.insert(make_pair(e0d_dof, sloc_level_set));
         }
         else
         {
            sloc_level_set = ite->second;
         }
         // add element in level set region
         if (debug) cout << " adding an element in region for level set for dof " << e0d_dof->getId() << endl;
         domain_for_sloc_level_set_e0d_dof.add(e);

         for (int i = 0; i < 3; ++i)
         {
            for (int j = 0; j < e->size(i); ++j)
            {
               mEntity* vv = e->get(i, j);
               domain_for_sloc_level_set_e0d_dof.add(vv);
            }
         }

         // set level set value
         for (int i = 0; i < e->size(0); ++i)
         {
            mEntity* v = e->get(0, i);
            std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator ite = e0d3d_e1d_dof.find(v);
            mEntity* ed = ite->second;
            std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator its = e0d3d_s_loc.find(v);
            double s_loc = its->second;
            double val;
            if (ed == ed1)
            {
               std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator itl = e1d_dof_length.find(ed);
               double length = itl->second;
               val = s_loc - length;
            }
            else
            {
               val = s_loc;
            }
            if (debug)
            {
               cout << " for node  " << v->getId() << endl;
               mVertex* vv = (mVertex*)v;
               cout << " at position " << vv->point() << endl;
               cout << " sloc  is " << s_loc << endl;
               cout << " level set  is " << val << endl;
            }
            (*sloc_level_set)(v) = val;
         }
      }
      else if (edges.size() > 2)
      {
         cerr << " edges.size() is " << edges.size() << " and thus > 2 : not coded yet but could be coded " << endl;
         assert(0);
         abort();
      }
   }

   if (debug)
   {
      cout << " the number of activated level sets is " << e0d_dof_sloc_level_set.size() << endl;
   }
}

void lCrack::orderEdges(mEdge*& ed1, mEdge*& ed2)
{
   const bool debug = false;
   mEntity* v0_ed1 = ed1->get(0, 0);
   mEntity* v1_ed1 = ed1->get(0, 1);
   mEntity* v0_ed2 = ed2->get(0, 0);
   mEntity* v1_ed2 = ed2->get(0, 1);
   if (v1_ed1 == v0_ed2)
      return;
   else if (v0_ed1 == v1_ed2)
   {
      mEdge* tmp = ed1;
      ed1 = ed2;
      ed2 = tmp;
   }
   else
   {
      cerr << "two non consecutive edges, too coarse mesh to compute sifs" << endl;
      cerr << "the edges are " << endl;
      if (debug) ed1->print();
      if (debug) ed2->print();
      cerr << "increase size of nb_internodes " << endl;
      assert(0);
      abort();
   }
}

void lCrack::set3DElementsToFrontRelation(
    const xRegion& domain_for_integral, const std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e0d_s_loc,
    const std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e1d_dof_length,

    const std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>& e0d_e1d_dof,
    std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e0d3d_s_loc,
    std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>& e0d3d_e1d_dof) const
{
   const bool debug = false;
   if (debug) cout << "inner 3D nodes informations " << endl;
   for (mEntity* pe : domain_for_integral.range(0))
   {
      mVertex* v = static_cast<mVertex*>(pe);
      mEntity* e1d;
      double percent_position;
      getClosestSegmentFromTo(v, mesh_crack_front, e1d, percent_position);
      mEntity* v1 = e1d->get(0, 0);
      mEntity* v2 = e1d->get(0, 1);

      if (debug)
      {
         cout << " for node " << endl;
         v->print();
         cout << " closest edge is " << endl;
         e1d->print();
         cout << " with nodes " << endl;
         v1->print();
         v2->print();
      }

      std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator it_v1 = e0d_s_loc.find(v1);
      assert(it_v1 != e0d_s_loc.end());
      double s_loc1 = it_v1->second;

      std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator it_v2 = e0d_s_loc.find(v2);
      assert(it_v2 != e0d_s_loc.end());
      double s_loc2 = it_v2->second;

      double s_loc;

      std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator ite_v1 = e0d_e1d_dof.find(v1);
      assert(ite_v1 != e0d_e1d_dof.end());
      mEntity* e1d_dof1 = ite_v1->second;

      std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>::const_iterator ite_v2 = e0d_e1d_dof.find(v2);
      assert(ite_v2 != e0d_e1d_dof.end());
      mEntity* e1d_dof2 = ite_v2->second;

      mEntity* e1d_dof;

      if ((s_loc1 == 0.) && (e1d_dof1 != e1d_dof2))
      {
         std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator it = e1d_dof_length.find(e1d_dof2);
         s_loc1 = it->second;
         e1d_dof = e1d_dof2;
      }
      else if ((s_loc2 == 0.) && (e1d_dof1 != e1d_dof2))
      {
         std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::const_iterator it = e1d_dof_length.find(e1d_dof1);
         s_loc2 = it->second;
         e1d_dof = e1d_dof1;
      }
      else
      {
         if (e1d_dof1 != e1d_dof2)
         {
            cerr << " information in closest edge for the node e0d " << v->getId() << endl;
            cerr << " located at " << v->point() << endl;
            cerr << " e1d is      [ " << e1d->get(0, 0)->getId() << " " << e1d->get(0, 1)->getId() << " ] " << endl;

            cerr << " s_loc1 " << s_loc1 << endl;
            cerr << " e1d_dof1 is [ " << e1d_dof1->get(0, 0)->getId() << " " << e1d_dof1->get(0, 1)->getId() << " ] " << endl;

            cerr << " s_loc2 " << s_loc2 << endl;
            cerr << " e1d_dof2 is [ " << e1d_dof2->get(0, 0)->getId() << " " << e1d_dof2->get(0, 1)->getId() << " ] " << endl;

            assert(0);
            abort();
         }
         e1d_dof = e1d_dof1;
      }

      s_loc = (1. - percent_position) * s_loc1 + percent_position * s_loc2;
      e0d3d_s_loc[v] = s_loc;
      e0d3d_e1d_dof[v] = e1d_dof;
      if (debug)
      {
         cout << "node " << v->getId() << " has s_loc " << s_loc << endl;
         cout << "node " << v->getId() << " is related to e1d_dof  [ " << e1d_dof->get(0, 0)->getId() << " "
              << e1d_dof->get(0, 1)->getId() << " ] " << endl;
      }
   }
}

void lCrack::extractMeshCrackFrontDofsForJint3DParks(
    xMesh* mesh_crack_front_dofs, std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e0d_s_loc,
    std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>& e1d_dof_length,
    std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>& e0d_dof_e0d,
    std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey>& e0d_e1d_dof, int nb_internodes) const
{
   AOMD::mMesh* mmesh_crack_front_dofs = &mesh_crack_front_dofs->getMesh();
   // code valable seulement pour un morceaui de fissure
   const bool debug = false;
   xIter iter = mesh_crack_front->begin(0);
   mVertex* e0d_start = (mVertex*)*iter;
   for (mEntity* pe : mesh_crack_front->range(0))
   {
      if (pe->size(1) == 1) e0d_start = static_cast<mVertex*>(pe);
   }

   if (debug) cout << endl << "e0d_start" << endl;
   if (debug) e0d_start->print();

   const int group_id = 1000;

   mVertex* start_dof = mmesh_crack_front_dofs->createVertex(e0d_start->point(), mmesh_crack_front_dofs->getGEntity(group_id, 0));

   e0d_dof_e0d.insert(make_pair(start_dof, e0d_start));

   mVertex* left_dof = start_dof;
   std::vector<mVertex*> vec_for_CL, vec_for_test, vec_for_map;
   std::vector<double> abs_loc_final, abs_for_CL1, abs_for_test;

   vec_for_CL.push_back(e0d_start);

   mEntity* e1d_start = e0d_start->get(1, 0);
   mVertex* v_current = e0d_start;
   mVertex* v_next;
   mEntity* ed_current = e1d_start;
   mEntity* ed_next;
   double abs_loc_CL = 0.0;
   double abs_loc_CL2 = 0.0;
   double abs_loc_test = 0.0;
   // double length_e1d_dof_CL = 0.0 ;
   double length_e1d_dof_test = 0.0;

   abs_for_CL1.push_back(0.);
   if (debug) cout << "debut premier remplissage" << endl;
   for (int i = 0; i < nb_internodes; i++)  // verifier les indices !!!!
   {
      v_next = (mVertex*)ed_current->get(0, 0);
      if (v_next == v_current)
      {
         v_next = (mVertex*)ed_current->get(0, 1);
      }
      if ((v_next == e0d_start) || (v_next->size(1) == 1))
      {
         cout << "not enough elements on front for this internode" << endl;
         assert(0);
         abort();
         return;
      }
      xtensor::xVector<> xvec(((mVertex*)ed_current->get(0, 0))->point(), ((mVertex*)ed_current->get(0, 1))->point());
      abs_loc_CL += xvec.mag();
      ed_next = v_next->get(1, 0);
      if (ed_next == ed_current)
      {
         ed_next = v_next->get(1, 1);
      }
      abs_for_CL1.push_back(abs_loc_CL);
      vec_for_CL.push_back(v_next);
      v_current = v_next;
      ed_current = ed_next;
   }
   for (int i = 0; i < nb_internodes; i++)
   {
      if (debug)
      {
         cout << "in e0d_s_loc.insert(make_pair(vec_for_CL[i], abs_for_CL1[i]));" << endl;
         vec_for_CL[i]->print();
      }
      e0d_s_loc.insert(make_pair(vec_for_CL[i], abs_for_CL1[i]));
   }
   length_e1d_dof_test = abs_loc_CL;
   if (debug) cout << "fin premier remplissage length_e1d_dof_test is " << length_e1d_dof_test << endl;

   bool test_sortie = false;
   while (!test_sortie)
   {
      abs_loc_CL2 = abs_loc_test;  // la taille de l avant dernier pour concatenation eventuelle
      abs_loc_test = 0.;
      vec_for_test.clear();
      abs_loc_final.clear();
      abs_for_test.clear();
      abs_for_test.push_back(abs_loc_test);
      vec_for_test.push_back(v_next);          // ou v_current
      for (int i = 0; i < nb_internodes; i++)  // verifier les indices !!!!
      {
         if (!test_sortie)
         {
            v_next = (mVertex*)ed_current->get(0, 0);
            if (v_next == v_current)
            {
               v_next = (mVertex*)ed_current->get(0, 1);
            }
         }
         if (test_sortie)
         {
            cout << "exiting mesh dof generation" << endl;
         }
         else if (((v_next == e0d_start) || (v_next->size(1) == 1)))
         {
            xtensor::xVector<> xvec(((mVertex*)ed_current->get(0, 0))->point(), ((mVertex*)ed_current->get(0, 1))->point());
            abs_loc_test += xvec.mag();
            abs_for_test.push_back(abs_loc_test);
            vec_for_test.push_back(v_next);

            if (vec_for_test.size() < vec_for_CL.size())
            {
               // attention abs_for_cl est peut etre plus court
               for (size_t ii = 0; ii < vec_for_test.size(); ii++)  // on ne prend pas le premier
               {
                  if (debug) cout << endl << endl << "vec_for_test.size() = " << vec_for_test.size() << endl;
                  if (ii > 0)
                  {  // en effet, on ne concatene pas le premier qui est identique
                     vec_for_CL.push_back(vec_for_test[ii]);
                     abs_for_CL1.push_back(abs_for_test[ii] + abs_loc_CL2);
                  }
                  for (size_t i = 0; i < vec_for_CL.size() - 1; i++)
                  {
                     if (debug)
                     {
                        cout << endl << "making pair : " << abs_for_CL1[i] << endl;
                        if (debug) vec_for_CL[i]->print();
                     }
                     e0d_s_loc.insert(make_pair(vec_for_CL[i], abs_for_CL1[i]));
                  }
                  if (v_next->size(1) == 1)
                  {
                     if (debug)
                     {
                        cout << "last node" << endl;
                        v_next->print();
                     }
                     e0d_s_loc.insert(make_pair(v_next, abs_loc_CL2 + abs_loc_test));
                     if (debug) cout << endl << "making pair : " << abs_loc_CL2 + abs_loc_test << endl << endl;
                  }
                  length_e1d_dof_test = abs_loc_CL2 + abs_loc_test;
               }
               mVertex* right_dof2;
               mEdge* e1d_dof;
               if (v_next->size(1) == 1)
               {
                  right_dof2 =
                      mmesh_crack_front_dofs->createVertex(v_next->point(), mmesh_crack_front_dofs->getGEntity(group_id, 0));

                  e0d_dof_e0d.insert(make_pair(right_dof2, v_next));
                  e1d_dof =
                      mmesh_crack_front_dofs->createEdge(left_dof, right_dof2, mmesh_crack_front_dofs->getGEntity(group_id, 1));
                  mEntity* eee1d_dof = (mEntity*)e1d_dof;
                  e0d_e1d_dof.insert(make_pair(v_next, eee1d_dof));
               }
               else
               {
                  right_dof2 = start_dof;
                  e1d_dof =
                      mmesh_crack_front_dofs->createEdge(left_dof, right_dof2, mmesh_crack_front_dofs->getGEntity(group_id, 1));
               }
               mEntity* ee1d_dof = (mEntity*)e1d_dof;
               if (debug)
               {
                  ee1d_dof->print();
                  left_dof->print();
                  right_dof2->print();
               }
               e1d_dof_length.insert(make_pair(ee1d_dof, length_e1d_dof_test));
               left_dof = right_dof2;  // utile?? non je ne crois pas
               for (size_t k = 0; k < (vec_for_CL.size() - 1); k++)
               {
                  e0d_e1d_dof.insert(make_pair(vec_for_CL[k], ee1d_dof));
               }
            }
            else  //(vec_for_test.size()=vec_for_CL.size())
            {
               for (size_t i = 0; i < vec_for_CL.size() - 1; i++)
               {
                  e0d_s_loc.insert(make_pair(vec_for_CL[i], abs_for_CL1[i]));
                  e0d_s_loc.insert(make_pair(vec_for_test[i], abs_for_test[i]));
               }
               if (v_next->size(1) == 1) e0d_s_loc.insert(make_pair(v_next, abs_loc_test));
               mVertex* right_dof;
               mEdge* e1d_dof;
               right_dof = mmesh_crack_front_dofs->createVertex(vec_for_test[0]->point(),
                                                                mmesh_crack_front_dofs->getGEntity(group_id, 0));
               e0d_dof_e0d.insert(make_pair(right_dof, vec_for_test[0]));
               e1d_dof = mmesh_crack_front_dofs->createEdge(left_dof, right_dof, mmesh_crack_front_dofs->getGEntity(group_id, 1));
               mEntity* ee1d_dof = (mEntity*)e1d_dof;
               if (debug)
               {
                  ee1d_dof->print();
                  left_dof->print();
                  right_dof->print();
               }
               e1d_dof_length.insert(make_pair(ee1d_dof, abs_for_CL1[nb_internodes]));
               for (size_t k = 0; k < (vec_for_CL.size() - 1); k++)
               {
                  e0d_e1d_dof.insert(make_pair(vec_for_CL[k], ee1d_dof));
               }
               left_dof = right_dof;  // on fait le e1d_dof sur vec_for_test
               if (v_next->size(1) == 1)
               {
                  right_dof =
                      mmesh_crack_front_dofs->createVertex(v_next->point(), mmesh_crack_front_dofs->getGEntity(group_id, 0));
                  e0d_dof_e0d.insert(make_pair(right_dof, v_next));
               }
               else
               {
                  right_dof = start_dof;
               }
               e1d_dof = mmesh_crack_front_dofs->createEdge(left_dof, right_dof, mmesh_crack_front_dofs->getGEntity(group_id, 1));
               left_dof = right_dof;  // utile?? non je ne crois pas
               ee1d_dof = (mEntity*)e1d_dof;
               if (debug)
               {
                  ee1d_dof->print();
                  left_dof->print();
                  right_dof->print();
               }
               e1d_dof_length.insert(make_pair(ee1d_dof, abs_loc_test));
               for (size_t k = 0; k < (vec_for_test.size() - 1); k++)
               {
                  e0d_e1d_dof.insert(make_pair(vec_for_test[k], ee1d_dof));
               }
               if (v_next->size(1) == 1)
               {
                  e0d_e1d_dof.insert(make_pair(v_next, ee1d_dof));
               }
            }
            test_sortie = true;
         }
         else  // on n est pas dans le cas de sortie
         {
            if (debug) cout << "on n est pas dans le cas de sortie" << endl;
            xtensor::xVector<> xvec(((mVertex*)ed_current->get(0, 0))->point(), ((mVertex*)ed_current->get(0, 1))->point());
            if (debug) cout << " xvec mag  " << xvec.mag() << endl;
            abs_loc_test += xvec.mag();
            if (debug) cout << " abs_loc_test  " << abs_loc_test << endl;
            abs_for_test.push_back(abs_loc_test);
            vec_for_test.push_back(v_next);
            if (vec_for_test.size() < vec_for_CL.size())
            {
               if (debug) cout << "vec_for_test.size() < vec_for_CL.size()" << endl;
            }

            else  //(vec_for_test.size() = vec_for_CL.size())
            {
               if (debug) cout << " on n est pas dans le cas de sortie, vec_for_test.size() = vec_for_CL.size()" << endl;
               // test
               length_e1d_dof_test = abs_for_CL1[nb_internodes];

               vec_for_map = vec_for_CL;
               vec_for_CL = vec_for_test;
               abs_for_CL1 = abs_for_test;
               vec_for_test.clear();
               abs_for_test.clear();

               if (debug)
               {
                  for (int i = 0; i < (nb_internodes); i++)
                  {
                     cout << "making pair in standard : " << abs_for_CL1[i] << endl;
                     vec_for_CL[i]->print();
                  }
               }

               for (int i = 0; i < (nb_internodes); i++)
               {
                  e0d_s_loc.insert(make_pair(vec_for_CL[i], abs_for_CL1[i]));
               }

               abs_for_test.push_back(0.);
               vec_for_test.push_back(vec_for_test[0]);
               mVertex* right_dof = mmesh_crack_front_dofs->createVertex(vec_for_test[0]->point(),
                                                                         mmesh_crack_front_dofs->getGEntity(group_id, 0));
               e0d_dof_e0d.insert(make_pair(right_dof, vec_for_test[0]));
               mEdge* e1d_dof =
                   mmesh_crack_front_dofs->createEdge(left_dof, right_dof, mmesh_crack_front_dofs->getGEntity(group_id, 1));
               mEntity* ee1d_dof = (mEntity*)e1d_dof;

               if (debug)
               {
                  ee1d_dof->print();
                  left_dof->print();
                  right_dof->print();
               }
               if (debug)
               {
                  cout << " on est pas dans le cas sortie " << endl;
                  cout << " e1d_dof [ " << ee1d_dof->get(0, 0)->getId() << " " << ee1d_dof->get(0, 1)->getId()
                       << " ] entre avec length " << abs_for_CL1[nb_internodes] << endl;
               }

               //		  e1d_dof_length.insert(make_pair(ee1d_dof, abs_for_CL1[nb_internodes]));
               // test
               e1d_dof_length.insert(make_pair(ee1d_dof, length_e1d_dof_test));

               if (debug)
               {
                  for (size_t k = 0; k < (vec_for_map.size() - 1); k++)
                  {
                     cout << "e0d_e1d_dof.insert(make_pair(vec_for_CL[k],ee1d_dof));" << endl;
                     vec_for_map[k]->print();
                     ee1d_dof->print();
                  }
               }

               for (size_t k = 0; k < (vec_for_map.size() - 1); k++)
               {
                  e0d_e1d_dof.insert(make_pair(vec_for_map[k], ee1d_dof));
               }

               left_dof = right_dof;
            }
         }
         ed_next = v_next->get(1, 0);
         if (ed_next == ed_current)
         {
            ed_next = v_next->get(1, 1);
         }
         v_current = v_next;
         ed_current = ed_next;  // utile?
      }

      if (debug)
      {
         if (v_current->getId() == 22)
         {
            cout << "having treated node 22" << endl;
            cout << "size(1) = " << v_current->size(1) << endl;
         }
      }
   }

   if (debug)
   {
      cout << "  node_e0d  s_loc  " << endl;
      for (std::pair<mEntity*, double> e_v : e0d_s_loc)
      {
         cout << e_v.first->getId() << " " << e_v.second << endl;
      }
      cout << "  e1d_dof  length  " << endl;
      for (std::pair<mEntity*, double> e_v : e1d_dof_length)
      {
         cout << " [ " << e_v.first->get(0, 0)->getId() << " " << e_v.first->get(0, 1)->getId() << " ] " << e_v.second << endl;
      }
      cout << "  node_e0d_dof  node_0d  " << endl;
      for (std::pair<mEntity*, mEntity*> e_e : e0d_dof_e0d)
      {
         cout << e_e.first->getId() << " " << e_e.second->getId() << endl;
      }
      cout << "  node_e0d  e1d_dof  " << endl;
      for (std::pair<mEntity*, mEntity*> e_e : e0d_e1d_dof)
      {
         cout << e_e.first->getId() << " [ " << e_e.second->get(0, 0)->getId() << " " << e_e.second->get(0, 1)->getId() << " ] "
              << endl;
      }
   }

   AOMD::classifyUnclassifiedVerices(mmesh_crack_front_dofs);
   mmesh_crack_front_dofs->modifyState(0, 1, true);

   if (debug)
   {
      cout << " mesh_crack_front_dofs infos " << endl;
      mesh_crack_front_dofs->getMesh().printAll();
   }
}

void xcSetParksRadialFunction::visit(xLevelSet& f, xRegion target)
{
   for (mEntity* e : target.range(0))
   {
      double n = lsn(e);
      double t = lst(e);
      double r = sqrt(n * n + t * t);
      if (r >= rho)
         f(e) = 0.0;
      else
         f(e) = 1.0;
      //      else f(e) = exp(-(r/(r-rho))*(r/(r-rho)));
   }
   return;
}

void lCrack::setQDirectionParks(xVectorField& q_dir) const
{
   const bool debug = false;
   const xRegion& target = q_dir.getSupport();
   for (mEntity* v : target.range(0))
   {
      q_dir(v) = lst.getGrad(v);
      if (debug) cout << " norm of the direction vector at point " << v->getId() << " : " << q_dir(v).mag() << endl;
   }
   // correction sur les bords

   //   for(xIter it = target.begin(2); it != target.end(2); ++it)
   //     {
   //       mEntity *f =  *it;
   //       assert(f->size(3) == 1 || f->size(3) == 2);
   //       if (f->size(3) == 1)
   // 	{
   // 	  if (debug) cout << " boundary face detected " << endl;
   // 	  mEntity* e = f->get(3,0);
   // 	  xtensor::xVector<> normal;
   // 	  xfem::xGeomElem geo(e);
   // 	  geo.normalVector(f, normal);
   // 	  for (int i = 0; i < f->size(0); ++i)
   // 	    {
   // 	      mEntity* v   = f->get(0,i);
   // 	      if (debug) cout << " q_dir before projection  " << q_dir(v) << endl;
   // 	      double proj  = q_dir(v) * normal;
   // 	      xtensor::xVector<> red  = normal * proj;
   // 	      q_dir(v) -= red;
   // 	      if (debug) cout << " q_dir  after projection  " << q_dir(v) << endl;
   // 	    }

   // 	}
   //     }
}

}  // namespace xcrack
