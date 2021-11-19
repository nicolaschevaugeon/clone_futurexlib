/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#include <fstream>
#include <cassert>

#include "main.h"
#include "xData.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"
#include "xLinearSystem.h"
#include "xLinearSystemSolverLU.h"
#include "xForm.h"
#include "xAssembler.h"
#include "xGmshAsciiExport.h"
#include "xEnv.h"
#include "xSpace.h"
#include "xGeomElem.h"
#include "xFemMatrix.h"
#include "xMaterialSensitivity.h"
#include "xIntegrationRule.h"
#include "xValue.h"
#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
#include "xField.h"
#include "xApproxFunction.h"
#include "xOperators.h"
#include "xSolVisitor.h"
#include "xMaterial.h"


using namespace xfem;


Mechanics_c::Mechanics_c() 
{ 
  xMaterialManagerSingleton::instance().registerMaterial("elastic", xMaterialCreator<xElastic>() );
}

Mechanics_c::~Mechanics_c() 
{ }


static int RECOVER_RESULTS = 0, FLAG_NAMEPROBLEM = 0, SAMPLE_WRITE = 10;
char mess[256];
int Get_Options(int narg, char  * arg[], char  * NameProblem) {
	
  int  i;
  
  strcpy(NameProblem, "") ; i = 1 ;
  
  while (i < narg){
    
    if (arg[i][0] == '-') {
      if      (!strcmp(arg[i]+1, "recover"  ))  { RECOVER_RESULTS = 1 ; i++ ; }
      else if (!strcmp(arg[i]+1, "sample"  ))   { i++ ; SAMPLE_WRITE = atoi(arg[i++]) ; }
      else {
	sprintf(mess, "Warning : Unknown option %s", arg[i++]) ;
      }		
    }
    else {
      FLAG_NAMEPROBLEM = 1 ; 
      strcpy(NameProblem, arg[i++]) ;
    }
    
  }
  return 0 ;
}

int main(int argc, char *argv[])  
{  

  xData data;
  char pname[256];
  Get_Options(argc,argv,pname);
  
  fprintf(stderr, "Starting reading the master data file\n");
  data.ReadInfo(pname);
  fprintf(stderr, "Done with reading the master data file\n");
  
  Mechanics_c Mechanics;
  Mechanics.TreatmentOfFormulation(&data);
  return 0; 
  
}


void Mechanics_c :: TreatmentOfFormulation (xData *data) {


  fprintf(stderr, "Starting reading the mesh file\n");
  data->ReadMesh();
  fprintf(stderr, "End      reading the mesh file\n");
  data->ReadZones();

  xRegion all(data->mesh);


#if 1
  xcut::xPhysSurf inclusion(data->mesh, 
		      xCylinder(Trellis_Util::mPoint(0.,0.,0.), xtensor::xVector(0.,0.,1.), 0.4),
		      xClassifyOn("inclusion", data->mesh), 
		      xClassifyOn("matrix", data->mesh));
  
#endif


  xSpaceLagrange lagx("DISPLACEMENT_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
  xSpaceLagrange lagy("DISPLACEMENT_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
  xSpaceLagrange lagz("DISPLACEMENT_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_ONE);
  xSpaceComposite  lagrange(lagx, lagy, lagz);


#if 1
  xValKeyExtend key_modifier("_INCLUSION");  
  xScalarFunctionDerivDiscXFEM enrichment_function(inclusion);
  xSpaceFiltered::filter_t filter(bind1st(mem_fun(&xcut::xPhysSurf::boundary_strict), &inclusion));
  xSpaceXFEM space_full(lagrange, enrichment_function, key_modifier);
  xSpaceFiltered enriched(space_full, filter);
  xField disp_l(&double_manager, lagrange, enriched);
#endif



  xValueCreator<xValueDouble>  creator;
  DeclareInterpolation(disp_l, creator, all.begin(), all.end());
  double_manager.PrintForDebug("dcl.dbg");


  printf("Hello 0 in thermic\n"); 

  TreatmentOfEssEnv(disp_l, data);
  double_manager.PrintForDebug("ess.dbg");

  printf("Hello 0.1 in thermic\n"); 

  xStateDofCreator<> snh(double_manager, "dofs");
  DeclareState(disp_l, snh, all.begin(), all.end());
//cout << "nbr dof after nothanging " << double_manager.size("dofs") << endl;
 

  double_manager.PrintForDebug("sym.dbg");
printf("Hello 2 in thermic nb dof is %d\n", double_manager.size("dofs"));
  xlinalg::xCSRVector b(double_manager.size("dofs"));
  xlinalg::xCSRVector sol(double_manager.size("dofs"));
  xlinalg::xCSRMatrix A(double_manager.size("dofs"));


printf("Hello 2.01 in thermic\n");
  xlinalg::xLinearSystemSolverLU solver;
printf("Hello 2.02 in thermic\n");
  xLinearSystem system(&A, &solver); // 

printf("Hello 2.03 in thermic\n");
  xAssemblerBasic<> assembler(A, b);
  xIntegrationRulePartition integration_rule_env(3);
  xIntegrationRulePartition integration_rule(3); //2 is the degree to integrate

printf("Hello 2.05 in thermic\n");
  TreatmentOfNatEnv(disp_l, assembler, integration_rule_env, data, data->allGroups);

printf("Hello 2.06 in thermic\n");
  xUniformMaterialSensitivity<xtensor::xTensor4> hooke("strain");
printf("Hello 2.06 in thermic\n");
  xFormBilinearWithLaw<xGradOperator<xtool::xIdentity<xtensor::xTensor2> >, 
                       xUniformMaterialSensitivity<xtensor::xTensor4>, 
                       xGradOperator<xtool::xIdentity<xtensor::xTensor2> > > diffusive(hooke);
  //DiffusiveVectorxFormBilinear diffusive(hooke);
printf("Hello 2.07 in thermic\n");
  Assemble(diffusive, assembler, integration_rule, disp_l, disp_l, all.begin(), all.end()); 
printf("Hello 2.08 in thermic\n");

  printf("Hello 2.3 in thermic\n");
//std::ofstream out("mat.m");
//A.OutputMatrixOctaveFormat(out);
//return;

  //system.Solve(double_manager, "dofs");
  system.Solve(b, sol);
  Visit(xWriteSolutionVisitor(sol.begin()), 
	double_manager.begin("dofs"), 
	double_manager.end("dofs"));


  //printf("Hello 2.6 in thermic\n");

  //DofManager.StoreResult(SOLUTION.GetArray());

  double_manager.PrintForDebug("res.dbg");

  printf("Hello 3 in thermic\n");

  xGmshAsciiExport pexport;

  printf("Hello 5 in thermic\n");
  xEvalField<xtool::xIdentity<xtensor::xVector> > eval_disp(disp_l);
  pexport.Export(eval_disp, "DISPLACEMENT", integration_rule, all.begin(), all.end());


//  xEvalGradField<xtensor::xSymmetrize> eval_strain(disp_l);
//  xUniformMaterialSensitivity<xtensor::xTensor4> hooke2("strain");
//  xEvalBinary< xtool::xMult<xtensor::xTensor4Isotropic, xtensor::xTensor2, xtensor::xTensor2> > stress(hooke, eval_strain);
//  pexport.Export(stress, "STRESS", integration_rule, all.begin(), all.end());


#if 0
///////////////////////////////////////////
//// Estimated Error computation ////////////////////
///////////////////////////////////////////
  xSpaceLagrange lagerr2x("ERROR_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_TWO); 
  xSpaceLagrange lagerr2y("ERROR_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_TWO);
  xSpaceLagrange lagerr2z("ERROR_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_TWO);
  xSpaceComposite lagerr2(lagerr2x, lagerr2y, lagerr2z); 

  xSpaceLagrange lagerr1x("ERROR_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE); 
  xSpaceLagrange lagerr1y("ERROR_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
  xSpaceLagrange lagerr1z("ERROR_Z", xSpace::VECTOR_Z, xSpaceLagrange::DEGREE_ONE);
  xSpaceComposite lagerr1(lagerr1x, lagerr1y, lagerr1z); 

  xSpaceDifference lagerr(lagerr2, lagerr1);

  xField err_l(&double_manager, lagerr);

  DeclareInterpolation(err_l, creator, all.begin(), all.end());
  double_manager.PrintForDebug("dcl_err.dbg");

  TreatmentOfEssEnv(err_l, data);  
  double_manager.PrintForDebug("ess_err.dbg");

  printf("Hello 3... in thermic\n"); 

  //declare temperature as fixed now
  printf("Hello 3.001 in thermic\n"); 
  DeleteState(disp_l, all.begin(), all.end());
  printf("Hello 3.002 in thermic\n"); 
  DeclareState(disp_l, xStateFixedCreator(), all.begin(), all.end());
  printf("Hello 3.003 in thermic_1\n"); 

 
  xStateDofCreator<> snh_err(double_manager, "dofs_err");
  DeclareState(err_l, snh_err, all.begin(), all.end());

  double_manager.PrintForDebug("sym.dbg");
printf("Hello 3 in thermic nb dof is %d\n", double_manager.size("dofs_err"));
  xlinalg::xCSRVector b_err(double_manager.size("dofs_err"));
  xlinalg::xCSRVector sol_err(double_manager.size("dofs_err"));
  xlinalg::xCSRMatrix A_err(double_manager.size("dofs_err"));

  xIntegrationRuleBasic integration_rule_env_err(3);
  xIntegrationRuleBasic integration_rule_err(3); //2 is the degree to integrate

printf("Hello 3.01 in thermic\n");
  xlinalg::xLinearSystemSolverLU solver_err;
printf("Hello 3.02 in thermic\n");
  xLinearSystem system_err(&A_err, &solver_err);

printf("Hello 3.03 in thermic\n");
  xAssemblerBasic<> assembler_err(A_err, b_err);


printf("Hello 3.05 in thermic\n");
  TreatmentOfNatEnv(err_l, assembler_err, integration_rule_env_err, data, data->allGroups);

printf("Hello 3.07 in thermic\n");
  Assemble(diffusive, assembler_err, integration_rule_err, err_l, err_l,  all.begin(), all.end()); 
printf("Hello 3.075 in thermic\n");
  Assemble(diffusive, assembler_err, integration_rule_err, err_l, disp_l, all.begin(), all.end()); 
printf("Hello 3.08 in thermic\n");

  //system_err.Solve(double_manager, "dofs_err");
  system_err.Solve(b_err, sol_err);
  Visit(xWriteSolutionVisitor(sol_err.begin()), 
	double_manager.begin("dofs_err"), 
	double_manager.end("dofs_err"));


  double_manager.PrintForDebug("res_err.dbg");


  //declare err as fixed now
  printf("Hello 3.001 in thermic\n"); 
  DeleteState(err_l, all.begin(), all.end());
  printf("Hello 3.002 in thermic\n"); 
  DeclareState(err_l, xStateFixedCreator(), all.begin(), all.end());
  printf("Hello 3.003 in thermic_2\n"); 


///////////////////////////////////////////////////////////////
//  Necessary for the error computation
////////////////////////////////////////////////////////////////

  xField err_info_l(&double_manager, xSpaceConstant("ERROR_INFO"));
  xValueCreator<xValueError>  creator_err;
  DeclareInterpolation(err_info_l, creator_err, all.begin(), all.end());

  //DiffusiveVectorZeroForm_c err_form(hooke, err_l);
  xEvalGradField<xtool::xIdentity<xtensor::xTensor2> > grad_err(err_l);
  xFormZeroEvalBilinearFormWithLaw<xEvalGradField<xtool::xIdentity<xtensor::xTensor2> >,
                                    xUniformMaterialSensitivity<xtensor::xTensor4>,
                                    xEvalGradField<xtool::xIdentity<xtensor::xTensor2> > > err_form(grad_err, hooke);

  xValueError::choice("ABS2");
  xFillFieldFromZeroForm fill_err(err_info_l, err_form);
  ApplyCommandOnIntegrationRule(fill_err, integration_rule_err, all.begin(), all.end());


  //DiffusiveVectorZeroForm_c eng_form(hooke, disp_l);
  xEvalGradField<xtool::xIdentity<xtensor::xTensor2> > grad_disp(disp_l);
  xFormZeroEvalBilinearFormWithLaw<xEvalGradField<xtool::xIdentity<xtensor::xTensor2> >,
                                    xUniformMaterialSensitivity<xtensor::xTensor4>,
                                    xEvalGradField<xtool::xIdentity<xtensor::xTensor2> > > eng_form(grad_disp, hooke);

  xValueError::choice("ENG");
  xFillFieldFromZeroForm fill_eng(err_info_l, eng_form);
  ApplyCommandOnIntegrationRule(fill_eng, integration_rule_err, all.begin(), all.end());  


  xFormZeroUnit unit_form;
  xValueError::choice("VOL");
  xFillFieldFromZeroForm fill_vol(err_info_l, unit_form);
  ApplyCommandOnIntegrationRule(fill_vol, integration_rule_err, all.begin(), all.end());
  //ZeroFormVisitor_c visit_vol(unit_form, integration_rule_err, err_info_l);
  //visit(visit_vol, all.begin(), all.end()); 

  std::ofstream out("error.txt");
  out << "energy  is  " << xValueError::total("ENG") << endl;
  out << "err abs is  " << xValueError::total("ABS") << endl;
  out << "err rel is  " << xValueError::total("REL") << endl;
  out << "vol     is  " << xValueError::total("VOL") << endl;

  xValueError::choice("REL");
  xEvalField<xtool::xIdentity<double> > val_err_info(err_info_l);
  pexport.Export(val_err_info, 
		 "ERROR_CONTRIBUTION_ESTIMATED", integration_rule, 
		 all.begin(), all.end());

  xValueError::choice("ENG");
  pexport.Export(val_err_info, 
		 "ENERGY", integration_rule, 
		 all.begin(), all.end());
  
  xValueError::choice("VOL");
  pexport.Export(val_err_info, 
		 "VOLUMES", integration_rule, 
		 all.begin(), all.end());

  xValueError::choice("DNS");
  pexport.Export(val_err_info, 
		 "DNS", integration_rule, 
		 all.begin(), all.end());

  
#endif


  return;
  
}


void Mechanics_c :: TreatmentOfEssEnv (const xField& listFunctionSpace, xData * data){

  for (xPhysicalEnv::const_iterator it = data->PhysicalEnv->begin(); it != data->PhysicalEnv->end(); ++it) {
    const xEnv& env = *it;
    string phys = env.Phys; 
    int type = env.Type; 
    int entity = env.Entity;
    //
    // FEM case
    //
    if (phys == "DISPLACEMENT_X" || phys == "DISPLACEMENT_Y" || phys == "DISPLACEMENT_Z") {       
      if (type == FIX) {
cout << "before Dirichmet entity " << entity << " dim " << env.getDimension() << endl;
        xClassRegion bc(data->mesh, entity, env.getDimension());
	DirichletBoundaryCondition (listFunctionSpace, phys, bc.begin(), bc.end());
      }
      else assert(1 == 0);      
    }
} // End loop over the environemnt info
return;
}



void Mechanics_c :: TreatmentOfNatEnv   (const xField& listFunctionSpace, 
				       xAssembler& assembler, xIntegrationRule& integration_rule,
				       xData * data, xBoundary& groups){
  

  for (xPhysicalEnv::const_iterator it = data->PhysicalEnv->begin(); it != data->PhysicalEnv->end(); ++it) {
    const xEnv& env = *it;
    if (env.Phys == "TRACTION_X" || env.Phys == "TRACTION_Y"  || env.Phys == "TRACTION_Z" ) {
      assert(env.Type == FIX);
      xtensor::xVector val;
      if (env.Phys == "TRACTION_X") val(0) =  env.getValue();
      if (env.Phys == "TRACTION_Y") val(1) =  env.getValue();
      if (env.Phys == "TRACTION_Z") val(2) =  env.getValue();
      cout << "val is " << val(0) << " " << val(1) << " " << val(2) << endl;
      xEvalConstant<xtensor::xVector>  flux(val);
      xFormLinearWithLoad<xValOperator<xtool::xIdentity<xtensor::xVector> >, xEvalConstant<xtensor::xVector> > lin(flux); 
      xClassRegion bc(data->mesh, env.Entity, env.getDimension());
      Assemble(lin, assembler, integration_rule, listFunctionSpace, bc.begin(), bc.end(), 
	       xUpperAdjacency()); 
    }
  } 
  return;
}



