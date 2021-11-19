/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#include "main.h"
#include "Data.h"
#include <cassert>
#include "CSR_Matrix.h"
#include "CSR_Vector.h"
#include "LinearSystem.h"
#include "LinearSystemSolverSparskit.h"

#include "FormsBase.h"
#include "AssemblerBasic.h"
#include "GMSHASCIIExport.h"
#include "PhysicalEnv.h"
#include "SpaceLagrange.h"
#include "SpaceFiltered.h"
#include "SpaceConstant.h"
#include "SpaceComposite.h"
#include "SpaceDifference.h"
#include "GeomElem.h"
#include "FemMatrix.h"
#include "UniformMaterialSensitivity.h"
#include "Load.h"
#include "mext.h"
#include "IntegratorBasic.h"
#include "IntegratorPartition.h"
#include "ValueDouble.h"
#include "ValueError.h"
#include "ValueCreators.h"
#include "StateCreators.h"
#include "lPhysSurf.h"
#include "lSimpleGeometry.h"
#include "mxIterator.h"
#include "Field.h"
#include "FunctionDerivDiscXFEM.h"
#include "SpaceXFEM.h"
#include "Operators.h"
#include "SolVisitor.h"
#include "Conductive.h"

#include <fstream>
using namespace xfem;

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

Mechanics_c::Mechanics_c() 
{ 
  MaterialManagerSingleton_c::Instance().registerMaterial("thermic", CreateMaterial_c<Conductive_c>() ); 

}

Mechanics_c::~Mechanics_c() 
{ 
}



int main(int argc, char *argv[])  
{  
  Data_c data;
  char pname[256];
  Get_Options(argc,argv,pname);
  
  fprintf(stderr, "Starting reading the master data file\n");
  data.ReadInfo(pname);
  fprintf(stderr, "Done with reading the master data file\n");


  
  Mechanics_c Mechanics;
  Mechanics.TreatmentOfFormulation(&data);
  return 0; 
  
}



void Mechanics_c :: TreatmentOfFormulation (Data_c *data) {


  //hole
//    lPhysSurf all(data->mesh, lCompl(lSphere(Trellis_Util::mPoint(0.,0.,0.), 0.7)));
//    SpaceFiltered_c::filter_t filter(&all, &lPhysSurf::covers);
//    Field_c temp_l(&DoubleManager,
//  		 SpaceFiltered_c(
//  		 SpaceLagrange_c("TEMPERATURE", SCALAR, DEGREE_ONE), filter));

  fprintf(stderr, "Starting reading the mesh file\n");
  data->ReadMesh();
  fprintf(stderr, "End      reading the mesh file\n");
  data->ReadZones();

  //Inclusion
  mxRegion all(data->mesh);
//		      lCylinder(Trellis_Util::mPoint(0.,0.,0.), xtensor::xVector(0.,0.,1.), 0.5),
//		      lPlane(Trellis_Util::mPoint(0.,0.,0.), xtensor::xVector(0.,1.,0.)),
//		      lPlane(Trellis_Util::mPoint(0.,0.,0.), xtensor::xVector(0.,1.,0.)),


  lPhysSurf inclusion(data->mesh, 
		      lCylinder(Trellis_Util::mPoint(-1.,-1.,0.), xtensor::xVector(0.,0.,1.), 1.0),
		      mxClassifyOn("inclusion", data->mesh), 
		      mxClassifyOn("matrix", data->mesh));
  SpaceLagrange_c lag("TEMPERATURE", SCALAR, DEGREE_ONE);
  SpaceComposite_c  lagrange(lag);

  ValKeyExtend_c key_modifier("_INCLUSION");  
  ScalarFunctionDerivDiscXFEM_c enrichment_function(inclusion);
  SpaceFiltered_c::filter_t filter(&inclusion, &lPhysSurf::boundary_strict);
  SpaceFiltered_c enriched(SpaceXFEM_c(lagrange, enrichment_function, key_modifier), filter);
  Field_c temp_l(&DoubleManager, lagrange, enriched);

  CreateValue_c<ValueDouble_c>  creator;
  DeclareInterpolation(temp_l, creator, all.begin(), all.end());
  DoubleManager.PrintForDebug("dcl.dbg");


  printf("Hello 0 in inclusion_xfem\n"); 




#if 1
  OuterEval_c dirichlet; 
  mxClassRegion bc110(data->mesh, 110, 1); //surface 110 of dimension 1
  mxClassRegion bc111(data->mesh, 111, 1); //surface 111 of dimension 1
  IntegratorBasic_c integrator_basic(2); 
  LinearSystemSolverSparskit_c solver_l2;
  L2Projection(temp_l, dirichlet, integrator_basic, solver_l2, bc110.begin(), bc110.end()); 
//  LinearSystemSolverSparskit_c solver_l2_bis;
//  L2Projection(temp_l, dirichlet, integrator_basic, solver_l2_bis, bc111.begin(), bc111.end()); 
#endif





  TreatmentOfEssEnv(temp_l, data);
  DoubleManager.PrintForDebug("ess.dbg");

  printf("Hello 0.1 in inclusion_xfem\n"); 

  CreateStateDof_c<> snh(DoubleManager, "dofs");
  DeclareState(temp_l, snh, all.begin(), all.end());
//cout << "nbr dof after nothanging " << DoubleManager.size("dofs") << endl;
 
  DoubleManager.PrintForDebug("sym.dbg");
printf("Hello 2 in inclusion_xfem nb dof is %d\n", DoubleManager.size("dofs"));
  CSR_Vector_c b(DoubleManager.size("dofs"));
  CSR_Vector_c sol(DoubleManager.size("dofs"));
  CSR_Matrix_c A(DoubleManager.size("dofs"));


printf("Hello 2.01 in inclusion_xfem\n");
  LinearSystemSolverSparskit_c solver;
printf("Hello 2.02 in inclusion_xfem\n");
  LinearSystem_c system(&A, &solver);

printf("Hello 2.03 in inclusion_xfem\n");
  AssemblerBasic_c assembler(A, b);
  IntegratorPartition_c integrator_env(3);
  IntegratorPartition_c integrator(3); //3 is the degree to integrate
  //IntegratorBasic_c integrator(2); //2 is the degree to integrate


printf("Hello 2.05 in inclusion_xfem\n");
  TreatmentOfNatEnv(temp_l, assembler, integrator_env, data, data->allGroups);

#if 1
  mxMesh* interface = inclusion.getMesh_bnd();
  ConstantEval_c<double> jump_in_flux(-0.432809);
  LinearFormWithLoad_c<ValOperator_c<xtool::xIdentity<double> >, ConstantEval_c<double> > lin(jump_in_flux); 
  IntegratorBasic_c integrator_bnd(2); //2 is the degree to integrate
  mxRegion bnd(interface);
  Assemble(lin, assembler, integrator_bnd, temp_l, bnd.begin(), bnd.end(), mxUpperCreator()); 
#endif


printf("Hello 2.06 in inclusion_xfem\n");
  UniformMaterialSensitivity_c<xtensor::xTensor2> fourier("temperature_gradient");
printf("Hello 2.06 in inclusion_xfem\n");
  //DiffusiveVectorBilinearForm_c diffusive(hooke);
  BilinearFormWithLaw_c<GradOperator_c<xtool::xIdentity<xtensor::xVector> >, 
                        UniformMaterialSensitivity_c<xtensor::xTensor2>, 
                        GradOperator_c<xtool::xIdentity<xtensor::xVector> > > diffusive(fourier);
printf("Hello 2.07 in inclusion_xfem\n");
  Assemble(diffusive, assembler, integrator, temp_l, temp_l, all.begin(), all.end()); 
printf("Hello 2.08 in inclusion_xfem\n");

  printf("Hello 2.3 in inclusion_xfem\n");
//std::ofstream out("mat.m");
//A.OutputMatrixOctaveFormat(out);
//return;

  //system.Solve(DoubleManager, "dofs");
  system.Solve(b, sol);
  visit(WriteSolutionVisitor_c(sol.begin()), 
	DoubleManager.begin("dofs"), 
	DoubleManager.end("dofs"));


  //printf("Hello 2.6 in inclusion_xfem\n");

  //DofManager.StoreResult(SOLUTION.GetArray());

  DoubleManager.PrintForDebug("res.dbg");

  printf("Hello 3 in inclusion_xfem\n");

  GMSHASCIIExport_c pexport;

  printf("Hello 5 in inclusion_xfem\n");
  GetValField_c<xtool::xIdentity<double> > eval_disp(temp_l);
  pexport.Export(eval_disp, "TEMPERATURE", integrator, all.begin(), all.end());


  GetGradField_c<xtool::xIdentity<xtensor::xVector> >  eval_grad(temp_l);
  EvalBinary_c< xmult<xtensor::xTensor2, xtensor::xVector, xtensor::xVector> > flux(fourier, eval_grad);
  pexport.Export(flux, "FLUX", integrator, all.begin(), all.end());



#if 0
///////////////////////////////////////////
//// Estimated Error computation ////////////////////
///////////////////////////////////////////
//  Field_c err_l(&DoubleManager,
//      SpaceFiltered_c(SpaceLagrange_c("ERROR", SCALAR, DEGREE_TWO_MINUS_ONE), filter));
  SpaceLagrange_c lag2("ERROR", SCALAR, DEGREE_TWO); 
  SpaceComposite_c lagerr2(lag2); 

  SpaceLagrange_c lag1("ERROR", SCALAR, DEGREE_ONE); 
  SpaceComposite_c lagerr1(lag1); 

  SpaceDifference_c lagerr(lagerr2, lagerr1);

  SpaceFiltered_c enriched_err(SpaceXFEM_c(lagerr, enrichment_function, key_modifier), filter);
  Field_c err_l(&DoubleManager, lagerr, enriched_err);

  DeclareInterpolation(err_l, creator, all.begin(), all.end());
  DoubleManager.PrintForDebug("dcl_err.dbg");

  TreatmentOfEssEnv(err_l, data);  
  DoubleManager.PrintForDebug("ess_err.dbg");

  printf("Hello 3... in inclusion_xfem\n"); 

  //declare temperature as fixed now
  printf("Hello 3.001 in inclusion_xfem\n"); 
  DeleteState(temp_l, all.begin(), all.end());
  printf("Hello 3.002 in inclusion_xfem\n"); 
  DeclareState(temp_l, CreateStateFixed_c(), all.begin(), all.end());
  printf("Hello 3.003 in inclusion_xfem\n"); 

 
  CreateStateDof_c<> snh_err(DoubleManager, "dofs_err");
  DeclareState(err_l, snh_err, all.begin(), all.end());


  DoubleManager.PrintForDebug("sym.dbg");
printf("Hello 3 in inclusion_xfem nb dof is %d\n", DoubleManager.size("dofs_err"));
  CSR_Vector_c b_err(DoubleManager.size("dofs_err"));
  CSR_Vector_c sol_err(DoubleManager.size("dofs_err"));
  CSR_Matrix_c A_err(DoubleManager.size("dofs_err"));

  IntegratorPartition_c integrator_env_err(4);
  IntegratorPartition_c integrator_err(4); //2 is the degree to integrate

printf("Hello 3.01 in inclusion_xfem\n");
  LinearSystemSolverSparskit_c solver_err;
printf("Hello 3.02 in inclusion_xfem\n");
  LinearSystem_c system_err(&A_err, &solver_err);

printf("Hello 3.03 in inclusion_xfem\n");
  AssemblerBasic_c assembler_err(A_err, b_err);


printf("Hello 3.05 in inclusion_xfem\n");
  TreatmentOfNatEnv(err_l, assembler_err, integrator_env_err, data, data->allGroups);

printf("Hello 3.07 in inclusion_xfem\n");
  Assemble(diffusive, assembler_err, integrator_err, err_l, err_l,  all.begin(), all.end()); 
printf("Hello 3.075 in inclusion_xfem\n");
  Assemble(diffusive, assembler_err, integrator_err, err_l, temp_l, all.begin(), all.end()); 
printf("Hello 3.08 in inclusion_xfem\n");

  //system_err.Solve(DoubleManager, "dofs_err");
  system_err.Solve(b_err, sol_err);
  visit(WriteSolutionVisitor_c(sol_err.begin()), 
	DoubleManager.begin("dofs_err"), 
	DoubleManager.end("dofs_err"));

  DoubleManager.PrintForDebug("res_err.dbg");

  //pexport.ExportXFEM(err_l, "ERROR", data, all.begin(), all.end(), mxAcceptAll());
  GetValField_c<xtool::xIdentity<xtensor::xVector> > val_err(err_l);
  pexport.Export(val_err, "ERROR", integrator, all.begin(), all.end());

  //declare err as fixed now
  printf("Hello 3.001 in inclusion_xfem\n"); 
  DeleteState(err_l, all.begin(), all.end());
  printf("Hello 3.002 in inclusion_xfem\n"); 
  DeclareState(err_l, CreateStateFixed_c(), all.begin(), all.end());
  printf("Hello 3.003 in inclusion_xfem\n"); 


///////////////////////////////////////////////////////////////
//  Necessary for the exact1 and estimated error computation
////////////////////////////////////////////////////////////////
  Field_c err_info_l(&DoubleManager, SpaceConstant_c("ERROR_INFO"));
  CreateValue_c<ValueError_c>  creator_err;
  DeclareInterpolation(err_info_l, creator_err, all.begin(), all.end());


  //DiffusiveVectorZeroForm_c err_form(hooke, err_l);
  GetGradField_c<xtool::xIdentity<xtensor::xTensor2> > grad_err(err_l);
  ZeroFormEvalBilinearFormWithLaw_c<GetGradField_c<xtool::xIdentity<xtensor::xTensor2> >,
                                    UniformMaterialSensitivity_c<xtensor::xTensor4>,
                                    GetGradField_c<xtool::xIdentity<xtensor::xTensor2> > > err_form(grad_err, hooke);


  ValueError_c::choice("ABS2");
  FillFieldFromZeroForm_c fill_err(err_info_l, err_form);
  CommandOnIntegratorPath(fill_err, integrator_err, all.begin(), all.end());
  //ZeroFormVisitor_c visit_err(err_form, integrator_err, err_info_l);
  //visit(visit_err, all.begin(), all.end()); 

  //DiffusiveVectorZeroForm_c eng_form(hooke, temp_l);
  GetGradField_c<xtool::xIdentity<xtensor::xTensor2> > grad_disp(temp_l);
  ZeroFormEvalBilinearFormWithLaw_c<GetGradField_c<xtool::xIdentity<xtensor::xTensor2> >,
                                    UniformMaterialSensitivity_c<xtensor::xTensor4>,
                                    GetGradField_c<xtool::xIdentity<xtensor::xTensor2> > > eng_form(grad_disp, hooke);

  ValueError_c::choice("ENG");
  FillFieldFromZeroForm_c fill_eng(err_info_l, eng_form);
  CommandOnIntegratorPath(fill_eng, integrator_err, all.begin(), all.end());
  //ZeroFormVisitor_c visit_eng(eng_form, integrator_err, err_info_l);
  //visit(visit_eng, all.begin(), all.end()); 

  UnitZeroForm_c unit_form;
  ValueError_c::choice("VOL");
  FillFieldFromZeroForm_c fill_vol(err_info_l, unit_form);
  CommandOnIntegratorPath(fill_vol, integrator_err, all.begin(), all.end());
  //ZeroFormVisitor_c visit_vol(unit_form, integrator_err, err_info_l);
  //visit(visit_vol, all.begin(), all.end()); 

  ofstream out("error.txt");
  out << "energy  is  " << ValueError_c::total("ENG") << endl;
  out << "err abs is  " << ValueError_c::total("ABS") << endl;
  out << "err rel is  " << ValueError_c::total("REL") << endl;
  
  ValueError_c::choice("REL");
  GetValField_c<xtool::xIdentity<double> > val_err_info(err_info_l);
  pexport.Export(val_err_info, 
		 "ERROR_CONTRIBUTION_ESTIMATED", integrator, 
		 all.begin(), all.end());


  ValueError_c::choice("ENG");
  pexport.Export(val_err_info, 
		 "ENERGY", integrator, 
		 all.begin(), all.end());

  ValueError_c::choice("VOL");
  pexport.Export(val_err_info, 
		 "VOLUMES", integrator, 
		 all.begin(), all.end());

  ValueError_c::choice("DNS");
  pexport.Export(val_err_info, 
		 "DNS", integrator, 
		 all.begin(), all.end());
  
#endif




  Field_c err_info_l(&DoubleManager, SpaceConstant_c("ERROR_INFO"));
  CreateValue_c<ValueError_c>  creator_err;
  DeclareInterpolation(err_info_l, creator_err, all.begin(), all.end());
  GetValField_c<xtool::xIdentity<double> > val_err_info(err_info_l);
  ofstream out("error.txt");
#if 1
///////////////////////////////////////////
//// Exact Error computation ////////////////////
///////////////////////////////////////////

  IntegratorPartition_c integrator_exa(4); 

  //DiffusiveVectorEnergyExactZeroForm_c eng_exa_form(hooke, exact);
  //GetReferenceStrain_c exact_strain(exact);
  ZeroFormEvalBilinearFormWithLaw_c<EvalExactFlux_c,
                                    UniformMaterialSensitivity_c<xtensor::xTensor2>,
                                    EvalExactFlux_c > eng_exa_form(exact_flux, fourier);
  ValueError_c::choice("ENG_EXA");
  FillFieldFromZeroForm_c fill_eng_exa(err_info_l, eng_exa_form);
  CommandOnIntegratorPath(fill_eng_exa, integrator_exa, all.begin(), all.end());

  EvalBinary_c< std::minus<xtensor::xVector> > diff_flux(eval_grad, exact_flux); 
  ZeroFormEvalBilinearFormWithLaw_c<EvalBinary_c< std::minus<xtensor::xVector> >,
                                    UniformMaterialSensitivity_c<xtensor::xTensor2>,
                                    EvalBinary_c< std::minus<xtensor::xVector> > > err_exa_form(diff_flux, fourier);
  ValueError_c::choice("ABS2_EXA");
  FillFieldFromZeroForm_c fill_err_exa(err_info_l, err_exa_form);
  CommandOnIntegratorPath(fill_err_exa, integrator_exa, all.begin(), all.end());

//   ValueError_c::choice("EFF_EXA");
//   pexport.Export(val_err_info,
// 		 "LOCAL_EFFECTIVITY", integrator, 
// 		 all.begin(), all.end());

  ValueError_c::choice("REL_EXA");
  pexport.Export(val_err_info,
		 "ERROR_CONTRIBUTION_EXACT", integrator, 
		 all.begin(), all.end());

  
  out << "energy exact is   " << ValueError_c::total("ENG_EXA") << endl;
  out << "err abs exact is  " << ValueError_c::total("ABS_EXA") << endl;
  out << "err rel exact is  " << ValueError_c::total("REL_EXA") << endl;
//  out << "indice effect is  " << ValueError_c::total("EFF_EXA") << endl;
  

#endif






  return;
  
}


void Mechanics_c :: TreatmentOfEssEnv (const Field_c& listFunctionSpace, Data_c * data){

  for (data->PhysicalEnv->First(); !data->PhysicalEnv->IsDone(); data->PhysicalEnv->Next()) {
    Env_c env = data->PhysicalEnv->CurrentItem();
    string phys = env.Phys; 
    int type = env.Type; 
    int entity = env.Entity;
    //
    // FEM case
    //
    if (phys == "TEMPERATURE") {       
      if (type == FIX) {
        mxClassRegion bc(data->mesh, entity, env.getDimension());
	DirichletBoundaryCondition (listFunctionSpace, phys, bc.begin(), bc.end());
      }
      else assert(1 == 0);      
    }
} // End loop over the environemnt info
return;
}



void Mechanics_c :: TreatmentOfNatEnv   (const Field_c& listFunctionSpace, 
				       Assembler_c& assembler, Integrator_c& integrator,
				       Data_c * data, Boundary_c& groups){
  

  for (data->PhysicalEnv->First(); !data->PhysicalEnv->IsDone(); data->PhysicalEnv->Next()) {
    Env_c env = data->PhysicalEnv->CurrentItem();
    if (env.Phys == "SURFACIC_HEAT_FLUX") {
      assert(env.Type == FIX);
      ConstantEval_c<double>  flux(env.GetValue(0.0));
      LinearFormWithLoad_c<ValOperator_c<xtool::xIdentity<double> >, 
                           ConstantEval_c<double> > lin(flux); 
      mxClassRegion bc(data->mesh, env.Entity, env.getDimension());
      Assemble(lin, assembler, integrator,
	       listFunctionSpace, bc.begin(), bc.end(), mxUpperAdjacency()); 
    }
  } // End loop over the environemnt info
  return;
}
