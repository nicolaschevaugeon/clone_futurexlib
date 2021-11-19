/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

/*
   INCOMPRESSIBILITE A PRESSION CONTINUE :
   Elements T2P1-c (aussi appele T6P3-c) (passe la condition INF-SUP (LBB)  )
   Elements T1P1-c (aussi appele T3P3-c) (echoue au INF-SUP (LBB)  )

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
#include "GMSHBINARYExport.h"
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
#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
#include "xIterator.h"
#include "Field.h"
#include "FunctionDerivDiscXFEM.h"
#include "SpaceXFEM.h"
#include "Operators.h"
#include "SolVisitor.h"

#include "NewtonianFluid.h"
#include "Eval.h"
#include "AssemblerTranspose.h"

#include <fstream>
using namespace xfem;




NewtonianFlow_c::NewtonianFlow_c() 
{ 
  MaterialManagerSingleton_c::Instance().registerMaterial("newtonian_fluid", CreateMaterial_c<NewtonianFluid_c>() );
}

NewtonianFlow_c::~NewtonianFlow_c() 
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

int main(int argc, char *argv[])  {


  Data_c data;
  char pname[256];
  Get_Options(argc,argv,pname);
  
  fprintf(stderr, "Starting reading the master data file\n");
  data.ReadInfo(pname);
  fprintf(stderr, "Done with reading the master data file\n");

  NewtonianFlow_c incompressible_mechanics;
  incompressible_mechanics.TreatmentOfFormulation(&data);

  return 0; 
  
}


void NewtonianFlow_c :: TreatmentOfFormulation (Data_c *data) {

  const bool debug = false;
  fprintf(stderr, "Starting reading the mesh file\n");
  data->ReadMesh();
  fprintf(stderr, "End      reading the mesh file\n");
  data->ReadZones();
  xRegion all(data->mesh);

  SpaceLagrange_c velox("VELOCITY_X", VECTOR_X, DEGREE_TWO);
  SpaceLagrange_c veloy("VELOCITY_Y", VECTOR_Y, DEGREE_TWO);
  SpaceLagrange_c press("PRESSURE", SCALAR, DEGREE_ONE);

  SpaceComposite_c  velo(velox, veloy);
  SpaceComposite_c  pressure(press);

  Field_c velo_l (&DoubleManager, velo);
  Field_c press_l(&DoubleManager, pressure);

  CreateValue_c<ValueDouble_c>  creator;
  DeclareInterpolation(velo_l,  creator, all.begin(), all.end());
  DeclareInterpolation(press_l, creator, all.begin(), all.end());
  if (debug) DoubleManager.PrintForDebug("dcl.dbg");

  TreatmentOfEssEnv(press_l, data);
  TreatmentOfEssEnv(velo_l, data);
  xClassRegion bc_left(data->mesh,  112, 1);
  xClassRegion bc_right(data->mesh, 110, 1);
  LinearSystemSolverSparskit_c solver_l2;
  IntegratorBasic_c integrator_l2(4);
  EvalPara parabole;
  EvalFctAtPoint_c<xtensor::xVector, EvalPara>  eval_para(parabole);
  L2Projection(velo_l, eval_para, integrator_l2, solver_l2, bc_left.begin(), bc_left.end());
  DeclareState(velo_l, CreateStateFixed_c(), bc_left.begin(), bc_left.end());
  L2Projection(velo_l, eval_para, integrator_l2, solver_l2, bc_right.begin(), bc_right.end());
  DeclareState(velo_l, CreateStateFixed_c(), bc_right.begin(), bc_right.end());


  if (debug) DoubleManager.PrintForDebug("ess.dbg");

  CreateStateDof_c<> snh(DoubleManager, "dofs");
  DeclareState(velo_l, snh, all.begin(), all.end());
  DeclareState(press_l, snh, all.begin(), all.end());

  if (debug) DoubleManager.PrintForDebug("sym.dbg");
  CSR_Vector_c b(DoubleManager.size("dofs"));
  CSR_Vector_c sol(DoubleManager.size("dofs"));
  CSR_Matrix_c A(DoubleManager.size("dofs"));


  LinearSystemSolverSparskit_c solver;
  LinearSystem_c system(&A, &solver);

  AssemblerBasic_c assembler(A, b);
  AssemblerTranspose_c assembler_transpose(A, b);
  

  IntegratorBasic_c integrator_uu(4); 
  IntegratorBasic_c integrator_up(3); 

  UniformMaterialSensitivity_c<xtensor::xTensor4> newton("strain_rate");
  
  BilinearFormWithLaw_c<GradOperator_c< xfem::symmetrize >, 
                        UniformMaterialSensitivity_c<xtensor::xTensor4>,
                        GradOperator_c< xfem::symmetrize > > kuu(newton);
  
  //DEF de KUP
  BilinearFormWithoutLaw_c<GradOperator_c< trace >, 
                           ValOperator_c<xtool::xIdentity<double> > > kup;

  Assemble(kuu, assembler, integrator_uu, velo_l, velo_l, all.begin(), all.end()); 
  assembler_transpose.setCoeff(-1.0);
  Assemble(kup, assembler_transpose, integrator_up, velo_l, press_l, all.begin(), all.end());


  if (debug) 
    {
      std::ofstream out1("mat.m");
      A.OutputMatrixOctaveFormat(out1);
      std::ofstream out2("vecB.m");
      b.OutputOctaveFormat(out2);
    }
  system.Solve(b, sol);
  visit(WriteSolutionVisitor_c(sol.begin()), 
	DoubleManager.begin("dofs"), 
	DoubleManager.end("dofs"));

  if (debug) 
    {
      std::ofstream out3("vecS.m");
      sol.OutputOctaveFormat(out3);
    }  

  if (debug) DoubleManager.PrintForDebug("res.dbg");


  GMSHASCIIExport_c  pexport;

  //pexport.setRecurSplit(4);
  GetValField_c<xtool::xIdentity<xtensor::xVector> > eval_velo(velo_l);
  pexport.Export(eval_velo, "VELOCITY", integrator_uu, all.begin(), all.end());
  //

  GetValField_c<xtool::xIdentity<double> > eval_press(press_l);
  pexport.Export(eval_press, "PRESSURE", integrator_up, all.begin(), all.end());



}

void NewtonianFlow_c :: TreatmentOfEssEnv (const Field_c& listFunctionSpace, Data_c * data){

  for (data->PhysicalEnv->First(); !data->PhysicalEnv->IsDone(); data->PhysicalEnv->Next()) {
    Env_c env = data->PhysicalEnv->CurrentItem();
    string phys = env.Phys; 
    int type = env.Type; 
    int entity = env.Entity;
    //
    // FEM case
    //
    if (phys == "VELOCITY_X" || phys == "VELOCITY_Y" || phys == "VELOCITY_Z" || phys == "PRESSURE" ) {       
      if (type == FIX) {
cout << "before Dirichmet entity " << entity << " dim " << env.getDimension() << endl;
        xClassRegion bc(data->mesh, entity, env.getDimension());
	DirichletBoundaryCondition (listFunctionSpace, phys, bc.begin(), bc.end(), env.GetValue());
      }
      else assert(1 == 0);      
    }
} // End loop over the environemnt info
return;
}



void NewtonianFlow_c :: TreatmentOfNatEnv   (const Field_c& listFunctionSpace, 
				       Assembler_c& assembler, Integrator_c& integrator,
				       Data_c * data, Boundary_c& groups){
  

  for (data->PhysicalEnv->First(); !data->PhysicalEnv->IsDone(); data->PhysicalEnv->Next()) {
    Env_c env = data->PhysicalEnv->CurrentItem();
    if (env.Phys == "TRACTION_X" || env.Phys == "TRACTION_Y"  || env.Phys == "TRACTION_Z" ) {
      assert(env.Type == FIX);
      xtensor::xVector val;
      if (env.Phys == "TRACTION_X") val(0) =  env.GetValue();
      if (env.Phys == "TRACTION_Y") val(1) =  env.GetValue();
      if (env.Phys == "TRACTION_Z") val(2) =  env.GetValue();
      cout << "val is " << val(0) << " " << val(1) << " " << val(2) << endl;
      ConstantEval_c<xtensor::xVector>  flux(val);
      LinearFormWithLoad_c<ValOperator_c<xtool::xIdentity<xtensor::xVector> >, ConstantEval_c<xtensor::xVector> > lin(flux); 
      xClassRegion bc(data->mesh, env.Entity, env.getDimension());
      Assemble(lin, assembler, integrator, listFunctionSpace, bc.begin(), bc.end(), 
	       xUpperAdjacency()); 
    }
  } 
  return;
}
