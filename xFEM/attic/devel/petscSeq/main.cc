/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <cassert>
#include "main.h"
#include "xAlgorithm.h"
#include "xData.h"
#include "xValueManager.h"
#include "xCSRMatrix.h"
#include "xCSRVector.h"
#include "xLinearSystem.h"
#include "xLinearSystemSolverLU.h"
#include "xForm.h"
#include "xAssembler.h"
#include "xExportGmsh.h"
#include "xEnv.h"
#include "xSpace.h"
#include "xGeomElem.h"
#include "xFemMatrix.h"
#include "xMaterial.h"
#include "xIntegrationRule.h"
#include "xValue.h"
#include "xPhysSurf.h"
#include "xSimpleGeometry.h"
#include "xField.h"
#include "xApproxFunction.h"
#include "xOperators.h"
#include "xSolVisitor.h"
#include "xEval.h"
#include "xMaterialSensitivity.h"

#include "xLinearSystemSolverPetscSeq.h"
#include "xLinearSystemPetsc.h"
#include "petscksp.h"


#include <fstream>
using namespace xfem;
using namespace AOMD;



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
//  PetscInitialize(&argc,&argv,(char *)0,(char *)0);
// PetscInitialize(0,(char)0,(char *)0,(char *)0);

  xData data;
  char pname[256];
  Get_Options(argc,argv,pname);
  
  fprintf(stderr, "Starting reading the master data file\n");
  data.ReadInfo(pname);
  fprintf(stderr, "Done with reading the master data file\n");

  
  Mechanics_c Mechanics;
  Mechanics.TreatmentOfFormulation(&data);

// PetscFinalize();
  return 0; 
  
}


class SplitCriteria : public mSplitCallbacks
{
public :
  /*<i> constructor </i>*/
  SplitCriteria(mMesh *m): theMesh(m), maxLevelOfRefinement(1) {}
  int maxLevelOfRefinement;
  mMesh *theMesh;

  /*<i> this operator will be used by AOMD as the
     criterion for refinement. If it returns -1
     mesh entity e is coarsened, if it is refined.
     If it is 0, then the element is unchanged     
  </i>*/
  int operator () (mEntity *e)
  {
    if(theMesh->getRefinementLevel(e) >= maxLevelOfRefinement)
      return 0;
    for (int i = 0; i < e->size(0); ++i) {
      if (e->get(0,i)->getId() == 7) return 1;
      //return 1;
    }
    return 0;
  }
  /*<i> These two operators are used when AOMD is coupled with 
     a solver. User can e.g. project the finite element solution
     from the coarse mesh to the refined or the opposite </i>*/
  virtual void splitCallback  (mEntity*){};
  virtual void unsplitCallback(mEntity*){};  
};


void Mechanics_c :: TreatmentOfFormulation (xData *data) {

// Necessary...
  PetscStuffInitialize();

  fprintf(stderr, "Starting reading the mesh file\n");
  data->ReadMesh();
  fprintf(stderr, "End      reading the mesh file\n");
  data->ReadZones();

  xRegion all(data->mesh);

  xSpaceLagrange lagx("DISPLACEMENT_X", xSpace::VECTOR_X, xSpaceLagrange::DEGREE_ONE);
  xSpaceLagrange lagy("DISPLACEMENT_Y", xSpace::VECTOR_Y, xSpaceLagrange::DEGREE_ONE);
  xSpaceComposite  lagrange(lagx, lagy);

  xField disp_l(&double_manager, lagrange);

  xValueCreator<xValueDouble>  creator;
  SetXfemDebugFlag(true);
  DeclareInterpolation(disp_l, creator, all.begin(), all.end());
  SetXfemDebugFlag(false);
  double_manager.PrintForDebug("dcl.dbg");


  printf("Hello 0 in thermic\n"); 

  //TreatmentOfEssEnv(disp_l, data);
  std::vector<string> vars_ess;
  vars_ess.push_back("DISPLACEMENT_X");
  vars_ess.push_back("DISPLACEMENT_Y");
  TreatmentOfEssentialEnv(disp_l, data->PhysicalEnv->begin(), data->PhysicalEnv->end(), 
			  vars_ess, data->mesh);


  double_manager.PrintForDebug("ess.dbg");

  printf("Hello 0.1 in thermic\n"); 

  xStateDofCreator<> snh(double_manager, "dofs");
  DeclareState(disp_l, snh, all.begin(), all.end());
//cout << "nbr dof after nothanging " << double_manager.size("dofs") << endl;
 

  double_manager.PrintForDebug("sym.dbg");
printf("Hello 2 in thermic nb dof is %d\n", double_manager.size("dofs"));

  xIntegrationRuleBasic integration_rule_env(3);
  xIntegrationRuleBasic integration_rule(3); //2 is the degree to integrate

printf("Hello 2.05 in thermic\n");
//TreatmentOfNatEnv(disp_l, assembler, integration_rule_env, data, data->allGroups);
  std::vector<string> vars_nat;
  vars_nat.push_back("TRACTION_X");
  vars_nat.push_back("TRACTION_Y");
  


printf("Hello 2.06 in thermic\n");
  xUniformMaterialSensitivity<xtensor::xTensor4> hooke("strain");
printf("Hello 2.06 in thermic\n");
//  DiffusiveVectorxFormBilinear diffusive(hooke);
   xFormBilinearWithLaw<xGradOperator<xtool::xIdentity<xtensor::xTensor2> >, 
                         xUniformMaterialSensitivity<xtensor::xTensor4>,
                         xGradOperator<xtool::xIdentity<xtensor::xTensor2> > > diffusive(hooke);
printf("Hello 2.07 in thermic\n");

#if 0
// // Classical CSR Version...-----------------------------------------------------
  xlinalg::xCSRVector b(double_manager.size("dofs"));
  xlinalg::xCSRVector sol(double_manager.size("dofs"));
  xlinalg::xCSRMatrix A(double_manager.size("dofs"));


printf("Hello 2.01 in thermic\n");
  xlinalg::xLinearSystemSolverLU solver;
printf("Hello 2.02 in thermic\n");
  xLinearSystem system(&A, &solver);


printf("Hello 2.03 in thermic\n");
  xAssemblerBasic<> assembler(A, b);


AssembleNaturalEnvVector(assembler, integration_rule_env, disp_l, 
			   data->PhysicalEnv->begin(), data->PhysicalEnv->end(), 
			   vars_nat, data->mesh, xUpperAdjacency());


Assemble(diffusive, assembler, integration_rule, disp_l, disp_l, all.begin(), all.end()); 
system.Solve(b, sol);
Visit(xWriteSolutionVisitor(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));
#endif


#if 0
// // Petsc Version 1 (dirty)...-----------------------------------------------------
// // Note that we directly use petsc objects (PetscVector and Matrix are typenames of petsc objects)


// Note that PetscStuffInitialize and Finalize are necessary... (see before and below...)

//Vecteurs...
  PetscVector bp,solp;
  int ierr;
  ierr = VecCreate(PETSC_COMM_WORLD,&bp);
  ierr = VecSetSizes(bp,PETSC_DECIDE,double_manager.size("dofs"));
  ierr = VecSetFromOptions(bp);
  ierr = VecDuplicate(bp,&solp);

//Matrice...
  PetscMatrix Ap;
 ierr = MatCreate(PETSC_COMM_WORLD,&Ap);
  ierr = MatSetSizes(Ap,PETSC_DECIDE,PETSC_DECIDE,double_manager.size("dofs"),double_manager.size("dofs"));
  ierr = MatSetFromOptions(Ap);

//Interfaces pour l'assembleur...
  pVectorInterface bpI(bp);
  pMatrixInterface ApI(Ap);
  xAssemblerBasic<pMatrixInterface, pVectorInterface,double> assemblerp(ApI, bpI);

AssembleNaturalEnvVector(assemblerp, integration_rule_env, disp_l, 
			   data->PhysicalEnv->begin(), data->PhysicalEnv->end(), 
			   vars_nat, data->mesh, xUpperAdjacency());


Assemble(diffusive, assemblerp, integration_rule, disp_l, disp_l, all.begin(), all.end()); 


  VecAssemblyBegin(bp);
  VecAssemblyEnd(bp);
  MatAssemblyBegin(Ap, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Ap, MAT_FINAL_ASSEMBLY);


  KSP            ksp;
  PC             pc;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
  ierr = KSPSetOperators(ksp,Ap,Ap,DIFFERENT_NONZERO_PATTERN);
  ierr = KSPGetPC(ksp,&pc);
  ierr = PCSetType(pc,PCJACOBI);
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  ierr = KSPSetFromOptions(ksp);
  ierr = KSPSolve(ksp,bp,solp);
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

std::vector<double> asol(double_manager.size("dofs"));
double    *asolD(&asol[0]);
ierr = VecGetArray(solp,&asolD);

  //--- On ecrit dans le DoubleManager...
 xWriteSolutionVisitorP visitor(asolD);

//  Visit(xWriteSolutionVisitor(asolD), double_manager.begin("dofs"), double_manager.end("dofs"));
 Visit(visitor,double_manager.begin("dofs"),double_manager.end("dofs"));
 VecRestoreArray(solp,&asolD);
 PetscFree(asolD);
#endif


#if 1
// // Petsc Version 2 (cleaner)...-----------------------------------------------------

// Note that PetscStuffInitialize and Finalize are necessary... (see before and below...)

//Vecteurs... 
 xPetscVector bx(double_manager.size("dofs")),solx;
 solx.setSize(double_manager.size("dofs"));
// Remark that solx and bx were created in a different way...


// Matrices...
// Default preallocation...
xPetscMatrix Ax(double_manager.size("dofs"));
// User defined Preallocation...
// petscMatrix Ax(double_manager.size("dofs"),ComputeNnz(disp_l,all.begin(), all.end()));


  xAssemblerBasic<xPetscMatrix,xPetscVector,double> assemblerx(Ax, bx);

AssembleNaturalEnvVector(assemblerx, integration_rule_env, disp_l, 
			   data->PhysicalEnv->begin(), data->PhysicalEnv->end(), 
			   vars_nat, data->mesh, xUpperAdjacency());



// On prealloue un peu plus que necessaire (20% : d'ou le 1.2) (car Ax a ete cree par preallocation automatique)
  Ax.preallocate(ComputeNnz(disp_l,all.begin(), all.end()),1.2);

  Assemble(diffusive, assemblerx, integration_rule, disp_l, disp_l, all.begin(), all.end()); 

  bx.finalAssembly();
  Ax.finalAssembly();





 xlinalg::xLinearSystemSolverPetscSeq solverx;
// Define PC and KSP : see different types in xlinalg::xLinearSystemSolverPetscSeq source.
solverx.setKSPType("KSPGMRES");
solverx.setPCType("PCJACOBI");
// Define tolerance
solverx.setTol(1e-14);



 solverx.solve(&Ax, bx, solx);


{

  //xWriteSolutionVisitorPetsc  visitorxx(solx);
xWriteSolutionVisitor<xPetscVector>  visitorxx(solx.begin());
Visit(visitorxx,double_manager.begin("dofs"),double_manager.end("dofs"));
}



#endif


#if 1
// // EXPORT :
  double_manager.PrintForDebug("res.dbg");

  printf("Hello 3 in thermic\n");

  xexport::xExportGmshAscii  pexport;
  xexport::xExportGmshBinary  pbinexport;

  printf("Hello 5 in thermic\n");
  xEvalField<xtool::xIdentity<xtensor::xVector> > eval_disp(disp_l);
  Export(eval_disp, pexport, "DISPLACEMENT", integration_rule, all.begin(), all.end());
  //

  xEvalGradField<xtensor::xSymmetrize> eval_strain(disp_l);
  xEvalBinary< xtool::xMult<xtensor::xTensor4, xtensor::xTensor2, xtensor::xTensor2> > stress(hooke, eval_strain);
  Export(stress, pbinexport, "STRESS", integration_rule, all.begin(), all.end());
#endif


// Necessary...
PetscStuffFinalize( );
  return;
  
}






