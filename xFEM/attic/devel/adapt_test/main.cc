/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "main.h"
#include <cassert>
#include "CSR_Matrix.h"
#include "CSR_Vector.h"
#include "LinearSystem.h"
#include "LinearSystemSolverITL.h"

#include "ZeroFormsDerived.h"
#include "BilinearFormsDerived.h"
#include "LinearFormsDerived.h"
#include "AssemblerBasic.h"
#include "ExportOwner.h"
#include "Export.h"
#include "PhysicalEnv.h"
#include "SpaceLagrange.h"
#include "SpaceFiltered.h"
#include "SpaceConstant.h"
#include "SpaceComposite.h"
#include "GeomElem.h"
#include "FemMatrix.h"
#include "Hooke.h"
#include "Bc.h"
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
#include "ReferenceSolutions.h"
#include "Visitor.h"

#include <fstream>

//for adaptation
#include "MeshAdapt.h"
#include "AdaptUtil.h"
#include "AOMD_Internals.h"
#include "MeshTools.h"



Mechanics_c::Mechanics_c() 
{ exact =  new InclusionSolution_c(1.0, 0.3,  1.0, 0.3, 0.4); }

Mechanics_c::~Mechanics_c() 
{ delete exact; }


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
  Data_c data;
  char pname[256];
  Get_Options(argc,argv,pname);
  
  fprintf(stderr, "Starting reading the master data file\n");
  data.ReadInfo(pname);
  fprintf(stderr, "Done with reading the master data file\n");
  fprintf(stderr, "Starting checking the master data file\n");
  data.CheckInfo();
  fprintf(stderr, "Passed checking the master data file\n");

  fprintf(stderr, "Starting reading the mesh file\n");
  data.ReadMesh();
  fprintf(stderr, "End      reading the mesh file\n");
  
  Mechanics_c Mechanics;
  Mechanics.TreatmentOfFormulation(&data);
  return 0; 
  
}

#if 1
int CB_count=0;
extern "C" void myCallback(pPList oldCavity, pPList newcavity,void *userdata)
{
  CB_count++;
  // PList_printx(oldCavity);
  // PList_printx(newcavity);
}
#endif

void Mechanics_c :: TreatmentOfFormulation (Data_c *data) {


  mxRegion all(data->mesh);


// mesh refinement
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
#if 1
  mxMesh* mesh = data->mesh;
  MeshAdapt rdr(mesh,0,0,1);    // snap off  
  rdr.setCallback(myCallback,0);
  int n,i;
  int ne=M_numEdges(mesh);
  int *ind=new int[ne];
  for( i=0; i<ne; i++ )  ind[i]=0;
    srand(1001);
  for( i=0; i<ne; i++ ) {
    n=rand()%ne;
    ind[n]=1;
  }
  ///////////
  EIter editer=M_edgeIter(mesh);
  i=0; n=0;
  pEdge edge;
  while( edge=EIter_next(editer) ) 
    if( ind[i++] ) {
      rdr.setAdaptLevel(edge,1);
        n++;
    }
  EIter_delete(editer);

  cout<< n <<" of "<<ne<<"mesh edges are marked for refinement"<<endl;
  delete[] ind;
    
  // do refinement
  cout<<"-------Begin refining and snapping-----"<<endl;
  rdr.run();
  cout << " refinement and snapping done !!!" << endl;
  cout << " call Callback function "<<CB_count<<" times"<<endl;


#endif
// END mesh refinement
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


    

  SpaceLagrange_c lagx("DISPLACEMENT_X", VECTOR_X, DEGREE_ONE);
  SpaceLagrange_c lagy("DISPLACEMENT_Y", VECTOR_Y, DEGREE_ONE);
//  SpaceLagrange_c lagz("DISPLACEMENT_Z", VECTOR_Z, DEGREE_TWO);

//  SpaceComposite_c  lagrange(lagx, lagy, lagz);
  SpaceComposite_c  lagrange(lagx, lagy);


  Field_c disp_l(&ValManager, lagrange);

  CreateValue_c<ValueDouble_c>  creator;
  DeclareInterpolation(disp_l, creator, all.begin(), all.end());
  ValManager.PrintForDebug("dcl.dbg");


  printf("Hello 0 in thermic\n"); 

  TreatmentOfEssEnv(disp_l, data);
  ValManager.PrintForDebug("ess.dbg");

  printf("Hello 0.1 in thermic\n"); 

  CreateStateDof_c<> snh(ValManager, "dofs");
  DeclareState(disp_l, snh, all.begin(), all.end());
//cout << "nbr dof after nothanging " << ValManager.size("dofs") << endl;
 

  ValManager.PrintForDebug("sym.dbg");
printf("Hello 2 in thermic nb dof is %d\n", ValManager.size("dofs"));
  CSR_Vector_c b(ValManager.size("dofs"));
  CSR_Matrix_c A(ValManager.size("dofs"));


printf("Hello 2.01 in thermic\n");
  LinearSystemSolverITL_c solver;
printf("Hello 2.02 in thermic\n");
  LinearSystem_c system(&A, &b, &solver);

printf("Hello 2.03 in thermic\n");
  AssemblerBasic_c assembler(A, b);
  IntegratorBasic_c integrator_env(3);
  IntegratorBasic_c integrator(3); //2 is the degree to integrate

printf("Hello 2.05 in thermic\n");
  TreatmentOfNatEnv(disp_l, assembler, integrator_env, data, data->allGroups);

printf("Hello 2.06 in thermic\n");
  Hooke_c hooke(data->zones);
printf("Hello 2.06 in thermic\n");
  DiffusiveVectorBilinearForm_c diffusive(hooke);
printf("Hello 2.07 in thermic\n");
  Assemble(diffusive, assembler, integrator, disp_l, disp_l, all.begin(), all.end()); 
printf("Hello 2.08 in thermic\n");

  printf("Hello 2.3 in thermic\n");
//std::ofstream out("mat.m");
//A.OutputMatrixOctaveFormat(out);
//return;

  system.Solve(ValManager, "dofs");


  //printf("Hello 2.6 in thermic\n");

  //DofManager.StoreResult(SOLUTION.GetArray());

  ValManager.PrintForDebug("res.dbg");

  printf("Hello 3 in thermic\n");

  Export_c * pexport = ExportOwner_c::Instance()->GetExport(&ValManager);


  printf("Hello 5 in thermic\n");
  pexport->ExportXFEM(disp_l, "DISPLACEMENT", data, all.begin(), all.end(), mxAcceptAll());
  //printf("Hello 6 in thermic\n");
  //pexport->ExportDualScalar(disp_l, "flux", data, all.begin(), all.end(), this);


#if 1
///////////////////////////////////////////
//// Estimated Error computation ////////////////////
///////////////////////////////////////////
  SpaceLagrange_c lagerrx("ERROR_X", VECTOR_X, DEGREE_TWO_MINUS_ONE); 
  SpaceLagrange_c lagerry("ERROR_Y", VECTOR_Y, DEGREE_TWO_MINUS_ONE);
  //SpaceLagrange_c lagerrz("ERROR_Z", VECTOR_Z, DEGREE_TWO_MINUS_ONE);

//  SpaceComposite_c lagerr(lagerrx, lagerry, lagerrz); 
  SpaceComposite_c lagerr(lagerrx, lagerry); 

  Field_c err_l(&ValManager, lagerr);

  DeclareInterpolation(err_l, creator, all.begin(), all.end());
  ValManager.PrintForDebug("dcl_err.dbg");

  TreatmentOfEssEnv(err_l, data);  
  ValManager.PrintForDebug("ess_err.dbg");

  printf("Hello 3... in thermic\n"); 

  //declare temperature as fixed now
  printf("Hello 3.001 in thermic\n"); 
  DeleteState(disp_l, all.begin(), all.end());
  printf("Hello 3.002 in thermic\n"); 
  DeclareState(disp_l, CreateStateFixed_c(), all.begin(), all.end());
  printf("Hello 3.003 in thermic\n"); 

 
  CreateStateDof_c<> snh_err(ValManager, "dofs_err");
  DeclareState(err_l, snh_err, all.begin(), all.end());

  ValManager.PrintForDebug("sym.dbg");
printf("Hello 3 in thermic nb dof is %d\n", ValManager.size("dofs_err"));
  CSR_Vector_c b_err(ValManager.size("dofs_err"));
  CSR_Matrix_c A_err(ValManager.size("dofs_err"));

  IntegratorBasic_c integrator_env_err(3);
  IntegratorBasic_c integrator_err(3); //2 is the degree to integrate

printf("Hello 3.01 in thermic\n");
  LinearSystemSolverITL_c solver_err;
printf("Hello 3.02 in thermic\n");
  LinearSystem_c system_err(&A_err, &b_err, &solver_err);

printf("Hello 3.03 in thermic\n");
  AssemblerBasic_c assembler_err(A_err, b_err);


printf("Hello 3.05 in thermic\n");
  TreatmentOfNatEnv(err_l, assembler_err, integrator_env_err, data, data->allGroups);

printf("Hello 3.07 in thermic\n");
  Assemble(diffusive, assembler_err, integrator_err, err_l, err_l,  all.begin(), all.end()); 
printf("Hello 3.075 in thermic\n");
  Assemble(diffusive, assembler_err, integrator_err, err_l, disp_l, all.begin(), all.end()); 
printf("Hello 3.08 in thermic\n");

  system_err.Solve(ValManager, "dofs_err");

  ValManager.PrintForDebug("res_err.dbg");

  //pexport->ExportXFEM(err_l, "ERROR", data, all.begin(), all.end(), mxAcceptAll());

  //declare err as fixed now
  printf("Hello 3.001 in thermic\n"); 
  DeleteState(err_l, all.begin(), all.end());
  printf("Hello 3.002 in thermic\n"); 
  DeclareState(err_l, CreateStateFixed_c(), all.begin(), all.end());
  printf("Hello 3.003 in thermic\n"); 


///////////////////////////////////////////////////////////////
//  Necessary for the error computation
////////////////////////////////////////////////////////////////
  Field_c err_info_l(&ValManager, SpaceConstant_c("ERROR_INFO"));
  CreateValue_c<ValueError_c>  creator_err;
  DeclareInterpolation(err_info_l, creator_err, all.begin(), all.end());


  DiffusiveVectorZeroForm_c err_form(hooke, err_l);
  ValueError_c::choice("ABS2");
  ZeroFormVisitor_c visit_err(err_form, integrator_err, err_info_l);
  Visit(visit_err, all.begin(), all.end()); 

  DiffusiveVectorZeroForm_c eng_form(hooke, disp_l);
  ValueError_c::choice("ENG");
  ZeroFormVisitor_c visit_eng(eng_form, integrator_err, err_info_l);
  Visit(visit_eng, all.begin(), all.end()); 

  UnitZeroForm_c unit_form;
  ValueError_c::choice("VOL");
  ZeroFormVisitor_c visit_vol(unit_form, integrator_err, err_info_l);
  Visit(visit_vol, all.begin(), all.end()); 

  std::ofstream out("error.txt");
  out << "energy  is  " << ValueError_c::total("ENG") << endl;
  out << "err abs is  " << ValueError_c::total("ABS") << endl;
  out << "err rel is  " << ValueError_c::total("REL") << endl;
  out << "vol     is  " << ValueError_c::total("VOL") << endl;
  
  ValueError_c::choice("REL");
  pexport->ExportXFEM(err_info_l, "ERROR_CONTRIBUTION_ESTIMATED", data, all.begin(), all.end(), 
		      mxAcceptAll());

  ValueError_c::choice("ENG");
  pexport->ExportXFEM(err_info_l, "ENERGY", data, all.begin(), all.end(), mxAcceptAll());

  ValueError_c::choice("VOL");
  pexport->ExportXFEM(err_info_l, "VOLUMES", data, all.begin(), all.end(), mxAcceptAll());

  ValueError_c::choice("DNS");
  pexport->ExportXFEM(err_info_l, "ERROR_DENSITY", data, all.begin(), all.end(), mxAcceptAll());
  
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
    if (phys == "DISPLACEMENT_X" || phys == "DISPLACEMENT_Y" || phys == "DISPLACEMENT_Z") {       
      if (type == FIX) {
cout << "before Dirichmet entity " << entity << " dim " << env.getDimension() << endl;
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
    if (env.Phys == "TRACTION_X" || env.Phys == "TRACTION_Y"  || env.Phys == "TRACTION_Z" ) {
      assert(env.Type == FIX);
      mVector val;
      if (env.Phys == "TRACTION_X") val(0) =  env.GetValue(0.0);
      if (env.Phys == "TRACTION_Y") val(1) =  env.GetValue(0.0);
      if (env.Phys == "TRACTION_Z") val(2) =  env.GetValue(0.0);
      cout << "val is " << val(0) << " " << val(1) << " " << val(2) << endl;
      ConstantVectorEnv_c flux(val);
      WorkVectorLinearForm_c lin(flux);
      mxClassRegion bc(data->mesh, env.Entity, env.getDimension());
      Assemble(lin, assembler, integrator, listFunctionSpace, bc.begin(), bc.end(), 
	       mxUpperAdjacency()); 
      //IntegratorBoundary_c bc_int(data->mesh, env.Entity, env.getDimension(), 2);
      //Assemble(lin, assembler, bc_int, listFunctionSpace,  bc_int.begin(), bc_int.end()); 
    }
    if (env.Phys == "STRESS") {
      assert(env.Type == FIX);
      NormalTensorEnv_c flux(exact);
      WorkVectorLinearForm_c lin(flux);
      mxClassRegion bc(data->mesh, env.Entity, env.getDimension());
      Assemble(lin, assembler, integrator, listFunctionSpace, bc.begin(), bc.end(), 
	       mxUpperAdjacency()); 

    }
  } 
  return;
}



