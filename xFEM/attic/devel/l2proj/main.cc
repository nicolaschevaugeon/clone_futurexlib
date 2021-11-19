/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "main.h" 

#include "xMesh.h"

using namespace xfem;
using namespace std;



int main(int argc, char *argv[])
{ 


//creation of a little mesh, two quadrilateral elements.
//
//
//  3        4        5
//  -------------------
//  |        |        |
//  |        |        |
//  |        |        |
//  |        |        |
//  |        |        |
//  -------------------
//  0        1        2
//
// The two elements are on the surface number 100.
//

//
//mesh creation
//

  xMesh mesh;
  
  mesh.createVertex(0, 0., 0., 0., 0);
  mesh.createVertex(1, 1., 0., 0., 0);
  mesh.createVertex(2, 2., 0., 0., 0);
  mesh.createVertex(3, 0., 1., 0., 0);
  mesh.createVertex(4, 1., 1., 0., 0);
  mesh.createVertex(5, 2., 1., 0., 0);
  
  
  mesh.createFaceWithVertices(0,1,4,3,mesh.getGEntity(100,2));
  mesh.createFaceWithVertices(1,2,5,4,mesh.getGEntity(100,2));
  AOMD::classifyUnclassifiedVerices(&mesh);
  
  mesh.exportGmsh("mesh.msh");  

 
//un exemple de projection L2
xValueManagerDist<double> vals;
xSpaceLagrange temp_space("TEMPERATURE", SCALAR, DEGREE_ONE);
xField         temp(&vals, temp_space);
xValueCreator<xValueDouble>  creator_double;

DeclareInterpolation(temp, creator_double, mesh.begin(2), mesh.end(2));

xStateDofCreator<> snh(double_manager, "dofs");
DeclareState(temp, snh, mesh.begin(2), matrix.end(2));         


xlinalg::xCSRVector b(vals.size("dofs"));
xlinalg::xCSRVector sol(vals.size("dofs"));
xlinalg::xCSRMatrix A(vals.size("dofs"));
xAssemblerBasic<> assembler(A, b);


xIntegrationRuleBasic integration_rule(2);



xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<double> >, 
                        xValOperator<xtool::xIdentity<double> > > bilin; 


Assemble(bilin, assembler, integration_rule, temp, temp, mesh.begin(2), mesh.end(2));


xFormLinearWithLoad<xValOperator<xtool::xIdentity<double> >, xEval<double> >  lin(t_exact); 

Assemble(lin, assembler, integration_rule, temp, mesh.begin(2), mesh.end(2));


 system.Solve(b, sol);
 Visit(xWriteSolutionVisitor(sol.begin()), 
       vals.begin("dofs"), 
       vals.end("dofs"));





return 0;
}
