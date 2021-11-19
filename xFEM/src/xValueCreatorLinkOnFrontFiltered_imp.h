/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _VALUECREATORLINKONFRONTFILTERED_IMP_H_
#define _VALUECREATORLINKONFRONTFILTERED_IMP_H_
// Xfem
#include "xAlgorithm.h"
#include "xSpacePolynomial.h"
#include "xValueCreatorLinkOnFrontFiltered.h"
#include "xAssembler.h"
#include "xGraphMatrix.h"
#include "xCSRMatrix.h"

namespace xfem
{
template <typename ITERFRONT>
xValueCreatorLinkOnFrontFiltered::xValueCreatorLinkOnFrontFiltered(const ITERFRONT& begin_, const ITERFRONT& end_,
                                                                   xValueManagerDist<double>* v_,
                                                                   const double tol_, bool link_isolated_,
                                                                   xEntityFilter filt_, xEntityToEntity entity_to_entity_)
 : value_creator(begin_, end_, v_, link_isolated_, filt_), double_manager(v_)
{
  // declare values and states
  xSpacePolynomialLagrange space("loff", xSpace::SCALAR, 1);
  xField<double> field(double_manager, space);
  xStateDofCreator<> state_dof_creator(*double_manager, "loff_dofs");
  DeclareValueField(field, value_creator, value_creator.beginVertexIter(), value_creator.endVertexIter());
  DeclareState(field, state_dof_creator, value_creator.beginVertexIter(), value_creator.endVertexIter());
  int size_dof = double_manager->size("loff_dofs");
  // create matrix graph
  xlinalg::xGraphMatrix graph(size_dof, size_dof, false);
  AssembleGraph(graph, field, begin_, end_, entity_to_entity_);
  // assemble lumped mass matrix
  xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<double> >, xValOperator<xtool::xIdentity<double> > > form;
  xlinalg::xCSRMatrix M(size_dof);
  xAssemblerLumped<> assembler(M);
  xIntegrationRuleBasic integration_rule(2);
  Assemble(form, assembler, integration_rule, field, field, begin_, end_, entity_to_entity_, entity_to_entity_);

  buildTable(size_dof, tol_, field, M, graph);
}
}

#endif
