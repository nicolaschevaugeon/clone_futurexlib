/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

// SolverBase
#include "xGraphMatrix.h"
#include "xCSRMatrix.h"
// Xfem
#include "xValueCreatorLinkOnFrontFiltered.h"
// std
#include <algorithm>
#include <iterator>
#include <utility>

using AOMD::mVertex;
using xlinalg::xGraphMatrix;
using xlinalg::xCSRMatrix;
using std::distance;
using std::make_pair;
using std::min_element;
using std::pair;
using std::vector;

namespace xfem
{
void xValueCreatorLinkOnFrontFiltered::buildTable(const int size_dof, const double tol, const xField<double>& field, const xlinalg::xCSRMatrix& M, xlinalg::xGraphMatrix& graph)
{
  // create numdof to vertices relation
  std::vector <AOMD::mVertex*> vertices(size_dof);
  for (std::vector<AOMD::mVertex*>::iterator it=value_creator.beginVertexIter(), end=value_creator.endVertexIter(); it!=end; ++it)
  {
    xFiniteElement FEM;
    FEM.setKeys(*it, field.begin(), field.end());
    std::vector<xValue<double>*> values;
    double_manager->getValPtr(FEM.beginKey(), FEM.endKey(), values);
    for (std::vector<xValue<double>*>::iterator itv=values.begin(), endv=values.end(); itv!=endv; ++itv)
    {
      xStateOfValue* state=(*itv)->getState();
      if (state->state==xStateOfValue::DOF)
        vertices[static_cast<xStateOfValueDof*>(state)->Numdof-1]=*it;
    }
  }
  // clean up double manager
  DeleteState(field, value_creator.beginVertexIter(), value_creator.endVertexIter());
  DeleteValueField(field, value_creator.beginVertexIter(), value_creator.endVertexIter());
  double_manager->clear_subset("loff_dofs");
  // copy mass into vector (just to facilitate life later)
  std::vector<double> masses;
  masses.reserve(size_dof);
  for (int i=0; i<size_dof; ++i)
    masses.push_back(M.GetMatrix(i+1,i+1)); // Fortran convention !!!
  // create table, modify the graph on-the-fly and dispatch mass to neighboors
  vector<double>::iterator first = masses.begin();
  vector<double>::iterator min_elt = min_element(masses.begin(), masses.end());
  double val = *min_elt;
  while (val<tol)
  {
    int j = distance(first, min_elt);
    int* col = graph.getCol(j);
    int n = graph.getSizeCol(j);
    if (n==2) // in 2D if DOF j is on boundary, we proceed with its neighboor
    {
      if (j==col[0]) j=col[1];
      else j=col[0];
      n = graph.getSizeCol(j);
      if (n==2) {
        std::cout << "Error in xValueCreatorLinkOnFrontFiltered" << std::endl;
        throw; // means that the front has an isolated edge
      }
      col = graph.getCol(j);
      val = masses[j];
    }
    val/=static_cast<double>(n-1);
    for (int i=0; i<n; ++i)
    {
      if (col[i]!=j)
      {
        table.insert(make_pair(vertices[j], vertices[col[i]]));
        masses[col[i]]+=val; // Fortran convention !!!
        graph.remove(j, col[i]);
        for (int k=0; k<n; ++k)
          if (col[k]!=j && col[k]!=col[i]) graph.add(col[k], col[i]);
      }
    }
    masses[j]=std::numeric_limits<double>::max();
    graph.clearCol(j);
    min_elt = std::min_element(masses.begin(), masses.end());
    val = *min_elt;
  }
}

xValueCreatorLinkOnFrontFiltered::~xValueCreatorLinkOnFrontFiltered() = default;

xValue<double>* xValueCreatorLinkOnFrontFiltered::operator()(const xValKey& key) const
{
  pair<Table::const_iterator, Table::const_iterator> equal_range = table.equal_range(static_cast<mVertex*>(key.getEnti()));
  int dist = distance(equal_range.first, equal_range.second);
  if (dist)
  {
    xValKey key_other = key;
    vector<xValue<double>*> values;
    values.reserve(dist);
    for (Table::const_iterator it=equal_range.first, end=equal_range.second; it!=end; ++it)
    {
      key_other.setEnti(it->second);
      xValue<double>* value_other = double_manager->find(key_other);
      if (!value_other) return nullptr;
      else values.push_back(value_other);
    }
    vector<double> coeffs(values.size(), 1./static_cast<double>(values.size()));
    return new xValueLinearCombination<double>(coeffs, values);
  }
  else return value_creator(key);
}

vector<mVertex*>::iterator xValueCreatorLinkOnFrontFiltered::beginVertexIter()
{
  return value_creator.beginVertexIter();
}

vector<mVertex*>::iterator xValueCreatorLinkOnFrontFiltered::endVertexIter()
{
  return value_creator.endVertexIter();
}
}
