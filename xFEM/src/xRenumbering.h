/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _RENUMBERING_HH
#define _RENUMBERING_HH

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------
#include <map>
#include <set>
#include <vector>

// -- Xfem --------------------------------------------------------------------
#include "xField.h"
#include "xFiniteElement.h"
#include "xRegion.h"
#include "xStateOfValue.h"
#include "xValue.h"

//#include "xRckVisitor.h"

// -- Boost -------------------------------------------------------------------
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

namespace xfem
{
// ----------------------------------------------------------------------------
// ReverseCuthillMcKeeNumbering
// ----------------------------------------------------------------------------
/*! \ingroup Utilities
    \brief Apply reverse Cuthill-McKee ordering to unknown nodal values.

    This algorithm assumes that the free nodal values have already been
    numbered. This means that these nodal values have a state pointer to
    a StateDof_c object, storing the actual number of the corresponding
    unknown (Numdof variable).

    The role of this function is to optimize the numbering of the unknowns in
    order to minimize the bandwidth of the linear system that will be build
    using this numbering.

    \param  field the field containing the unknowns.
    \param  start the first element of the assembling steps.
    \param  last   past the last element of the assembling steps.
    \return a pair containing the old and new bandwidth.
*/
template <typename VT>
pair<size_t, size_t> xReverseCutHillMcKeeNumbering(xField<VT>& field, xIter start, xIter last);
// ----------------------------------------------------------------------------
// xReverseCuthillMcKeeNumberingBoost
// ----------------------------------------------------------------------------
/*! \ingroup Utilities
    \brief Apply reverse Cuthill-McKee ordering to unknown nodal values.

    This algorithm assumes that the free nodal values have already been
    numbered. This means that these nodal values have a state pointer to
    a StateDof_c object, storing the actual number of the corresponding
    unknown (Numdof variable).

    The role of this function is to optimize the numbering of the unknowns in
    order to minimize the bandwidth of the linear system that will be build
    using this numbering.

    \param  field the field containing the unknowns.
    \param  start the first element of the assembling steps.
    \param  last   past the last element of the assembling steps.
    \return a pair containing the old and new bandwidth.

    This function uses boost graph library for graph and renumbering
    algorithm.


    \bug Problems occur when using this function with parallel test cases.
         Original measured bandwidth is a little bit lower than in the
         ReverseCuthillMcKeeRenumbering function and new bandwidth is often
         significantly lower than in the previous function (?).
         Work only on the degree of freedom in the group named "dofs" in the field.
*/
template <typename VT>
pair<int, int> xReverseCutHillMcKeeNumberingBoost(xField<VT>& field, xIter start, xIter last);

// ----------------------------------------------------------------------------
// CLASS xNumdofRenumberingVisitor
// ----------------------------------------------------------------------------
/*! \class NumdofRenumberingVisitor
    \ingroup Utilities
    \brief DoubleManager Visitor writer for changing dofs ordering.

    This visitor changes the Numdof variable for xStateOfValueDof values encountered
    in the range on which it is applied. New ordering is assumed to be stored
    in a random access container begining at index zero and new index also
    starts with zero. Thus an offset (which is equal to 1) is used by the
    visit method to map Numdof value with the new index.

    This visitor is relevant when unknowns numbering is modified by an
    ordering algorithm such as the reverse Cuthill-Mckee

*/
template <typename RANDOM_ITERATOR>
class xNumDofRenumberingVisitor
{
  public:
   // -- Constructors & destructor -------------------------------------------
   xNumDofRenumberingVisitor(RANDOM_ITERATOR s) : current(s) {}

   // -- Visitor method ------------------------------------------------------

   xValue<double>* operator()(xValue<double>* v)
   {
      size_t offset = 1;
      xStateOfValueDof* state = dynamic_cast<xStateOfValueDof*>(v->getState());

      if (state != nullptr)
      {
         state->Numdof = *(this->current + state->Numdof - offset) + offset;
      }
      return v;
   }

  private:
   RANDOM_ITERATOR current;
};

// ----------------------------------------------------------------------------
// CLASS WriteIndexedSolutionVisitor
// ----------------------------------------------------------------------------
/*! \class WriteIndexedSolutionVisitor
    \ingroup Utilities
    \brief DoubleManager Visitor writer for index based storage.

    The WriteSolutionVisitor of Xfem which assumes the same
    ordering for the input data and the nodal values storage.
    This visitor classes works with the index number stored in the xStateOfValueDof
    (Numdof variable). This index is used to retrieve the value in the input
    data, using iterator arithmetic and therefore this visitor could only
    work if the input data iterator is a random access iterator type.

    This visitor is relevant when unknowns numbering is modified by an
    ordering algorithm such as the reverse Cuthill-Mckee

*/
template <typename RANDOM_ITERATOR>
class xWriteIndexedSolutionVisitor
{
  public:
   // -- Constructors & destructor -------------------------------------------
   xWriteIndexedSolutionVisitor(RANDOM_ITERATOR s) : current(s) {}

   // -- Visitor method ------------------------------------------------------

   xValue<double>* operator()(xValue<double>* v)
   {
      size_t offset = 1;
      xStateOfValueDof* state = dynamic_cast<xStateOfValueDof*>(v->getState());

      if (state != nullptr)
      {
         size_t index = state->Numdof - offset;
         v->setVal(*(this->current + index));
         return v;
      }
      return nullptr;
   }

  private:
   RANDOM_ITERATOR current;
};

// ----------------------------------------------------------------------------
// CLASS xAddIndexedSolutionVisitor
// ----------------------------------------------------------------------------
/*!
  \ingroup Utilities
  \brief DoubleManager Visitor writer for index based storage.

*/
template <typename RANDOM_ITERATOR>
class xAddIndexedSolutionVisitor
{
  public:
   // -- Constructors & destructor -------------------------------------------
   xAddIndexedSolutionVisitor(RANDOM_ITERATOR s) : current(s) {}

   // -- Visitor method ------------------------------------------------------

   xValue<double>* operator()(xValue<double>* v)
   {
      size_t offset = 1;
      xStateOfValueDof* state = dynamic_cast<xStateOfValueDof*>(v->getState());
      if (state != nullptr)
      {
         size_t index = state->Numdof - offset;
         v->setVal(v->getVal() + *(this->current + index));
         return v;
      }
      return nullptr;
   }

  private:
   RANDOM_ITERATOR current;
};

// ----------------------------------------------------------------------------
// CLASS xAddIndexedScaledSolutionVisitor
// ----------------------------------------------------------------------------
/*!
Same as above xAddIndexedSolutionVisitor but instead of updating with the vector X it updates with fac*X
Usefull in a Newton Raphson resolution with a linear search.

*/

template <typename RANDOM_ITERATOR>
class xAddIndexedScaledSolutionVisitor
{
  public:
   // -- Constructors & destructor -------------------------------------------
   xAddIndexedScaledSolutionVisitor(RANDOM_ITERATOR s, double a) : current(s), fac(a) {}

   // -- Visitor method ------------------------------------------------------

   xValue<double>* operator()(xValue<double>* v)
   {
      size_t offset = 1;
      xStateOfValueDof* state = dynamic_cast<xStateOfValueDof*>(v->getState());
      if (state != nullptr)
      {
         size_t index = state->Numdof - offset;
         v->setVal(v->getVal() + *(this->current + index) * fac);
         return v;
      }
      return nullptr;
   }

  private:
   RANDOM_ITERATOR current;
   double fac;
};

#include "xRenumbering_imp.h"
}  // namespace xfem

#endif
// == END OF FILE =============================================================
