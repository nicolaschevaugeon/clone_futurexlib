/*
     This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef __INTEGRATORS_SINGULARCRACK2D__H
#define __INTEGRATORS_SINGULARCRACK2D__H

#include <string>

#include "lCrack.h"
#include "xEntityFilter.h"
#include "xIntegrationRule.h"

namespace xcrack
{
class IntegratorSingularCrack2D_c : public xfem::xIntegrationRule
{
  public:
   IntegratorSingularCrack2D_c(lCrack& crck, int dpart = 0, int dg = 2)
       : crack(crck), degree(dpart), degree_gen(dg), filter(xfem::xAcceptAll()), filter_front(xfem::xAcceptAll())
   {
   }

   template <class T>
   IntegratorSingularCrack2D_c(lCrack& crck, int dpart, int dg, const T& f2)
       : crack(crck), degree(dpart), degree_gen(dg), filter(xfem::xAcceptAll()), filter_front(f2)
   {
   }

   template <class T>
   IntegratorSingularCrack2D_c(lCrack& crck, const T& f, int dpart = 0, int dg = 2, const T& f2 = xfem::xAcceptAll())
       : crack(crck), degree(dpart), degree_gen(dg), filter(f), filter_front(f2)
   {
   }

   void accept(xfem::xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

   void getSingularNode(xfem::xGeomElem* geom_integ);
   void createSingularMapping(int);

  private:
   xtensor::xVector<> orientation(AOMD::mEntity* e) const;

   lCrack& crack;
   int degree;      // for crack front (tip) elements
   int degree_gen;  // for other elements in general case
   mutable int SingPoint;

   xfem::xEntityFilter filter;
   xfem::xEntityFilter filter_front;  // for elements close to the front
   // attention : filter_front MUST include info for nodes as well (support)
};

}  // namespace xcrack

#endif
