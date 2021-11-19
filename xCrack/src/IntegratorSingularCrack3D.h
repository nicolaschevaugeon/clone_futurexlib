/*
        This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef __INTEGRATORS_SINGULARCRACK3D__H
#define __INTEGRATORS_SINGULARCRACK3D__H

#include "xIntegrationRule.h"
//#include "IntegratorAbstract.h"
#include <string>

#include "lCrack.h"
#include "xEntityFilter.h"

namespace xcrack
{
class IntegratorSingularCrack3D_c : public xfem::xIntegrationRule
{
  public:
   IntegratorSingularCrack3D_c(const lCrack& crck, int dpart = 0, int dg = 2) : crack(crck), degree(dpart), degree_gen(dg)
   {  // filter(xAcceptAll()),filter_front(xAcceptAll()) {
      filter = xfem::xAcceptAll();
      filter_front = xfem::xAcceptAll();
   }

   template <class T>
   IntegratorSingularCrack3D_c(const lCrack& crck, int dpart, int dg, const T& f2)
       : crack(crck), degree(dpart), degree_gen(dg), degree_z(dg), filter(xfem::xAcceptAll()), filter_front(f2)
   {
   }

   template <class T>
   IntegratorSingularCrack3D_c(const lCrack& crck, int dpart, int dg, int dg_z, const T& f2)
       : crack(crck), degree(dpart), degree_gen(dg), degree_z(dg_z), filter(xfem::xAcceptAll()), filter_front(f2)
   {
   }

   template <class T>
   IntegratorSingularCrack3D_c(const lCrack& crck, const T& f, int dpart = 0, int dg = 2, const T& f2 = xfem::xAcceptAll())
       : crack(crck), degree(dpart), degree_gen(dg), filter(f), filter_front(f2)
   {
   }

   void accept(xfem::xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

   void getSingularNode(xfem::xGeomElem* geom_integ);
   void createSingularMapping(int);

  private:
   const lCrack& crack;
   int degree;      // for crack front (tip) elements
   int degree_gen;  // for other elements in general case
   int degree_z;    // along crack tip
   int SingPoint;

   xfem::xEntityFilter filter;
   xfem::xEntityFilter filter_front;  // for elements close to the front
   // attention : filter_front MUST include info for nodes as well (support)
};

}  // namespace xcrack
#endif
