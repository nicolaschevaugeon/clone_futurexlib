/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.

*/
  
#ifndef __INTEGRATORS_CRACK__H
#define __INTEGRATORS_CRACK__H

#include <string>
#include "mEntity.h"
#include "xIntegrationRule.h"
#include "xEntityFilter.h"


namespace xcrack
{

/// Integrator for crack problem. To degree for the integration are given.
/*!  One for most element, (dg) one for the element close to the crack tip, as defined by the user given filter. (f2) !*/
class IntegratorCrack_c : public xfem::xIntegrationRule {
    
public:
  
  IntegratorCrack_c(int dpart=0, int dg=2) 
    : degree(dpart),degree_gen(dg), filter(xfem::xAcceptAll()),filter_front(xfem::xAcceptAll()) {}
  
  template<class T>
  IntegratorCrack_c(int dpart, int dg, const T& f2) 
    : degree(dpart),degree_gen(dg), filter(xfem::xAcceptAll()),filter_front(f2) {}
  
  template<class T>
  IntegratorCrack_c(const T& f, int dpart=0, int dg=2, const T& f2=xfem::xAcceptAll()) 
    : degree(dpart), degree_gen(dg),filter(f),filter_front(f2) {}
  
  void accept(xfem::xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;
  
private:

  int degree; // for crack front (tip) elements
  int degree_gen; // for other elements in general case
  
  xfem::xEntityFilter filter;
  xfem::xEntityFilter filter_front; //for elements close to the front
  // attention : filter_front MUST include info for nodes as well (support)
  
};


}

#endif
