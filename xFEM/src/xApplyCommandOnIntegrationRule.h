#ifndef XAPPLYCOMMANDONINTEGRATIONRULE_H
#define XAPPLYCOMMANDONINTEGRATIONRULE_H

#include "xCommandOnGeomElem.h"
#include "xIntegrationRule.h"
#include "xDebug.h"

namespace xfem{
template <class ITER>
void ApplyCommandOnIntegrationRule(xCommandOnGeomElem& command,
                                   const xIntegrationRule& integration_rule,
                                   ITER it, ITER end,
                                   xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
  const bool debug = xdebug_flag;
  for ( ; it != end; ++it)
    {
      AOMD::mEntity* e_integ = *it;
      AOMD::mEntity* e_appro = integ2appro(e_integ);
      if (debug) std::cout<<"ApplyCommandOnIntegrationRule case #\n";
      if (debug) e_integ->print();
      if (debug) std::cout<<"E_appro:\n";
      if (debug) e_appro->print();
      xGeomElem geom_appro(e_appro);
      xGeomElem geom_integ(e_integ);
      command.setIntegElem(&geom_integ);
      command.openApproxElem(&geom_appro);
      integration_rule.accept(command, e_integ);
      command.closeApproxElem(&geom_appro);
    }
}


}

#endif // XAPPLYCOMMANDONINTEGRATIONRULE_H
