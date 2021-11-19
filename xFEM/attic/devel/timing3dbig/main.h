/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _MECHANICS_H
#define _MECHANICS_H 


#include "xValueManager.h"
#include "xData.h"
#include "xField.h"
#include "xAssembler.h"

using namespace xfem;

class Mechanics_c {
public :
  Mechanics_c ();
  ~Mechanics_c ();
  void TreatmentOfFormulation (xData *data);
  void TreatmentOfEssEnv   (const xField& listFunctionSpace, xData * data);
  void TreatmentOfNatEnv   (const xField& listFunctionSpace, 
			    xAssembler& assembler, xIntegrationRule& integration_rule,
			    xData * data, xBoundary& groups);

private:
  xValueManagerDist<double> double_manager;
};


 
#endif



