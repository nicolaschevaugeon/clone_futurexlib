/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _NEWTONIANFLOW_H
#define _NEWTONIANFLOW_H

class Integrator_c;
class ReferenceSolution_c;
class Data_c;
class Field_c;
class Boundary_c;
class Assembler_c;
class Integrator_c;
#include "DoubleManager.h"
#include "xVector.h"

class EvalPara 
{
public:
  xtensor::xVector operator()(const Trellis_Util::mPoint& p) 
    {
      return xtensor::xVector(1-p(1)*p(1), 0.,0.);
    }
};

class NewtonianFlow_c {
public :
  NewtonianFlow_c ();
  ~NewtonianFlow_c ();
  void TreatmentOfFormulation (Data_c *data);
  void TreatmentOfEssEnv   (const Field_c& listFunctionSpace, Data_c * data);
  void TreatmentOfNatEnv   (const Field_c& listFunctionSpace, 
			    Assembler_c& assembler, Integrator_c& integrator,
			    Data_c * data, Boundary_c& groups);

private:
  DoubleManager_c DoubleManager;
};


 
#endif
