/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _MECHANICS_H
#define _MECHANICS_H 

#include "Formulation.h"
#include "Material.h"
#include "Zone.h"
#include "BilinearFormsDerived.h"
#include "Fourier.h"
#include "ValueCreators.h"
#include <string>

class Assembler_c;
class Integrator_c;
class ReferenceSolution_c;

class ChangeKey_c 
{
public:
  ChangeKey_c(const string& e);
  void operator()(ValKey_c& key); 
private:
  string extension;
};




class Mechanics_c : public Formulation_c {
public :
  Mechanics_c ();
  ~Mechanics_c ();
  string  GetFormulationName(void)  const  { return "FORMULATION_MECHANICS"; }

  void TreatmentOfFormulation (Data_c *data);



  void TreatmentOfEssEnv   (const Field_c& listFunctionSpace, Data_c * data);
  void TreatmentOfNatEnv   (const Field_c& listFunctionSpace, 
			    Assembler_c& assembler, Integrator_c& integrator,
			    Data_c * data, Boundary_c& groups);

private:
  ReferenceSolution_c* exact;

};


// template <class R>
// class CreateRegularAndAverageValue_c {
//    public:
//      CreateRegularAndAverageValue_c(ValManager_c* v) : c_ave(v) {} 
//      Value_c* operator()(const ValKey_c& key) 
//        { 
//          if (Value_c* v = c_ave(key)) return v;
// 	 if (!c_ave.isAverage())      
// 	   { 
// 	     cout << "not average " << endl; 
// 	     return c_reg(key);
// 	   }
//          cout << "no val created " << endl;
// 	 return 0;
//        }
//    private:
//      CreateValueLinearCombination_c c_ave;
//      CreateValue_c<R>     c_reg;
// };
 

#endif
