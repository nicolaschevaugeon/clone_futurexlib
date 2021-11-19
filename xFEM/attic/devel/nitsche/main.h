/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _MECHANICS_H
#define _MECHANICS_H 

#include <iostream>
class Integrator_c;
class Data_c;
class Field_c;
class Boundary_c;
class Assembler_c;
class Integrator_c;
#include "DoubleManager.h"
#include "Eval.h"
#include "GeomElem.h"
#include "xTensor2.h"
#include "mPoint.h"
#include "mext.h"

using namespace xfem;
using namespace std;


class  EvalExactFlux_c : public Eval_c<xtensor::xVector>
{
public:
  EvalExactFlux_c() : center(-1.,-1., 0.), ui(5000.), u(20000.), ro(2.), ri(1.) 
    {
     const bool debug = true;
     c1 = u  - (u-ui)/log(ro/ri)*log(ro);
     c2 = (u-ui)/log(ro/ri);
     cout << " c1 " << c1 << " c2 " << c2 << endl; 
    }

 void operator()(GeomElem_c*  geo_appro, GeomElem_c* geo_integ, result_type& res) {
   const bool debug = true;
   xtensor::xVector pos(center, geo_integ->getXYZ());
   double r = sqrt(pos*pos);
   if (r <= ri) 
    { res(0) = 0.; res(1) = 0.; res(2) = 0.; return; }
   res(0) =  c2/r * pos(0)/r;
   res(1) =  c2/r * pos(1)/r;
   res(2) = 0.0; 
 }
private:
 Trellis_Util::mPoint center;
 double ui, u, ro, ri;
 double c1, c2;
 result_type val;
};



class OuterEval_c :  public Eval_c<double> {
public:
 OuterEval_c(void) : center(-1.,-1., 0.), ui(5000.), u(20000.), ro(2.), ri(1.) 
   {
     const bool debug = true;
     c1 = u  - (u-ui)/log(ro/ri)*log(ro);
     c2 = (u-ui)/log(ro/ri);
     cout << " c1 " << c1 << " c2 " << c2 << endl; 
   }
 void operator()(GeomElem_c*  geo_appro, GeomElem_c* geo_integ, result_type& res) {
   const bool debug = true;
   xtensor::xVector pos(center, geo_integ->getXYZ());
   double r = sqrt(pos*pos);
   if (debug) cout << " center " << center << " r " << r << endl;
   res = c1 + c2 * log(r);
 }
private:
 Trellis_Util::mPoint center;
 double ui, u, ro, ri;
 double c1, c2;
 result_type val;
};



class Mechanics_c {
public :
  Mechanics_c ();
  ~Mechanics_c ();
  void TreatmentOfFormulation (Data_c *data);
  void TreatmentOfEssEnv   (const Field_c& listFunctionSpace, Data_c * data);
  void TreatmentOfNatEnv   (const Field_c& listFunctionSpace, 
			    Assembler_c& assembler, Integrator_c& integrator,
			    Data_c * data, Boundary_c& groups);

private:
  DoubleManager_c DoubleManager;
  EvalExactFlux_c exact_flux;
  //EvalExactInclusionStress_c exact_stress;

};

 
#endif
