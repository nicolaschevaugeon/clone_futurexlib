/* 
     This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.

*/
  
#include "xcFormLinearEnergyRelease.h"
#include "xForm.h"
#include "xValKey.h"
#include "xEval.h"
#include "xGeomElem.h"
#include "xTensor2.h"
#include "xMaterial.h"
#include "xMaterialManager.h"
#include "CrackPostpro.h"
#include "xForm.h"
#include "xRegion.h"
#include "xMaterialSensitivity.h"
#include "xDebug.h"
#include <cmath>

using namespace std;
using namespace xfem;


namespace xcrack 
{




xcEvalEshelbyHPP::xcEvalEshelbyHPP(const xEval<xtensor::xTensor2<> >& g,
				   const xEval<xtensor::xTensor4<> >& p) 
  : eval_grad_disp(g), eval_phys(p) {}

void xcEvalEshelbyHPP::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
				  result_type& result) const
  {
    //const bool debug = false;
    eval_grad_disp(geo_appro, geo_integ, grad_disp);
    eval_phys(geo_appro, geo_integ, phys);
   
	
  
    xtensor::xTensor2<> stress = phys * grad_disp;
    
    //compute Eshelby tensor 
    xtensor::xTensor2<> stress_jk_grad_disp_ki(stress * grad_disp);
    
    double engh= 0.5 * stress_jk_grad_disp_ki.trace(); 
    
    for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) 
	  {
	    result(i,j) = - stress_jk_grad_disp_ki(j,i);
	  }
	result(i,i) += engh;
      } 
    return;
  }

xcEvalEshelbyHPPThermoElastic::xcEvalEshelbyHPPThermoElastic(const xEval<xtensor::xTensor2<> >& g,
				       const xEvalMaterialValues & s) 
  : eval_grad_disp(g), eval_meca_fields(s) {}

void xcEvalEshelbyHPPThermoElastic::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
				  result_type& result) const
  {
    //const bool debug = false;
    eval_grad_disp(geo_appro, geo_integ, grad_disp);
    xTensors meca_fields;
    eval_meca_fields(geo_appro, geo_integ, meca_fields);
      
    xtensor::xTensor2<> eps   = meca_fields.tensor2("elastic_strain"); 
    xtensor::xTensor2<> stress = meca_fields.tensor2("stress");   
    
    xtensor::xTensor2<> stress_jk_grad_disp_ki(stress * grad_disp);
  
    double engh  = 0.5*stress.contract(eps);
    for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) 
	  {
	    result(i,j) = - stress_jk_grad_disp_ki(j,i);
	  }
	result(i,i) += engh;
      } 
    return;
  }

  xcEvalEshelbyInteractionHPPThermoElastic::xcEvalEshelbyInteractionHPPThermoElastic(const xEval<xtensor::xTensor2<> >& g, const xEvalMaterialValues & _eval_meca_fields, const xcCrackBase & c , int _mode):eval_grad_disp(g), eval_meca_fields(_eval_meca_fields), mode(_mode), crack(c){} ;

  void xcEvalEshelbyInteractionHPPThermoElastic::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
						       result_type& result) const{
    xTensors mecafield;
    eval_meca_fields(geo_appro, geo_integ, mecafield);
    xtensor::xTensor2<> stress = mecafield.tensor2("stress");
    eval_grad_disp(geo_appro, geo_integ, grad_disp);
    crack.getAsymptoticFields(geo_appro, mode, stress_aux, grad_disp_aux); xtensor::xTensor2<> stress_jk_grad_disp_aux_ki(stress * grad_disp_aux);
    xtensor::xTensor2<> stress_aux_jk_grad_disp_ki(stress_aux * grad_disp);
    
    double work = stress_jk_grad_disp_aux_ki.trace(); 
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) 
	{
	  result(i,j) = - stress_jk_grad_disp_aux_ki(j,i)
	    - stress_aux_jk_grad_disp_ki(j,i);
	}
      result(i,i) += work;
    }
 
    return;
    
  }
  
  /*
xcEvalEshelbyHPPThermo::xcEvalEshelbyHPPThermo(const xEval<xtensor::xTensor2<> >& g,
					       const xEval<double >& t,
					       const xEval<xtensor::xTensor4<> >& p,
					       const xEval<double > &alpha,
					       double _T0)
  : eval_grad_disp(g), eval_temperature(t), eval_phys(p), eval_alpha(alpha), T0(_T0)  {}


void xcEvalEshelbyHPPThermo::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
				  result_type& result) const
{
  const bool debug = false;
  xtensor::xTensor2<> grad_disp;
  eval_grad_disp(geo_appro, geo_integ, grad_disp);
  xtensor::xTensor2<> eps(grad_disp );
  eps.symmetrize();
  xtensor::xTensor4<> phys;
  eval_phys(geo_appro, geo_integ, phys);
  double T;
  eval_temperature(geo_appro, geo_integ, T);
  double alpha;
  eval_alpha(geo_appro, geo_integ, alpha );

  xtensor::xTensor2<> epsth(eps);
  for (int i=0; i<3; ++i)
    epsth(i,i) -= alpha*(T-T0);
  
    
  xtensor::xTensor2<> stress = phys * epsth;
  double energy =0.5* stress.contract(epsth);
  
  //compute Eshelby tensor 
   xtensor::xTensor2<> tmp(stress * grad_disp);
   tmp.transpose();
   xtensor::xTensor2<> Eshelbyth = tmp *(-1.);

   
   for (int i=0; i<3; ++i){
     Eshelbyth(i,i) += energy;
   }
   result =Eshelbyth;
   return;
}
  */

xcEvalJintHPPThermoSourceTerm::xcEvalJintHPPThermoSourceTerm(const xEval<xTensors>& s,
							     const xEval<xtensor::xVector<> >& tg,			   
							     const xEval<double > &alpha, const xcCrackBase & _crack)
  : eval_meca(s), eval_grad_temperature(tg),  eval_alpha(alpha), crack(_crack) {}


void xcEvalJintHPPThermoSourceTerm::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
				  result_type& result) const
{
  //const bool debug = false;
  xTensors meca_val;
  eval_meca(geo_appro, geo_integ, meca_val);
  if (mode==0)
    stress = meca_val.tensor2("stress");
  else {
    crack.getAsymptoticFields(geo_appro, mode, stress_aux, grad_disp_aux);
    stress = stress_aux;
  }
  eval_grad_temperature(geo_appro, geo_integ, gT);
  eval_alpha(geo_appro, geo_integ, alpha );
  result =gT*alpha*stress.trace();
  // std::cout << "source " << result << std::endl;
  return;
}

  
xcEvalEshelbyHPPInteraction::xcEvalEshelbyHPPInteraction(const xEval<xtensor::xTensor2<> >& g,
							 const xEval<xtensor::xTensor4<> >& p,
							 const xcCrackBase& c, int m) 
  : eval_grad_disp(g), eval_phys(p), crack(c), mode(m)  {}

void xcEvalEshelbyHPPInteraction::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		result_type& result) const
  {
    const bool debug = false;
    eval_grad_disp(geo_appro, geo_integ, grad_disp);
    eval_phys(geo_appro, geo_appro, phys);
    xtensor::xTensor2<> stress = phys * grad_disp;
    crack.getAsymptoticFields(geo_appro, mode, stress_aux, grad_disp_aux);
       
    //compute Eshelby tensor 
    xtensor::xTensor2<> stress_jk_grad_disp_aux_ki(stress * grad_disp_aux);
    xtensor::xTensor2<> stress_aux_jk_grad_disp_ki(stress_aux * grad_disp);
    
    double work = stress_jk_grad_disp_aux_ki.trace(); 
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) 
	{
	  result(i,j) = - stress_jk_grad_disp_aux_ki(j,i)
	    - stress_aux_jk_grad_disp_ki(j,i);
	}
      result(i,i) += work;
    }
    
    xMaterial *mat = xMaterialManagerSingleton::instance().getMaterial(geo_appro);
    const xTensors* properties = mat->getProperties();
    double E   = properties->scalar("YOUNG_MODULUS");
    double nu  = properties->scalar("POISSON_RATIO");
  
    switch (mode){
    case 1:
    case 2:{
      double Estar = E /(1-nu*nu);
      result *= Estar/2.;
      break;
    }
    case 3:
      double mu = E/2/(1+nu);
      result *= mu;
      break;
    }
    
    if (debug) 
      {
        cout << " stress " << endl << stress << endl;
        cout << " grad_disp " << endl << grad_disp << endl;
	cout << " stress_jk_grad_disp_aux_ki " << endl << stress_jk_grad_disp_aux_ki << endl;
	cout << " stress_aux_jk_grad_disp_ki " << endl <<  stress_aux_jk_grad_disp_ki  << endl;
        cout << " work " << work << endl;
        cout << "eshelby interaction " << endl << result << endl;
      }

    return;
  }

xcEvalEshelbyHPPInteractionv2::xcEvalEshelbyHPPInteractionv2(const xEval<xtensor::xTensor2<> >& s,
							     const xEval<xtensor::xTensor2<> >& gd,
							     const xcCrackBase& c, int m) 
  : eval_stress(s), eval_grad_disp(gd), crack(c), mode(m)  {}

void xcEvalEshelbyHPPInteractionv2::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		result_type& result) const
  {
    //const bool debug = false;
    eval_grad_disp(geo_appro, geo_integ, grad_disp);
      
    eval_stress(geo_appro, geo_appro, stress); 
    crack.getAsymptoticFields(geo_appro, mode, stress_aux, grad_disp_aux);
    
    xtensor::xTensor2<> stress_jk_grad_disp_aux_ki(stress * grad_disp_aux);
    xtensor::xTensor2<> stress_aux_jk_grad_disp_ki(stress_aux * grad_disp);
    xtensor::xTensor2<> strain(grad_disp);
    strain.symmetrize();
    xtensor::xTensor2<> strain_aux(grad_disp_aux);
    strain_aux.symmetrize();
    
    double work = 0.5*(stress.contract(strain_aux) + stress_aux.contract(strain)); 
    

  
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) 
	{
	  result(i,j) = - stress_jk_grad_disp_aux_ki(j,i)
	    - stress_aux_jk_grad_disp_ki(j,i);
	}
      result(i,i) += work;
    }

    return;



  }

 
 
  /*
xcEvalEshelbyHPPInteractionThermoSource::xcEvalEshelbyHPPInteractionThermoSource ( const xEval<xtensor::xVector<> >& tg,
										   const xEval<double >& alpha,
										   const lCrack& c, int m):
  eval_temperature_gradient(tg), eval_alpha(alpha), crack(c), mode(m){};

void xcEvalEshelbyHPPInteractionThermoSource::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
						    result_type& result) const{
   getAuxFields(geo_appro, crack, mode, stress_aux, grad_disp_aux);

   xtensor::xTensor2<> ddilatthdT(0.);
   double alpha;
   eval_alpha(geo_appro, geo_integ, alpha ); 
   xtensor::xVector<> gT;
   eval_temperature_gradient(geo_appro, geo_integ, gT);
   for (int i=0; i<3 ; ++i)
     ddilatthdT(i,i) =alpha;
   result = gT *(stress_aux.contract(ddilatthdT));
   return;
}
  */
       
xcFormLinearEnergyRelease::xcFormLinearEnergyRelease (const xEval<xtensor::xVector<> >&  q,
			       const xEval<xtensor::xTensor2<> >& g, // a remplacer par un scalar_field
			       const xEval<xtensor::xTensor2<> >& e, 
			       const xEval<xtensor::xVector<> >&  s) 
  : eval_qdir(q), eval_grad_qdir(g), eval_eshelby(e), eval_source(s) 
{} 

void xcFormLinearEnergyRelease::accumulate_pnt(xfem::xGeomElem*  geo_integ)
{
  const bool debug = xdebug_flag;
  xtensor::xVector<> source;
  eval_source(geo_appro, geo_integ, source);
  //std::cout << source << std::endl;
  eval_qdir(geo_appro, geo_integ, qdir);
  eval_grad_qdir(geo_appro, geo_integ, grad_qdir);
  eval_eshelby(geo_appro, geo_integ, eshelby);
  
  double  before_val_test = xfem::xFormProducts::product(eshelby, grad_qdir); 
  
  xtensor::xVector<> before_grad_test = qdir * eshelby;
  
  if (debug) 
    {
      cout << " xyz " << geo_integ->getXYZ() << endl;
      cout << " qdir        " << endl << qdir << endl;
      cout << " grad_qdir   " << endl << grad_qdir << endl;
      cout << " eshelby         " << endl << eshelby << endl;
      cout << " before_val_test " << endl << before_val_test << endl;
      cout << " before_grad_test " << endl << before_grad_test << endl;
    }

  std::vector<xtensor::xVector<> > values_grad_test;
  GradTest.eval(test,geo_appro,geo_integ, values_grad_test);
  std::vector<double> values_val_test;
  ValTest.eval(test,geo_appro,geo_integ, values_val_test);
  
  double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();

  for(unsigned int i=0;i<test->size();++i){
    double val = - xfem::xFormProducts::product(before_grad_test, values_grad_test[i])   -  (before_val_test * values_val_test[i])  - (source*qdir*values_val_test[i]); 
    // val =  - (source*qdir*values_val_test[i]); 
    if (debug) 
      {
	cout << " test function   " << i << endl;
	cout << " values_grad_test[i] "   << values_grad_test[i] << endl;
	cout << " values_val_test[i]  "   << values_val_test[i]  << endl;
	cout << " before_grad_test" << before_grad_test << endl;
	cout << " before_val_test" << before_val_test << endl;
	cout << " - product(before_grad_test, values_grad_test[i])   -  (before_val_test * values_val_test[i]); " << val << endl;
	cout << "eshelby " << eshelby<< endl;
	cout << "source"  << source << endl;
	cout << "wdet " << wdet <<endl; 
      }
    val *= wdet;
    Vector.add(i, val);
  }
}


xcFormLinearEnergyReleaseBoundary::xcFormLinearEnergyReleaseBoundary (const xEval<xtensor::xVector<> >&  q,
							      const xEval<xtensor::xTensor2<> >& e) 
: eval_qdir(q), eval_eshelby(e) 
{} 

void xcFormLinearEnergyReleaseBoundary::accumulate_pnt(xfem::xGeomElem*  geo_integ)
{
  //const bool debug = false;
  eval_qdir(geo_appro, geo_integ, qdir);
  eval_eshelby(geo_appro, geo_integ, eshelby);
  eval_normal(geo_appro, geo_integ, normal);

  double  before_val_test = xfem::xFormProducts::product(qdir, eshelby, normal); 
  

  std::vector<double> values_val_test;
  ValTest.eval(test,geo_appro,geo_integ, values_val_test);
  
  double wdet = geo_integ->GetWeight() * geo_integ->GetDetJac();

  for(unsigned int i=0;i<test->size();++i){
    double val = before_val_test * values_val_test[i] ;         
    val *= wdet;
    Vector.add(i, val);
  }
}

}








