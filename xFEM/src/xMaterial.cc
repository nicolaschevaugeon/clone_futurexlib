/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <cstdio>
#include "xMaterial.h"
#include "xTensor2.h"
#include "xDebug.h"

using namespace std;

namespace xfem
{

// xUniformMaterialSensitivity::xUniformMaterialSensitivity(const std::string& variable_) 
//   : mat_manager(xMaterialManagerSingleton::Instance()), variable(variable_)   {}

// virtual xUniformMaterialSensitivity::~xUniformMaterialSensitivity() {}

// void xUniformMaterialSensitivity::operator()(xGeomElem*  geo_appro, 
// 					     xGeomElem* geo_integ, result_type& curr_law) 
//   {
//     const bool debug = xdebug_flag;
//     if(debug)
//       {
// 	cout<<"Appro=\n"; 
// 	geo_appro->getEntity()->print();
// 	cout<<"Integ=\n";
// 	geo_integ->getEntity()->print();
//       }
//     if (debug) cout << "in setcurrentelem zone_id is " << geo_integ->getZone() << endl;
//     xMaterial * mat = mat_manager.getMaterial(geo_integ);
//     if (debug) cout << "after mat_manager.getMaterial" << endl;
//     if (debug) assert(mat != 0);
//     mat->sensitivityTo(variable, curr_law);
//   }


xConductive::xConductive()  
{
  properties_signature.register_string("MATERIAL_CLASS");
  properties_signature.register_string("NAME");
  properties_signature.register_scalar("THERMIC_CONDUCTIVITY");
  properties.setSignature(&properties_signature);
  properties.astring("MATERIAL_CLASS") = "MATERIAL_CONDUCTIVE";
}

void xConductive::checkProperties()
{
  //check the properties
  double conductivity = properties.scalar("THERMIC_CONDUCTIVITY");
  assert(conductivity > 0.);
  
  //Set the derived informations
  SetThermicConductivityIsotropic (conductivity, thermic_conductivity);
}

void  xConductive::sensitivityTo(const std::string& phys_token, 
				  xtensor::xTensor2<>& sensitivity) {
  if (phys_token == "temperature_gradient") sensitivity = thermic_conductivity;
  else { 
    fprintf(stderr, "No sensitivity coded\n"); assert(0); }
  return;
}

void xConductive::SetThermicConductivityIsotropic(double k, xtensor::xTensor2<>& conductivity){
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      conductivity(i,j) = k * delta(i,j);
    }
  }
  return;
}



xNewtonianFluid::xNewtonianFluid() 
{
  variables_signature.register_tensor2("stress");
  variables_signature.register_tensor2("strain_rate");
 
  properties_signature.register_string("MATERIAL_CLASS");
  properties_signature.register_string("NAME");
  properties_signature.register_scalar("VISCOSITY");
  properties.setSignature(&properties_signature);
  properties.astring("MATERIAL_CLASS") = "MATERIAL_NEWTONIAN_FLUID";
  
}

void xNewtonianFluid::checkProperties() 
{
  //const bool debug = xdebug_flag;
  //const double viscosity = properties.scalar("VISCOSITY");
  assert(properties.scalar("VISCOSITY") >= 0);

  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
	for (int l = 0; l < 3; ++l) {
          // 0.5???
	  viscosity_tensor(i,j,k,l) = (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }


}


//protected member functions
void  xNewtonianFluid::sensitivityTo(const std::string& phys_token, 
			       xtensor::xTensor4<>& sensitivity) {
  if (phys_token == "strain_rate")       sensitivity = viscosity_tensor;
  else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
  return; 
}


xElastic::xElastic() 
{
  variables_signature.register_tensor2("stress");
  variables_signature.register_tensor2("strain");
 
  properties_signature.register_string("MATERIAL_CLASS");
  properties_signature.register_string("NAME");
  properties_signature.register_scalar("YOUNG_MODULUS");
  properties_signature.register_scalar("POISSON_RATIO");
  properties_signature.register_scalar("DENSITY");
  properties.setSignature(&properties_signature);
  properties.astring("MATERIAL_CLASS") = "MATERIAL_ELASTIC";
}




void xElastic::checkProperties() 
{
  const bool debug = xdebug_flag;
//check the properties
  double young_modulus = properties.scalar("YOUNG_MODULUS");
  double poisson_ratio = properties.scalar("POISSON_RATIO");
/// modification patrick
  //double density =
  properties.scalar("DENSITY");
///
  assert(young_modulus > 0);
  assert(poisson_ratio < 0.5  && poisson_ratio > -1.0);
/// modification patrick
//  assert(density > 0);
///

//Set the derived informations
  SetElasticStiffnessIsotropic   (young_modulus, poisson_ratio, elastic_stiffness);
  SetElasticComplianceIsotropic  (young_modulus, poisson_ratio, elastic_compliance);

//for xtensor::xTensor4Isotropic
  E = young_modulus;
  nu = poisson_ratio;  
  lam = nu*E/((1.+nu)*(1.-2.*nu));
  mu  = E/(2.*(1.+nu));

//Set the admissible formulations
  if (debug) cout << "Elastic material created done " << endl;  
}


//protected member functions
 void  xElastic::sensitivityTo(const std::string& phys_token, 
 			       xtensor::xTensor4<>& sensitivity) {
   const bool debug = xdebug_flag;
   if (debug) cout << "inside  xElastic::sensitivityTo" << endl;
   if (phys_token == "strain")       sensitivity = elastic_stiffness;
   else if (phys_token == "stress")  sensitivity = elastic_compliance;
   else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
   return; 
 }

void  xElastic::sensitivityTo(const std::string& phys_token, 
                  xtensor::xTensor4Isotropic& sensitivity) {
   const bool debug = xdebug_flag;
   if (debug) cout << "inside  xElastic::sensitivityTo" << endl;
   if (phys_token == "strain")       { sensitivity.lam = lam; sensitivity.mu = mu; }
   else if (phys_token == "stress")  { sensitivity.lam = -(nu/E);  sensitivity.mu =(1+nu)/2./E;  }
   else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
   return; 
 }



void  xElastic::sensitivityTo(const std::string& phys_token,
                  xtensor::xTensor4AnisoPlaneStrain& sensitivity) {
   const bool debug = xdebug_flag;
   if (debug) cout << "inside  xElastic::sensitivityTo xtensor::xTensor4AnisoPlaneStrain" << endl;
   if (phys_token == "strain")       { sensitivity.E = E; sensitivity.nu = nu; sensitivity.updateCoeffs(); }
   else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
   return;
 }


void  xElastic::sensitivityTo(const std::string& phys_token,
                  xtensor::xTensor4AnisoPlaneStress& sensitivity) {
   const bool debug = xdebug_flag;
   if (debug) cout << "inside  xElastic::sensitivityTo" << endl;
   if (phys_token == "strain")       { sensitivity.E = E; sensitivity.nu = nu; sensitivity.updateCoeffs();}
   else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
   return;
 }



void  xElastic::sensitivityTo(const std::string& phys_token, 
			       double& sensitivity) {
  if (phys_token == "acceleration") sensitivity = properties.scalar("DENSITY");
  else if (phys_token == "young_modulus") sensitivity = properties.scalar("YOUNG_MODULUS");
  else if (phys_token == "poisson_ratio") sensitivity = properties.scalar("POISSON_RATIO");
  else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
  return;
}

void xElastic::computeCurrentState()
{
  curr->tensor2("stress")  = elastic_stiffness * curr->tensor2("strain");
  return;
}


//Note the identity tensor in [3][3][3][3]
//0.5 * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k))
void xElastic::SetElasticStiffnessIsotropic(double E, double nu, 
					    xtensor::xTensor4<>& stiffness) {
  const double lam = nu*E/((1.+nu)*(1.-2.*nu));
  const double mu  = E/(2.*(1.+nu));
  int i, j, k, l;
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = lam * delta(i,j) * delta(k,l) 
	    + mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }
  return;
}
//Note the identity tensor in [3][3][3][3]
//0.5 * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k))
void xElastic::SetElasticStiffnessIsotropicPlate(double E, double nu, 
						 xtensor::xTensor4<>& stiffness) {
  const double lam = nu*E/((1.+nu)*(1.-2.*nu));
  const double mu  = E/(2.*(1.+nu));
  const double lamp = 2. * lam * mu / (lam + 2. *mu); 
  int i, j, k, l;
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = lamp * delta(i,j) * delta(k,l) 
	    + mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }
  return;
}

void xElastic::SetElasticStiffnessIsotropicFromLame(double lam, double mu, 
						      xtensor::xTensor4<>& stiffness) {
  int i, j, k, l;
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = lam * delta(i,j) * delta(k,l) 
	    + mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }
  return;
}
















void xElastic::SetElasticComplianceIsotropic(double E, double nu, 
					       xtensor::xTensor4<>& compliance){

  const double lam1 = -nu/E;
  const double mu1  = (1.+nu)/(2.*E);
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      for (int k = 0; k < 3; ++k){
	for (int l = 0; l < 3; ++l){
	  compliance(i,j,k,l) = lam1 * delta(i,j) * delta(k,l) 
	    + mu1 * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k)); 
	}
      }
    }
  }
  return;
}

void xElastic::SetElasticComplianceIsotropicPlate(double E, double nu, 
					       xtensor::xTensor4<>& compliance){

  
  const double lam = nu*E/((1.+nu)*(1.-2.*nu));
  const double mu  = E/(2.*(1.+nu));
  const double lamp = 2. * lam * mu / (lam + 2. *mu);
  const double Ep = mu*(3*lamp+2*mu)/(lamp+mu);
  const double nup= lamp/(2*(lamp+mu));

  const double lam1p = -nup/Ep;
  const double mu1p  = (1.+nup)/(2.*Ep);

//  const double lam1 = -nu/E;
//  const double mu1  = (1.+nu)/(2.*E);


  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      for (int k = 0; k < 3; ++k){
	for (int l = 0; l < 3; ++l){
	  compliance(i,j,k,l) = lam1p * delta(i,j) * delta(k,l) 
	    + mu1p * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k)); 
	}
      }
    }
  }
  return;
}


void xElastic::SetElasticStiffnessOrthotropic(double E1,   double E2,   double E3,
						double nu12, double nu23, double nu31,
						double mu12, double mu23, double mu31,
						xtensor::xTensor4<>& stiffness)
{
// permet de creer le tenseur d'elasticite pour un materiau orthotrope
// 9 coeff sont donnes dans les axes d'orthotropie
// les trois permier sont les modules de Young selon les axes 1, 2 et 3:  E1  E2    E3
// les tris suivant sont les coefficients de poisson dans l'ordre:       nu12 nu23 nu31
// les trois derniers sont les modules de cisaillement dans l'order       mu12 mu23 mu31

//On met d'abord a zero
  int i, j, k, l;
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = 0.0;
	}
      }
    }
  }

//
// les relations suivantes sont verifiees nu(ij)/E(i) = nu(ji)/E(j)
// Reference T.J. Chung Continuum Mechanics Prentice Hall 1988

  const double nu21  = E2 * nu12/E1; 
  const double nu32  = E3 * nu23/E2;
  const double nu13  = E1 * nu31/E3;
  const double J     = (1-nu12*nu21-nu23*nu32-nu31*nu13-2.0*nu21*nu32*nu13)/(E1*E2*E3);
 

//on rentre les valeurs non nulles du tableau ci-dessous (21 valeurs)
//indice 1,2,3 ATTENTION plus tard indice 0,1,2 
//  E1111  E1122  E1133  E1112  E1123  E1113
//         E2222  E2233  E2212  E2223  E2213
//                E3333  E3312  E3323  E3313
//                       E1212  E1223  E1213
//                              E2323  E2313
//                                     E1313


  stiffness(0,0,0,0) = (1.0 - nu23 * nu32) / (J * E2 * E3);
  stiffness(1,1,1,1) = (1.0 - nu13 * nu31) / (J * E1 * E3);
  stiffness(2,2,2,2) = (1.0 - nu12 * nu21) / (J * E1 * E2);

  stiffness(0,0,1,1) = (nu21 + nu31 * nu23) / (J * E2 * E3);
  stiffness(0,0,2,2) = (nu31 + nu21 * nu32) / (J * E2 * E3);
  stiffness(1,1,2,2) = (nu32 + nu12 * nu31) / (J * E1 * E3);

  stiffness(0,1,0,1) = mu12;
  stiffness(1,2,1,2) = mu23;
  stiffness(0,2,0,2) = mu31;

//On applique les symetries pour retrouver la 3x3x3x3 totale
//Les symmetries sont 
//  Eijkl = Eijlk = Ejikl = Eklij
//  ATTENTION c'est delicat

//   On fait d'abord Eklij = Eijkl avec (ij) dans 11 22 33 12 23 13 et
//                                      (kl) dans 11 22 33 12 23 13  
//   Cela donne 15 infos
  for (i = 0; i < 3; ++i){
    for (j = i+1; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = k+1; l < 3; ++l){
	  stiffness(k,l,i,j) = stiffness(i,j,k,l);
	}
      }
    }
  }
 
//   On fait Eijlk = Eijkl pout tout (ij) dans 11 22 33 12 23 13 et pour (kl) dans 12 23 13
//   Cela donne 18 infos
  for (i = 0; i < 3; ++i){
    for (j = i+1; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = k+1; l < 3; ++l){
	  stiffness(i,j,l,k) = stiffness(i,j,k,l);
	}
      }
    }
  }

//   On fait Ejikl = Eijkl pout tout (kl) dans 11 22 33 12 23 13 21 32 31 et pour (ij) dans 12 23 13
//   Cela donne 27 infos
  for (i = 0; i < 3; ++i){
    for (j = i+1; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(j,i,k,l) = stiffness(i,j,k,l);
	}
      }
    }
  }

//   au total on a les 21 de depart plus 15 + 18 + 27 ce qui fait bien 81 en tout

  return;
}

void xElastic::SetShearStiffnessIsotropic(double E, double nu, xtensor::xTensor4<>& stiffness){
  const double mu  = E/(2.*(1.+nu));
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      for (int k = 0; k < 3; ++k){
	for (int l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }
  return;
}


void xElastic::SetInvLameOneIsotropic(double E, double nu, double& val){
  val = ((1.+nu)*(1.-2.*nu))/(nu*E);
  return;
}


//
// For plates
//

// it is the plane stress matrix times (h^3/12)
void xElastic::SetPlateBendingStiffnessIsotropic(double E, double nu, double t, 
						   xtensor::xTensor4<>& stiffness)
{
  const double lam = nu*E/((1.+nu)*(1.-2.*nu));
  const double mu  = E/(2.*(1.+nu));
  const double lamp = 2. * lam * mu / (lam + 2. *mu);
  const double coeff = (t*t*t)/12.;
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      for (int k = 0; k < 3; ++k){
	for (int l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = coeff * (lamp * delta(i,j) * delta(k,l) 
		+ mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k)));
	}
      }
    }
  }
  return;
}


void xElastic::SetPlateShearStiffnessIsotropic(double E, double nu, double t, 
						 xtensor::xTensor2<>& stiffness)
{
  const double mu  = E/(2.*(1.+nu));

  for (int i = 0; i < 3; ++i){
    for (int j =  0; j < 3; ++j){
      stiffness(i,j) = mu*t*(5./6.)*delta(i,j);
    }
  }
  return;
}

void xElastic::SetElasticStiffnessTransverselyIsotropic(double E, double E3, double G3,
					      double nu, double nu3, xtensor::xTensor4<>& stiffness) {

//see http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_iso_transverse.cfm



// Z (3) is the fiber direction
// X and Y (1 and 2 ) are transverse direction

  int i, j, k, l;  
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = 0.0;
	}
      }
    }
  }

  const double Delta = (1.0+nu)*((1.0-nu)*E3-2.0*nu3*nu3*E);

  const double C11 = (E*E3-nu3*nu3*E*E)/Delta;
  const double C12 = (nu*E*E3-nu3*nu3*E*E)/Delta;
  const double C13 = (1.0+nu)*nu3*E*E3/Delta;
  const double C33 = (1.0-nu*nu)*E3*E3/Delta;
  const double C66 = (C11-C12);
  const double C44 = 2.0*G3;

  printf("Delta %8.4e, C11 %8.4e, C33 %8.4e\n",Delta, C11, C33);

  stiffness(0,0,0,0) = C11;
  stiffness(1,1,1,1) = C11;
  stiffness(2,2,2,2) = C33;

  stiffness(0,0,1,1) = C12;
  stiffness(1,1,0,0) = C12;

  stiffness(0,0,2,2) = C13;
  stiffness(2,2,0,0) = C13;

  stiffness(1,1,2,2) = C13;
  stiffness(2,2,1,1) = C13;

  stiffness(0,1,0,1) = C44;
  stiffness(1,0,1,0) = C44;

  stiffness(0,2,0,2) = C44;
  stiffness(2,0,2,0) = C44;

  stiffness(1,2,1,2) = C66;
  stiffness(2,1,2,1) = C66;


  return;
}
void xElastic::SetElasticComplianceTransverselyIsotropic(double E, double E3, double G3,
							  double nu, double nu3, 
							   xtensor::xTensor4<>& compliance) {

//see http://www.engr.utk.edu/~cmc/528/chapter1/sec4.html

// Z (3) is the fiber direction
// X and Y (1 and 2 ) are transverse direction


  int i, j, k, l;  
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  compliance(i,j,k,l) = 0.0;
	}
      }
    }
  }

//  const double nu31 = nu13*E3/E;
//  const double nu32 = nu13*E3/E;
  const double G  = E/2.0/(1.0+nu);

  compliance(0,0,0,0) = 1.0/E;
  compliance(1,1,1,1) = 1.0/E;
  compliance(2,2,2,2) = 1.0/E3;

  compliance(0,0,1,1) = -nu/E;
  compliance(1,1,0,0) = -nu/E;

  compliance(0,0,2,2) = -nu3/E3;
  compliance(2,2,0,0) = -nu3/E3;

  compliance(1,1,2,2) = -nu3/E3;
  compliance(2,2,1,1) = -nu3/E3;

  compliance(0,1,0,1) = 0.5/G3;
  compliance(1,0,1,0) = 0.5/G3;

  compliance(0,2,0,2) = 0.5/G3;
  compliance(2,0,2,0) = 0.5/G3;

  compliance(1,2,1,2) = 0.5/G;
  compliance(2,1,2,1) = 0.5/G;


  return;
}

  xEvalThermoElasticVariables:: xEvalThermoElasticVariables( const xEval<double > &_eval_temperature,  const xEval<xtensor::xTensor2<> > &_eval_strain)
 :eval_temperature(_eval_temperature), eval_strain(_eval_strain){}
 
  void xEvalThermoElasticVariables:: operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xTensors& tensors) const{
    xtensor::xTensor2<> strain;
    eval_strain( geo_appro,  geo_integ, strain);
    double T;
    eval_temperature( geo_appro,  geo_integ, T);
    tensors.scalar("temperature") = T;
    tensors.tensor2("strain") =strain ;
  }


  xThermoElastic::xThermoElastic(){
  values_signature.register_tensor2("stress");
  values_signature.register_tensor2("elastic_strain");
  
  values.setSignature(&values_signature);
  
  //variables_signature.register_tensor2("stress");
  variables_signature.register_tensor2("strain"); 
  variables_signature.register_scalar("temperature"); 
  variables.setSignature(&variables_signature);
  
  properties_signature.register_string("MATERIAL_CLASS");
  properties_signature.register_string("NAME");
  properties_signature.register_scalar("YOUNG_MODULUS");
  properties_signature.register_scalar("POISSON_RATIO");
  properties_signature.register_scalar("DENSITY");
  properties_signature.register_scalar("THERMIC_DILATATION");
  properties_signature.register_scalar("THERMIC_CONDUCTIVITY");

  properties_signature.register_scalar("T0");
  properties.setSignature(&properties_signature);

  //  std::cout << "set ////// " << std::endl;
  properties.astring("MATERIAL_CLASS") = "MATERIAL_THERMO_ELASTIC";
  return;   
} 
 
void xThermoElastic::checkProperties(){
  E =  properties.scalar("YOUNG_MODULUS");
  nu = properties.scalar("POISSON_RATIO");
  k =         properties.scalar("THERMIC_CONDUCTIVITY");
  setThermicConductivityIsotropic( k );
  lam = nu*E/((1.+nu)*(1.-2.*nu));
  mu  = E/(2.*(1.+nu));
  setElasticStiffnessIsotropicFromLame (elastic_stiffness);
  alpha  =    properties.scalar("THERMIC_DILATATION");
  T0 =    properties.scalar("T0");
  
  
  return;
}

void xThermoElastic::computeCurrentState(){
  if (!uptodate){
    // std::cout << " compoute  CS" << std::endl;
    xtensor::xTensor2<> strain (variables.tensor2("strain") );
    //    std::cout << strain << std::endl;
    const double T = variables.scalar("temperature") ;
    for (int i=0; i<3;++i){
       strain(i,i) -= alpha*(T-T0);
    }
    // std::cout << T <<" "  <<T0 <<" " << alpha << std::endl;
    
    xtensor::xTensor2<> stress( strain *(2.*mu));
    
    double lamtraceeps = lam*strain.trace();// + E*alpha/(1.-2.*nu)*(T-T0);
    for (int i=0; i< 3;++i) stress(i,i) +=  lamtraceeps;
    values.tensor2("stress") = stress;
    values.tensor2("elastic_strain") = strain;
    
    //std::cout << stress << std::endl;
    // std::cout << std::endl;
    uptodate =true;
  }
  return;
}

void xThermoElastic::setElasticStiffnessIsotropicFromLame( xtensor::xTensor4<>& stiffness) {
  int i, j, k, l;
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = lam * delta(i,j) * delta(k,l) 
	    + mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }
  return;
}

void xThermoElastic::SetElasticComplianceIsotropic(double E, double nu, 
					       xtensor::xTensor4<>& compliance) const {

  const double lam1 = -nu/E;
  const double mu1  = (1.+nu)/(2.*E);
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      for (int k = 0; k < 3; ++k){
	for (int l = 0; l < 3; ++l){
	  compliance(i,j,k,l) = lam1 * delta(i,j) * delta(k,l) 
	    + mu1 * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k)); 
	}
      }
    }
  }
  return;
}


void xThermoElastic::sensitivityTo(const std::string& phys_token, xtensor::xTensor4Isotropic& sensitivity){
  const bool debug = xdebug_flag;
  if (debug) cout << "inside  xThermoElastic::sensitivityTo" << endl;
  if (phys_token == "strain")      { sensitivity.lam = lam;  sensitivity.mu = mu; }
  //else if (phys_token == "stress")  sensitivity = elastic_compliance;
  else { 
    std::cerr << "No sensitivity coded for phys token "<< phys_token << std::endl; 
    assert(0); 
  }
  return; 
}

void xThermoElastic::sensitivityTo(const std::string& phys_token, xtensor::xTensor4<>& sensitivity){
  const bool debug = xdebug_flag;
  if (debug) cout << "inside  xThermoElastic::sensitivityTo" << endl;
  if (phys_token == "strain")      {sensitivity = elastic_stiffness; }
  
  else if (phys_token == "stress") {
    SetElasticComplianceIsotropic(E, nu, sensitivity);
  }
  else { 
    std::cerr << "No sensitivity coded for phys token "<< phys_token << std::endl; 
    assert(0); 
  }
  return; 
}


void xThermoElastic::sensitivityTo(const std::string& phys_token, double & sensitivity){
  const bool debug = xdebug_flag;
  if (debug) cout << "inside  xThermoElastic::sensitivityTo" << endl;
  if (phys_token == "temperature_gradient")      { sensitivity = k ;}
  else if (phys_token == "temperature_gradient_of_volume_forces")      { sensitivity = -alpha*E/(1-2*nu) ;}
  else if (phys_token == "temperature_jump_of_surface_tension")      { sensitivity = alpha*E/(1-2*nu) ;}
  else if (phys_token == "temperature_of_thermic_dilatation")      { sensitivity = alpha ;}
  else { 
    std::cerr << "No sensitivity coded for phys token "<< phys_token << std::endl; 
    assert(0); 
  }
  return; 
}

void  xThermoElastic::sensitivityTo(const std::string& phys_token, 
				  xtensor::xTensor2<>& sensitivity) {
  if (phys_token == "temperature_gradient") 
    {
      sensitivity = thermic_conductivity;
      return;
    }
  if (phys_token == "temperature_of_thermic_dilatation"){
    sensitivity = xtensor::xTensor2<>(0.0);
    for (int i=0; i<3; ++i) {sensitivity(i,i) = alpha;}
    return;
  } 
  else { 
    std::cerr << "No sensitivity coded for phys token "<< phys_token << std::endl; 
    assert(0); 
  }
  return;
}

void xThermoElastic::setThermicConductivityIsotropic(double k ){
  thermic_conductivity = xtensor::xTensor2<>(0.);
  for (int i = 0; i < 3; ++i){
      thermic_conductivity(i,i) = k;
  }
  return;
}

void xThermoElasticTemperatureDependant::sensitivityTo(const std::string& phys_token, xtensor::xTensor4Isotropic& sensitivity){
  const bool debug = xdebug_flag;
  if (debug) cout << "inside  xThermoElastic::sensitivityTo" << endl;
  if (phys_token == "strain")      { sensitivity.lam = lam;  sensitivity.mu = mu; }
  //else if (phys_token == "stress")  sensitivity = elastic_compliance;
  else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
  return; 
}

void xThermoElasticTemperatureDependant::sensitivityTo(const std::string& phys_token, xtensor::xTensor4<>& sensitivity){
  const bool debug = xdebug_flag;
  if (debug) cout << "inside  xThermoElastic::sensitivityTo" << endl;
  if (phys_token == "strain")      {sensitivity = elastic_stiffness; }
  
  //else if (phys_token == "stress")  sensitivity = elastic_compliance;
  else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
  return; 
}

void xThermoElasticTemperatureDependant::sensitivityTo(const std::string& phys_token, double & sensitivity){
  const bool debug = xdebug_flag;
  if (debug) cout << "inside  xThermoElastic::sensitivityTo" << endl;
  if (phys_token == "temperature_gradient")      {
    sensitivity = k ;
    // std::cout << sensitivity << std::endl;
  }
  else if (phys_token == "temperature_gradient_of_volume_forces")      { 
    sensitivity = -alpha*E/(1-2*nu) ;
    // std::cout << "tgv " << sensitivity << std::endl;
  }
  else if (phys_token == "temperature_jump_of_surface_tension")      {
    sensitivity = alpha*E/(1-2*nu) ;
    // std::cout << "tgs " << sensitivity << std::endl;
  }
  else  { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
  return; 
}

void  xThermoElasticTemperatureDependant::sensitivityTo(const std::string& phys_token, 
				  xtensor::xTensor2<>& sensitivity) {
  if (phys_token == "temperature_gradient") sensitivity = thermic_conductivity;
  else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
  // std::cout << sensitivity << std::endl;
  return;
}

void xThermoElasticTemperatureDependant::setThermicConductivityIsotropic(double k ){
  thermic_conductivity = xtensor::xTensor2<>(0.);
  for (int i = 0; i < 3; ++i){
      thermic_conductivity(i,i) = k;
  }
  return;
}

xThermoElasticTemperatureDependant::xThermoElasticTemperatureDependant (){
  //  std:: cout << "TD Const" << std::endl;
  variables_signature.register_tensor2("stress");
  variables_signature.register_tensor2("strain"); 
  
  properties_signature.register_string("MATERIAL_CLASS");
  properties_signature.register_string("NAME");
  properties_signature.register_piecewiselinear("YOUNG_MODULUS");
  properties_signature.register_piecewiselinear("POISSON_RATIO");
  //  properties_signature.register_piecewiselinear("DENSITY");
  properties_signature.register_piecewiselinear("THERMIC_DILATATION");
  properties_signature.register_scalar("THERMIC_CONDUCTIVITY");
 
 
 properties.setSignature(&properties_signature);
 properties.astring("MATERIAL_CLASS") = "MATERIAL_THERMO_ELASTIC_TEMP_DEP";
  // std::map<double, double > Ett;
  //Ett[0] = 3.;
  //Ett[50] = 2.;
  //Ett[100] = 1.;
  //Et = linearbypart(Ett);
}

void xThermoElasticTemperatureDependant::checkProperties(){
  Et = properties.piecewiselinear("YOUNG_MODULUS");
  nut = properties.piecewiselinear("POISSON_RATIO");
  alphat  =    properties.piecewiselinear("THERMIC_DILATATION");
  k =         properties.scalar("THERMIC_CONDUCTIVITY"); 
  setThermicConductivityIsotropic( k );
  //std::cout << "in check, k= " << k << std::endl;
  setTemperature(temperature);
  return;
}

void xThermoElasticTemperatureDependant::setTemperature(double T){
  if (temperature!=T){
    temperature = T;
    E =  Et(temperature);
    nu = nut(temperature);
    alpha  = alphat(temperature);
    //std::cout << "in set T =" << k << std::endl;
    setThermicConductivityIsotropic( k );
    setElasticStiffnessIsotropicFromLame(elastic_stiffness);
    lam = nu*E/((1.+nu)*(1.-2.*nu));
    mu  = E/(2.*(1.+nu)); 
  }   
}

void xThermoElasticTemperatureDependant::setElasticStiffnessIsotropicFromLame( xtensor::xTensor4<>& stiffness) {
  int i, j, k, l;
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = lam * delta(i,j) * delta(k,l) 
	    + mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }
  return;
}
  
} // end of namespace
