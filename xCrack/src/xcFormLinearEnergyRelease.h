/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.

*/
  
#ifndef _xcformelinearenergyrelease__H
#define _xcformelinearenergyrelease__H

#include "xForm.h"
#include "xValKey.h"
#include "xEval.h"
#include "xGeomElem.h"
#include "xTensor2.h"
#include "xOperators.h"
#include "xEnv.h"
#include "lCrack.h"
#include "xMaterialSensitivity.h"

using namespace std;
using namespace xfem;


namespace xcrack 
{

class xcCrackBase;

class xcEvalEshelbyHPP : public xEval<xtensor::xTensor2<> > 
  {
  public:
    xcEvalEshelbyHPP(const xEval<xtensor::xTensor2<> >& g,
		     const xEval<xtensor::xTensor4<> >& p);
    
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const override;
  private:
    const xEval<xtensor::xTensor2<> >& eval_grad_disp;
    const xEval<xtensor::xTensor4<> >& eval_phys;
    mutable xtensor::xTensor2<> grad_disp;
    mutable xtensor::xTensor4<> phys;
  };


  class xcEvalEshelbyHPPThermoElastic : public xEval<xtensor::xTensor2<> > 
  {
  public:
    xcEvalEshelbyHPPThermoElastic(const xEval<xtensor::xTensor2<> >& g, const xEvalMaterialValues & _eval_meca_fields);
    
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const override;
  private:
    const xEval<xtensor::xTensor2<> >& eval_grad_disp;
    const xEvalMaterialValues& eval_meca_fields;
    mutable xtensor::xTensor2<> grad_disp;
    mutable xtensor::xTensor2<> stress;
  };
 
 class xcEvalEshelbyInteractionHPPThermoElastic : public xEval<xtensor::xTensor2<> > 
  {
  public:
    xcEvalEshelbyInteractionHPPThermoElastic(const xEval<xtensor::xTensor2<> >& g, const xEvalMaterialValues & _eval_meca_fields, const xcCrackBase & c, int mode =0);
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const override;
    void setmode(int m) const    {mode =m;};

  private:
    const xEval<xtensor::xTensor2<> >& eval_grad_disp;
    const xEvalMaterialValues & eval_meca_fields;
    mutable int mode;
    const xcCrackBase& crack; 
    mutable xtensor::xTensor2<> grad_disp;
    mutable xtensor::xTensor2<> grad_disp_aux, stress_aux;

  };

/*
 class xcEvalEshelbyHPPThermo : public xEval<xtensor::xTensor2<> > 
  {
  public:
    xcEvalEshelbyHPPThermo(const xEval<xtensor::xTensor2<> >& g, const xEval<double >& temperature,
			   const xEval<xtensor::xTensor4<> >& p, const xEval<double > &alpha, double _T0);
    
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const;
  private:
    const xEval<xtensor::xTensor2<> >& eval_grad_disp;
    const xEval<double >& eval_temperature; 
    const xEval<xtensor::xTensor4<> >& eval_phys;
    const xEval<double >&  eval_alpha;
    double T0;
  };
*/
  
  class xcEvalJintSourceTerm : public xEval<xtensor::xVector<> >{
  public :
  xcEvalJintSourceTerm():mode (0){};
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
			    result_type& result) const override  {result = xtensor::xVector<>(0.); };
    void setAuxilliaryFieldMode(int _mode ) const {
      mode = _mode;
    }
    mutable int mode;
  };
  
  class xcEvalJintHPPThermoSourceTerm : public xcEvalJintSourceTerm
  {
  public:
    xcEvalJintHPPThermoSourceTerm (const xEval<xTensors>& eval_meca,
				   const xEval<xtensor::xVector<> >& temperature_gradient_eval,
				   const xEval<double  >& alpha,
				   const xcCrackBase &_crack
				   );
    
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const override;
  private:
    const xEval<xTensors>& eval_meca;
    const xEval<xtensor::xVector<> >&  eval_grad_temperature;
    const xEval<double >&  eval_alpha;
    const xcCrackBase& crack;
    mutable xtensor::xTensor2<> stress, stress_aux,  grad_disp_aux;
    mutable xtensor::xVector<> gT;
    mutable double alpha;
  };

  

class xcEvalEshelbyHPPInteraction : public xEval<xtensor::xTensor2<> >
  {
  public:
    xcEvalEshelbyHPPInteraction(const xEval<xtensor::xTensor2<> >& g,
				const xEval<xtensor::xTensor4<> >& p,
				const xcCrackBase& c, int m=1);
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const override;
    void setmode(int m) const    {mode =m;};

  private:
    const xEval<xtensor::xTensor2<> >& eval_grad_disp;
    const xEval<xtensor::xTensor4<> >& eval_phys;
    const xcCrackBase& crack; 
    mutable int mode;
    mutable xtensor::xTensor2<> grad_disp;
    mutable xtensor::xTensor4<> phys;
    mutable xtensor::xTensor2<> grad_disp_aux, stress_aux;
  };
  
class xcEvalEshelbyHPPInteractionv2 : public xEval< xtensor::xTensor2<> >
  {
  public:
    xcEvalEshelbyHPPInteractionv2(const xEval<xtensor::xTensor2<> >& s,
				const xEval<xtensor::xTensor2<> >& gd,
				const xcCrackBase& c, int m);
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const override;

       
  private:
    const xEval<xtensor::xTensor2<> >& eval_stress;
    const xEval<xtensor::xTensor2<> >& eval_grad_disp;
    const xcCrackBase& crack; 
    int mode;
    mutable xtensor::xTensor2<> grad_disp, stress;
    mutable xtensor::xTensor4<> phys;
    mutable xtensor::xTensor2<> grad_disp_aux, stress_aux;
  };

 class xcEvalJintTerms{
 public:
   xcEvalJintTerms():mode(1){};
   virtual ~xcEvalJintTerms()= default;;
   virtual const xEval<xtensor::xTensor2<> > &get_eshelby_evaluator() const =0;
   virtual const xEval<xtensor::xTensor2<> >& get_interaction_evaluator()const=0;
   virtual const xEval<xtensor::xVector<> >  &get_source_evaluator()const =0;
   virtual const xEval<xtensor::xVector<> > & get_source_interaction_evaluator()const= 0;
   void setmode(int _mode)  const {mode =_mode;};
   mutable int mode;
 };
 
 class xcEvalJintTermElasticity: public xcEvalJintTerms{
 public :
 xcEvalJintTermElasticity(const xEval<xtensor::xTensor2<> >& gd, const xEval<xtensor::xTensor4<> >& hooke, const xcCrackBase& crack ):eshelby_evaluator(gd, hooke), interaction_evaluator(gd, hooke, crack),  source_evaluator(), source_interaction_evaluator()
     { };
   
   const xEval<xtensor::xTensor2<> > &get_eshelby_evaluator() const override{
     return eshelby_evaluator;
   };
   
   const xEval<xtensor::xTensor2<> > &get_interaction_evaluator()const override {
     interaction_evaluator.setmode(mode);
     return interaction_evaluator;
   }
   const xEval<xtensor::xVector<> >  &get_source_evaluator() const override {return source_evaluator;} ;
   const xEval<xtensor::xVector<> >  &get_source_interaction_evaluator() const override  {return source_interaction_evaluator;};

   xcEvalEshelbyHPP eshelby_evaluator;
   xcEvalEshelbyHPPInteraction interaction_evaluator;
   xEvalZero <xtensor::xVector<> > source_evaluator ;
   xEvalZero <xtensor::xVector<> > source_interaction_evaluator;
 } ;

 class xcEvalJintTermThermoElasticity: public xcEvalJintTerms{
 public :
 xcEvalJintTermThermoElasticity(const xEval<xtensor::xTensor2<> >& gd, const xEval<xtensor::xVector<> >& gt , const xEval<double >& alpha ,   const xEvalMaterialValues & eval_meca , const xcCrackBase& crack ):eshelby_evaluator(gd, eval_meca),
     interaction_evaluator(gd, eval_meca ,  crack),
     source_evaluator(eval_meca, gt, alpha ,  crack)
     { };
   
   const xEval<xtensor::xTensor2<> > &get_eshelby_evaluator() const override{
     return eshelby_evaluator;
   };
   
   const xEval<xtensor::xTensor2<> > &get_interaction_evaluator()const override {
     interaction_evaluator.setmode(mode);
     return interaction_evaluator;
   }
   const xEval<xtensor::xVector<> >  &get_source_evaluator() const override {
     source_evaluator.setAuxilliaryFieldMode(0);
     return source_evaluator;
   } ;
   const xEval<xtensor::xVector<> >  &get_source_interaction_evaluator() const override  {
     source_evaluator.setAuxilliaryFieldMode(mode);
     return source_evaluator;
   };

   xcEvalEshelbyHPPThermoElastic            eshelby_evaluator;
   xcEvalEshelbyInteractionHPPThermoElastic interaction_evaluator;
   //xEvalConstant <xtensor::xVector<> > source_evaluator ;
   //xEvalConstant <xtensor::xVector<> > source_interaction_evaluator;
   xcEvalJintHPPThermoSourceTerm source_evaluator ;
 } ;
  
  
  
  /*class xcEvalEshelbyHPPInteractionThermoSource : public xEval<xtensor::xVector<> >
  {
  public:
    xcEvalEshelbyHPPInteractionThermoSource( const xEval<xtensor::xVector<> >& tg,
					     const xEval<double >& alpha,
				             const lCrack& c, int m);
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const;

       
  private: 
    const xEval<xtensor::xVector<> >& eval_temperature_gradient;
    const xEval<double >& eval_alpha;
    const lCrack& crack; 
    int mode;
    mutable xtensor::xTensor2<> grad_disp_aux, stress_aux;
  };
  */


  /*
  class xcEvalJintHPPThermoSourceTerm : public xEval<xtensor::xVector<> > 
  {
  public:
    xcEvalJintHPPThermoSourceTerm (const xEval<xtensor::xTensor2<> >& stress_eval,
				   const xEval<xtensor::xVector<> >& temperature_gradient_eval,
				   const xEval<double  >& alpha);
    
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const;
  private:
    const xEval<xtensor::xTensor2<> >& eval_stress;
    const xEval<xtensor::xVector<> >&  eval_grad_temperature;
    const xEval<double >&  eval_alpha;

  };
  */


  /*
class xcEvalEshelbyHPPThermoInteraction : public xEval<xtensor::xTensor2<> > 
  {
  public:
    xcEvalEshelbyHPPThermoInteraction(const xEval<xtensor::xTensor2<> >& g,
				const xEval<xtensor::xTensor4<> >& p,
				const lCrack& c, int m);
    
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, 
		    result_type& result) const;

    void getAuxFields(xfem::xGeomElem* geo_appro, 
		      const lCrack& crack, 
		      int mode, xtensor::xTensor2<>& stress_aux, 
		      xtensor::xTensor2<>& grad_disp_aux) const;
    
  private:
    const xEval<xtensor::xTensor2<> >& eval_grad_disp;
    const xEval<xtensor::xTensor4<> >& eval_phys;
    const lCrack& crack; 
    int mode;
    mutable xtensor::xTensor2<> grad_disp;
    mutable xtensor::xTensor4<> phys;
    mutable xtensor::xTensor2<> grad_disp_aux, stress_aux;
  };


  */


class xcFormLinearEnergyRelease : public  xFormLinear<double>
  {
  private:
    const xGradOperator<xtool::xIdentity<xtensor::xVector<> > > GradTest;
    const xValOperator<xtool::xIdentity<double> >    ValTest;
    
    const xEval<xtensor::xVector<> >&    eval_qdir;
    const xEval<xtensor::xTensor2<> >&   eval_grad_qdir;
    const xEval<xtensor::xTensor2<> >&   eval_eshelby;
    const xEval<xtensor::xVector<> >&    eval_source;
    
    xtensor::xTensor2<> eshelby;
    xtensor::xVector<>  qdir;
    xtensor::xTensor2<> grad_qdir;

  public:
    
    xcFormLinearEnergyRelease (const xEval<xtensor::xVector<> >&  q,
			       const xEval<xtensor::xTensor2<> >& g, // a remplacer par un scalar_field
			       const xEval<xtensor::xTensor2<> >& e,
			       const xEval<xtensor::xVector<> > & s
			       );
    
    void accumulate_pnt(xfem::xGeomElem*  geo_integ) override;
 
  };


class xcFormLinearEnergyReleaseBoundary : public  xFormLinear<double>
  {
  private:
    const xValOperator<xtool::xIdentity<double> >    ValTest;
    
    const xEval<xtensor::xVector<> >&    eval_qdir;
    const xEval<xtensor::xTensor2<> >&   eval_eshelby;
    xEvalNormal eval_normal;
    xtensor::xTensor2<> eshelby;
    xtensor::xVector<>  qdir, normal;
    
  public:
    
    xcFormLinearEnergyReleaseBoundary (const xEval<xtensor::xVector<> >&  q,
				       const xEval<xtensor::xTensor2<> >& e);
    
    void accumulate_pnt(xfem::xGeomElem*  geo_integ) override;
 
  };

}


#endif
