/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _XOPERATOR_H_H
#define _XOPERATOR_H_H

#include "mIterator.h"
#include "xLevelSet.h"
#include "xVectorField.h"
#include "xRegion.h"
#include "xElement.h"
#include "xEntityFilter.h"


namespace xfem
{

class xLevelSetModifier {
public:
  virtual void visit(xLevelSet& , xRegion target)  = 0;
  virtual ~xLevelSetModifier() = default;
};

class xLevelSetInspector {
public: 
  virtual void visit(const xLevelSet& , xRegion target) = 0;
  virtual ~xLevelSetInspector() = default;

};

/// This class visit a const xLevelSet and is made to fill a xLevelSet&
class xLevelSetCreator {
public: 
  virtual void visit(const xLevelSet&, xLevelSet&) = 0;
  virtual ~xLevelSetCreator() = default;
};

/// levelset inspector (inspector pattern) that find max and min value of the levelset.
 class xGetMinMaxLevelSetInspector : public xLevelSetInspector {
  public:
    void visit(const xLevelSet& ls, xRegion target) override;
    double getMin() const;
    double getMax() const;
  private:
    double min_ls, max_ls;
};
 
/// levelset modifier (visitor pattern) that shift (add to ) the levelset by the constructor given double value.
 class xShiftLevelSetModifier : public xLevelSetModifier {
 public:
  xShiftLevelSetModifier(const double& s );
  void visit(xLevelSet& ls, xRegion target) override;
 private:
  const double shift;
};

/// Created a "shifted" levelset. 
/*!
  Upon  visit, the created levelset get value of the input levelset + shift.
  !*/
 class  xShiftLevelSetCreator:public xLevelSetCreator {
  public: 
  xShiftLevelSetCreator(double shift);
  void visit(const xLevelSet& in, xLevelSet& out) override;
 private:
  const double shift;
};
 
 
/// Create a levelset which is the complement of the input one (out = -in)
 class  xComplementLevelSetCreator:public xLevelSetCreator {
 public: 
  xComplementLevelSetCreator();
  void visit(const xLevelSet& in, xLevelSet& out) override;
};
 


/// levelset modifier (visitor pattern) that scale (multiply )the levelset by the constructor given double value.
 class xScaleLevelSetModifier : public xLevelSetModifier {
 public:
  xScaleLevelSetModifier(const double& s );
  void visit(xLevelSet& ls, xRegion target) override;
 private:
  const double scale;
}; 

 

/// This class is exactly the same as xLevelSetCreator and is probably here for historical reasons (or to confuse peoples).
class xSpaceOperator : public xLevelSetCreator {
public: 
  void visit(const xLevelSet&, xLevelSet&) override = 0;
  ~xSpaceOperator() override = default;
};




//mets a zero les vitesse sur le maillage sauf sur
//un submesh ou les valeurs sont calculees sur base
//du vertex le plus proche sur le front
class xVelocityInitialization : public xLevelSetModifier {

public:
  xVelocityInitialization(const xLevelSet& vitfront, xRegion nodes_loc) : 
    velofront(vitfront), regloc(nodes_loc) {}
  void visit(xLevelSet& velo, xRegion target) override;

private:
  const xLevelSet& velofront;
  xRegion  regloc;

};
 
class xUniformValue : public xLevelSetModifier {

public:
  xUniformValue(const double& v) : val(v) {}
  void visit(xLevelSet& f, xRegion target) override;

private: 
  xUniformValue():val(0) {}
  const double val;

};

//for a custom function
class xNonUniformValue : public xLevelSetModifier {

public:
  xNonUniformValue(const double& v) /*: val(v)*/ {}
  void visit(xLevelSet& f, xRegion target) override;

private: 
  xNonUniformValue()/*:val(0)*/ {}
  //const double val;

};



  /// merge with the LS given in the constructor
  class xInterLevelSetModifier : public xLevelSetModifier 
  {
  public:
    xInterLevelSetModifier(const xLevelSet& ls1) : otherLS(ls1) {}
    void visit(xLevelSet& f, xRegion target) override{
      for(xIter it = target.begin(0); it != target.end(0); ++it) {
	AOMD::mVertex *v = (AOMD::mVertex*) *it;
        f(v) = std::min(f(v),otherLS(v));
      }
    }
  private: 
    const xLevelSet otherLS;
  };




//Visitors defined by means of a lambda functions
  
  //Here, the visit is done in the xLambdaModifierOnVertices (on vertices), and the modification is handled by the lambda function
class xLambdaModifierOnVertices : public xLevelSetModifier {
public:
  xLambdaModifierOnVertices(const std::function<double(double)> func_) : func(func_) {}
  void visit(xLevelSet& f, xRegion target) override;

private: 
  const std::function<double(double)> func;
};

//Here, the lambda function takes care of the visit (on the vertices, edges, elements...) AND ALSO of the modification of the value.
class xLambdaModifier : public xLevelSetModifier {
public:
  xLambdaModifier(const std::function<void(xLevelSet& f, xRegion& target)> func_) : func(func_) {}
  void visit(xLevelSet& f, xRegion target) override{
      func(f, target);
      return;
  }

private: 
  const std::function<void(xLevelSet& f, xRegion& target)> func;
};

  
  

class xPilotToStationary {

public:
  virtual bool   itersCompleted(const xLevelSet& f, xRegion target) = 0;
  virtual double getVirtualTimeStep() = 0;
  virtual void   initialize() = 0;
  virtual int    getIter() = 0;
  virtual double getLastError()=0 ;
  virtual ~xPilotToStationary() = default;
private: 

};
 class xComputeError;
class xPilotError : public xPilotToStationary {
public:
  xPilotError(xComputeError& c, const double& t, int i, const double& d) : 
   iter(0), epsrel(t), epsabs(t), itmax(i),  dt(d), compute_error(c) {}
  void   initialize() override {iter = 0;}
  bool   itersCompleted(const xLevelSet& f, xRegion target) override;
  double getVirtualTimeStep() override {return dt;};
  int    getIter() override {return iter;};
  double getLastError() override {return error;};
private:
  int    iter;
  double epsrel, epsabs;
  int    itmax;
  double dt;
  xComputeError& compute_error;
  double error;

};


class xComputeError : public xLevelSetInspector {

public: 
  void visit(const xLevelSet& current, xRegion target) override = 0;
  void loadReference(const xLevelSet& ref);
  std::pair<double, double > getError() const;
  ~xComputeError() override = default;
protected:
  xLevelSet reference;
  double errnorm;
  double refnorm;
};

class xL2norm : public xComputeError {

public: 
  void visit(const xLevelSet& current, xRegion target) override;
private:
};





class xEvolveToStationary : public xLevelSetModifier  {

public:

  xEvolveToStationary(xSpaceOperator& s, xPilotToStationary& p)
   : pilot(p) 
   {operators.push_back(&s); }
  xEvolveToStationary(std::list<xSpaceOperator*>& s, 
		      xPilotToStationary& p) : pilot(p), operators(s)  {}
  void visit(xLevelSet& f, xRegion target) override;

private:
  xPilotToStationary& pilot;
  std::list<xSpaceOperator*> operators;
  //xRegion target;
  
};

template <class TimeOp >
class xEvolveToStationaryT : public xLevelSetModifier  {

 public:

 xEvolveToStationaryT(xSpaceOperator& s, xPilotToStationary& p)
   : pilot(p) 
   {operators.push_back(&s); }
 xEvolveToStationaryT(std::list<xSpaceOperator*>& s, 
		      xPilotToStationary& p) : pilot(p), operators(s)  {}
   void visit(xLevelSet& f, xRegion target) override {
     pilot.initialize();
     while (!pilot.itersCompleted(f, target) ) {
       double dt = pilot.getVirtualTimeStep();
       std::list<xSpaceOperator*>::iterator it;
       for (it = operators.begin(); it !=operators.end(); it++) {
	 typename TimeOp::operator_type * ope =	(typename TimeOp::operator_type *) *it;
	 TimeOp step(dt, *ope);
	 f.accept(step,target);
       }
     }
     std::cout << "iterations " << pilot.getIter() << " error is " << pilot.getLastError() << std::endl;
   }
   
 private:
   xPilotToStationary& pilot;
   std::list<xSpaceOperator*> operators;
   //xRegion target;
   
 };




class xPilotTimeIntegration {

public:
  virtual bool   stepsCompleted() = 0;
  virtual double getTimeStep() = 0;
  virtual void   initialize() = 0;
  virtual ~xPilotTimeIntegration()= default;

};


class xPilotOneStep : public xPilotTimeIntegration {

public:
  xPilotOneStep(const double& d)  : dt(d)  {}
  void  initialize() override {step = 0;}
  bool  stepsCompleted() override;
  double getTimeStep() override {return dt;}
   
private:
  int step;
  double dt;

};


class xTimeIntegration : public xLevelSetModifier  {
public:
  xTimeIntegration(xSpaceOperator& s, xPilotTimeIntegration& p) 
    : pilot(p), soperator(s) {}
  void visit(xLevelSet& f, xRegion target) override;
private:
  xPilotTimeIntegration& pilot;
  xSpaceOperator&        soperator;
};



class xRungeKutta : public xLevelSetModifier {
  
 public:
 xRungeKutta(const double& t, xSpaceOperator& s) 
   : dt(t), space_operator(s) {}  
  void visit(xLevelSet& ls, xRegion target) override;
  typedef xSpaceOperator operator_type;
private: 
  xRungeKutta();
  double dt;
  xSpaceOperator& space_operator;
};

class xForwardEuler : public xLevelSetModifier {
  
 public:
 xForwardEuler(const double& t, xSpaceOperator& s) 
   : dt(t), space_operator(s) {}  
  void visit(xLevelSet& ls, xRegion target) override;
  typedef xSpaceOperator operator_type;
private: 
  xForwardEuler();
  double dt;
  xSpaceOperator& space_operator;
};
 
template <class T >
class xGaussSeidel : public xLevelSetModifier {
public:
  xGaussSeidel(const double& t, T& s) : dt(t), space_operator(s) {}  
  void visit(xLevelSet& ls, xRegion target) override{
    std::list <AOMD::mVertex * > listV;
    for(xIter it = target.begin(0); it != target.end(0); ++it) {
      AOMD::mVertex *v = (AOMD::mVertex*) *it;
      listV.push_back(v);
      double H= space_operator.evalH(ls, v);
      ls(v) += dt * H; 
    }
    for(std::list<AOMD::mVertex * >::reverse_iterator it = listV.rbegin(); it != listV.rend(); ++it) {
      AOMD::mVertex *v = (AOMD::mVertex*) *it;
      double H= space_operator.evalH(ls, v);
      ls(v) += dt * H; 
    }
  }
  typedef T operator_type;
private: 
  xGaussSeidel();
  double dt;
  T& space_operator;
};

// class xParallelLevelSetUnifier :public xLevelSetModifier{
//  public :
//   xParallelLevelSetUnifier(){};
//   void visit(xLevelSet& ls, xRegion target) override;
//   // Prerequist of DataExchange:
//   typedef double data_type;
//   void send(AOMD::mEntity *, data_type &) const;
//   void receive(AOMD::mEntity *, const std::vector<data_type > &) const;
//  private:
//   xLevelSet *ls;
// };

class xFitToVertices : public xLevelSetModifier
{
public: 
  xFitToVertices(const double& t, const double& f = 0.0) : TOL(t), fit_value(f) {}
  void visit(xLevelSet& ls, xRegion target) override; 

private: 
  const double TOL;
  const double fit_value;
  xFitToVertices();

  // private class for communication
  class xFitToVertexKeyManager
  {
      public:
        typedef const AOMD::mEntity * information_key_t;
        xFitToVertexKeyManager ( const partmanAOMD_t & _partman );
        information_key_t localObjectKey( const AOMD::mEntity  & o);
        std::set < int > getMessageRanks( const information_key_t & lk);
      private:
        const partmanAOMD_t &partman;
  };
  class xFitToVertexInfoManager
  {
      public:
        typedef xtool::homogeneous_data_style_trait data_style_trait;
        typedef xtool::send_only_keys_communication_trait communication_trait;
        typedef const AOMD::mEntity * information_key_t;
        typedef AOMD::mEntity * information_t;
        xFitToVertexInfoManager (xLevelSet& ls_, const partmanAOMD_t & partman_ ,double fit_value_);
        information_t getInfo(information_key_t key, int sendto);
        void setInfo(const std::vector < information_t > &info, int receivedfrom );
      private:
        const partmanAOMD_t &partman;
        xLevelSet& ls;
        const double fit_value;
  };
};



class xHJOperator // operateur donnant la valeur de F (vitesse) ou f (mb de droite) et du gradient
{
  public :
  virtual double operator()(AOMD::mVertex *v) { return 0.0;}
  virtual xtensor::xVector<> getGrad(AOMD::mEntity *e) { return xtensor::xVector<>(0.0,0.0,0.0);}
  virtual ~xHJOperator()= default;
};

class xCFL  {
public:
  static double getLength(xRegion);
  static double getDt(const xLevelSet& l);
  static double getDt(xRegion);
private:
};


class xLsetHamiltonJacobi : public xSpaceOperator 
{
  public:
  template <class F_,class f_> xLsetHamiltonJacobi(F_ &Fpar, f_ &fpar,bool r=true)
  {
    Fop=new F_(Fpar);
    fop=new f_(fpar);
    rentrant=r;
  }
  template <class F_> xLsetHamiltonJacobi(F_ &Fpar,bool r=true)
  {
    Fop=new F_(Fpar);
    fop=new xHJOperator;
    rentrant=r;
  }

  ~xLsetHamiltonJacobi() override {delete Fop;delete fop;}

  void visit(const xLevelSet& ls, xLevelSet& H) override;
  double evalH(const xLevelSet& ls, AOMD::mVertex *v);
  private :
  xHJOperator *Fop;
  xHJOperator *fop;
  bool rentrant;
  struct starvalue{
    double w, ww, l, ll;
  };
  void byelemcomputation(const std::vector<double > &lsv, 
	            const xtensor::xVector<> &gradlvlset,
		    const double  &F, const double &f, 
		    const xElementInfo & info ,
		    std::vector<double > & maxalphatild, 
		    double  &sommaxalphatild,
		    double &somdeltafi) const;
  
};


class xExtensionSpeed : public xHJOperator // pour l'extension
{
  public :
  xExtensionSpeed(const xLevelSet &v) : ls(v){}
  double operator()(AOMD::mVertex *v) override
  {
    double signe=ls(v);if (signe>0)signe=1.0;else if (signe<0) signe=-1.0;
    return signe;
  }

  xtensor::xVector<> getGrad(AOMD::mEntity *e) override
  {
    return ls.getGrad(e);
  }
  private :
  const xLevelSet& ls;
};


class xExtensionSpeedLs : public xHJOperator // pour l'extension d'une ls dans sa partie positive
{
  public :
  xExtensionSpeedLs(const xLevelSet &v,const xVectorField &vv) : ls(v),vec(vv) 
  { 
    xRegion s = ls.getSupport();
    lc = xCFL::getLength(s);
  }
  double operator()(AOMD::mVertex *v) override
  {
    return si(v);
  }

  xtensor::xVector<> getGrad(AOMD::mEntity *e) override
  {
    xtensor::xVector<> xx(0.0,0.0,0.0);
    for (int j = 0; j < e->size(0); j++) 
    {
      AOMD::mVertex* v = (AOMD::mVertex*) e->get(0,j);
      xx+=vec(v);
    }
    return xx.norm();
  }
  private :
  const xLevelSet& ls;
  const xVectorField& vec;  
  double lc;  
  double si(AOMD::mVertex* v)
  {
    double decal=3*lc;
    double lsv=ls(v);
    double s=(lsv-decal)/sqrt(lsv*lsv+4*decal*decal);
    if (s<0) s=0.0;
    return s;
  } 
};


class xExtensionSpeedRadial : public xHJOperator // pour l'extension (rayon)
{
  public :
  xExtensionSpeedRadial(const xLevelSet &v,const xLevelSet &w) : lsn(v),lst(w){}
  double operator()(AOMD::mVertex *v) override
  {
    return 1;
  }

  xtensor::xVector<> getGrad(AOMD::mEntity *e) override
  {
    xtensor::xVector<> vn=lsn.getGrad(e);
    xtensor::xVector<> vt=lst.getGrad(e);
    double n=0.0,t=0.0;
    for (int j = 0; j < e->size(0); j++) 
    {
      AOMD::mVertex* v = (AOMD::mVertex*) e->get(0,j);
      n+=lsn(v);
      t+=lst(v);
    }
    //t/e->size(0);
    //n/e->size(0);
    double r=sqrt(t*t+n*n);
    return ((vt*t+vn*n)/r);
  }
  private :
  const xLevelSet& lsn;
  const xLevelSet& lst; 
};

class xPropagSpeed : public xHJOperator // pour la propag
{
  public :
  xPropagSpeed(const xLevelSet& v,const xLevelSet& l) : velo(v),ls(l) {}
  double operator()(AOMD::mVertex *v) override
  {
    return velo(v);
  }  

  xtensor::xVector<> getGrad(AOMD::mEntity *e) override
  {
    return ls.getGrad(e);
  }

  private :
  const xLevelSet& velo;
  const xLevelSet& ls;
};


class xReReSpeed : public xHJOperator // pour reinit et reortho : idem extension mais fonction signe lissee
{
  public :
  xReReSpeed(const xLevelSet& l) : ls(l)
  {
    xRegion s = ls.getSupport();
    lc = xCFL::getLength(s);
  }
  
  double operator()(AOMD::mVertex *v) override
  {
    return si(v);
  }  

  xtensor::xVector<> getGrad(AOMD::mEntity *e) override
  {
    return ls.getGrad(e);
  }

  private :
  const xLevelSet& ls;
  double lc;
  double si(AOMD::mVertex* v)
  {
    double mag=ls.getGrad(v).mag();
 
    if (ls(v)==0)
      return 0.;
    return ls(v)/sqrt(ls(v)*ls(v)+16*lc*lc*mag*mag);
  }
};



// number 1
class xExtensionRadialOperator :  public xSpaceOperator
{

public:
  xExtensionRadialOperator(const xLevelSet& n,const xLevelSet& t) : lsn(n),lst(t) {}
  void visit(const xLevelSet& velo, xLevelSet& H) override
  {
//    xExtensionSpeed extsp(ls);
//    xReReSpeed extsp(ls);
    xExtensionSpeedRadial extsp(lsn,lst);
    xLsetHamiltonJacobi HJ(extsp);
    HJ.visit(velo,H);
  }
private:
  xExtensionRadialOperator();
  const xLevelSet& lsn;
  const xLevelSet& lst;
};

// number 1bis
class xExtensionOperator :  public xSpaceOperator
{

public:
  xExtensionOperator(const xLevelSet& ln) : ls(ln){}
  void visit(const xLevelSet& velo, xLevelSet& H) override
  {
    xExtensionSpeed extsp(ls);
    xLsetHamiltonJacobi HJ(extsp,false);
    HJ.visit(velo,H);
  }
  double evalH(const xLevelSet& velo, AOMD::mVertex *v){
    xExtensionSpeed extsp(ls);
    xLsetHamiltonJacobi HJ(extsp,false);
    return HJ.evalH(velo, v);
  }
private:
  xExtensionOperator();
  const xLevelSet& ls;
};

//number 2
class xPropagationOperator : public xSpaceOperator
{

public:
  xPropagationOperator(const xLevelSet& v) : velo(v) {}
  void visit(const xLevelSet& ls, xLevelSet& H) override
  {
    xPropagSpeed propagvelo(velo,ls);
    xLsetHamiltonJacobi HJ(propagvelo);
    HJ.visit(ls,H);
  }

private:
  xPropagationOperator();
  const xLevelSet& velo;
};

// number 3
class xReInitOperator : public xSpaceOperator
{

public:
  xReInitOperator() {} 
  void visit(const xLevelSet& ls, xLevelSet& H) override
  {
    xReReSpeed reinitvelo(ls);
    xLsetHamiltonJacobi HJ(reinitvelo,reinitvelo,false); // pas de valeurs rentrantes -> grad nul
    HJ.visit(ls,H);
  }
private:

};

//number 4
class xOrthoOperator : public xSpaceOperator
{

public:
 xOrthoOperator(const xLevelSet& _lother):reinitvelo(_lother), HJ(reinitvelo, false){}
    // lother(l){}
  void visit(const xLevelSet& ls, xLevelSet& H) override
  {    
    HJ.visit(ls,H);
  }
  double evalH(const xLevelSet& ls, AOMD::mVertex *v){ 
    return HJ.evalH(ls, v);
  }
private:
  xOrthoOperator();
  //  const xLevelSet& lother;
  xReReSpeed reinitvelo;
  xLsetHamiltonJacobi HJ;


};

class xLSExtensionOperator :  public xSpaceOperator
{

public:
  xLSExtensionOperator(const xLevelSet & lstt,const xVectorField& vect) : lst(lstt),vec(vect) {}
  void visit(const xLevelSet& ls, xLevelSet& H) override
  {
    xExtensionSpeedLs extsp(lst,vec);
    xLsetHamiltonJacobi HJ(extsp,false); // pas de valeurs rentrantes -> grad nul
    HJ.visit(ls,H);
  }
private:
  xLSExtensionOperator();
  const xLevelSet& lst;
  const xVectorField& vec;
};

#if 0
class QuadVec
{
public:
	QuadVec(std::vector<double> &X_, std::vector<double> &Y_, std::vector<double> &Z_, std::vector<double> &Velo_)
		:X(X_), Y(Y_), Z(Z_), Velo(Velo_)
	{	}

public:
	std::vector<double> &X,&Y,&Z,&Velo;

};
#endif

/// This operator smooth the input levelset by performing a least square optimization.
/*! member visit take an xLevelSet input lsin and output xLevelSet lsout which
   is a smoothed version of lsin using a leastsquare smoothing. 
   The lsout is identical to lsin where curvature of lsin is 0.
*/
 class xLevelSetSmoother : public xLevelSetCreator{
public: 
  void visit(const xLevelSet&, xLevelSet&) override;
  ~xLevelSetSmoother() override  = default;
};



/// levelset modifier (visitor pattern) that shift (add to ) the levelset by the constructor given xEval returning double value.
 class xEvalShiftLevelSetModifier : public xLevelSetModifier {
 public:
  xEvalShiftLevelSetModifier(const xEval<double>& eval, const xEntityFilter filter = xAcceptAll() );
  void visit(xLevelSet& ls, xRegion target) override;
  double getMaxShift();
 private:
  const xEval<double>& eval_shift;
  const xEntityFilter filter_upper;
  double max_shift;
};

 
} // end of namespace
	   






#endif
