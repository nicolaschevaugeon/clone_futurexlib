/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.

*/
  
#ifndef _FUNCTIONCRACKXFEM_H
#define _FUNCTIONCRACKXFEM_H

#include <cmath> 
#include "xApproxFunction.h"
#include "xEvalStorage.h"
#include "xEval.h"
#include "xVector.h"
#include "xGeomElem.h"

using namespace xfem;

namespace xcrack
{
  class lCrack;
  class xcCrackBase;
  class xcVectorFunctionCrack:public xApproxFunction, public xEval<xtensor::xVector<> >, private xEvalStorage {
  public: 
    xcVectorFunctionCrack(const xcCrackBase & crk, double nu );
    void getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
    void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
    void getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&, const xValKey &) const override;
    void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<>&, const xValKey &) const override;
    void getGradNum(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<>&, const xValKey & ) const;
    void operator()(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>& res) const override{ getVal( geo_appro, geo_integ, res); };
    static bool evaltrick;
  protected:
    virtual void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const=0;
    virtual void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res) const=0;
    const xcCrackBase& crack;
    double kappa;
    const double scalfac;
  };
  
  class xcVectorFunctionCrack1 : public xcVectorFunctionCrack {
  public:
    xcVectorFunctionCrack1(const xcCrackBase & crk, double nu );
  protected:
    void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
    void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res)const override;
  };

  class xcVectorFunctionCrack1r : public xcVectorFunctionCrack {
  public:
    xcVectorFunctionCrack1r(const xcCrackBase & crk, double nu );
  protected:
    void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
    void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res)const override;
  };


  class xcVectorFunctionCrack1r32 : public xcVectorFunctionCrack {
  public:
    xcVectorFunctionCrack1r32(const xcCrackBase & crk, double nu );
  protected:
    void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
    void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res)const override;
  };


  class xcVectorFunctionCrack2 : public xcVectorFunctionCrack {
  public:
    xcVectorFunctionCrack2(const xcCrackBase & crk, double nu );
  protected:
    void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
    void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res)const override;
  };

  class xcVectorFunctionCrack2r : public xcVectorFunctionCrack {
  public:
    xcVectorFunctionCrack2r(const xcCrackBase & crk, double nu );
  protected:
    void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
    void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res)const override;
  };

  class xcVectorFunctionCrack2r32 : public xcVectorFunctionCrack {
  public:
    xcVectorFunctionCrack2r32(const xcCrackBase & crk, double nu );
  protected:
    void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
    void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res)const override;
  };

  class xcVectorFunctionCrack3 : public xcVectorFunctionCrack {
  public:
    xcVectorFunctionCrack3(const xcCrackBase & crk, double nu );
  protected:
    void getValLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
    void getGradLocalAxis(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xTensor2<> &res)const override;
  };



class ScalarFunctionCrackXFEM_c : public xApproxFunction {
public:
  ScalarFunctionCrackXFEM_c(const xcCrackBase& crk);
  void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&)const override;
  void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
  void resetStorage() override{storage.reset();}
private:
  const xcCrackBase& crack;
  mutable xEvalStorage storage;
};


class ScalarFunctionCrackTip1XFEM_c : public xApproxFunction, public xEval<double> {
public:
  ScalarFunctionCrackTip1XFEM_c(const xcCrackBase& crk);
  void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&)const override;
  void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
  void operator()(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double& res) const override{ return getVal( geo_appro, geo_integ, res);}
  //    (const_cast<ScalarFunctionCrackTip1XFEM_c *>(this))->getVal( geo_appro, geo_integ, res); };
  void resetStorage() override{storage.reset();}
  static bool evaltrick;    
private:
  const xcCrackBase& crack;
  mutable xEvalStorage storage;
};


class ScalarFunctionCrackTip2XFEM_c : public xApproxFunction, public xEval<double> {
public:
  ScalarFunctionCrackTip2XFEM_c(const xcCrackBase& crk);
  void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&) const override;
  void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
  void operator()(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double& res) const override{ getVal( geo_appro, geo_integ, res); };
  void resetStorage() override{storage.reset();}
  static bool evaltrick;
private:
  const xcCrackBase& crack;
  mutable xEvalStorage storage;
};

class ScalarFunctionCrackTip3XFEM_c : public xApproxFunction, public xEval<double> {
public:
  ScalarFunctionCrackTip3XFEM_c(const xcCrackBase& crk);
  void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&) const override;
  void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
  void operator()(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double& res) const override{ getVal( geo_appro, geo_integ, res); };
  void resetStorage() override{storage.reset();}
  static bool evaltrick;
private:
  const xcCrackBase& crack;
  mutable xEvalStorage storage;
};

class ScalarFunctionCrackTip4XFEM_c : public xApproxFunction, public xEval<double> {
public:
  ScalarFunctionCrackTip4XFEM_c(const xcCrackBase& crk);
  void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&) const override;
  void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
  void operator()(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double& res) const override{ getVal( geo_appro, geo_integ, res); };
  void resetStorage() override{storage.reset();}
  static bool evaltrick;
private:
  const xcCrackBase& crack;
  mutable xEvalStorage storage;
};

class EnrichMod_c : public xApproxFunction {
public:
  EnrichMod_c(lCrack& crk,double r0p,double attp);
  void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&) const override;
  void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
  void resetStorage() override{storage.reset();}
private:
  lCrack* crack;
  double r0;
  double att;
  mutable xEvalStorage storage;
};

class axisChangeCurv{
public:
  axisChangeCurv(const xtensor::xVector<> &e0, const xtensor::xVector<> &e1, const xtensor::xVector<> &e2,
		 const xtensor::xTensor2<> &ge0, const xtensor::xTensor2<> & ge1, const xtensor::xTensor2<> & ge2);
  xtensor::xVector<>  localToGlobal(const xtensor::xVector<>&) const;
  xtensor::xTensor2<> localToGlobal(const xtensor::xVector<>&, const xtensor::xTensor2<>&) const;
private :
  const xtensor::xVector<> &e0;
  const xtensor::xVector<> &e1;
  const xtensor::xVector<> &e2;
  const xtensor::xTensor2<> &ge0;
  const xtensor::xTensor2<> &ge1;
  const xtensor::xTensor2<> &ge2;
  mutable xtensor::xTensor2<> Q;
  mutable xtensor::xTensor2<> Qinv;
};

class axisChange{
public:
  axisChange(const xtensor::xVector<> &e0, const xtensor::xVector<> &e1, const xtensor::xVector<> &e2);
  xtensor::xVector<>  localToGlobal(const xtensor::xVector<>&) const;
  xtensor::xTensor2<> localToGlobal(const xtensor::xTensor2<>&) const;
private :
  const xtensor::xVector<> &e0;
  const xtensor::xVector<> &e1;
  const xtensor::xVector<> &e2;
  mutable xtensor::xTensor2<> Q;
  mutable xtensor::xTensor2<> Qinv;
}; 

 class xcCrackHeavisideEnrichment: public xApproxFunction {
  public:
   xcCrackHeavisideEnrichment(const xcCrackBase& crack);
   void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&)const override;
   void  getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&)const override;
  private:
   const xcCrackBase& crack;
};


}

#endif




