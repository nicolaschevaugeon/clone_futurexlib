/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef _SHAPE_FUNCTIONSSS_
#define _SHAPE_FUNCTIONSSS_

#include <cassert>
#include <iostream>

#include "xApproxFctPtr.h"
#include "xEntityFilter.h"
#include "xEvalStorage.h"
#include "xValKey.h"
#include "xVector.h"

//
//  Warning :
//    Precondition
//         geo_appro and geo_integ must contain the current uvw point
//    Postcondition
//         geo_appro and geo_integ are not changed by the call

namespace xfem
{
class xLevelSet;
class xApproxFunction
{
  public:
   typedef approxFctPtr_t shapeFctPtr;
   virtual ~xApproxFunction() = default;
   ///< To reset the storage (only for enrichment functions)
   virtual void resetStorage(){};

  public:
   virtual std::string name() { return "unknown"; };
   /// return the value of a scalar shape function
   virtual void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const;
   /// return the gradient of a scalar shape function, with regard to global coordinates xyz
   virtual void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const;
   /// return the gradient of a  shape function with regard to local coordinates uvw (reference element)
   virtual void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const;
   /// return the hessian of a scalar shape function
   virtual void getHessian(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const;
   /// return the value of a vector shape function
   virtual void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const;
   /// return the gradient of a vector shape function with regard to global coordinates xyz
   virtual void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const;
   /// return the hessian of a vector shape function with regard to global coordinates xyz
   virtual void getHessian(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor4<double>&) const;
   /// return the value of a vector shape function
   virtual void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&, const xValKey&) const;
   /// return the gradient of a vector shape function with regard to global coordinates xyz
   virtual void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&, const xValKey&) const;
   /// return the gradient of a vector shape function with regard to local coordinates uvw (reference element)
   virtual void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const;
   virtual void getGradLocalAxpy(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>& y,
                                 const double& a) const;

  protected:
   xApproxFunction();
};

/// Specialization of xApproxFunction for constant shape functions.
class xApproxFunctionConstant : public xApproxFunction
{
  public:
   xApproxFunctionConstant() : xApproxFunction() {}
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
};

/// Create a vector shape function from a scalar one.
class xApproxFunctionVector : public xApproxFunction
{
  public:
   /// Input :  a scalar shape function and an integer representing the component of the vector shape function that is not zero.
   /*! The new shape function constructed is a vector shape function, with all the components of the vector equal to zero but
     the one selected by iComp which has the same value than the input shape function itself*/
   xApproxFunctionVector(shapeFctPtr sf, unsigned int comp);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   shapeFctPtr getScalarFunc() { return scalarFunction; }

  private:
   shapeFctPtr scalarFunction;
   const unsigned int component;
};

/// Create a shape function sum of a vector of shape function
class xApproxFunctionSummed : public xApproxFunction
{
  public:
   xApproxFunctionSummed(std::vector<shapeFctPtr>& approx);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   std::vector<shapeFctPtr> approx_functions;
};

class xApproxFunctionEnrichedXFEM : public xApproxFunction
{
  public:
   xApproxFunctionEnrichedXFEM(shapeFctPtr b, shapeFctPtr enr) : base(b), enrichment(enr) {}
   std::string name() override { return "xcApproxFunctionEnrichedXFEM"; }
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   void resetStorage() override { enrichment->resetStorage(); }

  protected:
   shapeFctPtr base;
   shapeFctPtr enrichment;
};

/// \class xApproxFunctionEnrichedXFEMShifted
/// \brief{Shifted approx function. Works only for nodes shifting.}
/// \param[in] e : Associated entity (vertex)
class xApproxFunctionEnrichedXFEMShifted : public xApproxFunction
{
  public:
   xApproxFunctionEnrichedXFEMShifted(shapeFctPtr b, shapeFctPtr enr, AOMD::mEntity* e_) : base(b), enrichment(enr), enti(e_)
   { /* std::cout<<"entity is"; enti->print();std::cout<<std::endl;*/
   }
   std::string name() override { return "xcApproxFunctionEnrichedXFEM SHIFTED"; }
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   void resetStorage() override { enrichment->resetStorage(); }

  protected:
   shapeFctPtr base;
   shapeFctPtr enrichment;
   AOMD::mEntity* enti;
};

/// Enriched shape function of vector type.
/*! The enrichment function is of type vector, while the classic function are scalar shape function */
class xApproxFunctionEnrichedVectorXFEM : public xApproxFunction
{
  public:
   /// Constructor : The enrichment function enr, need to be of tensor type vector, the base function need to be of tensor type
   /// scalar
   xApproxFunctionEnrichedVectorXFEM(shapeFctPtr b, shapeFctPtr enr) : base(b), enrichment(enr) {}
   std::string name() override { return "xcApproxFunctionEnrichedVectorXfem"; }
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   void resetStorage() override { enrichment->resetStorage(); }

  protected:
   shapeFctPtr base;
   shapeFctPtr enrichment;
};

/// Enriched shape function of vector type, that need to keep the enriched entity has a data of the object
/*! The enrichment function is of type vector, while the classic function are scalar shape function */
class xApproxFunctionEnrichedVector2XFEM : public xApproxFunction
{
  public:
   /// Constructor
   /*! The enrichment function enr, need to be of tensor type vector, the base function need to be of tensor type scalar.
     The entity e, is stored as a protected member of the class. Rational :  e is supposed to be the enriched entity.
     This class is used to define vector enrichment function for cracks.
   */
   xApproxFunctionEnrichedVector2XFEM(shapeFctPtr b, shapeFctPtr enr, const xValKey& _basekey)
       : base(b), enrichment(enr), basekey(_basekey)
   {
   }
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   void resetStorage() override { enrichment->resetStorage(); }
   // virtual void  getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&, AOMD::mEntity *);
   // virtual void  getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&, AOMD::mEntity *);
  protected:
   shapeFctPtr base;
   shapeFctPtr enrichment;
   xValKey basekey;
   // AOMD::mEntity* e;
};

/// This is the "historical" shape function in the xfem library.
/*! Despite the name (bound to change  in the future) this is not lagrange shape function, but a hierachical basis.
  HIERARCHICAL
*/
class xApproxFunctionScalarHierarchical : public xApproxFunction
{
  public:
   xApproxFunctionScalarHierarchical(int i) : ith(i) {}
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
   int ith;
};

class xApproxFunctionScalarPieceWiseHierarchical : public xApproxFunction
{
  public:
   xApproxFunctionScalarPieceWiseHierarchical(int i) : ith(i) {}
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
   int ith;
};

class xApproxFunctionVectorHierarchical : public xApproxFunction
{
  public:
   xApproxFunctionVectorHierarchical(int i, int c) : ith(i), comp(c) {}
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;
   void getGradLocalAxpy(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<>& y,
                         const double& a) const override;

  private:
   int ith;
   int comp;  // 0, 1 or 2
};

class xScalarFunctionDerivDiscXFEM : public xApproxFunction, private xEvalStorage
{
  public:
   xScalarFunctionDerivDiscXFEM(const xLevelSet& ls);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void resetStorage() override { storage.reset(); }

  private:
   const xLevelSet& field;
   mutable xEvalStorage storage;
};

class xScalarFunctionRegularizedHeaviside : public xApproxFunction, private xEvalStorage
{
  public:
   xScalarFunctionRegularizedHeaviside(const xEntityFilter& _f_in_isozero, const xEntityFilter& _f_subelement_regularized,
                                       const xEntityFilter& _f_value_one);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void resetStorage() override { storage.reset(); }

  private:
   xEntityFilter filter_in_iso_zero;
   xEntityFilter filter_subelement_regularized;
   xEntityFilter filter_value_one;
   mutable xEvalStorage storage;
};

class xScalarFunctionDiscXFEM : public xApproxFunction
{
  public:
   xScalarFunctionDiscXFEM(const xLevelSet& ls_);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void resetStorage() override { storage.reset(); }

  private:
   const xLevelSet& ls;
   mutable xEvalStorage storage;
};

//  inline bool sign( double & x) {
//    return x>0.;
//  };

double HierarchicalApproxFunction(int funct, const xGeomElem*);
xtensor::xVector<> GradHierarchicalApproxFunction(int funct, const xGeomElem*);

double HierarchicalApproxFunctionEdge(int ith, double u);
double HierarchicalApproxFunctionTriangle(int ith, double u, double v);
double HierarchicalApproxFunctionQuadrilateral(int ith, double u, double v);
double HierarchicalApproxFunctionTetrahedron(int ith, double u, double v, double w);
double HierarchicalApproxFunctionPrism(int ith, double u, double v, double w);
double HierarchicalApproxFunctionHexahedron(int ith, double u, double v, double w);

void GradHierarchicalApproxFunctionEdge(int ith, double u, double result[]);
void GradHierarchicalApproxFunctionTriangle(int ith, double u, double v, double result[]);
void GradHierarchicalApproxFunctionQuadrilateral(int ith, double u, double v, double result[]);
void GradHierarchicalApproxFunctionTetrahedron(int ith, double u, double v, double w, double result[]);
void GradHierarchicalApproxFunctionPrism(int ith, double u, double v, double w, double result[]);
void GradHierarchicalApproxFunctionHexahedron(int ith, double u, double v, double w, double result[]);

double SimplexApproxFunctionQuadrilateral(int ith, double u, double v);
double SimplexApproxFunctionHexahedron(int ith, double u, double v, double w);
void GradSimplexApproxFunctionQuadrilateral(int ith, double u, double v, double result[]);
void GradSimplexApproxFunctionHexahedron(int ith, double u, double v, double w, double result[]);

double Psi(double x, int n);
double Dpsi(double x, int n);
double Phi(double x, int n);
double Dphi(double x, int n);
double LegendrePol(double x, int n);

}  // namespace xfem

#endif
