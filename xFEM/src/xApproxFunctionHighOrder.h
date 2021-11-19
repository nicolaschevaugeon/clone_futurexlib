/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XAPPROXFUNCTIONHIGHORDER_
#define _XAPPROXFUNCTIONHIGHORDER_
#include "mEntity.h"
#include "xApproxFunction.h"

// macro for static allocation of trav member
#define XAPPROXFUNCTIONLAGRANGEMAXORDER 24

namespace xfem
{
/// The functions nbShapeEntitiestype return the
/*! number of shape function attached to the entity type for a given order
    This is not the total number of shape function that evaluate non-zero on that entity !
    for example nbShapeFunctionTri(3) return 1 :
    2
    |\
    | \
    7  6
    |   \
    |    \
    8  9  5
    |      \
    |       \
    0--3--4--1

    There is  a Total of 10 shape function for a triangle of order 3:
    10 = 3*1 + 3*2 + 1       = nbnode*nbshapenode +nbedge*nbshapeedge + nbshapetri
    nbShapeFunctionEdge(3) = 2;
    nbShapeFunctionTri(3)  = 1;

    for example nbShapeFunctionTri(4) return 3 :
    2
    |\
    | \
    9  8
    |   \
    |    \
   10 14  7
    |      \
    |       \
   11 12 13  6
    |         \
    |          \
    0--3--4--5--1

    There is  a Total of 15 shape function for a triangle of order 4:
    10 = 3*1 + 3*3 + 3       = nbnode*nbshapenode +nbedge*nbshapeedge + nbshapetri
    nbShapeFunctionEdge(4) = 3;
    nbShapeFunctionTri(4)  = 3;

 */
int nbShapeFunctionEdge(const int& order);
int nbShapeFunctionTri(const int& order);
int nbShapeFunctionTet(const int& order);
int nbShapeFunctionEdgeTotal(const int& order);
int nbShapeFunctionTriTotal(const int& order);
int nbShapeFunctionTetTotal(const int& order);

/// regular spaced lagrange. u  between 0 and 1.
double lag(int n, int i, double u);

/// compute recursivelly factorial of n
long int factorial(long int n);

/// Helper class to compute shape function in barycentric coordinate, attached to node.
class xSimplexVertex
{
  public:
   xSimplexVertex(int node);
   void getL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double& L) const;
   void getdL(const AOMD::mEntity::mType& type, double& dLdu, double& dLdv, double& dLdw) const;
   void getLdL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double& L, double& dLdu, double& dLdv,
               double& dLdw) const;

  protected:
   const int node;
};

/// Helper class to compute shape function in barycentric coordinate, attached to edge.
class xSimplexEdge
{
  public:
   xSimplexEdge(int edge);
   void getL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double* L) const;
   void getdL(const AOMD::mEntity::mType& type, double* dLdu, double* dLdv, double* dLdw) const;
   void getLdL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double* L, double* dLdu, double* dLdv,
               double* dLdw) const;

  protected:
   const int edge;
   mutable const xtensor::xPoint* lastpoint;
   mutable double storedL[2];
};

/// Helper class to compute shape function in barycentric coordinate, attached to a Triangle.
class xSimplexTri
{
  public:
   xSimplexTri(int face);
   void getL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double* L) const;
   void getdL(const AOMD::mEntity::mType& type, double* dLdu, double* dLdv, double* dLdw) const;
   void getLdL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double* L, double* dLdu, double* dLdv,
               double* dLdw) const;

  protected:
   const int face;
};

/// Helper class to compute shape function in barycentric coordinate, attached to Tet
class xSimplexTet
{
  public:
   xSimplexTet();
   void getL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double* L) const;
   void getdL(const AOMD::mEntity::mType& type, double* dLdu, double* dLdv, double* dLdw) const;
   void getLdL(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, double* L, double* dLdu, double* dLdv,
               double* dLdw) const;
   static const int Tfv[4][3];
   static const int Tev[6][2];
};

/// Approx function of the Bernstein Basis attached to a vertex.
/*! Note : For all shape function definined here, the entity type refer
    to the type of entity the shape function is attached to. Not the one it is evaluated on.
    For example on a tetrahedron, you have shape functions attached to node, edge, triangle and Tetrahedron.
 */
class xApproxFunctionScalarVertexBernstein : private xSimplexVertex, public xApproxFunction
{
  public:
   xApproxFunctionScalarVertexBernstein(int pL0, int node);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   int powL;
   mutable xEvalStorage storage;
};

class xApproxFunctionScalarEdgeBernstein : private xSimplexEdge, public xApproxFunction
{
  public:
   xApproxFunctionScalarEdgeBernstein(int pL0, int pL1, int edge);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   double getIntegrale() const;

  private:
   int powL[2];  // Allways casted as double, so make it double since the creation
   double factor;
   void getVal_(const double& u, double& res) const;
   void getVal_(const AOMD::mEntity* e, const xtensor::xPoint& uvw, double& res) const;
   mutable xEvalStorage storage;
};

class xApproxFunctionScalarTriBernstein : private xSimplexTri, public xApproxFunction
{
  public:
   xApproxFunctionScalarTriBernstein(int pL0, int pL1, int pL2, int face);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   int powL[3];
   double factor;
   mutable xEvalStorage storage;
};

class xApproxFunctionScalarTetBernstein : private xSimplexTet, public xApproxFunction
{
  public:
   xApproxFunctionScalarTetBernstein(int pL0, int pL1, int pL2, int pl3);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   int powL[4];
   double factor;
   mutable xEvalStorage storage;
};

class xApproxFunctionScalarVertexLagrange : public xApproxFunction, private xSimplexVertex
{
  public:
   xApproxFunctionScalarVertexLagrange(int order, int node, double coef = 1.);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   const int order;
   const double coef;
   const double coef_order;
   const int n;
};

class xApproxFunctionScalarEdgeLagrange : public xApproxFunction, private xSimplexEdge
{
  public:
   xApproxFunctionScalarEdgeLagrange(int order, int edge, int ionedge, double coef = 1.);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   static void getNodePos(int order, int ionedge, int* nodepos);
   void getNodeUvw(std::vector<double>& nodeuvw);

  private:
   const int order;
   const int ionedge;
   const double coef;
   const double coef_order;
   int n[2];
};

class xApproxFunctionScalarTriLagrange : public xApproxFunction, private xSimplexTri
{
  public:
   xApproxFunctionScalarTriLagrange(int order, int face, int ionface, double coef = 1.);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   static void getNodePos(int order, int ionface, int* nodepos);
   void getNodeUvw(std::vector<double>& nodeuvw);

  private:
   static std::map<int, std::vector<std::vector<int>>> noi;
   const int order;
   const int ionface;
   const double coef;
   const double coef_order;
   int n[3];
};

class xApproxFunctionScalarTetLagrange : public xApproxFunction, private xSimplexTet
{
  public:
   xApproxFunctionScalarTetLagrange(int order, int iontet, double coef = 1.);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   static void getNodePos(int order, int iontet, int* nodepos);

  private:
   static std::map<int, std::vector<std::vector<int>>> noi;
   const int order;
   const double coef;
   const double coef_order;
   int n[4];
};

}  // namespace xfem
#endif
