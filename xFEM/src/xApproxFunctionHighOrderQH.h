/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XAPPROXFUNCTIONHIGHORDER_QH_
#define _XAPPROXFUNCTIONHIGHORDER_QH_

#include "mEntity.h"
#include "xApproxFunction.h"

namespace xfem
{
class tensorSpaceDefinition
{
  public:
   tensorSpaceDefinition();

   static int nbShapeFunctionEdgeQH(const int& order);
   static int nbShapeFunctionQuad(const int& order);
   static int nbShapeFunctionHex(const int& order);
   static int nbShapeFunctionEdgeTotalQH(const int& order);
   static int nbShapeFunctionQuadTotal(const int& order);
   static int nbShapeFunctionHexTotal(const int& order);
   static bool acceptShapeFunc(const int& order, const int i, const int j, const int k = 0);  // Formally, not needed in practice
};

class trunkSpaceDefinition
{
  public:
   trunkSpaceDefinition();

   static int nbShapeFunctionEdgeQH(const int& order);
   static int nbShapeFunctionQuad(const int& order);
   static int nbShapeFunctionHex(const int& order);
   static int nbShapeFunctionEdgeTotalQH(const int& order);
   static int nbShapeFunctionQuadTotal(const int& order);
   static int nbShapeFunctionHexTotal(const int& order);
   static bool acceptShapeFunc(const int& order, const int i, const int j, const int k = 0);
};

////Swaptable
/////
///// \class xAOMDSwapTable
///// \brief{Compute the swap, aomd ordering}
// class xAOMDSwapTable{
// public :
//    xAOMDSwapTable();

//    static const int Fev[4][2];
//    static const int Hfv[6][4];
//    static const int Hev[12][2];
//};

class xApproxFunctionScalarVertexLagrangeQH : public xApproxFunction
{
  public:
   xApproxFunctionScalarVertexLagrangeQH(int order, int node);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   static void correctUvwForCurrentEntity(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, xtensor::xPoint& uvwCorr);

  private:
   const int order;
   const int inode;
   static const int Vindices[8][3];
};

class xApproxFunctionScalarEdgeLagrangeQH : public xApproxFunction
{
  public:
   xApproxFunctionScalarEdgeLagrangeQH(int order, int edge, int ionedge);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   const int order;
   const int ionedge;
   const int edge;
};

class xApproxFunctionScalarQuadLagrange : public xApproxFunction
{
  public:
   xApproxFunctionScalarQuadLagrange(int order, int face, int ionface);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionface, std::vector<std::vector<int>>& indices);
   const int order;
   const int ionface;
   const int face;
   std::vector<std::vector<int>> indices;
};

class xApproxFunctionScalarHexLagrange : public xApproxFunction
{
  public:
   xApproxFunctionScalarHexLagrange(int order, int ionhex);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionhex, std::vector<std::vector<int>>& indices);
   const int order;
   const int ionhex;
   std::vector<std::vector<int>> indices;
};

double lagPoly(int n, int i, double u);
double dLagPoly(int n, int i, double u);

/// \class xApproxFunctionScalarVertexBernsteinQH
/// \brief{Bernstein approx function for vertices (of quad/hex).}
class xApproxFunctionScalarVertexBernsteinQH : public xApproxFunction
{
  public:
   xApproxFunctionScalarVertexBernsteinQH(int order, int node);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   static void correctUvwForCurrentEntity(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, xtensor::xPoint& uvwCorr);

  private:
   const int order;
   const int inode;
   static const int Vindices[8][3];
};

/// \class xApproxFunctionScalarEdgeBernsteinQH
/// \brief{Bernstein approx function for edges (of quad/hex).}
class xApproxFunctionScalarEdgeBernsteinQH : public xApproxFunction
{
  public:
   xApproxFunctionScalarEdgeBernsteinQH(int order, int edge, int ionedge);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   const int order;
   const int ionedge;
   const int edge;
};

/// \class xApproxFunctionScalarQuadBernstein
/// \brief{Bernstein approx function for faces (of quad/hex).}
class xApproxFunctionScalarQuadBernstein : public xApproxFunction
{
  public:
   xApproxFunctionScalarQuadBernstein(int order, int face, int ionface);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionface, std::vector<std::vector<int>>& indices);
   const int order;
   const int ionface;
   const int face;
   std::vector<std::vector<int>> indices;
};

/// \class xApproxFunctionScalarHexBernstein
/// \brief{Bernstein approx function for Hex.}
class xApproxFunctionScalarHexBernstein : public xApproxFunction
{
  public:
   xApproxFunctionScalarHexBernstein(int order, int ionhex);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionhex, std::vector<std::vector<int>>& indices);
   const int order;
   const int ionhex;
   std::vector<std::vector<int>> indices;
};

inline double binomial(int n, int k);
inline double bernPoly(int v, int n, double u);
inline double dBernPoly(int v, int n, double u);

//------------------------- FF Hierarchiques Legendre ------------------------
// Swaptable
///
/// \class xSzaboSwapTable
/// \brief{Compute the swap, Szabo and Babuska ordering}
class xSzaboSwapTable
{
  public:
   xSzaboSwapTable();

   static const int Fev[4][2];
   static const int Hfv[6][4];
   static const int Hev[12][2];
};

/// \class xApproxFunctionScalarVertexHierarchicalLegendreQH
/// \brief{HierarchicalLegendre approx function for vertices (of quad/hex).}
class xApproxFunctionScalarVertexHierarchicalLegendreQH : public xApproxFunction
{
  public:
   xApproxFunctionScalarVertexHierarchicalLegendreQH(int order, int node);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;
   static void correctUvwForCurrentEntity(const AOMD::mEntity::mType& type, const xtensor::xPoint& uvw, xtensor::xPoint& uvwCorr);

  private:
   const int order;
   const int inode;
};

/// \class xApproxFunctionScalarEdgeHierarchicalLegendreQH
/// \brief{HierarchicalLegendre approx function for edges (of quad/hex).}
class xApproxFunctionScalarEdgeHierarchicalLegendreQH : public xApproxFunction
{
  public:
   xApproxFunctionScalarEdgeHierarchicalLegendreQH(int order, int edge, int ionedge);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   const int order;
   const int ionedge;
   const int edge;
};

/// \class xApproxFunctionScalarQuadHierarchicalLegendre
/// \brief{HierarchicalLegendre approx function for faces (of quad/hex).}

// template<typename SPACETYPE>
class xApproxFunctionScalarQuadHierarchicalLegendre : public xApproxFunction
{
  public:
   xApproxFunctionScalarQuadHierarchicalLegendre(int order, int face, int ionface, bool tensorSpace = true);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionface, std::vector<std::vector<int>>& indices, bool tensorSpace = true);
   const int order;
   const int ionface;
   const int face;
   std::vector<std::vector<int>> indices;
};

/// \class xApproxFunctionScalarHexHierarchicalLegendre
/// \brief{HierarchicalLegendre approx function for Hex.}
class xApproxFunctionScalarHexHierarchicalLegendre : public xApproxFunction
{
  public:
   xApproxFunctionScalarHexHierarchicalLegendre(int order, int ionhex, bool tensorSpace = true);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionhex, std::vector<std::vector<int>>& indices, bool tensorSpace = true);
   const int order;
   const int ionhex;
   std::vector<std::vector<int>> indices;
};

// inline double legendrePoly(int k, double u);
// inline double dLegendrePoly(int k, double u);

// Attention, k commence a 2 pour les phiFuncs
inline double phiFunc(int k, double u);
inline double dPhiFunc(int k, double u);

/// \class xApproxFunctionScalarEdgeHierarchicalLegendreQHSwap
/// \brief{HierarchicalLegendre approx function for edges (of quad/hex). SWAPPED VERSION}
class xApproxFunctionScalarEdgeHierarchicalLegendreQHSwap : public xApproxFunction
{
  public:
   xApproxFunctionScalarEdgeHierarchicalLegendreQHSwap(int order, int edge, int ionedge, int swap_ = 0);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   const int order;
   const int ionedge;
   const int edge;
   const int swap;
};

/// \class xApproxFunctionScalarQuadHierarchicalLegendreSwap
/// \brief{HierarchicalLegendre approx function for faces (of quad/hex). SWAPPED VERSION}

// template<typename SPACETYPE>
class xApproxFunctionScalarQuadHierarchicalLegendreSwap : public xApproxFunction
{
  public:
   xApproxFunctionScalarQuadHierarchicalLegendreSwap(int order, int face, int ionface, bool tensorSpace = true);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionface, std::vector<std::vector<int>>& indices, bool tensorSpace = true);
   const int order;
   const int ionface;
   const int face;
   std::vector<std::vector<int>> indices;
};

}  // namespace xfem

namespace xfem
{
//------------------------- FF Hierarchiques Legendre VERSION DEPRECIEE ------------------------

///// \class xApproxFunctionScalarVertexHierarchicalLegendreQH
///// \brief{HierarchicalLegendre approx function for vertices (of quad/hex).}
// class xApproxFunctionScalarVertexHierarchicalLegendreQH :  public xApproxFunction {
// public:
//    xApproxFunctionScalarVertexHierarchicalLegendreQH(int order, int node);
//    void  getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double &res) const override;
//    void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector& res) const override;
//    static void correctUvwForCurrentEntity(const AOMD::mEntity::mType& type, const xtensor::xPoint &uvw, xtensor::xPoint
//    &uvwCorr);
// private:
//    const int order;
//    const int inode;
//};

/// \class xApproxFunctionScalarEdgeHierarchicalLegendreQH
/// \brief{HierarchicalLegendre approx function for edges (of quad/hex).}
class xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated : public xApproxFunction
{
  public:
   xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated(int order, int edge, int ionedge, int swap_ = 0);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   const int order;
   const int ionedge;
   const int edge;
   const int swap;
};

/// \class xApproxFunctionScalarQuadHierarchicalLegendre
/// \brief{HierarchicalLegendre approx function for faces (of quad/hex).}

// template<typename SPACETYPE>
class xApproxFunctionScalarQuadHierarchicalLegendre_deprecated : public xApproxFunction
{
  public:
   xApproxFunctionScalarQuadHierarchicalLegendre_deprecated(int order, int face, int ionface, int swap_ = 0,
                                                            bool tensorspace = true);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& res) const override;

  private:
   static void setNodeIndices(int order, int ionface, std::vector<std::vector<int>>& indices, bool tensorSpace = true);
   static void applySwap(int* idx, double* signs, int swapcase);
   const int order;
   const int ionface;
   const int face;
   std::vector<std::vector<int>> indices;
   static const int swapvar[8][2];
   const int swap;
};

///// \class xApproxFunctionScalarHexHierarchicalLegendre
///// \brief{HierarchicalLegendre approx function for Hex.}
// class xApproxFunctionScalarHexHierarchicalLegendre :  public xApproxFunction {
// public:
//    xApproxFunctionScalarHexHierarchicalLegendre(int order, int ionhex, bool tensorSpace = true);
//    void  getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double &res) const override;
//    void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector& res) const override;
// private:
//    static void setNodeIndices(int order, int ionhex, std::vector<std::vector<int> > &indices, bool tensorSpace = true);
//    const int order;
//    const int ionhex;
//    std::vector<std::vector<int> >  indices;
//};

// inline double legendrePoly(int k, double u);
// inline double dLegendrePoly(int k, double u);

// Attention, k commence a 2 pour les phiFuncs
namespace deprecatedLegendre
{
inline double phiFunc(int k, double u);
inline double dPhiFunc(int k, double u);
}  // namespace deprecatedLegendre

}  // namespace xfem

#endif
