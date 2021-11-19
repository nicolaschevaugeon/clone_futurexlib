/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xSpacePolynomialQH.h"

#include <iostream>
#include <sstream>

#include "mEdge.h"
#include "mFace.h"
#include "mHex.h"
#include "xFiniteElement.h"
#include "xValKey.h"

using namespace xfem;
using AOMD::mEntity;
using std::cout;
using std::endl;
using xfem::xKeyInfo;

#define SWAPFUNCS 0

namespace xfem
{
xSpacePolynomialLagrangeQH::xSpacePolynomialLagrangeQH(const std::string &a, TensorialType_t c, const int p)
    : xSpacePolynomialQH<tensorSpaceDefinition>("Lagrange_HQ", a, c, p)
{
}

void xSpacePolynomialLagrangeQH::getFcts(const mEntity &e, femFcts &approscal)
{
   switch (e.getLevel())
   {
      case 0:
      {
         // pour le cas ou la fct n'est pas interpolante ? -> ca devrait tjrs valoir 1 non ?
         approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrangeQH(order, 0)));
         break;
      }
      case 1:
      {
         for (int node = 0; node < 2; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrangeQH(order, node)));
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrangeQH(order, 0, i)));
         }
         break;
      }
      case 2:
      {
         for (int node = 0; node < 4; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrangeQH(order, node)));
         }
         for (int edge = 0; edge < 4; ++edge)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrangeQH(order, edge, i)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadLagrange(order, 0, i)));
         }
         break;
      }
      case 3:
      {
         for (int node = 0; node < 8; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrangeQH(order, node)));
         }
         for (int edge = 0; edge < 12; ++edge)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrangeQH(order, edge, i)));
            }
         }
         for (int face = 0; face < 6; ++face)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadLagrange(order, face, i)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionHex(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarHexLagrange(order, i)));
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}

xSpacePolynomialBernsteinQH::xSpacePolynomialBernsteinQH(const std::string &a, TensorialType_t c, const int p)
    : xSpacePolynomialQH<tensorSpaceDefinition>("Bernstein_HQ", a, c, p)
{
}

void xSpacePolynomialBernsteinQH::getFcts(const mEntity &e, femFcts &approscal)
{
   switch (e.getLevel())
   {
      case 0:
      {
         // pour le cas ou la fct n'est pas interpolante ? -> ca devrait tjrs valoir 1 non ?
         approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernsteinQH(order, 0)));
         break;
      }
      case 1:
      {
         for (int node = 0; node < 2; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernsteinQH(order, node)));
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernsteinQH(order, 0, i)));
         }
         break;
      }
      case 2:
      {
         for (int node = 0; node < 4; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernsteinQH(order, node)));
         }
         for (int edge = 0; edge < 4; ++edge)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernsteinQH(order, edge, i)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadBernstein(order, 0, i)));
         }
         break;
      }
      case 3:
      {
         for (int node = 0; node < 8; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernsteinQH(order, node)));
         }
         for (int edge = 0; edge < 12; ++edge)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernsteinQH(order, edge, i)));
            }
         }
         for (int face = 0; face < 6; ++face)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadBernstein(order, face, i)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionHex(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarHexBernstein(order, i)));
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}

xSpacePolynomialHierarchicalLegendreTrunkSpaceQH::xSpacePolynomialHierarchicalLegendreTrunkSpaceQH(const std::string &a,
                                                                                                   TensorialType_t c, const int p)
    : xSpacePolynomialQH<trunkSpaceDefinition>("Hierarchical_Legendre_trunk_space_HQ", a, c, p)
{
}

void xSpacePolynomialHierarchicalLegendreTrunkSpaceQH::getFcts(const mEntity &e, femFcts &approscal)
{
   bool tensorSpace = false;
   switch (e.getLevel())
   {
      case 0:
      {
         // pour le cas ou la fct n'est pas interpolante ? -> ca devrait tjrs valoir 1 non ?
         approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, 0)));
         break;
      }
      case 1:
      {
         for (int node = 0; node < 2; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int i = 0; i < trunkSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH(order, 0, i)));
         }
         break;
      }
      case 2:
      {
         for (int node = 0; node < 4; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 4; ++edge)
         {
            for (int i = 0; i < trunkSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH(order, edge, i)));
            }
         }
         for (int i = 0; i < trunkSpaceDefinition::nbShapeFunctionQuad(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendre(order, 0, i, tensorSpace)));
         }
         break;
      }
      case 3:
      {
         for (int node = 0; node < 8; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 12; ++edge)
         {
            for (int i = 0; i < trunkSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH(order, edge, i)));
            }
         }
         for (int face = 0; face < 6; ++face)
         {
            for (int i = 0; i < trunkSpaceDefinition::nbShapeFunctionQuad(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendre(order, face, i, tensorSpace)));
            }
         }
         for (int i = 0; i < trunkSpaceDefinition::nbShapeFunctionHex(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarHexHierarchicalLegendre(order, i, tensorSpace)));
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}

#if !SWAPFUNCS
xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH::xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH(
    const std::string &a, TensorialType_t c, const int p)
    : xSpacePolynomialQH<tensorSpaceDefinition>("Hierarchical_Legendre_tensor_product_space_HQ", a, c, p)
{
}

void xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH::getFcts(const mEntity &e, femFcts &approscal)
{
   // e.print();

   switch (e.getLevel())
   {
      case 0:
      {
         // pour le cas ou la fct n'est pas interpolante ? -> ca devrait tjrs valoir 1 non ?
         approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, 0)));
         break;
      }
      case 1:
      {
         for (int node = 0; node < 2; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH(order, 0, i)));
         }
         break;
      }
      case 2:
      {
         for (int node = 0; node < 4; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 4; ++edge)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH(order, edge, i)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendre(order, 0, i)));
         }
         break;
      }
      case 3:
      {
         for (int node = 0; node < 8; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 12; ++edge)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH(order, edge, i)));
            }
         }
         for (int face = 0; face < 6; ++face)
         {
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendre(order, face, i)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionHex(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarHexHierarchicalLegendre(order, i)));
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}
#else

xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH::xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH(
    const std::string &a, TensorialType_t c, const int p)
    : xSpacePolynomialQH<tensorSpaceDefinition>("Hierarchical_Legendre_tensor_product_space_HQ", a, c, p)
{
}

void xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH::getFcts(const mEntity &e, femFcts &approscal)
{
   // e.print();

   switch (e.getLevel())
   {
      case 0:
      {
         // pour le cas ou la fct n'est pas interpolante ? -> ca devrait tjrs valoir 1 non ?
         approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, 0)));
         break;
      }
      case 1:
      {
         for (int node = 0; node < 2; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH(order, 0, i)));
         }
         break;
      }
      case 2:
      {
         for (int node = 0; node < 4; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 4; ++edge)
         {
            const int swap = computeSwapCaseQuadEdge<xSzaboSwapTable>(dynamic_cast<const AOMD::mFace &>(e), edge);
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQHSwap(order, edge, i, swap)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendre(order, 0, i)));
         }
         break;
      }
      case 3:
      {
         for (int node = 0; node < 8; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 12; ++edge)
         {
            const int swap = computeSwapCaseQuadEdge<xSzaboSwapTable>(dynamic_cast<const AOMD::mFace &>(*e), edge);
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQHSwap(order, edge, i, swap)));
            }
         }
         for (int face = 0; face < 6; ++face)
         {
            const int swap = computeSwapCaseHexFace<xSzaboSwapTable>(dynamic_cast<const AOMD::mHex &>(*e), face);
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendreSwap(order, face, i, swap)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionHex(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarHexHierarchicalLegendre(order, i)));
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}

#endif

// Specialization for tensor product space (more code could be factorized with trunk space)
template <>
void xSpacePolynomialQH<tensorSpaceDefinition>::create_geom_id_quad(const std::string &space_name, int order,
                                                                    std::array<std::vector<xValKey::ids_size_t>, 8> &quad_geom_id)
{
   const int nb_shape_quad = tensorSpaceDefinition::nbShapeFunctionQuad(order);
   for_each(quad_geom_id.begin(), quad_geom_id.end(), [&nb_shape_quad](std::vector<xValKey::ids_size_t> &quad_geom_id) {
      quad_geom_id.clear();
      quad_geom_id.resize(nb_shape_quad);
   });

   if (nb_shape_quad)
   {
      std::ostringstream name;
      std::string zero_string;
      int index = 0;
      for (int j = 0; j < order - 1; ++j)
      {
         for (int i = 0; i < order - 1; ++i)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TENSORPROD_" << i + 1 << "_" << j + 1;
            quad_geom_id[0][index] = xKeyInfo::getGeomId(name.str());
            ++index;
         }
      }

      const int nt = (order - 1);
      const std::array<std::function<int(int &, int &)>, 8> swapRenums = {
          [nt](int &i, int &j) { return i + nt * j; },
          [nt](int &i, int &j) { return j + nt * i; },
          [nt](int &i, int &j) { return (nt - 1 - j) + nt * i; },  // 2 et 6 inverses -> bug sur mon cahier ?
          [nt](int &i, int &j) { return (nt - 1 - i) + nt * j; },
          [nt](int &i, int &j) { return (nt - 1 - i) + nt * (nt - 1 - j); },
          [nt](int &i, int &j) { return (nt - 1 - j) + nt * (nt - 1 - i); },
          [nt](int &i, int &j) { return j + nt * (nt - 1 - i); },  // 2 et 6 inverses -> bug sur mon cahier ?
          [nt](int &i, int &j) { return i + nt * (nt - 1 - j); },
      };
      for (int swap = 1; swap < 8; ++swap)
      {
         std::vector<xValKey::ids_size_t> &quad_geom_id_s = quad_geom_id[swap];

         int idx = 0;
         for (int j = 0; j < order - 1; ++j)
         {
            for (int i = 0; i < order - 1; ++i)
            {
               quad_geom_id_s[idx] = quad_geom_id[0][swapRenums[swap](i, j)];
               //                    cout<<swap<<" "<<i<<" "<<j<<" "<<swapRenums[swap](i,j)<< " "<< i + nt * j<< " "<<j + nt * i
               //                    <<endl;
               ++idx;
            }
         }
      }

#if 0
        //Debug swaps
        for (int swap = 0; swap < 8; ++swap){
            auto keys = quad_geom_id[swap];
            cout<<"Generated keys for swap no "<<swap<<endl;
            for(int &key : keys){
                cout<<xKeyInfo::getGeomName(key)<<endl;
            }
        }
#endif
   }
}

// Specialization for trunk space (more code could be factorized with tensor product space)
// TODO: Do this properly rather than reconstructing the swaped keys by hand...
template <>
void xSpacePolynomialQH<trunkSpaceDefinition>::create_geom_id_quad(const std::string &space_name, int order,
                                                                   std::array<std::vector<xValKey::ids_size_t>, 8> &quad_geom_id)
{
   const int nb_shape_quad = trunkSpaceDefinition::nbShapeFunctionQuad(order);
   for_each(quad_geom_id.begin(), quad_geom_id.end(), [&nb_shape_quad](std::vector<xValKey::ids_size_t> &quad_geom_id) {
      quad_geom_id.clear();
      quad_geom_id.resize(nb_shape_quad);
   });

   if (nb_shape_quad)
   {
      std::ostringstream name;
      std::string zero_string;
      int index = 0;
      for (int j = 0; j < order - 3; ++j)
      {
         for (int i = j; i < order - 3; ++i)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[0][index] = xKeyInfo::getGeomId(name.str());
            ++index;
         }
      }

      // swap 1: swap i and j
      int idx = 0;
      for (int i = 0; i < order - 3; ++i)
      {
         for (int j = i; j < order - 3; ++j)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[1][idx] = xKeyInfo::getGeomId(name.str());
            ++idx;
         }
      }

      // swap 2: swap i and j && invert i
      idx = 0;
      for (int i = order - 4; i >= 0; --i)
      {
         for (int j = i; j < order - 3; ++j)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[2][idx] = xKeyInfo::getGeomId(name.str());
            ++idx;
         }
      }

      // swap 3: invert i, j unchanged
      idx = 0;
      for (int j = 0; j < order - 3; ++j)
      {
         for (int i = order - 4; i >= j; --i)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[3][idx] = xKeyInfo::getGeomId(name.str());
            ++idx;
         }
      }

      // swap 4:  invert i, invert j
      idx = 0;
      for (int j = order - 4; j >= 0; --j)
      {
         for (int i = order - 4; i >= j; --i)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[4][idx] = xKeyInfo::getGeomId(name.str());
            ++idx;
         }
      }

      // swap 5:  swap i and j && invert i and j
      idx = 0;
      for (int i = order - 4; i >= 0; --i)
      {
         for (int j = order - 4; j >= i; --j)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[5][idx] = xKeyInfo::getGeomId(name.str());
            ++idx;
         }
      }

      // swap 6:  swap i and j && invert j
      idx = 0;
      for (int i = 0; i < order - 3; ++i)
      {
         for (int j = order - 4; j >= i; --j)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[6][idx] = xKeyInfo::getGeomId(name.str());
            ++idx;
         }
      }

      // swap 7: invert j, i unchanged
      idx = 0;
      for (int j = order - 4; j >= 0; --j)
      {
         for (int i = j; i < order - 3; ++i)
         {
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_QUAD_TRUNK_" << i + 1 << "_" << j + 1;
            quad_geom_id[7][idx] = xKeyInfo::getGeomId(name.str());
            ++idx;
         }
      }

#if 0
        //Debug swaps
        for (int swap = 0; swap < 8; ++swap){
            auto keys = quad_geom_id[swap];
            cout<<"Generated keys for swap no "<<swap<<endl;
            for(int &key : keys){
                cout<<xKeyInfo::getGeomName(key)<<endl;
            }
        }
        //        throw;
#endif
   }
}

//---------------------------------------------------

xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH_deprecated::
    xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH_deprecated(const std::string &a, TensorialType_t c, const int p)
    : xSpacePolynomialQH<tensorSpaceDefinition>("Hierarchical_Legendre_tensor_product_space_HQ", a, c, p)
{
}

void xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH_deprecated::getKeys(AOMD::mEntity *e, femKeys *keys)
{
   switch (e->getLevel())
   {
      case 0:
      {
         keys->reserve(keys->size() + 1);
         getNodeKeys(dynamic_cast<AOMD::mVertex *>(e), keys);
         break;
      }
      case 1:
      {
         keys->reserve(keys->size() + tensorSpaceDefinition::nbShapeFunctionEdgeTotalQH(order));
         for (int node = 0; node < 2; ++node)
         {
            getNodeKeys(dynamic_cast<AOMD::mVertex *>(e->get(0, node)), keys);
         }
         getEdgeKeys((AOMD::mEdge *)(e), keys, 0);  // do not reorder keys
         break;
      }
      case 2:
      {
         if (e->getType() == AOMD::mEntity::QUAD)
         {
            keys->reserve(keys->size() + tensorSpaceDefinition::nbShapeFunctionQuadTotal(order));
            for (int node = 0; node < 4; ++node)
            {
               getNodeKeys(dynamic_cast<AOMD::mVertex *>(e->get(0, node)), keys);
            }
            for (int edge = 0; edge < 4; ++edge)
            {
               getEdgeKeys((AOMD::mEdge *)(e->get(1, edge)), keys, 0);  // do not reorder keys
            }
            getQuadKeys((AOMD::mFace *)(e), keys, 0);  // do not reorder keys
         }
         else
            throw;
         break;
      }
      case 3:
      {
         if (e->getType() == AOMD::mEntity::HEX)
         {
            keys->reserve(keys->size() + tensorSpaceDefinition::nbShapeFunctionHexTotal(order));
            for (int node = 0; node < 8; ++node)
            {
               getNodeKeys(dynamic_cast<AOMD::mVertex *>(e->get(0, node)), keys);
            }
            for (int edge = 0; edge < 12; ++edge)
            {
               getEdgeKeys(dynamic_cast<AOMD::mEdge *>(e->get(1, edge)), keys, 0);  // do not reorder keys
            }
            for (int face = 0; face < 6; ++face)
            {
               getQuadKeys(dynamic_cast<AOMD::mFace *>(e->get(2, face)), keys, 0);  // do not reorder keys
            }

            getHexKeys(dynamic_cast<AOMD::mHex *>(e), keys);
         }
         else
            throw;
         break;
      }
      default:
      {
         throw;
      }
   }
}

void xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH_deprecated::getFcts(const mEntity &e, femFcts &approscal)
{
   // e.print();

   switch (e.getLevel())
   {
      case 0:
      {
         // pour le cas ou la fct n'est pas interpolante ? -> ca devrait tjrs valoir 1 non ?
         approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, 0)));
         break;
      }
      case 1:
      {
         for (int node = 0; node < 2; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated(order, 0, i)));
         }
         break;
      }
      case 2:
      {
         for (int node = 0; node < 4; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 4; ++edge)
         {
            const int swap = computeSwapCaseQuadEdge<xSzaboSwapTable>(dynamic_cast<const AOMD::mFace &>(e), edge);
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(
                   shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated(order, edge, i, swap)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendre_deprecated(order, 0, i)));
         }
         break;
      }
      case 3:
      {
         for (int node = 0; node < 8; ++node)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexHierarchicalLegendreQH(order, node)));
         }
         for (int edge = 0; edge < 12; ++edge)
         {
            const int swap = computeSwapCaseHexEdge<xSzaboSwapTable>(dynamic_cast<const AOMD::mHex &>(e), edge);
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionEdgeQH(order); ++i)
            {
               approscal.push_back(
                   shapeFctPtr(new xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated(order, edge, i, swap)));
            }
         }
         for (int face = 0; face < 6; ++face)
         {
            const int swap = computeSwapCaseHexFace<xSzaboSwapTable>(dynamic_cast<const AOMD::mHex &>(e), face);
            for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionQuad(order); ++i)
            {
               approscal.push_back(
                   shapeFctPtr(new xApproxFunctionScalarQuadHierarchicalLegendre_deprecated(order, face, i, swap, true)));
            }
         }
         for (int i = 0; i < tensorSpaceDefinition::nbShapeFunctionHex(order); ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarHexHierarchicalLegendre(order, i)));
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}

//---------------------------------------------------

}  // end of namespace xfem
