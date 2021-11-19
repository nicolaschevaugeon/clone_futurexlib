/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XSPACEPOLYNOMIAL_QH_
#define _XSPACEPOLYNOMIAL_QH_
// stl include
#include <iosfwd>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Trellis include
#include "mEdge.h"
#include "mFace.h"
#include "mHex.h"

// xfem include
#include "xApproxFunctionHighOrder.h"
#include "xApproxFunctionHighOrderQH.h"
#include "xSpace.h"

namespace xfem
{
template <class SWAPTABLE>
int computeSwapCaseQuadEdge(const AOMD::mFace &quad, int edgenumber)
{
   return quad.get(1, edgenumber)->get(0, 0) != quad.get(0, SWAPTABLE::Fev[edgenumber][0]);
}

template <class SWAPTABLE>
int computeSwapCaseHexEdge(const AOMD::mHex &hex, int edgenumber)
{
   AOMD::mEntity *edge = hex.get(1, edgenumber);
   const int *hexnodeedge = SWAPTABLE::Hev[edgenumber];
   AOMD::mEntity *v0edge = edge->get(0, 0);
   AOMD::mEntity *v0hex = hex.get(0, hexnodeedge[0]);
   return (v0edge != v0hex);
}

// On fait l'hypothese que les quads ne sont pas retournes...
//  *--------------*      et pas ... *--------------*
//  |              +                    |         +
//  |              +                       |    +
//  |              +                         |
//  |              +                      +     |
//  |              +                   +           |
//  *--------------*                 *--------------*
//
//
template <class SWAPTABLE>
int computeSwapCaseHexFace(const AOMD::mHex &hex, int facenumber)
{
   AOMD::mEntity *face = hex.get(2, facenumber);
   const int *hexnodef = SWAPTABLE::Hfv[facenumber];
   AOMD::mEntity *v0hex = hex.get(0, hexnodef[0]);
   AOMD::mEntity *v1hex = hex.get(0, hexnodef[1]);

   AOMD::mEntity *v0face = face->get(0, 0);
   AOMD::mEntity *v1face = face->get(0, 1);
   AOMD::mEntity *v2face = face->get(0, 2);
   AOMD::mEntity *v3face = face->get(0, 3);
   //  std::cout<<"Face "<<v0face->getId()<<" "<<v1face->getId()<<" "<<v2face->getId()<<" "<<v3face->getId()<<std::endl;

   if (v0hex == v0face)
   {
      if (v1hex == v1face) return 0;
      return 1;  //=v4face
   }
   if (v0hex == v1face)
   {
      if (v1hex == v2face) return 2;
      return 3;  //=v5face
   }
   if (v0hex == v2face)
   {
      if (v1hex == v3face) return 4;
      return 5;  //=v6face
   }
   if (v0hex == v3face)
   {
      if (v1hex == v0face) return 6;
      return 7;  // v7face
   }
   // If none of the above was true
   return -1;
}

/// Class for high order shape function
/*! Only non-simplex element are considered here
      This is a pure virtual class.
      All derived class must implement getFcts, which is call by getKeysAndFunction
      and set the protected variable space_name at construction.
      SPACETYPE gives the definition of the space (tensor or trunk space)
   */

template <class SPACETYPE>
class xSpacePolynomialQH : public xSpaceRegular
{
  public:
   xSpacePolynomialQH(const std::string &space_name, const std::string &physical_name, TensorialType_t c, const int p);
   void getKeys(AOMD::mEntity *e, femKeys *keys) override;
   void getKeysAndFcts(AOMD::mEntity *e, femKeys *keys, femFcts *appro) override;

  protected:
   /// getFcts must push in approscal the scalar xApproxFunction that constitute the space, in the same order as the order of the
   /// keys, for all derived classes
   virtual void getFcts(const AOMD::mEntity &e, femFcts &approscal) = 0;
   void getNodeKeys(AOMD::mVertex *e, femKeys *keys);
   void getEdgeKeys(AOMD::mEdge *e, femKeys *keys, int swap);
   void getQuadKeys(AOMD::mFace *e, femKeys *keys, int swap);
   void getHexKeys(AOMD::mHex *e, femKeys *keys);
   const int order;
   const std::string space_name;
   const int nb_shape_edge;
   const int nb_shape_edge_tot;
   const int nb_shape_quad;
   const int nb_shape_quad_tot;
   const int nb_shape_hex;
   const int nb_shape_hex_tot;
   xValKey::ids_size_t node_geom_id;
   std::array<std::vector<xValKey::ids_size_t>, 2> edge_geom_id;
   std::array<std::vector<xValKey::ids_size_t>, 8> quad_geom_id;
   std::vector<xValKey::ids_size_t> hex_geom_id;

   // RHAAAA !
   std::array<std::vector<int>, 2> edge_orientation_swaps;
   std::array<std::vector<int>, 8> quad_orientation_swaps;

  private:
   void create_geom_id_node_qh(const std::string &space_name, int order, xValKey::ids_size_t &node_geom_id);
   void create_geom_id_edge_qh(const std::string &space_name, int order,
                               std::array<std::vector<xValKey::ids_size_t>, 2> &edge_geom_id);
   void create_geom_id_quad(const std::string &space_name, int order,
                            std::array<std::vector<xValKey::ids_size_t>, 8> &quad_geom_id);
   void create_geom_id_hex(const std::string &space_name, int order, std::vector<xValKey::ids_size_t> &hex_geom_id);

   void create_orientation_swaps_edge(int order);
   void create_orientation_swaps_quad(int order);
};

template <class SPACETYPE>
xSpacePolynomialQH<SPACETYPE>::xSpacePolynomialQH(const std::string &_space_name, const std::string &physical_name,
                                                  TensorialType_t c, const int p)
    : xSpaceRegular(physical_name, c),
      order(p),
      space_name(_space_name),
      nb_shape_edge(SPACETYPE::nbShapeFunctionEdgeQH(p)),
      nb_shape_edge_tot(SPACETYPE::nbShapeFunctionEdgeTotalQH(p)),
      nb_shape_quad(SPACETYPE::nbShapeFunctionQuad(p)),
      nb_shape_quad_tot(SPACETYPE::nbShapeFunctionQuadTotal(p)),
      nb_shape_hex(SPACETYPE::nbShapeFunctionHex(p)),
      nb_shape_hex_tot(SPACETYPE::nbShapeFunctionHexTotal(p))
{
   if (order > UCHAR_MAX)
   {
      std::cout << "In file " << __FILE__ << " Line " << __LINE__ << " oder greater then " << UCHAR_MAX
                << " !! Implementation have to be revised !" << std::endl;
   }

   create_geom_id_node_qh(space_name, order, node_geom_id);
   create_geom_id_edge_qh(space_name, order, edge_geom_id);
   create_geom_id_quad(space_name, order, quad_geom_id);
   create_geom_id_hex(space_name, order, hex_geom_id);
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::getNodeKeys(AOMD::mVertex *e, femKeys *keys)
{
   keys->push_back(xValKey(Phys, node_geom_id, (AOMD::mEntity *)e));
   return;
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::getEdgeKeys(AOMD::mEdge *e, femKeys *keys, int swap)
{
   const std::vector<xValKey::ids_size_t> &edge_geom_id_s = edge_geom_id[swap];
   for (int i = 0; i < nb_shape_edge; ++i)
   {
      keys->push_back(xValKey(Phys, edge_geom_id_s[i], (AOMD::mEntity *)(e)));
   }
   return;
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::getQuadKeys(AOMD::mFace *e, femKeys *keys, int swap)
{
   const std::vector<xValKey::ids_size_t> &quad_geom_id_s = quad_geom_id[swap];
   for (int i = 0; i < nb_shape_quad; ++i)
   {
      keys->push_back(xValKey(Phys, quad_geom_id_s[i], (AOMD::mEntity *)(e)));
   }
   return;
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::getHexKeys(AOMD::mHex *e, femKeys *keys)
{
   for (int i = 0; i < nb_shape_hex; ++i)
   {
      keys->push_back(xValKey(Phys, hex_geom_id[i], (AOMD::mEntity *)(e)));
   }
   return;
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::getKeys(AOMD::mEntity *e, femKeys *keys)
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
         keys->reserve(keys->size() + SPACETYPE::nbShapeFunctionEdgeTotalQH(order));
         for (int node = 0; node < 2; ++node)
         {
            getNodeKeys(dynamic_cast<AOMD::mVertex *>(e->get(0, node)), keys);
         }
         if (nb_shape_edge) getEdgeKeys((AOMD::mEdge *)(e), keys, 0);
         break;
      }
      case 2:
      {
         if (e->getType() == AOMD::mEntity::QUAD)
         {
            keys->reserve(keys->size() + nb_shape_quad_tot);
            for (int node = 0; node < 4; ++node)
            {
               getNodeKeys(dynamic_cast<AOMD::mVertex *>(e->get(0, node)), keys);
            }
            for (int edge = 0; edge < 4; ++edge)
            {
               const int swap = computeSwapCaseQuadEdge<xSzaboSwapTable>(dynamic_cast<const AOMD::mFace &>(*e), edge);
               //                e->get(1, edge)->print();
               //                std::cout<<"2D Swap is "<<swap<<std::endl;
               getEdgeKeys((AOMD::mEdge *)(e->get(1, edge)), keys, swap);
            }
            if (nb_shape_quad) getQuadKeys((AOMD::mFace *)(e), keys, 0);
         }
         else
            throw;
         break;
      }
      case 3:
      {
         if (e->getType() == AOMD::mEntity::HEX)
         {
            keys->reserve(keys->size() + nb_shape_hex_tot);
            for (int node = 0; node < 8; ++node)
            {
               getNodeKeys(dynamic_cast<AOMD::mVertex *>(e->get(0, node)), keys);
            }
            if (nb_shape_edge)
            {
               for (int edge = 0; edge < 12; ++edge)
               {
                  const int swap = computeSwapCaseHexEdge<xSzaboSwapTable>(dynamic_cast<const AOMD::mHex &>(*e), edge);
                  //                    e->get(1, edge)->print();
                  //                    std::cout<<"Edge Swap is "<<swap<<std::endl;
                  getEdgeKeys(dynamic_cast<AOMD::mEdge *>(e->get(1, edge)), keys, swap);
               }
            }

            if (nb_shape_quad)
            {
               for (int face = 0; face < 6; ++face)
               {
                  const int swap = computeSwapCaseHexFace<xSzaboSwapTable>(dynamic_cast<const AOMD::mHex &>(*e), face);
                  //                    e->get(2, face)->print();
                  //                    std::cout<<"Face Swap is "<<swap<<std::endl;
                  getQuadKeys(dynamic_cast<AOMD::mFace *>(e->get(2, face)), keys, swap);
               }
            }

            if (nb_shape_hex) getHexKeys(dynamic_cast<AOMD::mHex *>(e), keys);
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

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::getKeysAndFcts(AOMD::mEntity *e, femKeys *keys, femFcts *appro)
{
   getKeys(e, keys);
   int nfunc = keys->size();
   appro->reserve(appro->size() + nfunc);
   femFcts approscal;
   approscal.reserve(nfunc);
   getFcts(*e, approscal);
   switch (TensorialType)
   {
      case SCALAR:
         appro->insert(appro->end(), approscal.begin(), approscal.end());
         break;
      case VECTOR_X:
         for (femFcts::iterator it = approscal.begin(); it != approscal.end(); ++it)
         {
            appro->push_back(shapeFctPtr(new xApproxFunctionVector(*it, 0)));
         }
         break;
      case VECTOR_Y:
         for (femFcts::iterator it = approscal.begin(); it != approscal.end(); ++it)
         {
            appro->push_back(shapeFctPtr(new xApproxFunctionVector(*it, 1)));
         }
         break;
      case VECTOR_Z:
         for (femFcts::iterator it = approscal.begin(); it != approscal.end(); ++it)
         {
            appro->push_back(shapeFctPtr(new xApproxFunctionVector(*it, 2)));
         }
   }
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::create_geom_id_node_qh(const std::string &space_name, int order,
                                                           xValKey::ids_size_t &node_geom_id)
{
   std::ostringstream name;
   name << space_name << "_Order_" << order << "_Node";
   node_geom_id = xKeyInfo::getGeomId(name.str());
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::create_geom_id_edge_qh(const std::string &space_name, int order,
                                                           std::array<std::vector<xValKey::ids_size_t>, 2> &edge_geom_id)
{
   int nb_shape_edge = SPACETYPE::nbShapeFunctionEdgeQH(order);
   edge_geom_id[0].resize(nb_shape_edge);
   edge_geom_id[1].resize(nb_shape_edge);
   for (int i = 0; i < nb_shape_edge; ++i)
   {
      std::ostringstream name;
      name << space_name << "_Order_" << order << "_Edge_" << i + 1;
      edge_geom_id[0][i] = xKeyInfo::getGeomId(name.str());
   }
   std::copy(edge_geom_id[0].rbegin(), edge_geom_id[0].rend(), edge_geom_id[1].begin());
}

//// Specialization for tensor product space (more code could be factorized with trunk space)
// template<>
// void  xSpacePolynomialQH<tensorSpaceDefinition>::create_geom_id_quad(const std::string& space_name, int order, std::array <
// std::vector< xValKey::ids_size_t > , 8>  &quad_geom_id ){
//    const int nb_shape_quad = tensorSpaceDefinition::nbShapeFunctionQuad(order);
//    for_each(quad_geom_id.begin(), quad_geom_id.end(),
//             [&nb_shape_quad]( std::vector< int > & quad_geom_id ) {
//        quad_geom_id.clear();    quad_geom_id.resize(nb_shape_quad);
//    } );

//    if (nb_shape_quad){

//        std::ostringstream name;
//        std::string zero_string;
//        int index =0;
//        for (int j = 0; j < order-1; ++j){
//            for (int i = 0; i < order -1 ;  ++i ){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TENSORPROD_" << i+1 << "_" << j+1 ;
//                quad_geom_id[0][index]=xKeyInfo::getGeomId(name.str());
//                ++index;
//            }
//        }

//        const int nt = (order - 1) ;
//        const std::array<std::function<int(int &, int &)>, 8 > swapRenums = {
//                [nt](int &i, int &j){return i + nt * j;},
//        [nt](int &i, int &j){return j + nt * i;},
//        [nt](int &i, int &j){return (nt-1-j) + nt * i;},//2 et 6 inverses -> bug sur mon cahier ?
//        [nt](int &i, int &j){return (nt-1-i) + nt * j;},
//        [nt](int &i, int &j){return (nt-1-i) + nt * (nt-1-j);},
//        [nt](int &i, int &j){return (nt-1-j) + nt * (nt-1-i);},
//        [nt](int &i, int &j){return j + nt * (nt - 1 - i);},//2 et 6 inverses -> bug sur mon cahier ?
//        [nt](int &i, int &j){return i + nt * (nt-1-j);},
//    };
//        for (int swap = 1; swap < 8; ++swap){
//            std::vector<int> & quad_geom_id_s = quad_geom_id[swap];

//            int idx =0;
//            for (int j = 0; j < order-1; ++j){
//                for (int i = 0; i < order -1 ;  ++i ){
//                    quad_geom_id_s[idx] = quad_geom_id[0][swapRenums[swap](i,j)];
//                    //                    std::cout<<swap<<" "<<i<<" "<<j<<" "<<swapRenums[swap](i,j)<< " "<< i + nt * j<< "
//                    "<<j + nt * i  <<std::endl;
//                    ++idx;
//                }
//            }
//        }

//#if 0
//        //Debug swaps
//        for (int swap = 0; swap < 8; ++swap){
//            auto keys = quad_geom_id[swap];
//            std::cout<<"Generated keys for swap no "<<swap<<std::endl;
//            for(int &key : keys){
//                std::cout<<xKeyInfo::getGeomName(key)<<std::endl;
//            }
//        }
//#endif

//    }
//}

//// Specialization for trunk space (more code could be factorized with tensor product space)
//// TODO: Do this properly rather than reconstructing the swaped keys by hand...
// template<>
// void  xSpacePolynomialQH<trunkSpaceDefinition>::create_geom_id_quad(const std::string& space_name, int order, std::array <
// std::vector< xValKey::ids_size_t > , 8>  &quad_geom_id ){

//    const int nb_shape_quad = trunkSpaceDefinition::nbShapeFunctionQuad(order);
//    for_each(quad_geom_id.begin(), quad_geom_id.end(),
//             [&nb_shape_quad]( std::vector< int > & quad_geom_id ) {
//        quad_geom_id.clear();    quad_geom_id.resize(nb_shape_quad);
//    } );

//    if (nb_shape_quad){

//        std::ostringstream name;
//        std::string zero_string;
//        int index =0;
//        for (int j = 0; j < order - 3; ++j){
//            for (int i = j; i < order - 3; ++i){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[0][index]=xKeyInfo::getGeomId(name.str());
//                ++index;
//            }
//        }

//        //swap 1: swap i and j
//        int idx =0;
//        for (int i = 0; i < order - 3; ++i){
//            for (int j = i; j < order - 3; ++j){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[1][idx] = xKeyInfo::getGeomId(name.str());
//                ++idx;
//            }
//        }

//        //swap 2: swap i and j && invert i
//        idx =0;
//        for (int i = order -4; i >= 0; --i){
//            for (int j = i; j < order - 3; ++j){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[2][idx] = xKeyInfo::getGeomId(name.str());
//                ++idx;
//            }
//        }

//        //swap 3: invert i, j unchanged
//        idx =0;
//        for (int j = 0; j < order-3; ++j){
//            for (int i = order - 4; i >= j; --i){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[3][idx] = xKeyInfo::getGeomId(name.str());
//                ++idx;
//            }
//        }

//        //swap 4:  invert i, invert j
//        idx =0;
//        for (int j = order - 4; j >= 0; --j){
//            for (int i = order - 4; i >= j; --i){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[4][idx] = xKeyInfo::getGeomId(name.str());
//                ++idx;
//            }
//        }

//        //swap 5:  swap i and j && invert i and j
//        idx =0;
//        for (int i = order - 4; i >= 0; --i){
//            for (int j = order - 4; j >= i; --j){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[5][idx] = xKeyInfo::getGeomId(name.str());
//                ++idx;
//            }
//        }

//        //swap 6:  swap i and j && invert j
//        idx =0;
//        for (int i = 0; i < order - 3; ++i){
//            for (int j = order - 4; j >= i; --j){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[6][idx] = xKeyInfo::getGeomId(name.str());
//                ++idx;
//            }
//        }

//        //swap 7: invert j, i unchanged
//        idx =0;
//        for (int j = order -4; j >= 0; --j){
//            for (int i = j; i < order - 3; ++i){
//                name.str(zero_string);
//                name << space_name<< "_Order_"<< order <<"_QUAD_TRUNK_" << i+1 << "_" << j+1 ;
//                quad_geom_id[7][idx] = xKeyInfo::getGeomId(name.str());
//                ++idx;
//            }
//        }

//#if 0
//        //Debug swaps
//        for (int swap = 0; swap < 8; ++swap){
//            auto keys = quad_geom_id[swap];
//            std::cout<<"Generated keys for swap no "<<swap<<std::endl;
//            for(int &key : keys){
//                std::cout<<xKeyInfo::getGeomName(key)<<std::endl;
//            }
//        }
//        //        throw;
//#endif

//    }
//}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::create_geom_id_hex(const std::string &space_name, int order,
                                                       std::vector<xValKey::ids_size_t> &hex_geom_id)
{
   int nb_shape_hex = SPACETYPE::nbShapeFunctionHex(order);
   hex_geom_id.resize(nb_shape_hex);
   for (int i = 0; i < nb_shape_hex; ++i)
   {
      std::ostringstream name;
      name << space_name << "_Order_" << order << "_Hex_" << i;
      hex_geom_id[i] = xKeyInfo::getGeomId(name.str());
   }
}

template <class SPACETYPE>
void xSpacePolynomialQH<SPACETYPE>::create_orientation_swaps_edge(int order)
{
   const int nb_shape_edge = SPACETYPE::nbShapeFunctionEdgeQH(order);

   for_each(edge_orientation_swaps.begin(), edge_orientation_swaps.end(),
            [&nb_shape_edge](std::vector<int> &edge_orientation_swap) {
               edge_orientation_swap.clear();
               edge_orientation_swap.resize(nb_shape_edge);
               std::fill(edge_orientation_swap.begin(), edge_orientation_swap.end(), 1.);
            });

   // Swap 0 : no orientation change (edge_orientation_swaps[0] already filled with ones)

   // Swap 1 : change orientation only for odd shape functions
   for (int iShapeFunc = 1; iShapeFunc < nb_shape_edge; iShapeFunc += 2)
   {
      edge_orientation_swaps[1][iShapeFunc] = -1;
   }

   return;
}

/// \class xSpacePolynomialLagrangeQH
/// \brief{Lagrange polynomial tensor space for Hex/Quads.}
class xSpacePolynomialLagrangeQH : public xSpacePolynomialQH<tensorSpaceDefinition>
{
  public:
   xSpacePolynomialLagrangeQH(const std::string &a, TensorialType_t c, const int p);

  private:
   void getFcts(const AOMD::mEntity &e, femFcts &approscal) override;
};

/// \class xSpacePolynomialBernsteinQH
/// \brief{Bernstein polynomial tensor space for Hex/Quads.}
/// \bug{There may still have some troubles for non-zero Dirichlet B.C.}
class xSpacePolynomialBernsteinQH : public xSpacePolynomialQH<tensorSpaceDefinition>
{
  public:
   xSpacePolynomialBernsteinQH(const std::string &a, TensorialType_t c, const int p);

  private:
   void getFcts(const AOMD::mEntity &e, femFcts &approscal) override;
};

/// \class xSpacePolynomialHierarchicalLegendreTrunkSpaceQH
/// \brief{Hierarchical Legendre polynomial space for Hex/Quads (TRUNK SPACE).}
class xSpacePolynomialHierarchicalLegendreTrunkSpaceQH : public xSpacePolynomialQH<trunkSpaceDefinition>
{
  public:
   xSpacePolynomialHierarchicalLegendreTrunkSpaceQH(const std::string &a, TensorialType_t c, const int p);

  private:
   void getFcts(const AOMD::mEntity &e, femFcts &approscal) override;
};

/// \class xSpacePolynomialHierarchicalLegendreTrunkSpaceQH
/// \brief{Hierarchical Legendre polynomial space for Hex/Quads (TENSOR PRODUCT SPACE).}
class xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH : public xSpacePolynomialQH<tensorSpaceDefinition>
{
  public:
   xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH(const std::string &a, TensorialType_t c, const int p);

  private:
   void getFcts(const AOMD::mEntity &e, femFcts &approscal) override;
   //        void getKeys(AOMD::mEntity* e, femKeys* keys) override;//TEST
};

//    //TEST
//    void  xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH::getKeys(AOMD::mEntity* e, femKeys* keys) {
//        switch (e->getLevel()) {
//        case 0 : {
//            keys->reserve(keys->size()+1);
//            getNodeKeys(dynamic_cast< AOMD::mVertex *> (e), keys);
//            break;
//        }
//        case 1: {
//            keys->reserve(keys->size()+ tensorSpaceDefinition::nbShapeFunctionEdgeTotalQH(order) );
//            for (int node = 0; node <2 ; ++node) {
//                getNodeKeys(dynamic_cast<AOMD::mVertex *> (e->get(0, node)), keys);
//            }
//            if(nb_shape_edge) getEdgeKeys((AOMD::mEdge *) (e), keys, 0);
//            break;
//        }
//        case 2: {
//            if(e->getType() == AOMD::mEntity::QUAD) {
//                keys->reserve(keys->size()+ nb_shape_quad_tot );
//                for (int node=0; node < 4; ++node) {
//                    getNodeKeys(dynamic_cast<AOMD::mVertex *> (e->get(0, node)), keys);
//                }
//                e->print();
//                for (int edge=0; edge < 4; ++edge) {

//                    const int swap = computeSwapCaseQuadEdge<xSzaboSwapTable>(dynamic_cast < const AOMD::mFace & > ( *e ),
//                    edge); std::cout<<"edge "<<edge<<" "; e->get(1, edge)->print(); std::cout<<"--2D Swap is "<<swap<<std::endl;
//                    getEdgeKeys((AOMD::mEdge *) (e->get(1, edge)), keys, 0);
//                }
//                if(nb_shape_quad) getQuadKeys((AOMD::mFace *) (e), keys, 0);
//            }
//            else throw;
//            break;
//        }
//        case 3: {
//            if(e->getType() == AOMD::mEntity::HEX) {
//                keys->reserve(keys->size()+ nb_shape_hex_tot );
//                for (int node=0; node < 8; ++node) {
//                    getNodeKeys(dynamic_cast<AOMD::mVertex *> (e->get(0, node)), keys);
//                }
//                if(nb_shape_edge){
//                    for (int edge=0; edge < 12; ++edge) {
////                        const int swap = computeSwapCaseHexEdge<xSzaboSwapTable>( dynamic_cast < const AOMD::mHex & > ( *e ),
/// edge); /                        e->get(1, edge)->print(); /                        std::cout<<"--Edge Swap is
///"<<swap<<std::endl;

////                        e->print();
////                        printf("%i %i %i %i %i %i %i %i\n", e->get(0,0)->getId(), e->get(0,1)->getId(), e->get(0,2)->getId(),
/// e->get(0,3)->getId(), /                               e->get(0,4)->getId(), e->get(0,5)->getId(), e->get(0,6)->getId(),
/// e->get(0,7)->getId()); /                        std::cout<<"look for edge "<<edge<<" "; e->get(1, edge)->print(); /
/// std::cout<<"ids e "<<e->get(1, edge)->get(0,0)->getId()<<" "<<e->get(1, edge)->get(0,1)->getId()<<std::endl;
//                        getEdgeKeys(dynamic_cast<AOMD::mEdge *> (e->get(1, edge)), keys, 0);
//                    }
//                }

//                if(nb_shape_quad){
//                    for (int face=0; face < 6; ++face) {
//                        const int swap = computeSwapCaseHexFace<xSzaboSwapTable>( dynamic_cast < const AOMD::mHex & > ( *e ),
//                        face); e->get(2, face)->print(); std::cout<<"--Face Swap is "<<swap<<std::endl;
//                        getQuadKeys(dynamic_cast<AOMD::mFace *> (e->get(2, face)), keys, swap);
//                    }
//                }

//                if (nb_shape_hex) getHexKeys(dynamic_cast<AOMD::mHex *> (e), keys);
//            }
//            else throw;
//            break;
//        }
//        default: {
//            throw;
//        }
//        }
//    }

/// \class xSpacePolynomialHierarchicalLegendreTrunkSpaceQH
/// \brief{Hierarchical Legendre polynomial space for Hex/Quads (TRUNK SPACE).}
//    class xSpacePolynomialHierarchicalLegendreTrunkSpaceQH_deprecated : public xSpacePolynomialQH<trunkSpaceDefinition> {
//    public:
//        xSpacePolynomialHierarchicalLegendreTrunkSpaceQH_deprecated(const std::string& a, TensorialType_t c, const int  p);
//    private:
//        void getFcts(const AOMD::mEntity& e, femFcts& approscal) override;
//    };

/// \class xSpacePolynomialHierarchicalLegendreTrunkSpaceQH
/// \brief{Hierarchical Legendre polynomial space for Hex/Quads (TENSOR PRODUCT SPACE).}
class xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH_deprecated : public xSpacePolynomialQH<tensorSpaceDefinition>
{
  public:
   xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH_deprecated(const std::string &a, TensorialType_t c, const int p);

  private:
   void getFcts(const AOMD::mEntity &e, femFcts &approscal) override;
   void getKeys(AOMD::mEntity *e, femKeys *keys) override;  // TEST
};

}  // namespace xfem

#endif
