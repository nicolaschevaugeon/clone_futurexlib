/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xSpacePolynomial.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "mEdge.h"
#include "mFace.h"
#include "mTet.h"
#include "xFiniteElement.h"
#include "xNonLocalInfoForKeysAndFcts.h"
#include "xValKey.h"

using namespace xfem;
using AOMD::mEdge;
using AOMD::mEntity;
using AOMD::mFace;
using AOMD::mTet;
using AOMD::mVertex;
using std::cout;
using std::endl;
using xfem::xKeyInfo;

namespace xfem
{
int computeSwapCaseTriEdge(const mFace &tri, int edgenumber)
{
   return tri.get(1, edgenumber)->get(0, 0) != tri.get(0, edgenumber);
}

int computeSwapCaseTetEdge(const mTet &tet, int edgenumber)
{
   AOMD::mEntity *edge = tet.get(1, edgenumber);
   const int *tetnodeedge = xSimplexTet::Tev[edgenumber];
   mEntity *v0edge = edge->get(0, 0);
   mEntity *v0tet = tet.get(0, tetnodeedge[0]);
   return (v0edge != v0tet);
}

int computeSwapCaseTetFace(const mTet &tet, int facenumber)
{
   AOMD::mEntity *face = tet.get(2, facenumber);
   const int *tetnodef = xSimplexTet::Tfv[facenumber];
   AOMD::mEntity *v0tet = tet.get(0, tetnodef[0]);
   AOMD::mEntity *v1tet = tet.get(0, tetnodef[1]);
   AOMD::mEntity *v0face = face->get(0, 0);
   AOMD::mEntity *v1face = face->get(0, 1);
   AOMD::mEntity *v2face = face->get(0, 2);
   if (v0tet == v0face)
   {
      if (v1tet == v1face) return 0;
      return 3;
   }
   if (v0tet == v1face)
   {
      if (v1tet == v2face) return 1;
      return 5;
   }
   if (v0tet == v2face)
   {
      if (v1tet == v0face) return 2;
      return 4;
   }
   // If none of the above was true
   return -1;
}
void create_geom_id_node(const std::string &space_name, int order, xValKey::ids_size_t &node_geom_id)
{
   std::ostringstream name;
   name << space_name << "_Order_" << order << "_Node";
   node_geom_id = xKeyInfo::getGeomId(name.str());
}

void create_geom_id_edge(const std::string &space_name, int order, std::array<std::vector<xValKey::ids_size_t>, 2> &edge_geom_id)
{
   int nb_shape_edge = nbShapeFunctionEdge(order);
   edge_geom_id[0].resize(nb_shape_edge);
   edge_geom_id[1].resize(nb_shape_edge);
   for (int i = 0; i < nb_shape_edge; ++i)
   {
      std::ostringstream name;
      name << space_name << "_Order_" << order << "_Edge_" << i + 1 << "_" << order - i;
      edge_geom_id[0][i] = xKeyInfo::getGeomId(name.str());
   }
   std::copy(edge_geom_id[0].rbegin(), edge_geom_id[0].rend(), edge_geom_id[1].begin());
}

void create_geom_id_tri(const std::string &space_name, int order, std::array<std::vector<xValKey::ids_size_t>, 6> &tri_geom_id)
{
   const int nb_shape_tri = nbShapeFunctionTri(order);
   for_each(tri_geom_id.begin(), tri_geom_id.end(), [&nb_shape_tri](std::vector<xValKey::ids_size_t> &tri_geom_id_s) {
      tri_geom_id_s.clear();
      tri_geom_id_s.resize(nb_shape_tri);
   });
   if (nb_shape_tri)
   {
      std::vector<int> keytri_ijk[3] = {std::vector<int>(nb_shape_tri), std::vector<int>(nb_shape_tri),
                                        std::vector<int>(nb_shape_tri)};
      std::ostringstream name;
      std::string zero_string;
      int index = 0;
      for (int i = 0; i < order - 2; ++i)
      {
         for (int j = 0; j < order - 2 - i; ++j)
         {
            int k = order - 3 - i - j;
            keytri_ijk[0][index] = i;
            keytri_ijk[1][index] = j;
            keytri_ijk[2][index] = k;
            name.str(zero_string);
            name << space_name << "_Order_" << order << "_Tri_" << i + 1 << "_" << j + 1 << "_" << k + 1;
            tri_geom_id[0][index] = xKeyInfo::getGeomId(name.str());
            ++index;
         }
      }
      int cum = 0;
      std::vector<int> keytri(order - 2);
      for (int i = 0; i < order - 2; ++i)
      {
         keytri[i] = cum;
         cum += order - 2 - i;
      }
      const int swapcase[6][3] = {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
      for (int swap = 1; swap < 6; ++swap)
      {
         std::vector<xValKey::ids_size_t> &tri_geom_id_s = tri_geom_id[swap];
         const std::vector<int> &keytri_i = keytri_ijk[swapcase[swap][0]];
         const std::vector<int> &keytri_j = keytri_ijk[swapcase[swap][1]];
         // const  std::vector<int> &keytri_k = keytri_ijk [swapcase[swap][2] ];
         for (int kk = 0; kk < nb_shape_tri; ++kk)
         {
            const int i = keytri_i[kk];
            const int j = keytri_j[kk];
            int kk0 = keytri[i] + j;
            tri_geom_id_s[kk] = tri_geom_id[0][kk0];
         }
      }
   }
}

void create_geom_id_tet(const std::string &space_name, int order, std::vector<xValKey::ids_size_t> &tet_geom_id)
{
   int nb_shape_tet = nbShapeFunctionTet(order);
   tet_geom_id.resize(nb_shape_tet);
   for (int i = 0; i < nb_shape_tet; ++i)
   {
      std::ostringstream name;
      name << space_name << "_Order_" << order << "_Tet_" << i;
      tet_geom_id[i] = xKeyInfo::getGeomId(name.str());
   }
}

std::map<int, std::array<xSpace::femFcts, 4>> xSpacePolynomialBernstein::staticmap_order_approstored;
std::map<int, std::array<xSpace::femFcts, 4>> xSpacePolynomialLagrange::staticmap_order_approstored;

xSpacePolynomial::xSpacePolynomial(const std::string &_space_name, const std::string &physical_name, TensorialType_t c,
                                   const int p)
    : xSpaceRegular(physical_name, c),
      order(p),
      nb_shape_edge(nbShapeFunctionEdge(p)),
      nb_shape_edge_tot(nbShapeFunctionEdgeTotal(p)),
      nb_shape_tri(nbShapeFunctionTri(p)),
      nb_shape_tri_tot(nbShapeFunctionTriTotal(p)),
      nb_shape_tet(nbShapeFunctionTet(p)),
      nb_shape_tet_tot(nbShapeFunctionTetTotal(p)),
      space_name(_space_name)

{
   if (order > UCHAR_MAX)
   {
      cout << "In file " << __FILE__ << " Line " << __LINE__ << " order greater then " << UCHAR_MAX
           << " !! Implementation have to be revised !" << endl;
   }

   create_geom_id_node(space_name, order, node_geom_id);
   create_geom_id_edge(space_name, order, edge_geom_id);
   create_geom_id_tri(space_name, order, tri_geom_id);
   create_geom_id_tet(space_name, order, tet_geom_id);
}

void xSpacePolynomial::getNodeKeys(mVertex *e, femKeys *keys)
{
   keys->push_back(xValKey(Phys, node_geom_id, (mEntity *)e));
   return;
}

void xSpacePolynomial::getEdgeKeys(mEdge *e, femKeys *keys, int swap)
{
   const std::vector<xValKey::ids_size_t> &edge_geom_id_s = edge_geom_id[swap];
   for (int i = 0; i < nb_shape_edge; ++i)
   {
      keys->push_back(xValKey(Phys, edge_geom_id_s[i], (mEntity *)(e)));
   }
   return;
}

void xSpacePolynomial::getTriKeys(mFace *e, femKeys *keys, int swap)
{
   const std::vector<xValKey::ids_size_t> &tri_geom_id_s = tri_geom_id[swap];
   for (int i = 0; i < nb_shape_tri; ++i)
   {
      keys->push_back(xValKey(Phys, tri_geom_id_s[i], (mEntity *)(e)));
   }
   return;
}

void xSpacePolynomial::getTetKeys(mTet *e, femKeys *keys)
{
   for (int i = 0; i < nb_shape_tet; ++i)
   {
      keys->push_back(xValKey(Phys, tet_geom_id[i], (mEntity *)(e)));
   }
   return;
}

void xSpacePolynomial::getKeys(mEntity *e, femKeys *keys)
{
   switch (e->getLevel())
   {
      case 0:
      {
         keys->reserve(keys->size() + 1);
         getNodeKeys(static_cast<mVertex *>(e), keys);
         break;
      }
      case 1:
      {
         keys->reserve(keys->size() + nb_shape_edge_tot);
         for (int node = 0; node < 2; ++node)
         {
            getNodeKeys(static_cast<mVertex *>(e->get(0, node)), keys);
         }
         if (nb_shape_edge) getEdgeKeys(static_cast<mEdge *>(e), keys, 0);
         break;
      }
      case 2:
      {
         keys->reserve(keys->size() + nb_shape_tri_tot);
         if (e->getType() == mEntity::TRI)
         {
            for (int node = 0; node < 3; ++node)
            {
               getNodeKeys(static_cast<mVertex *>(e->get(0, node)), keys);
            }
            if (nb_shape_edge)
               for (int edge = 0; edge < 3; ++edge)
               {
                  const int swap = computeSwapCaseTriEdge(dynamic_cast<const mFace &>(*e), edge);
                  getEdgeKeys(static_cast<mEdge *>(e->get(1, edge)), keys, swap);
               }
            if (nb_shape_tri) getTriKeys(static_cast<mFace *>(e), keys, 0);
         }
         else
            throw;
         break;
      }
      case 3:
      {
         if (e->getType() == mEntity::TET)
         {
            keys->reserve(keys->size() + nb_shape_tet_tot);
            for (int node = 0; node < 4; ++node)
            {
               getNodeKeys(static_cast<mVertex *>(e->get(0, node)), keys);
            }
            if (nb_shape_edge)
               for (int edge = 0; edge < 6; ++edge)
               {
                  const int swap = computeSwapCaseTetEdge(dynamic_cast<const mTet &>(*e), edge);
                  getEdgeKeys(static_cast<mEdge *>(e->get(1, edge)), keys, swap);
               }
            if (nb_shape_tri)
               for (int face = 0; face < 4; ++face)
               {
                  const int swap = computeSwapCaseTetFace(dynamic_cast<const mTet &>(*e), face);
                  getTriKeys(static_cast<mFace *>(e->get(2, face)), keys, swap);
               }
            if (nb_shape_tet) getTetKeys(static_cast<mTet *>(e), keys);
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

void xSpacePolynomial::getKeysAndFcts(mEntity *e, femKeys *keys, femFcts *appro)
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

xSpacePolynomialLagrange::xSpacePolynomialLagrange(const std::string &a, TensorialType_t c, const int p)
    : xSpacePolynomial("Lagrange", a, c, p), approstored(xSpacePolynomialLagrange::staticmap_order_approstored[p])
{
   assert(p <= XAPPROXFUNCTIONLAGRANGEMAXORDER);
}

void xSpacePolynomialLagrange::getFcts(const mEntity &e, femFcts &approscal)
{
   switch (e.getLevel())
   {
      case 0:
      {
         if (approstored[0].empty())
         {
            approstored[0].push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, 0)));
         }
         approscal.insert(approscal.end(), approstored[0].begin(), approstored[0].end());
         break;
      }
      case 1:
      {
         if (approstored[1].empty())
         {
            for (int node = 0; node < 2; ++node)
            {
               approstored[1].push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, node)));
            }
            for (int i = 0; i < nb_shape_edge; ++i)
            {
               approstored[1].push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, 0, i)));
            }
         }
         approscal.insert(approscal.end(), approstored[1].begin(), approstored[1].end());
         break;
      }
      case 2:
      {
         if (e.getType() == mEntity::TRI)
         {
            if (approstored[2].empty())
            {
               for (int node = 0; node < 3; ++node)
               {
                  approstored[2].push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, node)));
               }
               for (int edge = 0; edge < 3; ++edge)
               {
                  for (int i = 0; i < nb_shape_edge; ++i)
                  {
                     approstored[2].push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, edge, i)));
                  }
               }
               for (int i = 0; i < nb_shape_tri; ++i)
               {
                  approstored[2].push_back(shapeFctPtr(new xApproxFunctionScalarTriLagrange(order, 0, i)));
               }
            }
            approscal.insert(approscal.end(), approstored[2].begin(), approstored[2].end());
         }
         else
            throw;
         break;
      }
      case 3:
      {
         if (e.getType() == mEntity::TET)
         {
            if (approstored[3].empty())
            {
               for (int node = 0; node < 4; ++node)
               {
                  approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, node)));
               }
               for (int edge = 0; edge < 6; ++edge)
               {
                  for (int i = 0; i < nb_shape_edge; ++i)
                  {
                     approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, edge, i)));
                  }
               }
               for (int face = 0; face < 4; ++face)
               {
                  for (int i = 0; i < nb_shape_tri; ++i)
                  {
                     approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarTriLagrange(order, face, i)));
                  }
               }

               for (int i = 0; i < nb_shape_tet; ++i)
               {
                  approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarTetLagrange(order, i)));
               }
            }
            approscal.insert(approscal.end(), approstored[3].begin(), approstored[3].end());
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

xSpacePolynomialBernstein::xSpacePolynomialBernstein(const std::string &a, TensorialType_t c, const int p)
    : xSpacePolynomial("Bernstein", a, c, p), approstored(xSpacePolynomialBernstein::staticmap_order_approstored[p])
{
}

void xSpacePolynomialBernstein::getFcts(const mEntity &e, femFcts &approscal)
{
   switch (e.getLevel())
   {
      case 0:
      {
         if (approstored[0].empty())
         {
            approstored[0].push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(order, 0)));
         }
         approscal.insert(approscal.end(), approstored[0].begin(), approstored[0].end());
         break;
      }
      case 1:
      {
         if (approstored[1].empty())
         {
            for (int node = 0; node < 2; ++node)
            {
               approstored[1].push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(order, node)));
               // approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(1 , node )));
            }
            for (int i = 1; i < order; ++i)
            {
               approstored[1].push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernstein(i, order - i, 0)));
            }
         }
         approscal.insert(approscal.end(), approstored[1].begin(), approstored[1].end());
         break;
      }
      case 2:
      {
         if (e.getType() == mEntity::TRI)
         {
            if (approstored[2].empty())
            {
               for (int node = 0; node < 3; ++node)
               {
                  approstored[2].push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(order, node)));
                  // approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(1 , node )));
               }
               for (int edge = 0; edge < 3; ++edge)
               {
                  for (int i = 1; i < order; ++i)
                  {
                     approstored[2].push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernstein(i, order - i, edge)));
                  }
               }
               for (int i = 1; i <= (order - 2); ++i)
               {
                  for (int j = 1; j <= (order - i - 1); ++j)
                  {
                     approstored[2].push_back(shapeFctPtr(new xApproxFunctionScalarTriBernstein(i, j, order - i - j, 0)));
                  }
               }
            }
            approscal.insert(approscal.end(), approstored[2].begin(), approstored[2].end());
         }
         else
            throw;
         break;
      }
      case 3:
      {
         if (e.getType() == mEntity::TET)
         {
            if (approstored[3].empty())
            {
               for (int node = 0; node < 4; ++node)
               {
                  approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(order, node)));
                  // approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(1 , node )));
               }
               for (int edge = 0; edge < 6; ++edge)
               {
                  for (int i = 1; i < order; ++i)
                  {
                     approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernstein(i, order - i, edge)));
                  }
               }
               for (int face = 0; face < 4; ++face)
               {
                  for (int i = 1; i <= (order - 2); ++i)
                  {
                     for (int j = 1; j <= (order - 1 - i); ++j)
                     {
                        approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarTriBernstein(i, j, order - i - j, face)));
                     }
                  }
               }
               for (int i = 1; i <= (order - 3); ++i)
               {
                  for (int j = 1; j <= (order - 2 - i); ++j)
                  {
                     for (int k = 1; k <= (order - 1 - i - j); ++k)
                     {
                        approstored[3].push_back(shapeFctPtr(new xApproxFunctionScalarTetBernstein(i, j, k, order - i - j - k)));
                     }
                  }
               }
            }
            approscal.insert(approscal.end(), approstored[3].begin(), approstored[3].end());
         }
         else
            throw;
         break;
      }
      default:
         throw;
   }
}

xSpacePolynomialDiscontinuous::xSpacePolynomialDiscontinuous(const std::string &a, TensorialType_t c, const int _meshDimension,
                                                             const int p)
    : xSpaceRegular(a, c), order(p), meshDimension(_meshDimension)
{
}

void xSpacePolynomialDiscontinuous::getEntityKeys(const int level, const int AOMD_id, mEntity *ref, femKeys *keys)
{
   // sets the N keys for entity "AOMD_id" of level "level" in reference entity "ref"
   // depending on order and entity's level
   // with N = 1, p-1 or (p-1)*(p-2)/2 depending on the "level"
   // key format:  "DiscontinuousLagrange_ID_type_id_i"
   // with type=0 Node, 1 Edge, 2 Face or 3 Vol and id the "AOMD_id" of the treated entity in "ref" adjancy
   if (order == 0)
   {
      if (level == ref->getLevel())
         keys->push_back(xValKey(Phys, xKeyInfo::getGeomId(string("DiscontinuousConstant")), (mEntity *)(ref)));
      return;
   }
   {
      std::string basename("DiscontinuousLagrange");
      int nbdofs = (level == 1 ? (order - 1) : 1);
      nbdofs = (level == 2 ? ((order - 1) * (order - 2) / 2) : nbdofs);
      nbdofs = (level == 3 ? ((order - 1) * (order - 2) * (order - 3) / 6) : nbdofs);
      if (level == 1)
      {
         bool swap = false;
         // if (AOMD_id > 0) swap = computeSwapCaseTriEdge(dynamic_cast < const mFace & > ( *ref ), AOMD_id);
         if (!swap)
         {
            for (int i = 0; i < nbdofs; ++i)
            {
               std::stringstream name;
               name << basename << "_ID_" << level << "_" << AOMD_id << "_" << i;
               keys->push_back(xValKey(Phys, xKeyInfo::getGeomId(name.str()), (mEntity *)(ref)));
            }
         }
         if (swap)
         {
            for (int i = nbdofs - 1; i >= 0; --i)
            {
               std::stringstream name;
               name << basename << "_ID_" << level << "_" << AOMD_id << "_" << i;
               keys->push_back(xValKey(Phys, xKeyInfo::getGeomId(name.str()), (mEntity *)(ref)));
            }
         }
         return;
      }
      for (int i = 0; i < nbdofs; i++)
      {
         std::stringstream name;
         name << basename << "_ID_" << level << "_" << AOMD_id << "_" << i;
         keys->push_back(xValKey(Phys, xKeyInfo::getGeomId(name.str()), (mEntity *)(ref)));
      }
   }
   return;
}
// special function to get entities from key definition. Put here to keep constitancy with above definition in getEntityKeys
void xSpacePolynomialDiscontinuous::getEntitiesFromKey(mEntity **e, mEntity **ref, int *n, xValKey &key)
{
   int level, id;
   *ref = key.getEnti();
   string name(xKeyInfo::getGeomName(key.getGeom()));
   sscanf(name.c_str(), "DiscontinuousLagrange_ID_%d_%d_%d", &level, &id, n);
   if ((*ref)->getLevel() != level)
      *e = (*ref)->get(level, id);
   else
      *e = *ref;
   return;
}

void xSpacePolynomialDiscontinuous::getKeys(mEntity *e, femKeys *keys)
{
   switch (e->getLevel())
   {
      case 0:
      {
         cout << "WARNING: xSpaceDiscontinuousLag(getKeys) : cannot attribute a key to a single vertex, needs the element ! No "
                 "key returned !"
              << endl;
         throw;
         break;
      }
      case 1:
      {
         if (meshDimension != 1)
         {
            cout << "WARNING: xSpaceDiscontinuousLag(getKeys) : cannot attribute a key to a single edge for a 2D problem, needs "
                    "the upper face ! No key returned !"
                 << endl;
         }
         else
         {
            for (int node = 0; node < 2; ++node) getEntityKeys(0, node, e, keys);
            getEntityKeys(1, -1, e, keys);
         }
         break;
      }
      case 2:
      {
         if (meshDimension != 2)
         {
            cout << "WARNING: xSpaceDiscontinuousLag(getKeys) : cannot attribute a key to a single face for a 3D problem, needs "
                    "the upper element ! No key returned !"
                 << endl;
         }
         else
         {
            if (e->getType() == mEntity::TRI)
            {
               for (int node = 0; node < 3; ++node) getEntityKeys(0, node, e, keys);
               for (int edge = 0; edge < 3; ++edge) getEntityKeys(1, edge, e, keys);
               getEntityKeys(2, -1, e, keys);
            }
            else
            {
               cout << "xSpaceDiscontinuousLag(getKeys) 2D not yet implemented for element type mEntity::enum mType::"
                    << e->getType() << endl;
               throw;
            }
         }
         break;
      }
      case 3:
      {
         if (e->getType() == mEntity::TET)
         {
            for (int node = 0; node < 4; ++node) getEntityKeys(0, node, e, keys);
            for (int edge = 0; edge < 6; ++edge) getEntityKeys(1, edge, e, keys);
            for (int face = 0; face < 4; ++face) getEntityKeys(2, face, e, keys);
            getEntityKeys(3, -1, e, keys);
         }
         else
         {
            cout << "xSpaceDiscontinuousLag(getKeys) 3D not yet implemented for element type mEntity::enum mType::"
                 << e->getType() << endl;
            throw;
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}

xSpacePolynomialDiscontinuousLagrange::xSpacePolynomialDiscontinuousLagrange(const std::string &a, TensorialType_t c,
                                                                             const int _meshDimension, const int p)
    : xSpacePolynomialDiscontinuous(a, c, _meshDimension, p)
{
   assert(p <= XAPPROXFUNCTIONLAGRANGEMAXORDER);
}

void xSpacePolynomialDiscontinuousLagrange::getKeysAndFcts(mEntity *e, femKeys *keys, femFcts *appro)
{
   getKeys(e, keys);
   int nfunc = keys->size();
   appro->reserve(appro->size() + nfunc);
   femFcts approscal;
   approscal.reserve(nfunc);

   if (order == 0)
      approscal.push_back(shapeFctPtr(new xApproxFunctionConstant()));
   else
   {
      switch (e->getLevel())
      {
         case 0:
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, 0)));
            break;
         }
         case 1:
         {
            for (int node = 0; node < 2; ++node)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, node)));
            }
            for (int i = 0; i < (order - 1); ++i)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, 0, 0, i)));
            }
            break;
         }
         case 2:
         {
            if (e->getType() == mEntity::TRI)
            {
               for (int node = 0; node < 3; ++node)
               {
                  approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, node)));
               }
               if (order > 1)
               {
                  for (int edge = 0; edge < 3; ++edge)
                  {
                     bool swap = false;
                     if (order > 2) swap = computeSwapCaseTriEdge(dynamic_cast<const mFace &>(*e), edge);
                     if (!swap)
                     {
                        for (int i = 0; i < (order - 1); ++i)
                        {
                           approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, edge, i)));
                        }
                     }
                     else
                     {
                        for (int i = order - 2; i >= 0; --i)
                        {
                           approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, edge, i)));
                        }
                     }
                  }
               }
               if (order > 2)
               {
                  for (int i = 0; i < ((order - 1) * (order - 2) / 2); ++i)
                  {
                     approscal.push_back(shapeFctPtr(new xApproxFunctionScalarTriLagrange(order, 0, 0, i)));
                  }
               }
            }
            else
            {
               throw;
            }
            break;
         }
         case 3:
         {
            if (e->getType() == mEntity::TET)
            {
               for (int node = 0; node < 4; ++node)
               {
                  approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexLagrange(order, node)));
               }
               if (order > 1)
               {
                  for (int edge = 0; edge < 6; ++edge)
                  {
                     bool swap = false;
                     if (order > 2) swap = computeSwapCaseTetEdge(dynamic_cast<const mTet &>(*e), edge);
                     if (!swap)
                     {
                        for (int i = 0; i < (order - 1); ++i)
                        {
                           approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, edge, i)));
                        }
                     }
                     else
                     {
                        for (int i = order - 2; i >= 0; --i)
                        {
                           approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeLagrange(order, edge, i)));
                        }
                     }
                  }
               }
               if (order > 2)
               {
                  for (int face = 0; face < 4; ++face)
                  {
                     femFcts approtri;
                     int nb_shape_tri = nbShapeFunctionTri(order);
                     approtri.reserve(nb_shape_tri);

                     std::vector<int> keytri_ijk[3] = {std::vector<int>(nb_shape_tri), std::vector<int>(nb_shape_tri),
                                                       std::vector<int>(nb_shape_tri)};
                     int indexs = 0;
                     for (int i = 0; i < order - 2; ++i)
                     {
                        for (int j = 0; j < order - 2 - i; ++j)
                        {
                           int k = order - 3 - i - j;
                           keytri_ijk[0][indexs] = i;
                           keytri_ijk[1][indexs] = j;
                           keytri_ijk[2][indexs] = k;
                           approtri.push_back(shapeFctPtr(new xApproxFunctionScalarTriLagrange(order, face, indexs)));
                           ++indexs;
                        }
                     }
                     int cum = 0;
                     std::vector<int> keytri(order - 2);
                     for (int i = 0; i < order - 2; ++i)
                     {
                        keytri[i] = cum;
                        cum += order - 2 - i;
                     }
                     const int swap = computeSwapCaseTetFace(dynamic_cast<const mTet &>(*e), face);
                     const int swapcase[6][3] = {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
                     const std::vector<int> &keytri_i = keytri_ijk[swapcase[swap][0]];
                     const std::vector<int> &keytri_j = keytri_ijk[swapcase[swap][1]];
                     // const  std::vector<int> &keytri_k = keytri_ijk [swapcase[swap][2] ];
                     for (int kk = 0; kk < nb_shape_tri; ++kk)
                     {
                        const int i = keytri_i[kk];
                        const int j = keytri_j[kk];
                        int kk0 = keytri[i] + j;
                        approscal.push_back(approtri[kk0]);
                     }
                  }
                  /*
                  for (int i = 0; i < nbShapeFunctionTri(order); ++i)
                    {
                      approscal.push_back(shapeFctPtr(new xApproxFunctionScalarTriLagrange(order, face, i)));
                      }*/
               }

               for (int i = 0; i < nbShapeFunctionTet(order); ++i)
               {
                  approscal.push_back(shapeFctPtr(new xApproxFunctionScalarTetLagrange(order, i)));
               }
            }
            else
            {
               throw;
            }
            break;
         }
         default:
         {
            throw;
         }
      }
   }

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

xSpacePolynomialOctree::xSpacePolynomialOctree(const std::string &_space_name, const std::string &physical_name,
                                               TensorialType_t c, const int p)
    : xSpacePolynomial(_space_name, physical_name, c, p), nb_free(0), cur_nli(nullptr)
{
}

void xSpacePolynomialOctree::getKeys(mEntity *e, femKeys *keys)
{
   femKeys keysl;
   xSpacePolynomial::getKeys(e, &keysl);
   nb_free = keysl.size();
   cur_nli = getNonLocalInfo(*e);
   if (cur_nli)
   {
      other_master_keys.resize(0);
      nb_free = cur_nli->getAllKeysFromSlavesAndFree(Phys, &keysl, &other_master_keys);
   }
   keys->insert(keys->end(), keysl.begin(), keysl.end());
}

void xSpacePolynomialOctree::getKeysAndFcts(mEntity *e, femKeys *keys, femFcts *appro)
{
   femKeys keysl;
   getKeys(e, &keysl);
   int nfunc = keysl.size();
   appro->reserve(appro->size() + nfunc);
   femFcts approscal;
   approscal.reserve(nfunc);
   getFcts(*e, keysl, approscal);
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
   keys->insert(keys->end(), keysl.begin(), keysl.end());
}

xSpacePolynomialOctreeLagrange::xSpacePolynomialOctreeLagrange(const std::string &a, TensorialType_t c, const int p,
                                                               xNonLocalInfoGeneratorForKeysAndFcts *_extended_generator)
    : xSpacePolynomialOctree("Lagrange", a, c, p), extended_generator(_extended_generator)
{
   assert(p <= XAPPROXFUNCTIONLAGRANGEMAXORDER);
   std_space.reset(new xSpacePolynomialLagrange(a, SCALAR, p));
   _extended_generator->generateNonLocalInfoForKeysAndFcts(std_space, p);
   std_space.reset(new xSpacePolynomialLagrange(a, c, p));
}
const xNonLocalInfoForKeysAndFcts *xSpacePolynomialOctreeLagrange::getNonLocalInfo(const AOMD::mEntity &e) const
{
   return extended_generator->getNonLocalInfo(e);
}
void xSpacePolynomialOctreeLagrange::getFcts(mEntity &e, const femKeys &keys, femFcts &approscal)
{
   int nb_keys = keys.size();
   if (!cur_nli)
   {
      (std::static_pointer_cast<xSpacePolynomialLagrange>(std_space))->getFcts(e, approscal);
      return;
   }
   else
   {
      switch (e.getLevel())
      {
         case 0:
         case 1:
            break;
         case 2:
         {
            if (e.getType() == mEntity::TRI)
            {
               it_begin_node = e.begin(0);
               it_end_node = e.end(0);
               if (nb_shape_edge)
               {
                  it_begin_edge = e.begin(1);
                  it_end_edge = e.end(1);
               }
            }
            else
            {
               cout << "In file " << __FILE__ << " unknow 2D element type" << endl;
               throw;
            }
            break;
         }
         case 3:
         {
            if (e.getType() == mEntity::TET)
            {
               it_begin_node = e.begin(0);
               it_end_node = e.end(0);
               if (nb_shape_edge)
               {
                  it_begin_edge = e.begin(1);
                  it_end_edge = e.end(1);
                  if (nb_shape_tri)
                  {
                     it_begin_face = e.begin(2);
                     it_end_face = e.end(2);
                  }
               }
            }
            else
            {
               cout << "In file " << __FILE__ << " unknow 3D element type" << endl;
               throw;
            }
            break;
         }
         default:
         {
            cout << "In file " << __FILE__ << " unknow element type" << endl;
            throw;
         }
      }
      // Loop on keys
      int k = 0;
      femKeys::const_iterator it = keys.begin();
      femKeys::const_iterator itend = keys.end();
      femKeys::const_iterator itombeg = other_master_keys.begin();
      femKeys::const_iterator itomend = other_master_keys.end();
      xApproxFunction *approx_func;
      for (; it != itend; ++it)
      {
         std::vector<shapeFctPtr> vshape;

         // free and potentialy master key
         if (k < nb_free)
         {
            if (find(itombeg, itomend, *it) != itomend)
            {
               vshape.reserve(nb_keys);
               // this test is useless as it is a free key
               // just in case ...
               if ((approx_func = getFctFromKey(e, *it, 1.))) vshape.emplace_back(approx_func);
               const xNonLocalInfoForKeysAndFcts::femRelKeys &slaves = cur_nli->getCoefAndKeysForKey(*it);
               xNonLocalInfoForKeysAndFcts::femRelKeys::const_iterator it_pair = slaves.begin();
               xNonLocalInfoForKeysAndFcts::femRelKeys::const_iterator it_pair_end = slaves.end();
               for (; it_pair != it_pair_end; ++it_pair)
               {
                  if ((approx_func = getFctFromKey(e, it_pair->first, it_pair->second))) vshape.emplace_back(approx_func);
               }
               // at least their should be the free key but as a test is present ->assert
               assert(!vshape.empty());
               approscal.push_back(shapeFctPtr(new xApproxFunctionSummed(vshape)));
            }
            // single function : not test on returned value as it is a free key so it
            // exist in element
            else
            {
               // this test is useless as it is a free key
               // just in case ...
               if ((approx_func = getFctFromKey(e, *it, 1.))) approscal.emplace_back(approx_func);
            }
            ++k;
         }
         // master key outside entity
         else
         {
            vshape.reserve(nb_keys);
            const xNonLocalInfoForKeysAndFcts::femRelKeys &slaves = cur_nli->getCoefAndKeysForKey(*it);
            xNonLocalInfoForKeysAndFcts::femRelKeys::const_iterator it_pair = slaves.begin();
            xNonLocalInfoForKeysAndFcts::femRelKeys::const_iterator it_pair_end = slaves.end();
            for (; it_pair != it_pair_end; ++it_pair)
            {
               if ((approx_func = getFctFromKey(e, it_pair->first, it_pair->second))) vshape.emplace_back(approx_func);
            }
            // normaly as it have been identifyed as a master key by somme slave of the entity it should'nt be empty
            assert(!vshape.empty());
            approscal.push_back(shapeFctPtr(new xApproxFunctionSummed(vshape)));
         }
      }
   }
}

xApproxFunction *xSpacePolynomialOctreeLagrange::getFctFromKey(const mEntity &e, const xValKey &key, const double coef)
{
   const mEntity *e_it = key.getEnti();
   const xValKey::ids_size_t geomid = key.getGeom();
   switch (e.getLevel())
   {
      case 0:
      {
         if (&e == e_it)
            return (new xApproxFunctionScalarVertexLagrange(order, 0, coef));
         else
            return nullptr;
         break;
      }
      case 1:
      {
         if (geomid == node_geom_id)
         {
            if (e.get(0, 0) == e_it)
               return (new xApproxFunctionScalarVertexLagrange(order, 0, coef));
            else if (e.get(0, 1) == e_it)
               return (new xApproxFunctionScalarVertexLagrange(order, 1, coef));
            else
               return nullptr;
         }
         else
         {
            if (&e == e_it)
               return (new xApproxFunctionScalarEdgeLagrange(order, 0, getShapeIdFromGeom(geomid, edge_geom_id[0]), coef));
            else
               return nullptr;
         }
         break;
      }
      case 2:
      {
         if (e.getType() == mEntity::TRI)
         {
            switch (e_it->getLevel())
            {
               case 0:
               {
                  AOMD::mAdjacencyContainer::iter it_find = find(it_begin_node, it_end_node, e_it);
                  if (it_find != it_end_node)
                     return (new xApproxFunctionScalarVertexLagrange(order, (it_find - it_begin_node), coef));
                  else
                     return nullptr;
                  break;
               }
               case 1:
               {
                  AOMD::mAdjacencyContainer::iter it_find = find(it_begin_edge, it_end_edge, e_it);
                  if (it_find != it_end_edge)
                  {
                     int i_edge = it_find - it_begin_edge;
                     return (
                         new xApproxFunctionScalarEdgeLagrange(order, i_edge, getShapeIdFromGeom(geomid, edge_geom_id[0]), coef));
                  }
                  else
                     return nullptr;
                  break;
               }
               case 2:
               {
                  if (&e == e_it)
                  {
                     return (new xApproxFunctionScalarTriLagrange(order, 0, getShapeIdFromGeom(geomid, tri_geom_id[0])));
                  }
                  else
                     return nullptr;
                  break;
               }
               default:
                  throw;
            }
         }
         else
            throw;
         break;
      }
      case 3:
      {
         if (e.getType() == mEntity::TET)
         {
            switch (e_it->getLevel())
            {
               case 0:
               {
                  AOMD::mAdjacencyContainer::iter it_find = find(it_begin_node, it_end_node, e_it);
                  if (it_find != it_end_node)
                     return (new xApproxFunctionScalarVertexLagrange(order, (it_find - it_begin_node), coef));
                  else
                     return nullptr;
                  break;
               }
               case 1:
               {
                  AOMD::mAdjacencyContainer::iter it_find = find(it_begin_edge, it_end_edge, e_it);
                  if (it_find != it_end_edge)
                  {
                     int i_edge = it_find - it_begin_edge;
                     return (
                         new xApproxFunctionScalarEdgeLagrange(order, i_edge, getShapeIdFromGeom(geomid, edge_geom_id[0]), coef));
                  }
                  else
                     return nullptr;
                  break;
               }
               case 2:
               {
                  AOMD::mAdjacencyContainer::iter it_find = find(it_begin_face, it_end_face, e_it);
                  if (it_find != it_end_face)
                  {
                     int i_face = it_find - it_begin_face;
                     return (
                         new xApproxFunctionScalarTriLagrange(order, i_face, getShapeIdFromGeom(geomid, tri_geom_id[0]), coef));
                  }
                  else
                     return nullptr;
                  break;
               }
               case 3:
               {
                  if (&e == e_it)
                     return (new xApproxFunctionScalarTetLagrange(order, getShapeIdFromGeom(geomid, tet_geom_id)));
                  else
                     return nullptr;
                  break;
               }
               default:
                  throw;
            }
         }
         else
            throw;
         break;
      }
      default:
         throw;
   }
}

int xSpacePolynomialOctree::getShapeIdFromGeom(xValKey::ids_size_t geomid, const std::vector<xValKey::ids_size_t> &ids)
{
   return std::distance(ids.begin(), lower_bound(ids.begin(), ids.end(), geomid));
}

xSpacePolynomialVariableP::xSpacePolynomialVariableP(const std::string &_space_name, const std::string &physical_name,
                                                     TensorialType_t c, int pmax,
                                                     std::function<int(const mEntity *e)> &orderPredicate_)
    : xSpaceRegular(physical_name, c), ordermax(pmax), space_name(_space_name), orderPredicate(orderPredicate_)
{
   if (ordermax > UCHAR_MAX)
   {
      cout << "In file " << __FILE__ << " Line " << __LINE__ << " order greater then " << UCHAR_MAX
           << " !! Implementation have to be revised !" << endl;
   }
}

inline int xSpacePolynomialVariableP::getEntityDegree(const mEntity *e) const { return orderPredicate(e); }

inline int xSpacePolynomialVariableP::get_nb_shape_edge(int p) { return nbShapeFunctionEdge(p); }
inline int xSpacePolynomialVariableP::get_nb_shape_edge_tot(int p) { return nbShapeFunctionEdgeTotal(p); }
inline int xSpacePolynomialVariableP::get_nb_shape_tri(int p) { return nbShapeFunctionTri(p); }
inline int xSpacePolynomialVariableP::get_nb_shape_tri_tot(int p) { return nbShapeFunctionTriTotal(p); }
inline int xSpacePolynomialVariableP::get_nb_shape_tet(int p) { return nbShapeFunctionTet(p); }
inline int xSpacePolynomialVariableP::get_nb_shape_tet_tot(int p) { return nbShapeFunctionTetTotal(p); }

void xSpacePolynomialVariableP::getNodeKeys(mVertex *e, femKeys *keys, int order)
{
   if (o_node_geom_id.find(order) == o_node_geom_id.end())
   {
      xValKey::ids_size_t &node_geom_id = o_node_geom_id[order];
      create_geom_id_node(space_name, order, node_geom_id);
   }
   const int &node_geom_id = o_node_geom_id[order];
   keys->push_back(xValKey(Phys, node_geom_id, (mEntity *)e));
   return;
}

void xSpacePolynomialVariableP::getEdgeKeys(mEdge *e, femKeys *keys, int order, int swap)
{
   if (o_edge_geom_id.find(order) == o_edge_geom_id.end())
   {
      std::array<std::vector<xValKey::ids_size_t>, 2> &edge_geom_id = o_edge_geom_id[order];
      create_geom_id_edge(space_name, order, edge_geom_id);
   }
   std::array<std::vector<xValKey::ids_size_t>, 2> &edge_geom_id = o_edge_geom_id[order];
   for (int i = 0; i < get_nb_shape_edge(order); ++i)
   {
      keys->push_back(xValKey(Phys, edge_geom_id[swap][i], (mEntity *)(e)));
   }
   return;
}

void xSpacePolynomialVariableP::getTriKeys(mFace *e, femKeys *keys, int order, int swap)
{
   if (o_tri_geom_id.find(order) == o_tri_geom_id.end())
   {
      std::array<std::vector<xValKey::ids_size_t>, 6> &tri_geom_id = o_tri_geom_id[order];
      create_geom_id_tri(space_name, order, tri_geom_id);
   }
   std::array<std::vector<xValKey::ids_size_t>, 6> &tri_geom_id = o_tri_geom_id[order];
   for (int i = 0; i < get_nb_shape_tri(order); ++i)
   {
      keys->push_back(xValKey(Phys, tri_geom_id[swap][i], (mEntity *)(e)));
   }
   return;
}

void xSpacePolynomialVariableP::getTetKeys(mTet *e, femKeys *keys, int order)
{
   if (o_tet_geom_id.find(order) == o_tet_geom_id.end())
   {
      std::vector<xValKey::ids_size_t> &tet_geom_id = o_tet_geom_id[order];
      create_geom_id_tet(space_name, order, tet_geom_id);
   }
   std::vector<xValKey::ids_size_t> &tet_geom_id = o_tet_geom_id[order];

   for (int i = 0; i < get_nb_shape_tet(order); ++i)
   {
      keys->push_back(xValKey(Phys, tet_geom_id[i], (mEntity *)(e)));
   }
   return;
}

void xSpacePolynomialVariableP::getKeys(mEntity *e, femKeys *keys)
{
   switch (e->getLevel())
   {
      case 0:
      {
         const int order = getEntityDegree(e);
         keys->reserve(keys->size() + 1);
         getNodeKeys(static_cast<mVertex *>(e), keys, order);
         break;
      }
      case 1:
      {
         // keys->reserve(keys->size() + get_nb_shape_edge_tot(order));
         for (int node = 0; node < 2; ++node)
         {
            mVertex *v = static_cast<mVertex *>(e->get(0, node));
            const int order = getEntityDegree(v);
            getNodeKeys(v, keys, order);
         }
         const int order = getEntityDegree(e);
         getEdgeKeys(static_cast<mEdge *>(e), keys, order, 0);
         break;
      }
      case 2:
      {
         //    keys->reserve(keys->size()+ get_nb_shape_tri_tot(order) );
         if (e->getType() == mEntity::TRI)
         {
            mFace *f = static_cast<mFace *>(e);
            for (int node = 0; node < 3; ++node)
            {
               mVertex *v = static_cast<mVertex *>(f->get(0, node));
               const int order = getEntityDegree(v);
               getNodeKeys(v, keys, order);
            }
            for (int edge = 0; edge < 3; ++edge)
            {
               mEdge *ed = static_cast<mEdge *>(f->get(1, edge));
               const int order = getEntityDegree(ed);
               const int swap = computeSwapCaseTriEdge(*f, edge);
               getEdgeKeys(ed, keys, order, swap);
            }
            const int order = getEntityDegree(f);
            getTriKeys(f, keys, order, 0);
         }
         else
            throw;
         break;
      }
      case 3:
      {
         throw; /*
          if (e->getType() == mEntity::TET)
            {
              keys->reserve(keys->size()+ get_nb_shape_tet_tot(ordermax) );
              for (int node = 0; node < 4; ++node) {
                getNodeKeys(dynamic_cast < mVertex * > ( e->get(0, node)), keys);
              }
              for (int edge = 0; edge < 6; ++edge) {
                getEdgeKeys(dynamic_cast < mEdge * > ( e->get(1, edge)), keys);
              }
              for (int face = 0; face < 4; ++face) {
                getTriKeys(dynamic_cast < mFace * > ( e->get(2, face)), keys);
              }

              getTetKeys(dynamic_cast < mTet * > ( e ), keys);
            }
          else
          throw;*/
         break;
      }
      default:
      {
         throw;
      }
   }
}

void xSpacePolynomialVariableP::getKeysAndFcts(mEntity *e, femKeys *keys, femFcts *appro)
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

xSpacePolynomialBernsteinVariableP::xSpacePolynomialBernsteinVariableP(const std::string &a, TensorialType_t c, const int p,
                                                                       std::function<int(const mEntity *e)> &orderPredicate_)
    : xSpacePolynomialVariableP("BernsteinVariableP", a, c, p, orderPredicate_)
{
}

void xSpacePolynomialBernsteinVariableP::getFcts(const mEntity &e, femFcts &approscal)
{
   switch (e.getLevel())
   {
      case 0:
      {
         approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(1, 0)));
         // approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(1 , 0 )));
         break;
      }
      case 1:
      {
         for (int node = 0; node < 2; ++node)
         {
            const int pn = getEntityDegree(e.get(0, node));
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(pn, node)));
         }
         const int pe = getEntityDegree(&e);
         for (int i = 1; i < pe; ++i)
         {
            approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernstein(i, pe - i, 0)));
            //                 approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernstein( 1, ordermax-i,  0, 0)));//
            //                 on prend la premiere a chaque fois
         }
         break;
      }
      case 2:
      {
         if (e.getType() == mEntity::TRI)
         {
            for (int node = 0; node < 3; ++node)
            {
               const int pn = getEntityDegree(e.get(0, node));
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(pn, node)));
            }
            for (int edge = 0; edge < 3; ++edge)
            {
               const int pe = getEntityDegree(e.get(1, edge));
               for (int i = 1; i < pe; ++i)
               {
                  approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernstein(i, pe - i, edge)));
               }
            }
            const int pf = getEntityDegree(&e);
            for (int i = 1; i <= (pf - 2); ++i)
            {
               for (int j = 1; j <= (pf - i - 1); ++j)
               {
                  approscal.push_back(shapeFctPtr(new xApproxFunctionScalarTriBernstein(i, j, pf - i - j, 0)));
               }
            }
         }
         else
            throw;
         break;
      }
      case 3:
      {
         if (e.getType() == mEntity::TET)
         {
            cout << "xSpacePolynomialBernsteinVariableP not implmented in 3D...\n";
            throw;

            for (int node = 0; node < 4; ++node)
            {
               approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(1, node)));
               // approscal.push_back(shapeFctPtr(new xApproxFunctionScalarVertexBernstein(1 , node )));
            }
            for (int edge = 0; edge < 6; ++edge)
            {
               // const int swap = computeSwapCaseTetEdge( dynamic_cast < const mTet & > ( e ), edge);
               //	cout << swap << endl;
               for (int i = 1; i < ordermax; ++i)
               {
                  approscal.push_back(shapeFctPtr(new xApproxFunctionScalarEdgeBernstein(1, i, edge)));
               }
            }
            for (int face = 0; face < 4; ++face)
            {
               //    const int swap = computeSwapCaseTetFace( dynamic_cast < const mTet & > ( e ), face);
               for (int i = 1; i <= (ordermax - 2); ++i)
               {
                  for (int j = 1; j <= (ordermax - 1 - i); ++j)
                  {
                     approscal.push_back(shapeFctPtr(new xApproxFunctionScalarTriBernstein(i, j, ordermax - i - j, face)));
                  }
               }
            }
            for (int i = 1; i <= (ordermax - 3); ++i)
            {
               for (int j = 1; j <= (ordermax - 2 - i); ++j)
               {
                  for (int k = 1; k <= (ordermax - 1 - i - j); ++k)
                  {
                     approscal.push_back(shapeFctPtr(new xApproxFunctionScalarTetBernstein(i, j, k, ordermax - i - j - k)));
                  }
               }
            }
         }
         else
         {
            throw;
         }
         break;
      }
      default:
      {
         throw;
      }
   }
}

}  // end of namespace xfem
