/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XSPACEPOLYNOMIAL_
#define _XSPACEPOLYNOMIAL_
#include <map>
#include <vector>

#include "xApproxFunctionHighOrder.h"
#include "xSpace.h"
#include "xSpacePtr.h"

namespace AOMD
{
class mEdge;
class mFace;
class mTet;
}  // namespace AOMD

namespace xfem
{
class xGeomElem;
class xNonLocalInfoForKeysAndFcts;
class xNonLocalInfoGeneratorForKeysAndFcts;
class xSpacePolynomialOctreeLagrange;

int computeSwapCaseTriEdge(const AOMD::mFace& tri, int edgenumber);
int computeSwapCaseTetEdge(const AOMD::mTet& tet, int edgenumber);
int computeSwapCaseTetFace(const AOMD::mTet& tet, int facenumber);

void create_geom_id_node(const std::string& space_name, int order, xValKey::ids_size_t& node_geom_id);
void create_geom_id_edge(const std::string& space_name, int order, std::array<std::vector<xValKey::ids_size_t>, 2>& edge_geom_id);
void create_geom_id_tri(const std::string& space_name, int order, std::array<std::vector<xValKey::ids_size_t>, 6>& tri_geom_id);
void create_geom_id_tet(const std::string& space_name, int order, std::vector<xValKey::ids_size_t>& tet_geom_id);

/// Base class for high order shape function
/*! Only simplex element are considered here
    This is a pure virtual class.
    All derived class must implement getFcts, which is call by getKeysAndFunction
    and set the protected variable space_name at construction.
 */
class xSpacePolynomial : public xSpaceRegular
{
  public:
   xSpacePolynomial(const std::string& space_name, const std::string& physical_name, TensorialType_t c, const int p);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  protected:
   /// getFcts must push in approscal the scalar xApproxFunction that constitute the space, in the same order as the order of the
   /// keys, for all derived classes
   virtual void getFcts(const AOMD::mEntity& e, femFcts& approscal) = 0;
   void getNodeKeys(AOMD::mVertex* e, femKeys* keys);
   void getEdgeKeys(AOMD::mEdge* e, femKeys* keys, int swap);
   void getTriKeys(AOMD::mFace* e, femKeys* keys, int swap);
   void getTetKeys(AOMD::mTet* e, femKeys* keys);
   const int order;
   const int nb_shape_edge;
   const int nb_shape_edge_tot;
   const int nb_shape_tri;
   const int nb_shape_tri_tot;
   const int nb_shape_tet;
   const int nb_shape_tet_tot;
   const std::string space_name;
   xValKey::ids_size_t node_geom_id;
   std::array<std::vector<xValKey::ids_size_t>, 2> edge_geom_id;
   std::array<std::vector<xValKey::ids_size_t>, 6> tri_geom_id;
   std::vector<xValKey::ids_size_t> tet_geom_id;
};

class xSpacePolynomialLagrange : public xSpacePolynomial
{
  public:
   xSpacePolynomialLagrange(const std::string& a, TensorialType_t c, const int p);

  private:
   friend class xSpacePolynomialOctreeLagrange;
   void getFcts(const AOMD::mEntity& e, femFcts& approscal) override;
   std::array<femFcts, 4>& approstored;
   static std::map<int, std::array<femFcts, 4>> staticmap_order_approstored;
};

class xSpacePolynomialBernstein : public xSpacePolynomial
{
  public:
   xSpacePolynomialBernstein(const std::string& a, TensorialType_t c, const int p);

  private:
   void getFcts(const AOMD::mEntity& e, femFcts& approscal) override;
   std::array<femFcts, 4>& approstored;
   static std::map<int, std::array<femFcts, 4>> staticmap_order_approstored;
};

/// Base class for discontinuous high order shape function
/*! Only simplex element are considered here */
class xSpacePolynomialDiscontinuous : public xSpaceRegular
{
  public:
   xSpacePolynomialDiscontinuous(const std::string& a, TensorialType_t c, const int _meshDimension, const int p);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   static void getEntitiesFromKey(AOMD::mEntity** e, AOMD::mEntity** ref, int* n, xValKey& key);

  protected:
   void getEntityKeys(const int level, const int AOMD_id, AOMD::mEntity* ref, femKeys* keys);
   const int order, meshDimension;
};

class xSpacePolynomialDiscontinuousLagrange : public xSpacePolynomialDiscontinuous
{
  public:
   xSpacePolynomialDiscontinuousLagrange(const std::string& a, TensorialType_t c, const int _meshDimension, const int p);
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
};

/// Base for space Octree
/*! this space use non local informations for keys and function generations */
class xSpacePolynomialOctree : public xSpacePolynomial
{
  public:
   xSpacePolynomialOctree(const std::string& space_name, const std::string& physical_name, TensorialType_t c, const int p);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  protected:
   int nb_free;
   const xNonLocalInfoForKeysAndFcts* cur_nli;
   int getShapeIdFromGeom(xValKey::ids_size_t geomid, const std::vector<xValKey::ids_size_t>& ids);
   femKeys other_master_keys;
   AOMD::mAdjacencyContainer::iter it_begin_node;
   AOMD::mAdjacencyContainer::iter it_end_node;
   AOMD::mAdjacencyContainer::iter it_begin_edge;
   AOMD::mAdjacencyContainer::iter it_end_edge;
   AOMD::mAdjacencyContainer::iter it_begin_face;
   AOMD::mAdjacencyContainer::iter it_end_face;
   virtual void getFcts(AOMD::mEntity& e, const femKeys& keys, femFcts& approscal) = 0;
   void getFcts(const AOMD::mEntity& e, femFcts& approscal) override {}
   virtual const xNonLocalInfoForKeysAndFcts* getNonLocalInfo(const AOMD::mEntity& e) const = 0;
};

/// Derive class of octree space for Lagrange function
/*! Shape function is the only part wich differs from base classe */
/*! Berstein derived class have to be done */
class xSpacePolynomialOctreeLagrange : public xSpacePolynomialOctree
{
  public:
   xSpacePolynomialOctreeLagrange(const std::string& a, TensorialType_t c, const int p,
                                  xNonLocalInfoGeneratorForKeysAndFcts* extended_generator);

  private:
   void getFcts(AOMD::mEntity& e, const femKeys& keys, femFcts& approscal) override;
   xApproxFunction* getFctFromKey(const AOMD::mEntity& e, const xValKey& key, const double coef);
   const xNonLocalInfoForKeysAndFcts* getNonLocalInfo(const AOMD::mEntity& e) const override;
   const xNonLocalInfoGeneratorForKeysAndFcts* extended_generator;
   spacePtr_t std_space;
};

/// Base class for h-p Hierarchicalal high order shape functions. Non homogeneous polynomial order is allowed on the mesh.
/*! Only simplex element are considered here
    This is a pure virtual class.
    All derived class must implement getFcts, which is call by getKeysAndFunction
    and set the protected variable space_name at construction.
 */

class xSpacePolynomialVariableP : public xSpaceRegular
{
  public:
   xSpacePolynomialVariableP(const std::string& space_name, const std::string& physical_name, TensorialType_t c, int pmax,
                             std::function<int(const AOMD::mEntity* e)>& orderPredicate_);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  protected:
   /// getFcts must push in approscal the scalar xApproxFunction that constitute the space, in the same order as the order of the
   /// keys, for all derived classes
   virtual void getFcts(const AOMD::mEntity& e, femFcts& approscal) = 0;
   void getNodeKeys(AOMD::mVertex* e, femKeys* keys, int order);
   void getEdgeKeys(AOMD::mEdge* e, femKeys* keys, int order, int swap);
   void getTriKeys(AOMD::mFace* e, femKeys* keys, int order, int swap);
   void getTetKeys(AOMD::mTet* e, femKeys* keys, int order);
   const int ordermax;

   const std::string space_name;
   std::function<int(const AOMD::mEntity* e)>& orderPredicate;

   std::map<int, xValKey::ids_size_t> o_node_geom_id;
   std::map<int, std::array<std::vector<xValKey::ids_size_t>, 2>> o_edge_geom_id;
   std::map<int, std::array<std::vector<xValKey::ids_size_t>, 6>> o_tri_geom_id;
   std::map<int, std::vector<xValKey::ids_size_t>> o_tet_geom_id;

   inline int getEntityDegree(const AOMD::mEntity* e) const;
   inline int get_nb_shape_edge(int order);
   inline int get_nb_shape_edge_tot(int order);
   inline int get_nb_shape_tri(int order);
   inline int get_nb_shape_tri_tot(int order);
   inline int get_nb_shape_tet(int order);
   inline int get_nb_shape_tet_tot(int order);
};

/// Derived class for p  Bernstein shape functions
class xSpacePolynomialBernsteinVariableP : public xSpacePolynomialVariableP
{
  public:
   xSpacePolynomialBernsteinVariableP(const std::string& a, TensorialType_t c, const int pmax,
                                      std::function<int(const AOMD::mEntity* e)>& orderPredicate_);

  private:
   void getFcts(const AOMD::mEntity& e, femFcts& approscal) override;
};

}  // namespace xfem
#endif
