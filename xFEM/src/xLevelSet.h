/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _LSET_H
#define _LSET_H

#include <iosfwd>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
// Trellis
#include "AOMDfwd.h"
#include "mEntity.h"
#include "mEntityContainer.h"
#include "mIterator.h"
#include "mVertex.h"
// xfem
#include "xGeomElem.h"
#include "xMesh.h"
#include "xPointToDouble.h"
#include "xRegion.h"
#include "xTensor2.h"
#include "xVector.h"

// eXlibris_tools
#include "xDataExchangerTraits.h"

namespace xtool
{
template <typename T>
class xKeyContainerSendAndRecv;
class xMpiInputBuffer;
class xMpiOutputBuffer;

}  // namespace xtool
namespace xfem
{
template <typename T>
class xEval;
template <class T>
/// class to store values at node. Template on the value to store
class xVertexValuesStorage
{
  public:
   xVertexValuesStorage(const xRegion& s, const T& val = T())
   {
      for (xIter it = s.begin(0); it != s.end(0); ++it)
      {
         AOMD::mVertex* v = (AOMD::mVertex*)*it;
         values[v] = val;
      }
   }
   /// get Functions
   const T& operator()(AOMD::mVertex* v) const
   {
      assert(values.find(v) != values.end());
      return values.find(v)->second;
   }
   /// set vals at vertex
   T& operator()(AOMD::mVertex* v) { return values[v]; }

  private:
   typedef std::unordered_map<AOMD::mEntity*, T, AOMD::EntityHashKey, AOMD::EntityEqualKey> storagetype;
   storagetype values;
};

class xLevelSetModifier;
class xLevelSetInspector;
class xLevelSetCreator;

template <class UnaryOperator>
class xEvalLevelSet;

class xLevelSet
{
  public:
   virtual ~xLevelSet();
   xLevelSet();
   xLevelSet(const xLevelSet& in);
   xLevelSet(const xRegion& s, const double& val = 0.0);  // attention; implicit cast avec xMesh * le cas echeant
   xLevelSet(const xRegion& s, const xPointToDouble& func);
   void load(const xPointToDouble& func);
   void loadOnElement(const xPointToDouble& func, AOMD::mEntity* e);
   xLevelSet& operator=(const xLevelSet& other);
   /// Transform level set to ls = -ls
   void complement();
   /// set Function ... obsolete ??
   //! This function could add a new value in the internal map and break the support it should be removed.
   double& operator()(AOMD::mEntity* v);  /// 03-07-09 :  ajout virtual
   /// This function return a reference to the ls value associated to v if it already exist. otherwise throw;
   double& at(const AOMD::mVertex& v);
   /// This function return a const reference to the ls value associated to v if it already exist. otherwise throw;
   const double& at(const AOMD::mVertex& v) const;

   /// get value associated to vertex. throw if nno value is associated to v.
   double operator()(const AOMD::mEntity* v) const;
   double operator()(const AOMD::mVertex& v) const;
   // accept
   void accept(xLevelSetModifier& visitor);
   void accept(xLevelSetModifier& visitor, const xRegion& target);
   void accept(xLevelSetInspector& visitor) const;
   void accept(xLevelSetInspector& visitor, const xRegion& target) const;
   void accept(xLevelSetCreator& visitor, xLevelSet& f) const;
   void accept(xLevelSetCreator& visitor, xLevelSet& f, const xRegion& target) const;

   /// get Vals at the nodes of entity e (e may be a node, a edge , ...)
   virtual std::vector<double> getVals(const AOMD::mEntity* e) const;
   virtual std::vector<double> getVals(const AOMD::mEntity& e) const;

   /// get Vals at the nodes given by argument
   virtual std::vector<double> getVals(std::vector<AOMD::mEntity*>& v) const;
   /// get grad at a node or on an element (On a node, its the mean of the neighboring element)
   /*! (On an element its the (constant) gradiant of the level set)*/
   xtensor::xVector<> getGrad(AOMD::mEntity* e) const;
   /// get val on an element interpolated at uvw from the nodes of e. The AOMD::mEntity is not know by the level set, an exception
   /// is thrown
   double getVal(AOMD::mEntity* e, const xtensor::xPoint& uvw) const;
   /// get val on an element interpolated at uvw from the nodes of e. return 0 if the element is not known by the levelset ...
   /*! return 1 if the element is known by the levelSet, and val is updated to the value of the levelset at this point.!*/
   bool getVal(AOMD::mEntity* e, const xtensor::xPoint& uvw, double& val) const;
   /// get grad on an element interpolated at uvw from the nodes of e. The AOMD::mEntity is not know by the level set, an
   /// exception is thrown
   xtensor::xVector<> getGrad(AOMD::mEntity* e, const xtensor::xPoint& uvw) const;
   /// get grad on an element interpolated at uvw from the nodes of e. return 0 if the element is not known by the levelset ...
   /*! return 1 if the element is known by the levelSet, and vals is updated to the gradient value of the levelset at this
    * point.!*/
   bool getGrad(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& vals) const;
   /// second gradient of the level set
   xtensor::xTensor2<> getCurv(AOMD::mEntity*) const;
   /// gradient of the normalized gradient of the level set
   /*!  the trace of this quantity is what is called the curvature */
   xtensor::xTensor2<> getTrueCurv(AOMD::mEntity*) const;
   /// get mesh
   const xRegion& getSupport() const;
   void setSupport(const xRegion& m, const double& val = 0.0);
   void reduceSupport(const xRegion& m);
   void clear();
   void printDebug() const;
   void exportMatlab(std::ostream& fout, const std::string& field_name, int level = 0) const;
   /// Add to the xLevelSet trace, the values of the xLevelSet object at the nodes of the xMesh mesh
   /*! For this function to work, nodes of the input mesh must know there "parents " in the support mesh of
     the xLevetSet object, via a "was_created_by_tag" attached entitie
     !*/

   void takeTraceOn(const xMesh& target_mesh, const xMesh::datamanager_t<AOMD::mEntity*>& was_created_by,
                    const xMesh::datamanager_t<double>& r_on_edge, xLevelSet& trace) const;
   // void takeTraceOn(const xMesh& mesh, xLevelSet& trace) const;

   void interpolateTo(const xMesh& mesh_new, const xMesh::datamanager_t<AOMD::mEntity*>& was_duplicated_from,
                      const xMesh::datamanager_t<AOMD::mEntity*>& was_created_by, xLevelSet& lsnew) const;
   /// change the support of the levelset. new_support is the new support.
   /*!
     Each vertex that where not known by the old support are set to value init_val.
     All vertices that were on the old support but not on the new are removed from the ls map.
     All vertices that were on both support are maintain in the datasttucture and keep there old values.
     !*/
   void restrictTo(xRegion new_support, double init_val = 0.0);
   bool isDefinedAt(AOMD::mEntity*) const;
   bool isDefinedOnElement(AOMD::mEntity* e) const;
   int side_of(const xGeomElem* geo_appro, const xGeomElem* geo_integ) const;
   void forceGradComputation();

  private:
   unsigned int was_created_by_tag;
   unsigned int r_on_edge_tag;
   unsigned int was_duplicated_from_tag;
   mutable bool key_up_to_date;
   mutable bool grad_up_to_date;
   mutable bool curv_up_to_date;
   mutable bool true_curv_up_to_date;
   mutable xRegion support;
   template <typename T>
   using containers_t = std::unordered_map<AOMD::mEntity*, T, AOMD::EntityHashKey>;
   typedef containers_t<double> lstype;
   lstype ls;
   mutable containers_t<xtensor::xVector<>> grad;
   mutable containers_t<xtensor::xTensor2<>> curv;
   mutable containers_t<xtensor::xTensor2<>> true_curv;
   // contains the gradient at the vertex (obtained by smoothing)
   // and the gradient inside each element
   void create_key() const;
   void compute_grad() const;
   void compute_curv() const;
   void compute_true_curv() const;
   class xLevelSetKeyManager
   {
     public:
      typedef const AOMD::mEntity* information_key_t;
      // Types related to data that are used to get information key
      // With this class we will iterate on lsr entries to set keys. Data type is then the pair mEntity*/double
      typedef std::iterator_traits<lstype::iterator>::value_type data_t;

      xLevelSetKeyManager(const partmanAOMD_t& _partman);

      information_key_t localObjectKey(const data_t& o) const;
      information_key_t remoteObjectKey(const xtool::xRemoteObject<AOMD::mEntity>& ro, const data_t& lo) const;
      xtool::xConstPartitionObject<AOMD::mEntity> getConstPartitionObject(const data_t& o) const;

     private:
      const partmanAOMD_t& partman;
   };
   mutable xtool::xKeyContainerSendAndRecv<xLevelSetKeyManager::information_key_t>* keys_container;
   // info manager for gradiant computation
   class xLevelSetGradInfoManager
   {
     public:
      typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
      typedef xtool::send_and_recv_keys_communication_trait communication_trait;
      typedef xLevelSetKeyManager::information_key_t information_key_t;

      xLevelSetGradInfoManager(containers_t<xtensor::xVector<>>& grad_, containers_t<int>& nbconnect_);
      void getInfo(information_key_t key, xtool::xMpiInputBuffer& buff, int sendto);
      void setInfo(information_key_t key, const xtool::xMpiOutputBuffer& buff, int receivedfrom);
      size_t getApproxDataSize();

     private:
      const size_t s;
      containers_t<xtensor::xVector<>>& grad;
      containers_t<int>& nbconnect;
   };
};

template <class UnaryOperator>
class xEvalLevelSet : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalLevelSet(const xLevelSet& ls_) : ls(ls_) {}
   xEvalLevelSet(const xLevelSet& ls_, const UnaryOperator& _funct) : ls(ls_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v;
      v = ls.getVal(appro->getEntity(), appro->getUVW());
      result = funct(v);
   }

  private:
   const xLevelSet& ls;
   UnaryOperator funct;
};

template <class UnaryOperator>
class xEvalGradLevelSet : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalGradLevelSet(const xLevelSet& ls_) : ls(ls_) {}
   xEvalGradLevelSet(const xLevelSet& ls_, const UnaryOperator& _funct) : ls(ls_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v;
      v = ls.getGrad(appro->getEntity());
      result = funct(v);
   }

  private:
   const xLevelSet& ls;
   UnaryOperator funct;
};

template <class UnaryOperator>
class xEvalGradSmoothLevelSet : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalGradSmoothLevelSet(const xLevelSet& ls_) : ls(ls_) {}
   xEvalGradSmoothLevelSet(const xLevelSet& ls_, const UnaryOperator& _funct) : ls(ls_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v;
      v = ls.getGrad(appro->getEntity(), appro->getUVW());
      result = funct(v);
   }

  private:
   const xLevelSet& ls;
   UnaryOperator funct;
};

// second gradient of the level set function
template <class UnaryOperator>
class xEvalHessianLevelSet : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalHessianLevelSet(const xLevelSet& ls_) : ls(ls_) {}
   xEvalHessianLevelSet(const xLevelSet& ls_, const UnaryOperator& _funct) : ls(ls_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const
   {
      typename UnaryOperator::argument_type v;
      v = ls.getCurv(appro->getEntity());
      result = funct(v);
   }

  private:
   const xLevelSet& ls;
   UnaryOperator funct;
};

// gradient of the normalized gradient of the level set.
template <class UnaryOperator>
class xEvalCurvatureLevelSet : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalCurvatureLevelSet(const xLevelSet& ls_) : ls(ls_) {}
   xEvalCurvatureLevelSet(const xLevelSet& ls_, const UnaryOperator& _funct) : ls(ls_), funct(_funct) {}
   // non optimized version, an optimized version exists when the UnaryOperator is Identity
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const
   {
      typename UnaryOperator::argument_type v;
      v = ls.getTrueCurv(appro->getEntity());
      result = funct(v);
   }

  private:
   const xLevelSet& ls;
   UnaryOperator funct;
};

}  // namespace xfem

#endif
