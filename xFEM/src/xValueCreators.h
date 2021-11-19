#ifndef _VALUE_CREATORS_H
#define _VALUE_CREATORS_H
// boost
#include <boost/graph/adjacency_list.hpp>
// xtensor
#include "xVector.h"
// xfem
#include "xApproxFunction.h"
#include "xField.h"
#include "xLinkOnFrontLinkGenerator.h"
#include "xValue.h"
#include "xValueLinearCombination.h"
#include "xValueManager.h"

namespace xfem
{
// cree une valeur si la clef n'est pas deja dans ValManager
//
// cree une valeur seulement pour la clef concernee et pas
//  pour d'autres clefs i.e. pas d'effets de bord
// Il se peut que la valeur ne soit pas cree auquel cqs le retour est 0

template <class T>
class xValueCreator
{
  public:
   T* operator()(const xValKey& k) const
   {
      const bool debug = xdebug_flag;
      if (debug) std::cout << " inside value creator " << std::endl;
      return new T;
   }
};

template <class T>
class xNValueCreator
{
  public:
   xNValueCreator(size_t _n) : n(_n) {}
   xNValue<T>* operator()(const xValKey& k) const
   {
      const bool debug = xdebug_flag;
      if (debug) std::cout << " inside value creator " << std::endl;
      return new xNValue<T>(n);
   }

  private:
   size_t n;
};

template <class T>
class xValueCreatorForKey
{
  public:
   T* operator()(const xValKey& k) const { return new T(k); }
};

class xTensorsValueCreator
{
  public:
   xTensorsValueCreator() = default;
   virtual ~xTensorsValueCreator() = default;
   xTensorsValueCreator(xTensorsSignature& sig) : signature(&sig) {}
   virtual xValue<tensorsPtr_t>* operator()(const xValKey& key) const;
   void setSignature(xTensorsSignature* s) { signature = s; }

  protected:
   xTensorsSignature* signature;
};

template <typename VT>
class xMirrorCreator
{
   //    typedef xfem::xValueManagerDist<VT>  xValueManagerDist<double>;
  public:
   xMirrorCreator(xValueManagerDist<VT>* v, std::function<std::list<AOMD::mEntity*>(AOMD::mEntity*)> _getmirrors);
   xValue<VT>* operator()(const xValKey& key) const;

  private:
   xValueManagerDist<VT>* double_manager;
   std::function<std::list<AOMD::mEntity*>(AOMD::mEntity*)> getmirrors;
};

// mesh->lookupForMirrorEntities(key.getEnti(), mirrors)
template <typename VT>
xMirrorCreator<VT>::xMirrorCreator(xValueManagerDist<VT>* v, std::function<std::list<AOMD::mEntity*>(AOMD::mEntity*)> _getmirrors)
    : double_manager(v), getmirrors(_getmirrors)
{
}

// create a mirror if a value exist within the mirrors
template <typename VT>
xValue<VT>* xMirrorCreator<VT>::operator()(const xValKey& key) const
{
   const bool debug = xdebug_flag;
   if (debug) std::cout << " inside xValue<double>* xMirrorCreator::operator()(const xValKey& key) " << std::endl;
   std::list<AOMD::mEntity*> mirrors = getmirrors(key.getEnti());
   if (mirrors.empty()) return nullptr;

   xValKey keyother = key;
   if (debug) double_manager->PrintForDebug("in_create_mirror.dbg");
   if (debug) std::cout << "In xMirrorCreator key  is " << key << std::endl;
   std::list<AOMD::mEntity*>::const_iterator it;
   for (it = mirrors.begin(); it != mirrors.end(); ++it)
   {
      keyother.setEnti(*it);
      if (debug) std::cout << "In xMirrorCreator key of mirror is " << keyother << std::endl;
      // check if the other value exists
      if (xValue<double>* other = double_manager->find(keyother))
      {
         if (debug) std::cout << "In CreateMirror_c, mirror is created " << std::endl;
         if (debug) other->getValPtr();  // to check the sanity of other
         if (debug) std::cout << "after sanity check " << std::endl;
         return new xValueLinearCombination<VT>(1.0, other->getValPtr());
      }
   }
   return nullptr;
}

template <template <class> class DATAMANAGER, typename VT = double>
class xValueCreatorRegularAndHanging
{
   //    typedef xfem::xValueManagerDist<double>  xValueManagerDist<double>;
  public:
   xValueCreatorRegularAndHanging(xValueManagerDist<VT>* v, int degree, const DATAMANAGER<AOMD::mEntity*>& isHangingOn);
   xValue<VT>* operator()(const xValKey& key) const;
   xValue<VT>* degreeOne(const xValKey& key) const;
   xValue<VT>* degreeTwo(const xValKey& key) const;

  private:
   xValueManagerDist<VT>* double_manager;
   xValueCreator<xSingleValue<VT>> creator_reg;
   int degree;
   const DATAMANAGER<AOMD::mEntity*>& isHangingOn;
};

template <template <class> class DATAMANAGER, typename VT>
xValueCreatorRegularAndHanging<DATAMANAGER, VT>::xValueCreatorRegularAndHanging(xValueManagerDist<VT>* v, int degree_,
                                                                                const DATAMANAGER<AOMD::mEntity*>& isHangingOn_)
    : double_manager(v), degree(degree_), isHangingOn(isHangingOn_)
{
}

template <template <class> class DATAMANAGER, typename VT>
xValue<VT>* xValueCreatorRegularAndHanging<DATAMANAGER, VT>::degreeOne(const xValKey& key) const
{
   const bool debug = xdebug_flag;

   AOMD::mEntity* e = key.getEnti();
   if (debug)
   {
      std::cout << "in hanging entity is " << std::endl;
      e->print();
   }
   AOMD::mEntity* const* pph = isHangingOn.getData(*e);
   if (!pph)
   {
      xValue<VT>* val = creator_reg(key);
      if (debug)
      {
         std::cout << " regular value created " << std::endl;
      }
      return val;
   }
   const AOMD::mEntity& h = **pph;
   if (debug)
   {
      std::cout << "hanging is " << std::endl;
      h.print();
   }
   xValKey keyother = key;
   typename xValueLinearCombination<VT>::coeffs_t coeffs;
   typename xValueLinearCombination<VT>::xvalues_t values;
   for (int i = 0; i < h.size(0); ++i)
   {
      AOMD::mEntity* v = h.get(0, i);
      keyother.setEnti(v);
      xValue<VT>* other;
      if (!(other = double_manager->find(keyother))) return nullptr;
      coeffs.push_back(1. / h.size(0));
      values.push_back(other);
   }
   if (debug) std::cout << "creating a value linear combination" << std::endl;
   return new xValueLinearCombination<VT>(coeffs, values);
}

template <template <class> class DATAMANAGER, typename VT>
xValue<VT>* xValueCreatorRegularAndHanging<DATAMANAGER, VT>::degreeTwo(const xValKey& key) const
{
   const bool debug = xdebug_flag;

   AOMD::mEntity* e = key.getEnti();
   assert(e->size(3) == 0);  // on vérifie que l'on est en 2D, le 3D n'est pas codé.

   if (debug)
   {
      std::cout << " in hanging entity is " << std::endl;
      e->print();
   }
   AOMD::mEntity* const* pph = isHangingOn.getData(*e);
   if (!pph)
   {
      xValue<VT>* val = creator_reg(key);
      if (debug)
      {
         std::cout << " regular value created " << std::endl;
      }
      return val;
   }
   const AOMD::mEntity& h = **pph;
   if (debug)
   {
      std::cout << "hanging is " << std::endl;
      h.print();
   }
   xValKey keyother = key;
   typename xValueLinearCombination<VT>::coeffs_t coeffs;
   typename xValueLinearCombination<VT>::xvalues_t values;
   xValue<VT>* other;

   if (e->getLevel() == 0)
   {
      for (int i = 0; i < h.size(0); ++i)
      {
         AOMD::mEntity* v = h.get(0, i);
         keyother.setEnti(v);
         if (!(other = double_manager->find(keyother))) return nullptr;
         coeffs.push_back(1. / h.size(0));
         values.push_back(other);
      }
      keyother.setEnti(const_cast<AOMD::mEntity*>(&h));
      if (!(other = double_manager->find(keyother))) return nullptr;
      values.push_back(other);
      coeffs.push_back(HierarchicalApproxFunctionEdge(3, 0.));
      if (debug) std::cout << "creating a value linear combination" << std::endl;
      return new xValueLinearCombination<VT>(coeffs, values);
   }
   else
   {
      for (int i = 0; i < h.size(0); ++i)
      {
         AOMD::mEntity* v = h.get(0, i);
         keyother.setEnti(v);
         if (!(other = double_manager->find(keyother))) return nullptr;
         values.push_back(other);
      }
      keyother.setEnti(const_cast<AOMD::mEntity*>(&h));
      if (!(other = double_manager->find(keyother))) return nullptr;
      values.push_back(other);
      // getting the right coeffs
      if (debug)
      {
         std::cout << " small edge " << std::endl;
         e->print();
         std::cout << " big edge   " << std::endl;
         h.print();
      }
      xGeomElem geo_small(e);
      geo_small.setUVW({0., 0., 0.});
      xGeomElem geo_big(const_cast<AOMD::mEntity*>(&h));
      geo_big.setUVWForXYZ(geo_small.getXYZ());
      xtensor::xPoint uvw = geo_big.getUVW();
      double u = uvw(0);
      if (debug) std::cout << " In hanging degree 2 u is " << u << std::endl;
      if (u < 0.)
      {
         coeffs.push_back((HierarchicalApproxFunctionEdge(1, u) -
                           HierarchicalApproxFunctionEdge(1, 0.) * HierarchicalApproxFunctionEdge(2, 0.) -
                           HierarchicalApproxFunctionEdge(1, 0.)) /
                          (HierarchicalApproxFunctionEdge(3, 0.)));
         coeffs.push_back((HierarchicalApproxFunctionEdge(2, u) -
                           HierarchicalApproxFunctionEdge(2, 0.) * HierarchicalApproxFunctionEdge(2, 0.)) /
                          (HierarchicalApproxFunctionEdge(3, 0.)));
         coeffs.push_back((HierarchicalApproxFunctionEdge(3, u) -
                           HierarchicalApproxFunctionEdge(3, 0.) * HierarchicalApproxFunctionEdge(2, 0.)) /
                          (HierarchicalApproxFunctionEdge(3, 0.)));
      }
      else
      {
         coeffs.push_back((HierarchicalApproxFunctionEdge(1, u) -
                           HierarchicalApproxFunctionEdge(1, 0.) * HierarchicalApproxFunctionEdge(1, 0.)) /
                          (HierarchicalApproxFunctionEdge(3, 0.)));
         coeffs.push_back((HierarchicalApproxFunctionEdge(2, u) -
                           HierarchicalApproxFunctionEdge(2, 0.) * HierarchicalApproxFunctionEdge(1, 0.) -
                           HierarchicalApproxFunctionEdge(2, 0.)) /
                          (HierarchicalApproxFunctionEdge(3, 0.)));
         coeffs.push_back((HierarchicalApproxFunctionEdge(3, u) -
                           HierarchicalApproxFunctionEdge(3, 0.) * HierarchicalApproxFunctionEdge(1, 0.)) /
                          (HierarchicalApproxFunctionEdge(3, 0.)));
      }
      return new xValueLinearCombination<VT>(coeffs, values);
   }
}

template <template <class> class DATAMANAGER, typename VT>
xValue<VT>* xValueCreatorRegularAndHanging<DATAMANAGER, VT>::operator()(const xValKey& key) const
{
   // const bool debug = xdebug_flag;
   if (degree == 1)
      return degreeOne(key);
   else if (degree == 2)
      return degreeTwo(key);
   else
      assert(0);  // to be coded.
   return nullptr;
}

// the following class creates and links dofs around a front
// on the front support according to the paper entitled
// "A stable Lagrange multiplier space for stiff interface conditions within the extended finite element method", IJNME,
// Bechet, Moes and Wohlmuth DOI: 10.1002/nme.2515
// The algorithm for linking dofs of elements cut by the interface is given in xLinkOnFrontLinkGenerator
class xValueCreatorLinkOnFront
{
   //    typedef xfem::xValueManagerDist<double>  xValueManagerDist<double>;
   // typedef map<mEdge*,double> orthoinfo;// the higher is the double, the more orthogonal to the front is the edge... double
   // in [0,1]
  public:
   // front is given by iterators on edges (1D front) or faces (2D front)
   // this constructor is more generic since one can filter the front
   template <typename ITERFRONT>
   xValueCreatorLinkOnFront(ITERFRONT it, ITERFRONT end, xValueManagerDist<double>* v, bool link_isolated = true,
                            xEntityFilter filter = xAcceptAll(), xEntityToEntity interf2appro = xCreator());

   // front is given by the mesh
   xValueCreatorLinkOnFront(xValueManagerDist<double>* v, const xMesh* m, bool link_isolated = true,
                            xEntityFilter filt = xAcceptAll(), xEntityToEntity interf2appro = xCreator());
   xValue<double>* operator()(const xValKey& key) const;
   // gives iterators on vertices around the front
   vector<AOMD::mVertex*>::iterator beginVertexIter();
   vector<AOMD::mVertex*>::iterator endVertexIter();

  private:
   // link_generator indicates what to do during the dofs creation
   xValueManagerDist<double>* double_manager;
   xLinkOnFrontLinkGenerator link_generator;
};

//----------------------------------------------------------------------------------------------------------------

template <typename ITERFRONT>
xValueCreatorLinkOnFront::xValueCreatorLinkOnFront(ITERFRONT it, ITERFRONT end, xValueManagerDist<double>* v, bool link_isolated,
                                                   xEntityFilter filter, xEntityToEntity interf2appro)
    : double_manager(v), link_generator(it, end, link_isolated, filter, interf2appro)
{
}
// The following class creates Values Linear Combinations in vertex_begin - vertex_end region
// to the two dofs on the front, related to the closest edge
class xValueCreatorLinkAwayFromFront
{
   //    typedef xfem::xValueManagerDist<double>  xValueManagerDist<double>;
   typedef std::map<AOMD::mEntity*, std::pair<AOMD::mEntity*, xtensor::xPoint>> linkinfo;

  public:
   template <class ITER>
   xValueCreatorLinkAwayFromFront(xValueManagerDist<double>* v, const xMesh* m, ITER vertex_begin, ITER vertex_end);
   xValue<double>* operator()(const xValKey& key) const;

  private:
   void getIntersect(xtensor::xPoint first, xtensor::xPoint second, xtensor::xPoint origin, xtensor::xPoint& intersect,
                     double& distance, bool& valid, double& d1, double& d2);
   xValueManagerDist<double>* double_manager;
   linkinfo tab;
   const xMesh* frontmesh;
};

// The following class creates Values Linear Combinations in vertex_begin - vertex_end region
// to the closer dofs on the front, according to the tolerance parameter "radius_tol"
class xValueCreatorLinkAwayFromFrontInsideRadius
{
   //    typedef xfem::xValueManagerDist<double>  xValueManagerDist<double>;
   typedef std::map<AOMD::mEntity*, std::pair<xtensor::xPoint, double>> coeffsinfo;
   typedef std::map<AOMD::mEntity*, coeffsinfo> linkinfo;

  public:
   template <class ITER>
   xValueCreatorLinkAwayFromFrontInsideRadius(double radius_tol, xValueManagerDist<double>* v, const xMesh* m, ITER vertex_begin,
                                              ITER vertex_end);
   xValue<double>* operator()(const xValKey& key) const;

  private:
   xValueManagerDist<double>* double_manager;
   linkinfo tab;
   const xMesh* frontmesh;
};

// links discontinuous keys on both sides of a given interface
// typically used to properly declare a ridge enrichment on the interface
class xValueLinkDiscontinuousSpaceOnFrontBase
{
  public:
   xValueLinkDiscontinuousSpaceOnFrontBase(const xMesh* m, const xLevelSet& ls, const xField<double>& _field);
   ~xValueLinkDiscontinuousSpaceOnFrontBase();

  protected:
   //    typedef xDoubleManager  xValueManagerDist<double>;
   typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
   typedef std::multimap<int, long> groupe_type;

   xValue<double>* basicCreator(const xValKey& key) const;
   xValueManagerDist<double>* double_manager;
   const xMesh* frontmesh;
   const int highest_dim;

   vector<xValKey> keys;
   groupe_type sets_of_keys;
   mutable vector<bool> is_treated;
   mutable vector<long> master_per_group;
   vector<int> component;
   Graph G;
};

class xValueCreatorLinkDiscontinuousSpaceOnFront : public xValueLinkDiscontinuousSpaceOnFrontBase
{
  public:
   xValueCreatorLinkDiscontinuousSpaceOnFront(const xMesh* m, const xLevelSet& ls, const xField<double>& _field);
   ~xValueCreatorLinkDiscontinuousSpaceOnFront();
   xValue<double>* operator()(const xValKey& key) const;
};

class xValueUpdatorLinkDiscontinuousSpaceOnFront : public xValueLinkDiscontinuousSpaceOnFrontBase
{
  public:
   xValueUpdatorLinkDiscontinuousSpaceOnFront(const xMesh* m, const xLevelSet& ls, const xField<double>& _field,
                                              const std::string& sub, const xValKey::ids_size_t, const xValKey::ids_size_t,
                                              const xValKey::ids_size_t);
   ~xValueUpdatorLinkDiscontinuousSpaceOnFront();
   xValue<double>* operator()(const xValKey& key) const;
   void postUpdate();

  private:
   std::set<xValKey> old_disco_keys;
   mutable std::set<xValKey> new_disco_keys;
   const std::string& sub;
};

//#ifdef WITH_XLEGACYSIMPLECUT
//  // same as xValueCreatorLinkDiscontinuousSpaceOnFront, but for "ramp heaviside" enrichment:
//  // we have N colors (N groups of dofs around one entity) plus the classical dofs (classical shape function) of this entity
//  // -> we have one dof in surplus
//  // -> only one of these "sets" of nodes (with the same integer in the "component" vector) has to be "slave", i.e. linked to
//  one other "master" group
//  class xValueCreatorRampHeavisideBase {
//    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
//    typedef std::pair<AOMD::mEntity*,AOMD::mEntity*> entity_id;

//  public:
//     xValueCreatorRampHeavisideBase(xMesh* m, const xcut::xPhysSurfByTagging &front, const xField<double> &field);
//     ~xValueCreatorRampHeavisideBase();

//  protected:
////    typedef xfem::xValueManagerDist<double>  xValueManagerDist<double>;
//    xValue<double>* basicCreator(const xValKey& key);
//    void add_to_graph(const vector<entity_id>& temp);
//    xValueManagerDist<double>* double_manager;
//    xMesh* frontmesh;
//    const int highest_dim;

//    // variables used for the operator()
//    // the graph is composed of xValKeys, the groups are called components
//    vector<xValKey> keys;
//    vector<int> component;
//    map<int,xValKey> master_key_per_component;// regrouping information of masters and slaves...
//    map<int,double> coefficient_per_component;// either 1 or -1, regrouping information of masters and slaves...

//    // variables used only for the graph
//    // graph is composed of discontinuous entities, the groups are called colors
//    vector<entity_id> dentities_tobecolored;
//    vector<int> graph_color;
//    multimap<int,entity_id> color_to_dentities;
//    map<int,int> master_color_per_color;
//    map<AOMD::mEntity*,int> number_of_colors_per_entity;

//    Graph G;
//  };

//  class xValueCreatorRampHeaviside : public xValueCreatorRampHeavisideBase {
//  public:
//    xValueCreatorRampHeaviside( xMesh* m, const xcut::xPhysSurfByTagging &front, const xField<double> &_field);
//    ~xValueCreatorRampHeaviside();
//    xValue<double>* operator()(const xValKey& key);
//  };

//  class xValueUpdatorRampHeaviside  : public xValueCreatorRampHeavisideBase
//  {
//  public:
//    xValueUpdatorRampHeaviside( xMesh* m, const xcut::xPhysSurfByTagging &front, const xField<double> &_field,const
//    std::string& sub,const xValKey::ids_size_t,const xValKey::ids_size_t,const xValKey::ids_size_t );
//    ~xValueUpdatorRampHeaviside();
//    xValue<double>* operator()(const xValKey& key);
//    void postUpdate();
//  private:
//    set<xValKey> old_disco_keys;
//    set<xValKey> new_disco_keys;
//    const std::string& sub;

//  };

//  //does not work and I do not know why !!
//  //  class xFirstIfNotSecondValueCreator {

//  //  public:
//  //    xFirstIfNotSecondValueCreator(const std::function<xValue<double>* (const xValKey&)>& f,
//  //  			         const std::function<xValue<double>* (const xValKey&)>& s) :
//  //                                  cfirst(f), csecond(s)  {}
//  //    xValue<double>* operator()(const xValKey& key)
//  //      {
//  //        const bool debug = false;
//  //        if (debug) std::cout << " key is " << key << std::endl;

//  //        if (xValue<double>* ret =  cfirst(key))
//  //  	{
//  //  	  if (debug)
//  //  	    {
//  //  	      std::cout << "cfirst did create a value " << std::endl;
//  //  	      ret->print(std::cout);
//  //  	    }
//  //  	  return ret;
//  //  	}
//  //        xValue<double>* ret2 = csecond(key);
//  //        if (debug)
//  // 	 {
//  // 	   if (!ret2) std::cout << "csecond did create a value " << std::endl;

//  // 	 }
//  //        return ret2;
//  //      }

//  //    private:
//  //      const std::function<xValue<double>* (const xValKey&)>&            cfirst;
//  //      const std::function<xValue<double>* (const xValKey&)>&           csecond;
//  //    };
//#endif // WITH_XLEGACYSIMPLECUT
}  // namespace xfem

#endif
