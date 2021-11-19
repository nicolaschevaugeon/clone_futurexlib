/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _VALUE_H
#define _VALUE_H

#include <iosfwd>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "xTensorsPtr.h"
#include "xValKey.h"
#include "xValManager.h"
#include "xDebug.h"
#include "xVector.h"
#include "xNearestNeighborInterface.h"
#include "xEntityFilter.h"




namespace xfem
{
  class xLevelSet;
class xField;
class xcut::xPhysSurfByTagging;
class xMesh;
class xTensors;
class xTensorsSignature;

using namespace std;


//! State of Value base class
//! 
class xStateOfValue {
  public:
  enum state_type { DOF, NONE, FIXED, LINEARCOMBINATION, UNDEFINED};
 xStateOfValue():state(UNDEFINED){};
 xStateOfValue(const state_type &_state ):state(_state){};
  virtual ~xStateOfValue(void) {}
  virtual std::ostream& print(std::ostream& o) const = 0;
  state_type state;
};

//! Value base class
//! 
template <class T>
class xValue 
{
public:
  xValue(void);
  virtual ~xValue();
  virtual T getVal() const = 0;
  virtual void setVal(T in) = 0;
  virtual std::ostream& print(std::ostream& o) const;
  virtual std::ostream& printVal(std::ostream& o) const = 0;
  virtual std::ostream& printState(std::ostream& o) const;
  virtual const xStateOfValue* getState() const;
  virtual       xStateOfValue* getState();
  virtual void setState(xStateOfValue* s);
  virtual void delState();
  virtual xValue<T>* getValPtr();
protected:
  mutable xStateOfValue* state;
private:
};


template <class T>
xValue<T>::xValue () : state(0)  {}         
template <class T>
xValue<T>::~xValue()   {delete state;}
template <class T>
std::ostream& xValue<T>::print(std::ostream& o) const 
{ printVal(o); o << std::endl; printState(o); return o; }
template <class T>
std::ostream& xValue<T>::printState(std::ostream& o) const 
{ if (state) return state->print(o);  return o;  }
template <class T>
const xStateOfValue* xValue<T>::getState() const  {return state;}
template <class T>
xStateOfValue* xValue<T>::getState() {return state;}
template <class T>
void xValue<T>::setState(xStateOfValue* s)  {
  if (state) delState();
  state = s;
}
template <class T>
void xValue<T>::delState() {if (state) {delete state; state = 0;} }
template <class T>
xValue<T>* xValue<T>::getValPtr()             {return this;}


class xValueDouble : public xValue<double> {
 public:
  xValueDouble();
  virtual double   getVal() const           {return value;}
  virtual   void   setVal(double in) {value = in;}
  virtual std::ostream& printVal(std::ostream & o) const {
    o << value;
    return o;
  }
private:
  double  value;
};


template <class T>
class xNValue : public xValue< T >{
 public:
 xNValue(size_t N):value(N, 0.){}
  static void set_current(size_t i){
    current = i;
  }
  virtual T   getVal() const           {return value[current];}
  virtual   void   setVal(double in) {value[current] = in;}
  virtual std::ostream& printVal(std::ostream & o) const {
    o << value[current];
    return o;
  };
 private:
  std::vector <T>  value;
  static size_t current;
};
 
  


// a value is tied to other values and a constant
// v = sum_i c_i v_i + b
// c_i are the coefficients, v_i are the other value
// b is a constant which by default is zero.

// in the implementation, b is an extra coefficient c_last 
// multiplying a fixed
// value v_last of value 1.0

class xStateOfValueLinearCombination;

class xValueLinearCombination : public xValue<double> {
  friend class xStateOfValueLinearCombination;

public:
  typedef std::vector<double> coeffs_t;
  typedef std::vector<xValue<double>*> values_t;

  xValueLinearCombination(const coeffs_t& c, const values_t& v, double coeff_last = 0.0);
  xValueLinearCombination(double c, xValue<double>* v, double coeff_last = 0.0);
  virtual double getVal() const;           
  virtual   void   setVal(double in);
  virtual const xStateOfValue* getState()  const;
  virtual       xStateOfValue* getState();
  void setState(xStateOfValue* s);
  virtual std::ostream& printVal(std::ostream& o) const;

private:
  typedef coeffs_t::const_iterator iterator_coeff;
  coeffs_t coeffs;
  typedef values_t::const_iterator iterator_value;
  values_t values;
  xValueDouble * v_last;
};

 
class xValueError : public xValue<double> {

public:
  xValueError();
  virtual double   getVal() const;
  virtual   void   setVal(double in);
  virtual std::ostream& printVal(std::ostream & o) const { return o; }
  static void choice(const std::string& s);
  static double total(const std::string& s);
  static void   clear();

private:
  typedef enum {ENG,     ABS2,     ABS,     REL, 
                ENG_EXA, ABS2_EXA, ABS_EXA, REL_EXA, EFF_EXA,
                EFF, DNS, VOL} choice_t;
  static choice_t _choice;
  double abs2_loc;
  double abs2_exa_loc;
  double eng_loc;
  double eng_exa_loc;
  double vol_loc;
  static double abs2_tot;
  static double abs2_exa_tot;
  static double eng_tot;
  static double eng_exa_tot;
  static double vol_tot;

};


//cr� une valeur si la clef n'est pas d��dans ValManager
//
//cr� une valeyr seulement pour la clef concern� et pas
//  pour d'autres clefs i.e. pas d'effets de bord
//Il se peut que la valeur ne soit pas cr�e auquel cqs le retour est 0


template <class T>
class xValueCreator {
   public:
   T* operator()(const xValKey& k) 
  { 
    const bool debug = xdebug_flag;
    if (debug) std::cout << " inside value creator " << std::endl;
    return new T; 
  }
};


template <class T>
class xNValueCreator {
 public:
 xNValueCreator( size_t _n):n(_n){};
  
  xNValue< T > * operator()(const xValKey& k) { 
    const bool debug = xdebug_flag;
    if (debug) std::cout << " inside value creator " << std::endl;
    return new xNValue<T>(n); 
  }
 private :
  size_t n;
 };
 

template <class T>
class xValueCreatorForKey {
    public:
    T* operator()(const xValKey& k) { return new T(k); }
};
 

class xMirrorCreator {
  typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
public:
  xMirrorCreator(xValueManagerDist<double>* v, xMesh* m);
  xValue<double>* operator()(const xValKey& key);
private:
  xValueManagerDist<double>* double_manager;
  xMesh* mesh;
};




class xValueCreatorRegularAndHanging {
  typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
public:
  xValueCreatorRegularAndHanging(xValueManagerDist<double>* v, xMesh* m, int degree);
  xValue<double>* operator()(const xValKey& key);
  xValue<double>* degreeOne(const xValKey& key);
  xValue<double>* degreeTwo(const xValKey& key);
private:
  xValueManagerDist<double>* double_manager;
  xMesh* mesh; 
  xValueCreator<xValueDouble> creator_reg;
  int degree;
};

// the following class creates and links dofs around a front
// on the front support according to the paper entitled
// "A stable Lagrange multiplier space for stiff interface conditions within the extended finite element method", IJNME, Bechet, Moes and Wohlmuth DOI: 10.1002/nme.2515

class xValueCreatorLinkOnFront {
	typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
	//typedef map<mEdge*,double> orthoinfo;// the higher is the double, the more orthogonal to the front is the edge... double in [0,1]

	// the following map indicates what to do during the dofs creation 
	// if the second AOMD::mVertex* is not NULL, the two vertices must be linked (linear combination state for the first vertex)
	// if not, the first vertex has to be declared as new dof ("double" state)
	typedef map <AOMD::mVertex*, AOMD::mVertex* > creationinfo;



	public:

        // front is given by iterators on edges (1D front) or faces (2D front)
        // this constructor is more generic since one can filter the front
        template <typename ITERFRONT>
        xValueCreatorLinkOnFront(ITERFRONT it, ITERFRONT end, xValueManagerDist<double>* v, bool link_isolated=true, xEntityFilter filter=xAcceptAll());

        // front is given by the mesh
	xValueCreatorLinkOnFront(xValueManagerDist<double>* v, const xMesh* m, bool link_isolated=true, xEntityFilter filt = xAcceptAll() );
	~xValueCreatorLinkOnFront();
	xValue<double>* operator()(const xValKey& key);
	// gives iterators on vertices around the front
	// should use iterators directly on the map, instead of creating a copy in a vector... ?
	vector<AOMD::mVertex*>::iterator beginVertexIter();
	vector<AOMD::mVertex*>::iterator endVertexIter();

	private:
	void buildVertexVector();
        template <typename ITERFRONT>
        void buildNormalAtNodes(ITERFRONT it, ITERFRONT end);
        void buildTab(bool link_isolated, xEntityFilter filter);
	//double compute_orthogonality(mEdge *e, xtensor::xVector vector1, AOMD::mVertex *v);
	//orthoinfo orthotab;
	xValueManagerDist<double>* double_manager;
	vector<AOMD::mVertex*> vertexvector;
	creationinfo tab;
    typedef std::map < AOMD::mEntity *,xtensor::xVector> normal_at_nodes_t;
    typedef normal_at_nodes_t::iterator normal_at_nodes_it_t;
    normal_at_nodes_t normal_at_nodes;
    
    // private class to sort AOMD::mEntity viewed as AOMD::mVertex with criteria on edge.normale value and coordinate acending order
    class nodeSortingCriteria 
    {
        public:
            nodeSortingCriteria(normal_at_nodes_t * norms_);
            bool operator()(AOMD::mEntity* ent1, AOMD::mEntity* ent2) const;
        private:
            normal_at_nodes_t * norms;
    };

    typedef set<AOMD::mEntity *,nodeSortingCriteria> node_set_t;

    // normal computation
    // nota :
    //    normForElem2D is dedicated to edge
    //    normForElem3D is dedicated to triangle
    //
    std::function<void (AOMD::mEntity *,xtensor::xVector&)> normForElem;
    void normForElem2D(AOMD::mEntity *elem, xtensor::xVector& norm);
    void normForElem3D(AOMD::mEntity *elem, xtensor::xVector& norm);
};

// The following class creates Values Linear Combinations in vertex_begin - vertex_end region
// to the two dofs on the front, related to the closest edge 
class xValueCreatorLinkAwayFromFront {

	typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
	typedef map<AOMD::mEntity*, pair<AOMD::mEntity*, Trellis_Util::mPoint> > linkinfo;

	public:
	template <class ITER>
		xValueCreatorLinkAwayFromFront(xValueManagerDist<double>* v, const xMesh* m, ITER vertex_begin, ITER vertex_end);
	xValue<double>* operator()(const xValKey& key);

	private:
	void getIntersect(Trellis_Util::mPoint first,Trellis_Util::mPoint second,Trellis_Util::mPoint origin,Trellis_Util::mPoint &intersect,double &distance, bool &valid, double &d1, double &d2);
	xValueManagerDist<double>* double_manager;
	linkinfo tab;
	const xMesh* frontmesh; 
};

// The following class creates Values Linear Combinations in vertex_begin - vertex_end region
// to the closer dofs on the front, according to the tolerance parameter "radius_tol"
class xValueCreatorLinkAwayFromFrontInsideRadius {

	typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
	typedef map<AOMD::mEntity*, pair<Trellis_Util::mPoint,double> > coeffsinfo;
	typedef map<AOMD::mEntity*,  coeffsinfo> linkinfo;

	public:
	template <class ITER>
		xValueCreatorLinkAwayFromFrontInsideRadius(double radius_tol, xValueManagerDist<double>* v, const xMesh* m, ITER vertex_begin, ITER vertex_end);
	xValue<double>* operator()(const xValKey& key);

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
    xValueLinkDiscontinuousSpaceOnFrontBase(const xMesh* m, const xLevelSet &ls, const xField &_field);
    ~xValueLinkDiscontinuousSpaceOnFrontBase();

  protected:
    typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
    typedef multimap <int,long > groupe_type;

    xValue<double>* basicCreator(const xValKey& key);
    xValueManagerDist<double>* double_manager;
    const xMesh* frontmesh;
    const int highest_dim;

    vector<xValKey> keys;
    groupe_type sets_of_keys;
    vector<bool> is_treated;
    vector<long> master_per_group;
    vector<int> component;
    Graph G;

};

class xValueCreatorLinkDiscontinuousSpaceOnFront : public xValueLinkDiscontinuousSpaceOnFrontBase 
{

  public:
    xValueCreatorLinkDiscontinuousSpaceOnFront(const xMesh* m, const xLevelSet &ls, const xField &_field);
    ~xValueCreatorLinkDiscontinuousSpaceOnFront();
    xValue<double>* operator()(const xValKey& key);

};

class xValueUpdatorLinkDiscontinuousSpaceOnFront  : public xValueLinkDiscontinuousSpaceOnFrontBase 
{
  public:
    xValueUpdatorLinkDiscontinuousSpaceOnFront(const xMesh* m, const xLevelSet &ls, const xField &_field,const std::string& sub,const int,const int,const int );
    ~xValueUpdatorLinkDiscontinuousSpaceOnFront();
    xValue<double>* operator()(const xValKey& key);
    void postUpdate(void);
  private:
    set<xValKey> old_disco_keys;
    set<xValKey> new_disco_keys;
    const std::string& sub;

};


// same as xValueCreatorLinkDiscontinuousSpaceOnFront, but for "ramp heaviside" enrichment:
// we have N colors (N groups of dofs around one entity) plus the classical dofs (classical shape function) of this entity
// -> we have one dof in surplus
// -> only one of these "sets" of nodes (with the same integer in the "component" vector) has to be "slave", i.e. linked to one other "master" group
class xValueCreatorRampHeavisideBase {
  typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
  typedef std::pair<AOMD::mEntity*,AOMD::mEntity*> entity_id;

  public:
  xValueCreatorRampHeavisideBase(xMesh* m, const xcut::xPhysSurfByTagging &front, const xField &field);
  ~xValueCreatorRampHeavisideBase();

  protected:
  typedef xValManager<xValKey, xValue<double>, xHashValKey, xEqualValKey> xValueManagerDist<double>;
  xValue<double>* basicCreator(const xValKey& key);
  void add_to_graph(const vector<entity_id>& temp);
  xValueManagerDist<double>* double_manager;
  xMesh* frontmesh;
  const int highest_dim;

  // variables used for the operator()
  // the graph is composed of xValKeys, the groups are called components
  vector<xValKey> keys;
  vector<int> component;
  map<int,xValKey> master_key_per_component;// regrouping information of masters and slaves...
  map<int,double> coefficient_per_component;// either 1 or -1, regrouping information of masters and slaves...

  // variables used only for the graph
  // graph is composed of discontinuous entities, the groups are called colors
  vector<entity_id> dentities_tobecolored;
  vector<int> graph_color;
  multimap<int,entity_id> color_to_dentities;
  map<int,int> master_color_per_color;
  map<AOMD::mEntity*,int> number_of_colors_per_entity;

  Graph G;
};

class xValueCreatorRampHeaviside : public xValueCreatorRampHeavisideBase {
  public:
    xValueCreatorRampHeaviside( xMesh* m, const xcut::xPhysSurfByTagging &front, const xField &_field);
    ~xValueCreatorRampHeaviside();
    xValue<double>* operator()(const xValKey& key);
};

class xValueUpdatorRampHeaviside  : public xValueCreatorRampHeavisideBase 
{
  public:
    xValueUpdatorRampHeaviside( xMesh* m, const xcut::xPhysSurfByTagging &front, const xField &_field,const std::string& sub,const int,const int,const int );
    ~xValueUpdatorRampHeaviside();
    xValue<double>* operator()(const xValKey& key);
    void postUpdate(void);
  private:
    set<xValKey> old_disco_keys;
    set<xValKey> new_disco_keys;
    const std::string& sub;

};


//does not work and I do not know why !!
//  class xFirstIfNotSecondValueCreator {

//  public:
//    xFirstIfNotSecondValueCreator(const std::function<xValue<double>* (const xValKey&)>& f, 
//  			         const std::function<xValue<double>* (const xValKey&)>& s) : 
//                                  cfirst(f), csecond(s)  {} 
//    xValue<double>* operator()(const xValKey& key)
//      {
//        const bool debug = false;
//        if (debug) std::cout << " key is " << key << std::endl;

//        if (xValue<double>* ret =  cfirst(key)) 
//  	{
//  	  if (debug) 
//  	    {
//  	      std::cout << "cfirst did create a value " << std::endl;
//  	      ret->print(std::cout);
//  	    }
//  	  return ret;
//  	}
//        xValue<double>* ret2 = csecond(key);
//        if (debug) 
// 	 {
// 	   if (!ret2) std::cout << "csecond did create a value " << std::endl;

// 	 }
//        return ret2;
//      }

//    private:
//      const std::function<xValue<double>* (const xValKey&)>&            cfirst;
//      const std::function<xValue<double>* (const xValKey&)>&           csecond;
//    };


class xTensorsValueCreator 
{
  public:
    xTensorsValueCreator() {}
    virtual ~xTensorsValueCreator() {}
    xTensorsValueCreator(xTensorsSignature &sig):signature(&sig) {}
    virtual xValue< tensorsPtr_t >* operator()(const xValKey& key);
    void setSignature(xTensorsSignature* s) { signature =  s;}
  protected:
    xTensorsSignature* signature;
};

} // end of namespace

#include "xValue_imp.h"
#endif
