/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include "xValue.h"
#include "xSpacePolynomial.h"
#include "xStateOfValue.h"
#include "xMesh.h"
#include "xTensors.h"
#include "xPhysSurfByTagging.h"
#include "xTensorsPtr.h"
#include "xApproxFunction.h"
#include "xField.h"
#include "xLevelSet.h"

namespace xfem
{
using std::endl;
using AOMD::mEdge;
using AOMD::mEntity;
using Trellis_Util::mPoint;

xValueDouble::xValueDouble() : value(0.0) {}



xValueLinearCombination::xValueLinearCombination(const coeffs_t& c, const values_t& v, double coeff_last) 
          : coeffs(c), values(v)
{ 
  const bool debug = xdebug_flag;
  if (debug) std::cout << "creating value average with " << values.size() << " values " << std::endl;
  v_last = new xValueDouble;
  v_last->setVal(1.0);
  v_last->setState(new xStateOfValueFixed(v_last));
  coeffs.push_back(coeff_last);
  values.push_back(v_last);
}
xValueLinearCombination::xValueLinearCombination(double c, xValue<double>* v, double coeff_last)
{ 
  coeffs.push_back(c); 
  values.push_back(v); 
  v_last = new xValueDouble;
  v_last->setVal(1.0);
  v_last->setState(new xStateOfValueFixed(v_last));
  coeffs.push_back(coeff_last);
  values.push_back(v_last);
}
//   xValueLinearCombination(double c1, xValue<double>* v1, double c2, xValue<double>* v2) 
//     { coeffs.push_back(c1); values.push_back(v1); 
//       coeffs.push_back(c2); values.push_back(v2); }

double xValueLinearCombination::getVal() const           
{ 
  double v = 0.0;
  iterator_value itv = values.begin();
  iterator_coeff itc = coeffs.begin();
  for (; itv != values.end(); ++itv, ++itc) 
    v += *itc * (*itv)->getVal();
  return v;
}
void   xValueLinearCombination::setVal(double in) {}
const xStateOfValue* xValueLinearCombination::getState()  const
{ 
  if (!state) { state = new xStateOfValueLinearCombination(this); }
  return state;
}
xStateOfValue* xValueLinearCombination::getState() 
{ 
  if (!state) { state = new xStateOfValueLinearCombination(this); }
  return state;
}
void xValueLinearCombination::setState(xStateOfValue* s) {}
std::ostream& xValueLinearCombination::printVal(std::ostream& o) const 
{ 
  o << "linear combination of the " << values.size() << " values below " << endl;
  o << "the coefficients are "     << endl;
  for (iterator_coeff it = coeffs.begin(); it != coeffs.end(); ++it)  o << *it            << endl;
  o << "and the values are   "     << endl;
  for (iterator_value it = values.begin(); it != values.end(); ++it) (*it)->print(o) << endl;
      return o;
}


xValueError::xValueError() : abs2_loc(0.0), eng_loc(0.0), vol_loc(0.0) {}


double xValueError::vol_tot =  0;
double xValueError::abs2_tot = 0;
double xValueError::eng_tot =  0;
double xValueError::abs2_exa_tot = 0;
double xValueError::eng_exa_tot =  0;
xValueError::choice_t xValueError::_choice = xValueError::ABS2;


void xValueError::clear()
{
  vol_tot =  0.;
  abs2_tot  = 0.;
  eng_tot =  0.;
  abs2_exa_tot =  0.;
  eng_exa_tot = 0.;  
}


double xValueError::total(const std::string& s)
{
  if (s == "VOL")          return vol_tot;
  else if (s == "ABS")     return std::sqrt(abs2_tot);
  else if (s == "ABS_EXA") return std::sqrt(abs2_exa_tot);
  else if (s == "ENG")     return eng_tot;
  else if (s == "ENG_EXA") return eng_exa_tot;
  else if (s == "REL")     return std::sqrt(abs2_tot/eng_tot);
  else if (s == "REL_EXA") return std::sqrt(abs2_exa_tot/eng_exa_tot);
  else if (s == "EFF_EXA") return std::sqrt(abs2_tot/abs2_exa_tot);
  else { assert(0), throw;};
}

void xValueError::choice(const std::string& s)
{
  if (s == "ABS2")     _choice = ABS2;
  else if (s == "ABS2_EXA") _choice = ABS2_EXA;
  else if (s == "ABS") _choice = ABS;
  else if (s == "ABS_EXA") _choice = ABS_EXA;
  else if (s == "VOL") _choice = VOL;
  else if (s == "REL") _choice = REL;
  else if (s == "REL_EXA") _choice = REL_EXA;
  else if (s == "DNS") _choice = DNS;
  else if (s == "ENG") _choice = ENG;
  else if (s == "ENG_EXA") _choice = ENG_EXA;
  else if (s == "EFF_EXA") _choice = EFF_EXA;
  else assert(0);
}

void   xValueError::setVal(double in)
{
  switch(_choice)
    {
    case ABS2:     abs2_loc = in;      abs2_tot     +=  in;           break;
    case VOL:      vol_loc = in;       vol_tot      +=  in;           break;
    case ENG:      eng_loc = in;       eng_tot      +=  in;           break;
    case ABS2_EXA: abs2_exa_loc = in;  abs2_exa_tot +=  in;           break;
    case ENG_EXA:  eng_exa_loc  = in;  eng_exa_tot  +=  in;           break;
    default: assert(0); break;
    }
  
}

double   xValueError::getVal() const
{
  switch(_choice)
    {
    case ABS2:     return abs2_loc;            break;
    case ABS2_EXA: return abs2_exa_loc;            break;
    case ABS:      return std::sqrt(abs2_loc); break;
    case ABS_EXA:  return std::sqrt(abs2_exa_loc); break;
    case VOL:      return vol_loc;             break;
    case ENG:      return eng_loc;             break;
    case ENG_EXA:  return eng_exa_loc;             break;
    case REL:      return std::sqrt(abs2_loc/eng_tot); break;
    case REL_EXA:  return std::sqrt(abs2_exa_loc/eng_exa_tot); break;
    case DNS:      return std::sqrt(abs2_loc*vol_tot/(eng_tot*vol_loc)); break;
    case EFF_EXA:  return std::sqrt(abs2_loc/abs2_exa_loc); break;
    default : assert(0); throw; break;
    }
}


xMirrorCreator::xMirrorCreator(xValueManagerDist<double>* v, xMesh* m) : double_manager(v), mesh(m) {} 

//create a mirror if a value exist within the mirrors
xValue<double>* xMirrorCreator::operator()(const xValKey& key) 
{ 
  const bool debug = xdebug_flag;
  if (debug) cout << " inside xValue<double>* xMirrorCreator::operator()(const xValKey& key) " << endl;
  std::list<AOMD::mEntity*> mirrors;
  if (!mesh->lookupForMirrorEntities(key.getEnti(), mirrors)) return 0;
  xValKey keyother = key; 
  if (debug) double_manager->PrintForDebug("in_create_mirror.dbg");
  if (debug) cout << "In xMirrorCreator key  is " << key << endl;
  std::list<mEntity*>::const_iterator it;
  for (it = mirrors.begin() ; it != mirrors.end(); ++it)
    {
      keyother.setEnti(*it);
      if (debug) cout << "In xMirrorCreator key of mirror is " << keyother << endl;
      //check if the other value exists
      if (xValue<double>* other = double_manager->find(keyother)) {
	if (debug) cout << "In CreateMirror_c, mirror is created " << endl;
	if (debug) other->getValPtr(); //to check the sanity of other
	if (debug) cout << "after sanity check " << endl;
	return new xValueLinearCombination(1.0, other->getValPtr());	    
      }
    }
  return 0;
}
xValueCreatorRegularAndHanging::xValueCreatorRegularAndHanging(xValueManagerDist<double>* v, xMesh* m, int degree_) 
  : double_manager(v), mesh(m), degree(degree_) {} 


xValue<double>* xValueCreatorRegularAndHanging::degreeOne(const xValKey& key)
{
  const bool debug = xdebug_flag;

  mEntity* e = key.getEnti();
  if (debug) { cout << "in hanging entity is " << endl; e->print(); }
  mEntity* h = e->getAttachedEntity(xMesh::get_is_hanging_on_tag());
  if (!h) 
    {
      xValue<double>* val = creator_reg(key);
      if (debug) { cout << " regular value created "  << endl; }
      return val;
    }

  if (debug && h) { cout << "hanging is " << endl; h->print(); }
  xValKey keyother = key; 
  xValueLinearCombination::coeffs_t coeffs;
  xValueLinearCombination::values_t values;
  for (int i = 0; i < h->size(0); ++i)
    {
      mEntity* v = h->get(0,i);
      keyother.setEnti(v);
      xValue<double>* other;
      if (!(other = double_manager->find(keyother))) return 0;
      coeffs.push_back(1/(double)h->size(0));
      values.push_back(other);
    }
  if (debug) cout << "creating a value linear combination" << endl;
  return new xValueLinearCombination(coeffs, values);
}


xValue<double>* xValueCreatorRegularAndHanging::degreeTwo(const xValKey& key)
{
  const bool debug = xdebug_flag;

  mEntity* e = key.getEnti();
  assert(e->size(3) == 0); // on vérifie que l'on est en 2D, le 3D n'est pas codé. 

  if (debug) { cout << " in hanging entity is " << endl; e->print(); }
  mEntity* h = e->getAttachedEntity(xMesh::get_is_hanging_on_tag());
  if (!h) 
    {
      xValue<double>* val = creator_reg(key);
      if (debug) { cout << " regular value created "  << endl; }
      return val;
    }

  if (debug && h) { cout << "hanging is " << endl; h->print(); }
  xValKey keyother = key; 
  xValueLinearCombination::coeffs_t coeffs;
  xValueLinearCombination::values_t values;
  xValue<double>* other;
  

  if (e->getLevel() == 0)
    {
      for (int i = 0; i < h->size(0); ++i)
	{
	  mEntity* v = h->get(0,i);
	  keyother.setEnti(v);
	  if (!(other = double_manager->find(keyother))) return 0;
	  coeffs.push_back(1/(double)h->size(0));
	  values.push_back(other);
	} 
      keyother.setEnti(h);
      if (!(other = double_manager->find(keyother))) return 0;
      values.push_back(other);
      coeffs.push_back(HierarchicalApproxFunctionEdge(3, 0.));
      if (debug) cout << "creating a value linear combination" << endl;
      return new xValueLinearCombination(coeffs, values);
    }
  else 
    {
      for (int i = 0; i < h->size(0); ++i)
	{
	  mEntity* v = h->get(0,i);
	  keyother.setEnti(v);
	  if (!(other = double_manager->find(keyother))) return 0;
	  values.push_back(other);
	} 
      keyother.setEnti(h);
      if (!(other = double_manager->find(keyother))) return 0;
      values.push_back(other);
      //getting the right coeffs
      if (debug) 
	{
	  cout << " small edge " << endl;
	  e->print();
	  cout << " big edge   " << endl;
	  h->print();
	}
      xGeomElem geo_small(e);
      geo_small.setUVW(Trellis_Util::mPoint(0.,0.,0.));
      xGeomElem geo_big(h);
      geo_big.setUVWForXYZ(geo_small.getXYZ());
      Trellis_Util::mPoint uvw = geo_big.getUVW();
      double u = uvw(0);
      if (debug) cout << " In hanging degree 2 u is " << u << endl;
      if (u < 0.)
	{
	  coeffs.push_back((HierarchicalApproxFunctionEdge(1,u) - HierarchicalApproxFunctionEdge(1,0.) * HierarchicalApproxFunctionEdge(2,0.) - HierarchicalApproxFunctionEdge(1,0.))/(HierarchicalApproxFunctionEdge(3,0.)));
	  coeffs.push_back((HierarchicalApproxFunctionEdge(2,u) - HierarchicalApproxFunctionEdge(2,0.) * HierarchicalApproxFunctionEdge(2,0.))/(HierarchicalApproxFunctionEdge(3,0.)));
	  coeffs.push_back((HierarchicalApproxFunctionEdge(3,u) - HierarchicalApproxFunctionEdge(3,0.) * HierarchicalApproxFunctionEdge(2,0.))/(HierarchicalApproxFunctionEdge(3,0.)));
	}
      else
	{
	  coeffs.push_back((HierarchicalApproxFunctionEdge(1,u) - HierarchicalApproxFunctionEdge(1,0.) * HierarchicalApproxFunctionEdge(1,0.))/(HierarchicalApproxFunctionEdge(3,0.)));
	  coeffs.push_back((HierarchicalApproxFunctionEdge(2,u) - HierarchicalApproxFunctionEdge(2,0.) * HierarchicalApproxFunctionEdge(1,0.)- HierarchicalApproxFunctionEdge(2,0.))/(HierarchicalApproxFunctionEdge(3,0.)));
	  coeffs.push_back((HierarchicalApproxFunctionEdge(3,u) - HierarchicalApproxFunctionEdge(3,0.) * HierarchicalApproxFunctionEdge(1,0.))/(HierarchicalApproxFunctionEdge(3,0.)));
	}
      return new xValueLinearCombination(coeffs, values);
    }
}






xValue<double>* xValueCreatorRegularAndHanging::operator()(const xValKey& key) 
{ 
  //const bool debug = xdebug_flag;
  if (degree == 1) return  degreeOne(key);
  else if (degree == 2) return  degreeTwo(key);
  else assert(0); //to be coded.
  return 0;
}

//----------------------------------------------------------------------------------------------------------------

xValueCreatorLinkOnFront::~xValueCreatorLinkOnFront(){
	vertexvector.clear();
}

//----------------------------------------------------------------------------------------------------------------

void xValueCreatorLinkOnFront::buildVertexVector(){
	creationinfo::iterator it = tab.begin();
	//cout << "building vector vextex " << endl;
	for (;it!=tab.end();it++){
		vertexvector.push_back(it->first);
	/*	cout << it->first->getId() << " " << it->first->point() << " ";
        if (it->second)
		    cout << it->second->getId() << " " << it->second->point() << endl;
        else
            cout << "NULL"<<endl;
    */
	}
	return;
}

//----------------------------------------------------------------------------------------------------------------

vector<AOMD::mVertex*>::iterator xValueCreatorLinkOnFront::beginVertexIter(){
			if (!vertexvector.size())
				buildVertexVector();
			return vertexvector.begin();
		}

//----------------------------------------------------------------------------------------------------------------

vector<AOMD::mVertex*>::iterator xValueCreatorLinkOnFront::endVertexIter(){
			if (!vertexvector.size())
				buildVertexVector();
			return vertexvector.end();
		}

//----------------------------------------------------------------------------------------------------------------

xValueCreatorLinkOnFront::nodeSortingCriteria::nodeSortingCriteria(normal_at_nodes_t * norms_)
    : norms(norms_)
{}

//----------------------------------------------------------------------------------------------------------------
void xValueCreatorLinkOnFront::normForElem2D(mEntity *e, xtensor::xVector& norm)
{
    xtensor::xVector edge_vect((static_cast<AOMD::mVertex *>(e->get(0,0)))->point(), (static_cast<AOMD::mVertex *>(e->get(0,1)))->point());
    xtensor::xVector plan_norm(0.,0.,1.);
    norm=plan_norm%edge_vect;
    norm.norm();
    return;
}
//----------------------------------------------------------------------------------------------------------------
void xValueCreatorLinkOnFront::normForElem3D(mEntity *e, xtensor::xVector& norm)
{
    Trellis_Util::mPoint p0 = (static_cast<AOMD::mVertex *>(e->get(0,0)))->point();
    xtensor::xVector edge_1(p0, (static_cast<AOMD::mVertex *>(e->get(0,1)))->point());
    xtensor::xVector edge_2(p0, (static_cast<AOMD::mVertex *>(e->get(0,2)))->point());
    norm=edge_1%edge_2;
    norm.norm();
    return;
}
//----------------------------------------------------------------------------------------------------------------
bool xValueCreatorLinkOnFront::nodeSortingCriteria::operator()(mEntity* ent1, mEntity* ent2) const
{

  normal_at_nodes_it_t f1 =  norms->find(ent1);
  assert( f1 != norms->end() );
  normal_at_nodes_it_t f2 =  norms->find(ent2);
  assert( f2 != norms->end() );

  if ((*f1).second(0)>(*f2).second(0)) return true;
  if ((*f1).second(0)<(*f2).second(0)) return false;
  return ((AOMD::mVertex *)ent1)->point().lexicographicLessThan (((AOMD::mVertex *)ent2)->point(),0.);
}



//----------------------------------------------------------------------------------------------------------------

xValueCreatorLinkOnFront::xValueCreatorLinkOnFront(xValueManagerDist<double>* v, const xMesh* m, bool link_isolated, xEntityFilter filter) :
  double_manager(v)
{
  if (m->size(0)==0)
    return;

  // nota : 
  //  sorting is base on 2 criteria :
  //    * node of front having related edge localy perpendicular to front mesh are treated first
  //    * node of front having same edge orthogonality criteria are sorted by theire coordinates
  //
  buildNormalAtNodes(m->begin(m->dim()), m->end(m->dim()));
  buildTab(link_isolated, filter);
}

void xValueCreatorLinkOnFront::buildTab(bool link_isolated, xEntityFilter filter)
{
  xEntityToEntity interf2appro = xCreator();
  const double zero=0.;
  double prod;

  // second, compute edge normal product at the nodes of front
  // and store it according to sorting criteria
  // =========================================================

  // generate container to store ordered nodes
  nodeSortingCriteria sc(&normal_at_nodes);
  node_set_t front_nodes(sc);

  // loop on container to avoid finding operation
  normal_at_nodes_it_t itvn = normal_at_nodes.begin();
  normal_at_nodes_it_t itvne = normal_at_nodes.end();
  for (; itvn != itvne; ++itvn )
  {
    mEntity *node = ( *itvn ).first;
    mEntity * e = interf2appro(node);
    if (e && filter(e))
    {
      ( *itvn ).second.norm();
      if(e->getLevel() == 1)
      {
        xtensor::xVector edge_vect((static_cast<AOMD::mVertex *>(e->get(0,0)))->point(), (static_cast<AOMD::mVertex *>(e->get(0,1)))->point());
        edge_vect.norm();
        prod =  (( *itvn ).second) * edge_vect ;
        prod = fabs( prod );
      }
      else
        prod=zero;

      // to minimize storage consumption use of xtensor::xVector first component to store product
      (*itvn).second(0)=prod;

      // this is inserting and sorting 
      front_nodes.insert(node);
    }
  }

  // map is no more used => clear it
  normal_at_nodes.clear();


  //set of nodes already visited  
  AOMD::mMeshEntityContainer visited_nodes;

  // third loop to select the vital edge
  // ====================================
  node_set_t::iterator beg_nodes=front_nodes.begin();
  node_set_t::iterator end_nodes=front_nodes.end();
  node_set_t::iterator it_nodes;
  for(it_nodes = beg_nodes; it_nodes != end_nodes; ++it_nodes){
    mEntity * e_bnd = *it_nodes;
    mEntity * e = interf2appro(e_bnd);

    if (e->getLevel() == 0)
    {
      // regular dof will be created
      tab.insert(pair<AOMD::mVertex*, AOMD::mVertex*>((AOMD::mVertex*)e,(AOMD::mVertex *)NULL));
      visited_nodes.add(e);
    }
    else if(e->getLevel() == 1)
    {
      //check if 0, 1 or 2 nodes have dofs
      //if 0, will have to create two dofs and link them
      mEntity* v1 = e->get(0,0);
      mEntity* v2 = e->get(0,1);
      if (!visited_nodes.find(v1) && !visited_nodes.find(v2)){
        tab.insert(pair<AOMD::mVertex*, AOMD::mVertex*>((AOMD::mVertex*)v1,(AOMD::mVertex *)NULL));
        tab.insert(pair<AOMD::mVertex*, AOMD::mVertex*>((AOMD::mVertex*)v2,(AOMD::mVertex*)v1));
        visited_nodes.add(v1);
        visited_nodes.add(v2);
      }
    }
  }
  // fourth loop for remaining vertices
  // ==================================
  for(it_nodes = beg_nodes; it_nodes != end_nodes; ++it_nodes)
  {
    mEntity * e_bnd = *it_nodes;
    mEntity * e = interf2appro(e_bnd);
    if(e->getLevel() == 1)
    {
      //check if 0, 1 or 2 nodes have dofs
      //if 1, create the other dof and link it
      mEntity* v1 = e->get(0,0);
      mEntity* v2 = e->get(0,1);
      bool b1 = (visited_nodes.find(v1));
      bool b2 = (visited_nodes.find(v2));
      assert(b1 || b2);
      if (!b1 || !b2)	{
        if (b2) swap(v1, v2);
        if (link_isolated)
          tab.insert(pair<AOMD::mVertex*, AOMD::mVertex*>((AOMD::mVertex*)v2,(AOMD::mVertex*)v1));
        else
          tab.insert(pair<AOMD::mVertex*, AOMD::mVertex*>((AOMD::mVertex*)v2,(AOMD::mVertex*)NULL));
        visited_nodes.add(v2);
      }
    }
  }
} 

//----------------------------------------------------------------------------------------------------------------

/*double xValueCreatorLinkOnFront::compute_orthogonality(mEdge *e, xtensor::xVector vector1, AOMD::mVertex *v){
  AOMD::mVertex *previous = (AOMD::mVertex*)(e->get(0,0));
	if (previous==v) previous = (AOMD::mVertex*)(e->get(0,1));
	xtensor::xVector vector2(previous->point(),v->point());
	vector2.normValue();// normalize...
	return fabs(vector1*vector2);
}*/

//----------------------------------------------------------------------------------------------------------------

xValue<double>* xValueCreatorLinkOnFront::operator()(const xValKey& key) { 
	bool debug=false;

	// Remark: xValueCreator MUST return a xValue*
	// It cannot just "do nothing"
	// Therefore, for the "Lag Mult Space", DeclareInterpolation must receive an iterator on the nodes, and not on the elements around the front
	// which means that the iterator given to DeclareInterpolation must be striclty equal to the vertex tab of the creator...

	mEntity* e = key.getEnti();

	// check if mEntity *e is known in the table
	creationinfo::const_iterator it = tab.find((AOMD::mVertex*)e);
	if (it==tab.end())
		throw;// should not happen...cf remark above...

    AOMD::mVertex *vertex = it->first;
    AOMD::mVertex *vertex2 = it->second;

	if (vertex2==NULL){// create one simple value
		xValueCreator<xValueDouble> creator_reg;
		if (debug) cout << "creating regular value for node " << vertex << endl;
		return creator_reg(key);
	}
	else{// create with linear combination state to the second
		xValKey keyother = key; 
		keyother.setEnti(vertex2);
		xValue<double>* other;

		if (!(other = double_manager->find(keyother))){
			if (debug) cout << "value for node " << vertex2 << " not declared yet, no value yet for node " << vertex << endl;
			return 0; // if the second value does not exist, wait until it's created...
		}
		else{
			if (debug) cout << "linking " << vertex << " to " << vertex2 << endl;
			return new xValueLinearCombination(1., other);
		}
	}
}

//----------------------------------------------------------------------------------------------------------------

void xValueCreatorLinkAwayFromFront::getIntersect(Trellis_Util::mPoint first,Trellis_Util::mPoint second,Trellis_Util::mPoint origin,Trellis_Util::mPoint &intersect,double &distance, bool &valid, double &d1, double &d2){
	// computes distance between point origin and the line first-second
	// point intersect is on the segment and minimizes the distance origin-intersect
	double xa = second(0)-first(0);
	double ya = second(1)-first(1);
	double xb = first(0)-origin(0);
	double yb = first(1)-origin(1);
	double lengtha = sqrt(xa*xa+ya*ya);
	double oriented_distance = (xa*yb-ya*xb)/lengtha;
	distance = fabs(oriented_distance);
	// vector v is OP normalized
	double vx = -ya;
	double vy = xa;
	double normv = sqrt(vx*vx+vy*vy);
	vx = vx/normv;
	vy = vy/normv;

	intersect(0) = origin(0) + oriented_distance*vx;
	intersect(1) = origin(1) + oriented_distance*vy;

	// relative distance intersection-first and intersection-second in [0,1]
	double xi = intersect(0)-first(0);
	double yi = intersect(1)-first(1);
	d1 = (xi*xa+yi*ya)/lengtha/lengtha;
	d2 = 1-d1;

	// now, if intersection is outside the edge, need to choose one edge's vertex
	if ((d1<0)||(d2<0)){
		d1 = max(min(d1,1.),0.);
		d2 = max(min(d2,1.),0.);
		valid=false;
		if (d1<0)
			for (int k=0;k<3;k++)
				intersect(k) = first(k);
		else
			for (int k=0;k<3;k++)
				intersect(k) = second(k);
	}
	else
		valid=true;// intersect is already on the edge
}

//----------------------------------------------------------------------------------------------------------------

xValue<double>* xValueCreatorLinkAwayFromFront::operator()(const xValKey& key) { 

	mEntity* e = key.getEnti();

	// check that mEntity *e is known in the table
	linkinfo::const_iterator it = tab.find(e);
	if (it==tab.end())
		throw;// should not happen...
		//return NULL;


	mEntity *element = it->second.first;
	Trellis_Util::mPoint uvw(it->second.second);

	xValKey keyother = key; 
	xValueLinearCombination::coeffs_t coeffs;
	xValueLinearCombination::values_t values;

	// computing coefficients for linear combination
	switch(element->getLevel()){
		case (0):
			coeffs.push_back(1.);
			break;
		case (1):
			coeffs.push_back(0.5*(1-uvw(0)));
			coeffs.push_back(0.5*(1+uvw(0)));
			break;
		case (2):
			coeffs.push_back(1-uvw(0)-uvw(1));
			coeffs.push_back(uvw(0));
			coeffs.push_back(uvw(1));
			break;
		case (3):
			coeffs.push_back(1-uvw(0)-uvw(1)-uvw(2));
			coeffs.push_back(uvw(0));
			coeffs.push_back(uvw(1));
			coeffs.push_back(uvw(2));
		default:
			throw;
	}
	int number_of_values_to_deal_with = std::max(element->size(0),1);
	mEntity* v;
	for (int i=0;i<number_of_values_to_deal_with;i++){
		if (element->getLevel()!=0) v = element->get(0,i);
		else v = element;
		keyother.setEnti(v);
		xValue<double>* other;
		if (!(other = double_manager->find(keyother))) throw;// the "lagrange" space must already exist along the front  //return 0;
		values.push_back(other);
	}

	return new xValueLinearCombination(coeffs, values);
}

//----------------------------------------------------------------------------------------------------------------

xValue<double>* xValueCreatorLinkAwayFromFrontInsideRadius::operator()(const xValKey& key) {
	const bool debug = xdebug_flag;


	mEntity* e = key.getEnti();

	if (debug){
		std::string filename = "testANN_operator.txt";
		std::ofstream out(filename.c_str(),ios_base::app);
		out << "operator() for Entity " << e->getId() << " of lvl " << e->getLevel() << endl;
		out.close();
	}

	// check that mEntity *e is known in the table
	linkinfo::const_iterator it = tab.find(e);
	if (it==tab.end())
		throw;// should not happen...
		//return NULL;

	const coeffsinfo *ci = &(it->second);

	xValueLinearCombination::coeffs_t coeffs;
	xValueLinearCombination::values_t values;
	xValKey keyother = key; 

	for (coeffsinfo::const_iterator itci=ci->begin();itci!=ci->end();itci++){
	
		mEntity *element = itci->first;
		Trellis_Util::mPoint uvw(itci->second.first);
		double frontnodecoef = itci->second.second;

		// computing coefficients for linear combination
		switch(element->getLevel()){
			case (0):
				coeffs.push_back(frontnodecoef);
				break;
			case (1):
				coeffs.push_back(0.5*(1-uvw(0))*frontnodecoef);
				coeffs.push_back(0.5*(1+uvw(0))*frontnodecoef);
				break;
			case (2):
				coeffs.push_back((1-uvw(0)-uvw(1))*frontnodecoef);
				coeffs.push_back(uvw(0)*frontnodecoef);
				coeffs.push_back(uvw(1)*frontnodecoef);
				break;
			case (3):
				coeffs.push_back((1-uvw(0)-uvw(1)-uvw(2))*frontnodecoef);
				coeffs.push_back(uvw(0)*frontnodecoef);
				coeffs.push_back(uvw(1)*frontnodecoef);
				coeffs.push_back(uvw(2)*frontnodecoef);
			default:
				throw;
		}
		int number_of_values_to_deal_with = std::max(element->size(0),1);
		mEntity* v;
		for (int i=0;i<number_of_values_to_deal_with;i++){
			if (element->getLevel()!=0) v = element->get(0,i);
			else v = element;
			keyother.setEnti(v);
			xValue<double>* other;
			if (!(other = double_manager->find(keyother))) throw;// the "lagrange" space must already exist along the front  //return 0;
			values.push_back(other);
		}
	}

	return new xValueLinearCombination(coeffs, values);
}

//----------------------------------------------------------------------------------------------------------------
xValueLinkDiscontinuousSpaceOnFrontBase::~xValueLinkDiscontinuousSpaceOnFrontBase(){
}

//----------------------------------------------------------------------------------------------------------------
xValueCreatorLinkDiscontinuousSpaceOnFront::~xValueCreatorLinkDiscontinuousSpaceOnFront(){
}

//----------------------------------------------------------------------------------------------------------------
xValueUpdatorLinkDiscontinuousSpaceOnFront::~xValueUpdatorLinkDiscontinuousSpaceOnFront(){
}

//----------------------------------------------------------------------------------------------------------------
xValueLinkDiscontinuousSpaceOnFrontBase::xValueLinkDiscontinuousSpaceOnFrontBase(const xMesh* m, const xLevelSet &ls, const xField &field) :
										double_manager((xValueLinkDiscontinuousSpaceOnFrontBase::xValueManagerDist<double> *)field.getDoubleManager()),frontmesh(m), highest_dim(m->getDim()+1)
{

	if (frontmesh->size(0)==0) return;

	xEntityToEntity interf2appro = xCreator();
    AOMD::mMeshEntityContainer visited_nodes, visited_faces;
	int nb_connected_edges,nb_connected_faces,nb_connected_volumes;

	for(xIter it = frontmesh->begin(0); it != frontmesh->end(0); ++it){
		mEntity * front_vertex = *it;
		mEntity * cut_edge = interf2appro(front_vertex);// either a vertex or an edge

		// if "cut_edge" is a vertex, nothing to be done !
		if (cut_edge->getLevel() == 1){
			for (int inode=0;inode<cut_edge->size(0);inode++){// for the two vertices of edge, treating the dofs of order 1 (on vertices)
                AOMD::mVertex *v = (AOMD::mVertex*)cut_edge->get(0,inode);

				// if vertex has not been treated yet
				if (!visited_nodes.find(v)){
					visited_nodes.add(v);
					nb_connected_edges = v->size(1);
					for (int iedge=0;iedge<nb_connected_edges;iedge++){// for all connected edges
						mEdge* ce = (mEdge*)v->get(1,iedge);
						nb_connected_faces = ce->size(highest_dim);
						if (nb_connected_faces!=1){
//							mFace *face_right = (mFace*)ce->get(2,0);
//							mFace *face_left = (mFace*)ce->get(2,1);
							double ls_val_1 = ls(ce->get(0,0));
							double ls_val_2 = ls(ce->get(0,1));
							if (ls_val_1*ls_val_2<0.){// if edge strictly cut

								// recovering all the keys of a couple < face - vertex >
								vector<xValKey> thekeys;

								for (int iface=0;iface<ce->size(highest_dim);iface++){
//									cout << "working on vertex " << v << "(" << v->getId() << ") and face " << v->get(highest_dim,iface) << endl;
									xField::const_iterator itspace = field.begin();
									xField::const_iterator itenspace = field.end();

									for (;itspace!=itenspace;itspace++){

										vector<xValKey> keyvec;
										(*itspace)->getKeys(ce->get(highest_dim,iface),&keyvec);

										// find the right vertex
//										cout << "size = " << keyvec.size() << endl;
										vector<xValKey>::iterator itvec=keyvec.begin();

										// export
//																			cout << "all keys:" << endl;
//																				for (;itvec!=keyvec.end();itvec++)
//																				cout << (*itvec) << endl;
//																				itvec=keyvec.begin();
										for (;itvec!=keyvec.end();itvec++)
											if ((*itvec).getEnti()==v){
												thekeys.push_back((*itvec));
//												cout << "considering key " << *itvec << endl;
											}

//																			// export
//																				vector<xValKey>::iterator itkeys=thekeys.begin();
//																				cout << "found the key(s) " << endl;
//																				for (;itkeys!=thekeys.end();itkeys++)
//																				cout << (*itkeys) << endl;
									}
								}

								// insert keys in the vector if not inserted yet
								vector<xValKey>::iterator thekeysit = thekeys.begin();
								vector<int> position;
								for (;thekeysit!=thekeys.end();thekeysit++){
									vector<xValKey>::iterator itfind = find(keys.begin(),keys.end(),(*thekeysit));
									if (itfind==keys.end()){
										keys.push_back((*thekeysit));
										position.push_back(keys.size()-1);
									}
									else
										position.push_back(itfind-keys.begin());// key position...
								}
//								for (int i=0;i<position.size();i++)
//									cout << position[i] << endl;
//								cout << "treating " << position.size() << " keys" << endl;
//								cout << "total " << keys.size() << " keys" << endl;

								// associating the keys in the graph
								int nbphys = position.size()/nb_connected_faces;
								int ilast = position.size();
                                assert(nbphys*nb_connected_faces==ilast);
								for (int ipos=nbphys;ipos<ilast;++ipos){
									add_edge(position[ipos%nbphys], position[ipos], G);
/*
 *  bugged version in 3D 
                                int mid = position.size()/2;
								for (int ipos =0;ipos<mid;ipos++){
									add_edge(position[ipos], position[mid+ipos], G);
*/
//									cout << "creating couple " << keys[position[ipos]] << "------" << keys[position[mid+ipos]] << endl;
								}
							}
						}
					}
				}
			}

			// now, still have to treat the high order keys on both sides of the edge... valid for 2D and 3D
			nb_connected_faces = cut_edge->size(highest_dim);
			if (nb_connected_faces!=1){

				vector<xValKey> thekeys;
				for (int iface=0;iface<cut_edge->size(highest_dim);iface++){
					//cout << "working on edge " << cut_edge << "(" << cut_edge->getId() << ") with nodes " << cut_edge->get(0,0)->getId() << ", " << cut_edge->get(0,1)->getId() <<  " and face " << cut_edge->get(2,iface) << "(" << cut_edge->get(2,iface)->getId() << ")" << endl;
					xField::const_iterator itspace = field.begin();
					xField::const_iterator itenspace = field.end();

					for (;itspace!=itenspace;itspace++){
						vector<xValKey> keyvec;
						(*itspace)->getKeys(cut_edge->get(highest_dim,iface),&keyvec);

						// find the right edge 
						vector<xValKey>::iterator itvec=keyvec.begin();

						for (;itvec!=keyvec.end();itvec++)
							if ((*itvec).getEnti()==cut_edge){
								thekeys.push_back((*itvec));
							}
					}
				}
				// insert keys in the vector
				vector<xValKey>::iterator thekeysit = thekeys.begin();
				vector<int> position;
				for (;thekeysit!=thekeys.end();thekeysit++){
						keys.push_back((*thekeysit));
						position.push_back(keys.size()-1);
				}
				// associating the keys in the graph
/*
 * bugged version in 3D
				int mid = position.size()/2;
				for (int ipos =0;ipos<mid;ipos++){
					add_edge(position[ipos], position[mid+ipos], G);
*/
			    int nbphys = position.size()/nb_connected_faces;
				int ilast = position.size();
                assert(nbphys*nb_connected_faces==ilast);
			    for (int ipos=nbphys;ipos<ilast;++ipos){
					add_edge(position[ipos%nbphys], position[ipos], G);
					//cout << "creating couple " << keys[position[ipos]] << "------" << keys[position[mid+ipos]] << endl;
				}
			}


			// and the faces... valid in 3D only
			if (highest_dim==3){
				nb_connected_faces = cut_edge->size(2);
				if (nb_connected_faces>2){
					for (int iface=0;iface<cut_edge->size(2);iface++){// link keys of each face surrounding cut_edge... if not done yet !
						mEntity *current_face = cut_edge->get(2,iface);

						if (!visited_faces.find(current_face)){// if face not visited yet
							visited_faces.add(current_face);

							nb_connected_volumes = current_face->size(3);
							if (nb_connected_volumes!=1){
								vector<xValKey> thekeys;
								for (int ivol=0;ivol<nb_connected_volumes;ivol++){
									xField::const_iterator itspace = field.begin();
									xField::const_iterator itenspace = field.end();

									for (;itspace!=itenspace;itspace++){
										vector<xValKey> keyvec;
										(*itspace)->getKeys(current_face->get(3,ivol),&keyvec);

										// find the right face 
										vector<xValKey>::iterator itvec=keyvec.begin();

										for (;itvec!=keyvec.end();itvec++)
											if ((*itvec).getEnti()==current_face){
												thekeys.push_back((*itvec));
											}
									}
								}
								// insert keys in the vector
								vector<xValKey>::iterator thekeysit = thekeys.begin();
								vector<int> position;
								for (;thekeysit!=thekeys.end();thekeysit++){
									keys.push_back((*thekeysit));
									position.push_back(keys.size()-1);
								}
								// associating the keys in the graph
								int mid = position.size()/2;
								for (int ipos =0;ipos<mid;ipos++){
									add_edge(position[ipos], position[mid+ipos], G);
									//cout << "creating couple " << keys[position[ipos]] << "------" << keys[position[mid+ipos]] << endl;
								}
							}
						}
					}
				}
			}
		}
	}

	// connectivity of the dofs computed with boost graph library
	std::vector<int> component_temp(boost::num_vertices(G));
	// this copy shouldn't be necessary... ?!
	component = component_temp;
	int num = connected_components(G, &component[0]);

	std::vector<int>::size_type i;
        const std::vector<int>::size_type ic = component.size();
	//cout << "Total number of components: " << num << endl;
	// creation of the multimap
	for (i = 0; i < ic; ++i){
		//cout << "Vertex " << i << " = key " << keys[i] <<" is in component " << component[i] << endl;
                sets_of_keys.insert(make_pair(component[i], (long)i));
	}
	//cout << endl;
	
        is_treated.resize(ic,0);
	master_per_group.resize(num,-1);
}
//----------------------------------------------------------------------------------------------------------------
xValue<double>* xValueLinkDiscontinuousSpaceOnFrontBase::basicCreator(const xValKey& key)
{ 

	// first check if the key entity is a node, vertex, face or volume
	// if highest dimension, then creates a value
	if (key.getEnti()->getLevel()==highest_dim){
		xValueCreator<xValueDouble> creator_reg;
		return creator_reg(key);
	}


	// find the key in the list
	vector<xValKey>::iterator itfind = find(keys.begin(),keys.end(),key);

	if (itfind==keys.end()){
		// creates a single value
		xValueCreator<xValueDouble> creator_reg;
		return creator_reg(key);
	}

	// recover its set number
	int position = itfind-keys.begin();

	// recover all keys matching this set "position" (same group of keys to link)
	pair<groupe_type::iterator, groupe_type::iterator> pos;
	pos = sets_of_keys.equal_range(component[position]);

	// see if one of these keys has been created yet
	bool already_treated=false;
	xValue<double> *other;
	groupe_type::iterator it=pos.first;
        groupe_type::iterator itend=pos.second;
	for (;it!=itend;it++)
        {
                     // get position in keys of this key of the set 
                     int itposition = it->second;

                     // if already treated key, check if every things are OK and 
                     // set to "other" the master per groupe
                     if (is_treated[itposition])
                     {
                          // check that there is not a bug with is_treated
                          if (itposition!=position)
                          {
                            already_treated=true;
                            if ( (other = double_manager->find(keys[master_per_group[component[position]]])) ){
                                      break;
                            }
                            else
                            {
                               cout<<"already treated but not in the double manager ????"<<endl;
                               throw 1;
                            }
                          }
                          else
                          {
                               cout<<"itposition=position shouldn't be possible as is_treated[itposition] is true and position is curently in treatement and not flaged as is_treated"<<endl;
                               throw 1;
                          }
                     }
             }


             // Wathever a other key of the set have been treated or not, this key is
             // now considere as treated
             // Put the corect flag 
             is_treated[position]=1;

             // if one of the key of the set have been already treated
             // linking to the existing value 
             if (already_treated){
                return new xValueLinearCombination(1., other);
             }
             // if it's the first key of the set to be treated
             // create a new value 
             // and set it as the master for the groupe
             else{
                master_per_group[component[position]]=position;
                xValueCreator<xValueDouble> creator_reg;
                return creator_reg(key);
             }


}
//----------------------------------------------------------------------------------------------------------------
xValueCreatorLinkDiscontinuousSpaceOnFront::xValueCreatorLinkDiscontinuousSpaceOnFront(const xMesh* m, const xLevelSet &ls, const xField &field) :
          xValueLinkDiscontinuousSpaceOnFrontBase(m,ls,field) 
{
}

//----------------------------------------------------------------------------------------------------------------

xValue<double>* xValueCreatorLinkDiscontinuousSpaceOnFront::operator()(const xValKey& key) 
{ 
        return basicCreator(key);
}

//----------------------------------------------------------------------------------------------------------------
xValueUpdatorLinkDiscontinuousSpaceOnFront::xValueUpdatorLinkDiscontinuousSpaceOnFront(const xMesh* m, const xLevelSet &ls, const xField &field,
                                                                                       const std::string& sub_,const int physid1, const int physid2, const int physid3) :
                                                                                         xValueLinkDiscontinuousSpaceOnFrontBase(m,ls,field),sub(sub_) 
{
                     // loop on double_manager to :
                     //    fill the set of old dicontinus key
                     //    remove from subset all old discontinus keys
                     //    replace all lineare combination by independante value (with correct value)
                     //    remove state for values associeted to any kind of keys
                     //
                     xValueManagerDist<double>::map_iterator it = double_manager->begin();
                     xValueManagerDist<double>::map_iterator itend = double_manager->end();
                     for (;it!=itend;++it)
                     {
                          
                          
                          const xValKey key = it->first;
                          const int phys = key.getPhys();
                          // if it is a discontinuous key
                          if (phys==physid1||phys==physid2||phys==physid3) 
                          {
			       // add key to the groupe of old key
                               old_disco_keys.insert(key);

                               // clean subset
                               double_manager->erase(it->second,sub); 

                               // replace linear combination by it's value
                               xValueLinearCombination * linear_combi = dynamic_cast<xValueLinearCombination *> (it->second);
                               if (linear_combi)
                               {
                                 double val = linear_combi->getVal();
                                 delete linear_combi;
                                 it->second = new xValueDouble;
                                 xValueDouble * new_val =  dynamic_cast<xValueDouble *> (it->second);
                                 if (new_val)
                                     new_val->setVal(val);
                                 else throw 1;
                               }
                               else
                               {
                                  // reset state for dofs
                                  xValueDouble * dof =  dynamic_cast<xValueDouble *> (it->second);
                                  if (dof) dof->delState();
                               }


                          }
                          // other field
                          // here in future there may be a reset strategie of state of 
                          // other field. Why ?
                          // This is for the case of mixed dof : continus and discontinus
                          // are presentent in subset in random order. When removing discontinus dof from
                          // subset in above loop it change dof numbering as subset size is 
                          // decremented. The risk is to have then 2 value with the same dof number at creation.
                          // To remove this problem  reset here the state of dof to force
                          // renumbering of the wall subset is a solution.
                          // The other solution is to prevent this mixed order and create first state value for other
                          // filds and at last creates states for discontinus field. Now watever append above the
                          // numbering of the subset will be ok and state of other field can be enchanged.
                     }

}

//----------------------------------------------------------------------------------------------------------------
xValue<double>* xValueUpdatorLinkDiscontinuousSpaceOnFront::operator()(const xValKey& key) 
{ 

	// store all keys of the new discontinus field
        new_disco_keys.insert(key);
        
        // first check if the key existe already in the double manager
        xValue<double> *exist = double_manager->find(key);

        // if it exist update procedure
        if (exist)
        {
              // first check if the key entity is a node, vertex, face or volume
             // if highest dimension, then normal procedur creates a value
             // here we return the existing one
             if (key.getEnti()->getLevel()==highest_dim){
                return exist;
             }

             // find the key in the list
             vector<xValKey>::iterator itfind = find(keys.begin(),keys.end(),key);

             // if not find in the list of keys the normal procedur creates a single value
             // here we return the existing one
             if (itfind==keys.end()){
                return exist;
             }
        
             // recover its set number
             int position = itfind-keys.begin();

             // recover all keys matching this set "position" (same group of keys to link)
             pair<groupe_type::iterator, groupe_type::iterator> pos;
             pos = sets_of_keys.equal_range(component[position]);

             // see if one of these keys has been treated yet
             bool already_treated=false;
             xValue<double> *other;
             groupe_type::iterator it=pos.first;
             groupe_type::iterator itend=pos.second;
             for (;it!=itend;it++)
             {
                     int itposition = (it->second);
                     // already treated
                     if (is_treated[itposition])
                     {
                          if (itposition!=position)
                          {
                            already_treated=true;
                            if ( (other = double_manager->find(keys[master_per_group[component[position]]])) ){
                                      break;
                            }
                            else
                            {
                                cout<<"already treated but not in the double manager ????"<<endl;
                             throw 1;
                            }
                          }
                          else
                          {
                             cout<<"itposition=position shouldn't be possible as is_treated[itposition] is true and position exist"<<endl;
                             throw 1;
                          }
                     }
             }


             // Wathever a other key of the set have been treated or not, this key is
             // now considere as treated
             // Put the corect flag 
             is_treated[position]=1;

             // if one of the key of the set have been already treated
             // linking to the existing value 
             // update master value if null (new disco dof)
             if (already_treated){
 
                
                // if master value is null, as it is a new disco dof, set it 
                // with the actual value of exist. 
                // nota the first existing key treated in the groupe win's : it's value is used for the whole groupe
                // if no test it's the last 
                // leave the test for now to save a getVal/setVal
                double val = other->getVal();
                if (val==(double)0.)
                {
                   other->setVal(exist->getVal());
                }

                // remove existing value from the double manager (avoid memory leeks)
                //delete exist; // comented as done in the update of xValManager, better ? to be thinking
                return new xValueLinearCombination(1., other);
             }
             // if it's the first key of the set to be treated
             // normal procedure create a new value
             // here we return the existing one
             // and set it as the master for the groupe
             else{
                master_per_group[component[position]]=position;
                return exist;
             }
        }
        // it doesn't exist : back to standard creation procedure 
        else
            return basicCreator(key);
}

//----------------------------------------------------------------------------------------------------------------
void xValueUpdatorLinkDiscontinuousSpaceOnFront::postUpdate(void)
{

    // find old keys not used
    vector<xValKey> diff_disco_keys(old_disco_keys.size());
    vector<xValKey>::iterator it    = diff_disco_keys.begin();
    vector<xValKey>::iterator itend = set_difference (old_disco_keys.begin(), old_disco_keys.end(),
                                                      new_disco_keys.begin(), new_disco_keys.end(),
                                                      it);
    //destroy old keys not used from double manager
    for (;it!=itend;++it)
        double_manager->erase(*it); 
}

//----------------------------------------------------------------------------------------------------------------
xValue< tensorsPtr_t >* xTensorsValueCreator::operator()(const xValKey& key)
{
	xTensorsValuePtr* v = new xTensorsValuePtr;
	tensorsPtr_t t = v->getVal();
	t->setSignature(signature);
  return v;
}

  template <>
  size_t xNValue<double > ::current = 0;


//----------------------------------------------------------------------------------------------------------------
xValueCreatorRampHeavisideBase::~xValueCreatorRampHeavisideBase()
{
}
xValueCreatorRampHeaviside::~xValueCreatorRampHeaviside()
{
}
xValueUpdatorRampHeaviside::~xValueUpdatorRampHeaviside()
{
}

//----------------------------------------------------------------------------------------------------------------

xValueCreatorRampHeavisideBase::xValueCreatorRampHeavisideBase(xMesh* fm, const xcut::xPhysSurfByTagging &front, const xField &field) : double_manager(( xValueCreatorRampHeavisideBase::xValueManagerDist<double> * )field.getDoubleManager()), frontmesh(fm), highest_dim(fm->getDim()+1)
{

    const bool debug = false;

    if (frontmesh->size(0) == 0)
        return;

    const xMesh & mesh = front.getMesh();
    int nb_connected_element;
    int dim,highest_dim_m1 = highest_dim-1;

    // loop on dimension of sub entities to be treated
    for (dim = 0; dim < highest_dim_m1; dim++)
    {


        // loop on mesh to find  sub entities to be treated
        xIter itend = mesh.end(dim);
        for ( xIter it = mesh.begin(dim); it != itend; ++it)
        {
            mEntity *e = ( *it );

            // if entity is to be treated
            if (front.supportCutStrictlyEltWise(e) && front.strictOut(e))
            {
                // init local graph and container
                Graph local_G;
                vector < entity_id > local_dentities_tobecolored;

                // loop on descriptive adjancy entity of this entitity
                // this descriptive adjancy entitie is either a edge in 2D or a face in 3D
                // it give associated element
                // in following case it add link to local graph :
                //     - if all of it's associated element are cut
                //     - if one of it's associated element is cut and the other is touching
                //     - if all of it's associated element are touching
                // A special case is when  all of it's associated element are touching. With this
                // first pass we create something in the local graph but a verification of component must
                // be done after local graph creation. Because this may generated component without a single
                // element cut
                const int nb_adjancy = e->size(highest_dim_m1);
                for (int iadj = 0; iadj < nb_adjancy; iadj++)
                {
                    // if this descriptive adjancy entity is cut by or touch iso-zero inspect
                    // associated element
                    mEntity * adj = e->get(highest_dim_m1,iadj);
                    if ( front.cutStrictly(adj) || front.touchingOut(adj) )
                    {
                        bool add_graph = false;
                        vector < entity_id > temp;
                        nb_connected_element = adj->size(highest_dim);
                        // nb_connected_element should be <= 2 => verif as remaining implementation
                        // rely on that assertion
                        if (nb_connected_element  > 2)
                        {
                            cout<<"File "<<__FILE__<<" line "<<__LINE__<<" nb_connected_element="<<nb_connected_element<<" but should be 2 or 1 !\n";
                            throw;
                        }
                        // if only one element conected -> connect itself if cut
                        else if (nb_connected_element < 2)
                        {
                            mEntity *el1 = adj->get(highest_dim,0);

                            // if element is stricly cut
                            if (front.cutStrictly(el1))
                            {
                                temp.push_back(make_pair(e,el1));
                                temp.push_back(make_pair(e,el1));
                                add_graph = true;
                            }
                        }
                        // otherwise 2 elements
                        else
                        {
                            mEntity *el1 = adj->get(highest_dim,0);
                            mEntity *el2 = adj->get(highest_dim,1);

                            // if first element is stricly cut
                            if (front.cutStrictly(el1))
                            {
                                // if second element is cut or touching add to graph
                                if (front.cutStrictly(el2) || front.touchingOut(el2) )
                                    add_graph = true;
                            }
                            // if first element is touching
                            else if (front.touchingOut(el1))
                            {
                                // if second element is cut or touching add to graph
                                if (front.cutStrictly(el2) || front.touchingOut(el2) )
                                    add_graph = true;
                            }
                            if (add_graph)
                            {
                                temp.push_back(make_pair(e,el1));
                                temp.push_back(make_pair(e,el2));
                            }
                        }     // end else 2 elements

                        // add to local graph
                        if (add_graph)
                        {

                            // insert the elements of temp in the local vector if not inserted yet
                            vector < entity_id >::const_iterator tempit = temp.begin();
                            vector < int > position;
                            for (; tempit != temp.end(); tempit++)
                            {
                                vector < entity_id >::iterator itfind = find(local_dentities_tobecolored.begin(),local_dentities_tobecolored.end(),( *tempit ));
                                if (itfind == local_dentities_tobecolored.end())
                                {
                                    local_dentities_tobecolored.push_back(( *tempit ));
                                    position.push_back(local_dentities_tobecolored.size()-1);
                                }
                                else
                                    position.push_back(itfind-local_dentities_tobecolored.begin());            // key position...
                            }

                            const int n = position.size();
                            for (int ipos = 1; ipos < n; ipos++)
                            {
                                boost::add_edge(position[0], position[ipos], local_G);
                            }
                        }     // end add to local graph

                    }     // end if adjancy is touching or cut

                }     // end loop on adjency

                // now by inspecting local graph we check that component have at least one element cut and if yes add this component
                // to global graph
                //
                // first create componnent
                int local_num_vertice = boost::num_vertices(local_G);
                int local_container_size = local_dentities_tobecolored.size();
                assert(local_num_vertice == local_container_size);
                vector < int > local_graph_color(local_num_vertice);
                int nb_local_color = connected_components(local_G, &local_graph_color[0]);

                // second check and add
                for (int color = 0; color < nb_local_color; ++color)
                {
                    bool is_cut = false;
                    vector < entity_id > temp;
                    for (int j = 0; j < local_container_size; ++j)
                    {
                        if (local_graph_color[j] == color)
                        {
                            temp.push_back(local_dentities_tobecolored[j]);
                            if (front.cutStrictly(local_dentities_tobecolored[j].second))
                                is_cut = true;
                        }
                    }
                    if (is_cut)
                        add_to_graph(temp);
                }

            } // end if entitity is to be treated

        } // end loop on entitity of dimention dim

    } // end loop on dimention

    // treating entity of dimention highest_dim_m1 in a different way has theire adjancy is themself
    // 2D adjancy edge of a edge is itself
    // 3D adjancy face of a face is itself
    // loop on mesh to find  sub entities to be treated
    xIter itend = mesh.end(highest_dim_m1);
    for ( xIter it = mesh.begin(highest_dim_m1); it != itend; ++it)
    {
        mEntity *e = ( *it );


        // if entity is to be treated :
        //   - it's support is Cut with at least on element strictly cut
        //   - it is strictly in
        if (front.supportCutStrictlyEltWise(e) && front.noTouchingOut(e))
        {
            nb_connected_element = e->size(highest_dim);
            // nb_connected_element should be <= 2 => verif as remaining implementation
            // rely on that assertion
            if (nb_connected_element  > 2)
            {
                cout<<"File "<<__FILE__<<" line "<<__LINE__<<" nb_connected_element="<<nb_connected_element<<" but should be 2 or 1 !\n";
                throw;
            }
            // if only one element conected -> nothing to do
            // otherwise 2 elements
            else if (nb_connected_element > 1)
            {
                mEntity *el1 = e->get(highest_dim,0);
                mEntity *el2 = e->get(highest_dim,1);

                // if first and second element strictly cut
                if ( front.cutStrictly(el1) && front.cutStrictly(el2)  )
                {
                    vector < entity_id > temp1;
                    temp1.push_back(make_pair(e,el1));
                    temp1.push_back(make_pair(e,el1));
                    vector < entity_id > temp2;
                    temp2.push_back(make_pair(e,el2));
                    temp2.push_back(make_pair(e,el2));
                    add_to_graph(temp1);
                    add_to_graph(temp2);
                }
            }     // end else 2 elements

        } // end if entitity is to be treated

    } // end loop on entitity of dimention highest_dim_m1

    // computing connectivity of the dofs with boost graph library
    std::vector < int > component_temp(boost::num_vertices(G));
    // this copy shouldn't be necessary... ?!
    graph_color = component_temp;
    int nb_dentity_components = connected_components(G, &graph_color[0]);

    //cout << "Total number of components: " << num << endl;
    for (size_t i = 0; i != graph_color.size(); ++i)
    {
        //cout << "Vertex " << i << " = key " << keys[i] <<" is in component " << graph_color[i] << endl;
        color_to_dentities.insert(make_pair(graph_color[i], dentities_tobecolored[i]));
    }
    //cout << endl;

    // initialization
    for (vector < entity_id >::iterator it = dentities_tobecolored.begin(); it != dentities_tobecolored.end(); it++)
    {
        number_of_colors_per_entity.insert(make_pair(it->first,0));
    }
    // the number of compoennts for each (colored) entity
    multimap < mEntity *, int > entities_to_components;
    for (int icomp = 0; icomp < nb_dentity_components; icomp++)
    {
        const int color = icomp;
        pair < multimap < int,entity_id >::iterator, multimap < int,entity_id >::iterator > pos
	  = color_to_dentities.equal_range(color);
        mEntity *e = ( pos.first )->second.first;
        number_of_colors_per_entity[e] += 1;
        entities_to_components.insert(make_pair(e,( pos.first )->first));
    }
    for (vector < entity_id >::iterator it = dentities_tobecolored.begin(); it != dentities_tobecolored.end(); it++)
    {
        mEntity* e = it->first;
        pair < multimap < mEntity *,int >::iterator, multimap < mEntity *, int >::iterator > pos;
        pos = entities_to_components.equal_range(e);
        int master_component = ( pos.first )->second;
        for (multimap < mEntity *,int >::iterator it = pos.first; it != pos.second; it++)
            master_color_per_color[it->second] = master_component;
    }

    // ----------------------------------
    //          graph completed
    // ----------------------------------
    // export
    if (debug)
    {
        cout << endl << "--------------------------------------------------------------" << endl;
        cout << "--------------------------- GRAPH INFOS ----------------------" << endl;
        cout << "--------------------------------------------------------------" << endl;
        cout << " keys --- component --- master component" << endl;
        vector < int >::iterator itcolor = graph_color.begin();
        for (vector < entity_id >::iterator it = dentities_tobecolored.begin(); it != dentities_tobecolored.end(); it++,itcolor++)
        {
            mEntity *face = it->second;
            cout << "entity " << it->first->getId() << " of level " << it->first->getLevel() << " made of (";
            if (it->first->getLevel() != 0)
                for (int k = 0; k < it->first->size(0); k++)
                    cout << it->first->get(0,k)->getId() << ",";
            else
                cout << it->first->getId();
            cout << ") on element " << face->getId() << " made of (";
            for (int k = 0; k < face->size(0); k++)
                cout << face->get(0,k)->getId() << ",";
            cout << ") --- color " << ( *itcolor ) << " --- master color " << master_color_per_color[( *itcolor )] << endl;
        }
        cout << "--------------------------------------------------------------" << endl << endl;
    }


    // ----------------------------------
    // now, creating the keys components, including master keys and slave keys, for all spaces
    // the following becomes specific for ramp heaviside
    // a key is considered if and only if:
    // - its entity is entirely in the positive ls zone
    // - the number of attached components is >= 2
    // (this means its support is cut at least twice and is in the positive ls zone...)
    // ----------------------------------

    map < int,int > master_component_per_component;
    map < int,int > slave_component_per_component;

    map < int,xValKey > master_key_per_component_temp;
    for (int icomp = 0; icomp < nb_dentity_components; icomp++)
    {
        const int color = icomp;
        if (debug)
            cout << endl << "working on color " << color << "/" << nb_dentity_components << endl;
        pair < multimap < int,entity_id >::iterator, multimap < int,entity_id >::iterator > pos
	  = color_to_dentities.equal_range(color);

        mEntity *e = pos.first->second.first;
        if (debug)
        {
            cout << "working on mEntity " << e->getId() << " of level " << e->getLevel() << " made of ";
            if (e->getLevel() != 0)
                for (int k = 0; k < e->size(0); k++)
                    cout << e->get(0,k)->getId() << " ";
            else
                cout << e->getId();
            cout << endl;

            cout << " number of components attched to this mEntity: " << number_of_colors_per_entity[e] << endl;
        }

        // if strictly in positive ls zone
        // then, enrichement needed
        if (front.noTouchingOut(e))
        {
            if (number_of_colors_per_entity[e] >= 2)
            {

                xField::const_iterator itspace = field.begin();
                xField::const_iterator itenspace = field.end();
                int spacenumber = 0;
                for (; itspace != itenspace; itspace++,spacenumber++)
                {
                    for (multimap < int,entity_id >::iterator it = pos.first; it != pos.second; it++)
                    {
                        vector < xValKey > thekeys;
                        mEntity *element = ( *it ).second.second;

                        if (debug)
                        {
                            cout << "working on element " << element->getId() << " of level " << element->getLevel() << " made of ";
                            for (int k = 0; k < element->size(0); k++)
                                cout << element->get(0,k)->getId() << " ";
                            cout << endl;
                        }

                        vector < xValKey > keyvec;
                        ( *itspace )->getKeys(element,&keyvec);
                        vector < xValKey >::iterator itvec = keyvec.begin();
                        for (; itvec != keyvec.end(); itvec++)
                        {
                            mEntity * ref, *elm;
                            int n;
                            xSpacePolynomialDiscontinuous::getEntitiesFromKey(&elm,&ref,&n,( *itvec ));
                            if (elm == e)
                            {
                                thekeys.push_back(( *itvec ));
                                if (debug)
                                    cout << "storing key " << thekeys.back() << endl;
                            }
                        }

                        // storing the keys and giving a component number
                        int subspacenumber = 0;
                        for (vector < xValKey >::iterator itvec = thekeys.begin(); itvec != thekeys.end(); itvec++,subspacenumber++)
                        {
                            keys.push_back(( *itvec ));
                            component.push_back(color + nb_dentity_components*subspacenumber + nb_dentity_components*thekeys.size()*spacenumber);
                            master_component_per_component[color + nb_dentity_components*subspacenumber + nb_dentity_components*thekeys.size()*spacenumber] = master_color_per_color[color] +  nb_dentity_components*subspacenumber + nb_dentity_components *thekeys.size()*spacenumber;
                            master_key_per_component_temp[component.back()] = keys.back(); // the last inserted key is by default the master of its component...
                        }
                    }

                }                // end space iterator
            }            // end enrichment needed
        }
    }


    // each group of components has a master... and also a slave:
    int mastercompo,compo;
    multimap < int,int > master_component_to_all_surrounding_components_temp;
    for (map < int,int >::iterator it = master_component_per_component.begin(); it != master_component_per_component.end(); it++) // for all components
    {
        compo = it->first;
        mastercompo = it->second;
        master_component_to_all_surrounding_components_temp.insert(make_pair(mastercompo,compo));
    }
    for (map < int,int >::iterator it = master_component_per_component.begin(); it != master_component_per_component.end(); it++)
    {
        compo = it->first;
        mastercompo = it->second;
        const pair < multimap < int,int >::iterator, multimap < int,int >::iterator > range_of_slaves = master_component_to_all_surrounding_components_temp.equal_range(mastercompo);
        for (multimap < int,int >::iterator itslaves = range_of_slaves.first; itslaves != range_of_slaves.second; itslaves++)  // the first surrounding component which is not the master is declared slave
        {
            int currentcompo = itslaves->second;
            if (currentcompo != mastercompo)
            {
                slave_component_per_component.insert(make_pair(compo,currentcompo));
                break;
            }
        }
    }


    // export
    if (debug)
    {
        cout << endl << "--------------------------------------------------------------" << endl;
        cout << "--------------------------- LINK INFO -----------------------------" << endl;
        cout << "--------------------------------------------------------------" << endl;
        cout << " keys --- component --- master " << endl;
        vector < int >::iterator itc = component.begin();
        for (vector < xValKey >::iterator it = keys.begin(); it != keys.end(); it++,itc++)
        {
            mEntity * face, *ref;
            int n;
            xSpacePolynomialDiscontinuous::getEntitiesFromKey(&face,&ref,&n,( *it ));
            cout << ( *it ) <<  " made of ";
            for (int k = 0; k < ref->size(0); k++)
                cout << ref->get(0,k)->getId() << ",";
                cout << " target entity "<<face->getId()<<" --- component " << ( *itc ) << " --- master component " << master_component_per_component[( *itc )] << " --- slave component " << slave_component_per_component[( *itc )]  << endl;
        }
        cout << "--------------------------------------------------------------" << endl << endl;
    }

    // finally, regrouping final informations and computing link coefficients...
    for (vector < int >::iterator it = component.begin(); it != component.end(); it++)
    {
        compo = *it;
        if (compo == slave_component_per_component[compo])       // the component is the slave of its group: link to master of group
        {
            master_key_per_component.insert(make_pair(compo, master_key_per_component_temp[master_component_per_component[compo]]));
            coefficient_per_component.insert(make_pair(compo,-1.));
        }
        else         // just link to its component master...
        {
            master_key_per_component.insert(make_pair(compo, master_key_per_component_temp[compo]));
            coefficient_per_component.insert(make_pair(compo,1.));
        }
    }
}

//----------------------------------------------------------------------------------------------------------------

xValue < double > * xValueCreatorRampHeavisideBase::basicCreator(const xValKey & key)
{
    const bool debug = false;

    // first check if the key exists in the list
    // if not, nothing has to be done. But since the operator() has to do someting (return an xvalue),
    // an xValue is created with zero value...
    vector < xValKey >::iterator itfind = find(keys.begin(),keys.end(),key);
    if (itfind == keys.end())
    {
        xValueDouble *val = new xValueDouble;
        val->setVal(0.0);
        val->setState(new xStateOfValueFixed(val));
        return val;
    }

    // debug /// debug
#if 0
    //if ( key.getEnti()->getId()!=42 || key.getPhys()!=4)
    //if ( key.getEnti()->getId()!=170 || key.getPhys()!=4)
    //if ( key.getEnti()->getId()!=42 )
    if (  ( key.getEnti()->getLevel() > 0 && ( key.getEnti()->get(0,0)->getId() != 79  || key.getEnti()->get(0,1)->getId() != 35 ) )
          || ( key.getEnti()->getLevel() < 1 )
          )
    {
        xValueDouble *val = new xValueDouble;
        val->setVal(0.0);
        val->setState(new xStateOfValueFixed(val));
        return val;
    }
#endif
    // debug /// debug

    // recover its component number
    int compo = component[itfind-keys.begin()];

    if (key == master_key_per_component[compo])     // this key is a master key: create value
    {
        if (debug)
        {
            cout << "creating new value for master key " << key << endl;
            mEntity * face, *ref;
            int n;
            xSpacePolynomialDiscontinuous::getEntitiesFromKey(&face,&ref,&n,( *const_cast<xValKey*>(&key) ));
            cout << "    on Element " << ref->getId() << " made of ";
            for (int k = 0; k < ref->size(0); k++)
                cout << ref->get(0,k)->getId() << ",";
            cout << endl;
        }
        xValueCreator < xValueDouble > creator_reg;
        return creator_reg(key);
    }

    xValue < double > *master;
    if ( ( master = double_manager->find(master_key_per_component[compo])) )    // master already created: link value to master with appropriate coef
    {
        if (debug)
        {
            cout << "linking the key " << key << " to its master key " << master_key_per_component[compo] << " with coef " << coefficient_per_component[compo] << endl;
            mEntity * face, *ref;
            int n;
            xSpacePolynomialDiscontinuous::getEntitiesFromKey(&face,&ref,&n,( *const_cast<xValKey*>(&key) ));
            cout << "    on Element " << ref->getId() << " made of ";
            for (int k = 0; k < ref->size(0); k++)
                cout << ref->get(0,k)->getId() << ",";
            cout << endl;
        }
        return new xValueLinearCombination(coefficient_per_component[compo], master);
    }

    // else: wait
    return 0;

}

//----------------------------------------------------------------------------------------------------------------

void xValueCreatorRampHeavisideBase::add_to_graph(const vector < entity_id > & temp)
{
    // insert the elements of temp in the main vector if not inserted yet
    vector < entity_id >::const_iterator tempit = temp.begin();
    vector < int > position;
    for (; tempit != temp.end(); tempit++)
    {
        vector < entity_id >::iterator itfind = find(dentities_tobecolored.begin(),dentities_tobecolored.end(),( *tempit ));
        if (itfind == dentities_tobecolored.end())
        {
            dentities_tobecolored.push_back(( *tempit ));
            position.push_back(dentities_tobecolored.size()-1);
        }
        else
            position.push_back(itfind-dentities_tobecolored.begin());            // key position...
    }

    const int n = position.size();
    for (int ipos = 1; ipos < n; ipos++)
    {
        add_edge(position[0], position[ipos], G);
    }
}

//----------------------------------------------------------------------------------------------------------------

xValueCreatorRampHeaviside::xValueCreatorRampHeaviside(xMesh* m, const xcut::xPhysSurfByTagging &front, const xField &field) : xValueCreatorRampHeavisideBase(m,front,field)
{
}

//----------------------------------------------------------------------------------------------------------------

xValue < double > * xValueCreatorRampHeaviside::operator()(const xValKey & key)
{
    return basicCreator(key);
}

//----------------------------------------------------------------------------------------------------------------

xValueUpdatorRampHeaviside::xValueUpdatorRampHeaviside( xMesh* m, const xcut::xPhysSurfByTagging &front, const xField &field, const std::string & sub_,const int physid1, const int physid2, const int physid3) : xValueCreatorRampHeavisideBase(m,front,field),sub(sub_)
{
    // loop on double_manager to :
    //    fill the set of old dicontinus key
    //    remove from subset all old discontinus keys
    //    replace all lineare combination by independante value (with correct value)
    //    remove state for values associeted to any kind of keys
    //
    xValueManagerDist<double>::map_iterator it = double_manager->begin();
    xValueManagerDist<double>::map_iterator itend = double_manager->end();
    for (; it != itend; ++it)
    {
        const xValKey key = it->first;
        const int phys = key.getPhys();
        // if it is a discontinuous key
        if (phys == physid1 || phys == physid2 || phys == physid3)
        {
            // add key to the groupe of old key
            old_disco_keys.insert(key);

            // clean subset
            double_manager->erase(it->second,sub);

            // replace linear combination by it's value
            xValueLinearCombination * linear_combi = dynamic_cast < xValueLinearCombination * > ( it->second );
            if (linear_combi)
            {
                double val = linear_combi->getVal();
                delete linear_combi;
                it->second = new xValueDouble;
                xValueDouble * new_val = dynamic_cast < xValueDouble * > ( it->second );
                if (new_val)
                    new_val->setVal(val);
                else
                    throw 1;
            }
            else
            {
                // reset state for dofs
                xValueDouble * dof = dynamic_cast < xValueDouble * > ( it->second );
                if (dof)
                    dof->delState();
            }
        }
        // other field
        // here in future there may be a reset strategie of state of
        // other field. Why ?
        // This is for the case of mixed dof : continus and discontinus
        // are presentent in subset in random order. When removing discontinus dof from
        // subset in above loop it change dof numbering as subset size is
        // decremented. The risk is to have then 2 value with the same dof number at creation.
        // To remove this problem  reset here the state of dof to force
        // renumbering of the wall subset is a solution.
        // The other solution is to prevent this mixed order and create first state value for other
        // filds and at last creates states for discontinus field. Now watever append above the
        // numbering of the subset will be ok and state of other field can be enchanged.
    }

}

//----------------------------------------------------------------------------------------------------------------
xValue < double > * xValueUpdatorRampHeaviside::operator()(const xValKey & key)
{
    // store all keys of the new discontinuous field
    new_disco_keys.insert(key);

    // first check if the key exists already in the double manager
    xValue < double > *exist = double_manager->find(key);

    // if it exists: update procedure
    if (exist)
    {
        // first check if the key exists in the list of keys to enrich
        // if not, nothing has to be done. But since the operator() has to do someting (return an xvalue),
        // an xValue is created with zero value in the standard procedure
        // here: returning the existing value
        vector < xValKey >::iterator itfind = find(keys.begin(),keys.end(),key);
        if (itfind == keys.end())
        {
            exist->setVal(0.0);
            exist->setState(new xStateOfValueFixed(exist));
            return exist;
        }

        // recover its component number
        int compo = component[itfind-keys.begin()];

        if (key == master_key_per_component[compo])
        {
            // this key is a master key: create value in standard procedure
            // here: return existing...
            return exist;
        }

        xValue < double > *master;
        if ( ( master = double_manager->find(master_key_per_component[compo])) )
        {
            // master already created: link value to master with appropriate coef
            // ... and check if master has the right value... or if it's new !
            double masterval = master->getVal();
            if (masterval == (double) 0.)
            {
                master->setVal(exist->getVal()*coefficient_per_component[compo]); // this might be the wrong coef, since the coef might have changed between the two fields... whatever...
            }
            return new xValueLinearCombination(coefficient_per_component[compo], master);
        }

        // else: wait
        return 0;
    }
    // it doesn't exist : back to standard creation procedure
    else
        return basicCreator(key);
}

//----------------------------------------------------------------------------------------------------------------

void xValueUpdatorRampHeaviside::postUpdate(void)
{

    // find old keys not used
    vector < xValKey > diff_disco_keys(old_disco_keys.size());
    vector < xValKey >::iterator it = diff_disco_keys.begin();
    vector < xValKey >::iterator itend = set_difference (old_disco_keys.begin(), old_disco_keys.end(),
                                                         new_disco_keys.begin(), new_disco_keys.end(),
                                                         it);
    //destroy old keys not used from double manager
    for (; it != itend; ++it)
    {
        double_manager->erase(*it);
    }
}



} // end of namespace










