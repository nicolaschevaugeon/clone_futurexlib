/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <cstring>
#include <cstdio>
#include <cmath>		
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>

#include "mIterator.h"
#include "mEdge.h"
#include "mFace.h"
#include "xMesh.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mTet.h"
#include "mAttachableDataContainer.h"

#include "xLevelSet.h"
#include "xLevelSetOperators.h"
#include "xElement.h"
#include "xRegion.h"
#include "xDebug.h"
#include "xIntegrationRule.h"
#include "xExportEnsight.h"

// to be removed when fully parallel
#include "workInProgress.h"

#ifdef PARALLEL
#include "ParUtil.h"
using std::ostringstream;
using AOMD::ParUtil;
#endif

using AOMD::mEntity;
using AOMD::mVertex;

namespace xfem
{

xLevelSet::xLevelSet(void) : grad_up_to_date(false), curv_up_to_date(false)  { }


xLevelSet::xLevelSet(const xRegion& s, const double& val) : 
      grad_up_to_date(false), curv_up_to_date(false), true_curv_up_to_date(false), support(s) 
{
  for(xIter it = support.begin(0); it != support.end(0); ++it) {
    mVertex *v = (mVertex *) *it;
    ls[v] = val;
  }
}

xLevelSet::xLevelSet(const xRegion& s, const xPointToDouble&  func) : 
      grad_up_to_date(false), curv_up_to_date(false), true_curv_up_to_date(false), support(s) 
{
  load(func);
}

  void xLevelSet::load(const xPointToDouble&  func)
  {
    for(xIter it = support.begin(0); it != support.end(0); ++it ) {
      AOMD::mVertex *v = (AOMD::mVertex *) *it;
      ls[v] = func(v->point());
    }
  }

  void  xLevelSet::loadOnElement(const xPointToDouble& func, AOMD::mEntity* e)
  {
    for(int i=0; i<e->size(0); ++i)
      {
	AOMD::mVertex* v = (AOMD::mVertex*) e->get(0,i);
	ls[v] = func(v->point());
      }
  }

  
xLevelSet & xLevelSet::operator = (const xLevelSet &other){
    support=other.support;

    for(xIter it = support.begin(0); it != support.end(0); ++it) {
      mVertex *v = (mVertex *) *it;
     ls[v] = other(v);
    }
    
    grad_up_to_date=false;
    curv_up_to_date=false;
    true_curv_up_to_date=false;
    return *this;
}


void xLevelSet::complement()
{
  for(xIter it = support.begin(0); it != support.end(0); ++it ) {
    AOMD::mVertex *v = (AOMD::mVertex *) *it;
    ls[v] *= -1.0 ;
  }
}

double& xLevelSet::operator () (mEntity* v) {
  // cout << "arg" ; 
  grad_up_to_date =  false; 
  curv_up_to_date = false; 
  true_curv_up_to_date = false; 
  return ls[v];
}


bool xLevelSet::isDefinedAt(mEntity* v) const { return (ls.find(v) != ls.end());}


bool xLevelSet::isDefinedOnElement(mEntity* e) const 
  { 
    for (int i = 0; i < e->size(0); ++i)
      {
	mEntity* v = e->get(0,i);
	if (ls.find(v) == ls.end()) return false;
      }
    return true;
  }

int xLevelSet::side_of(const xGeomElem* geo_appro, const xGeomElem* geo_integ) const{
  const bool debug=false;
  unsigned int side_tag;
  side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
  int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);
  
  if (orientation==0)
    {
      if(debug) cout<<"in sideof : in domain\n";
      mEntity* e = geo_appro->getEntity();
      //mEntity* e_sub; 
      xElement elem(e);
      elem.setUvw(geo_appro->getUVW());
      std::vector<double> vals = getVals(e); 
      double f = elem.getInterpoSca(vals);
      if (fabs(f)<1e-3)
        {
          if (geo_integ!=NULL)
            {
              //e_sub = geo_integ->getEntity();
              elem.xyz2uvw(geo_integ->getCDGxyz());
              
            }
          else {
            //e_sub = geo_appro->getEntity();
            elem.setUvw(geo_appro->getCDGuvw());
          }
          
          std::vector<double> vals = getVals(e); 
          f = elem.getInterpoSca(vals);
        }
      
      if (f>=0) 
        return (1);
      else 
        return (-1); 
    }
  else
    {
      if(debug)cout<<"in sideof : on slice\n";
      if(debug)
        {
          cout << "Entity in sideof :";
          geo_integ->getEntity()->print();
          cout<<endl;
        }
      if(debug) cout<<"Orientation in sideof\n"<<orientation<<endl;
      if (orientation>=0) 
        return (1);
      else 
        return (-1); 
    }
}


const double& xLevelSet::operator () (mEntity* v)  const {
  if (ls.find(v) == ls.end()) throw;
  assert(ls.find(v)!=ls.end());
  return ls.find(v)->second;
}


const xRegion& xLevelSet::getSupport() const {return support;}


void xLevelSet::setSupport(const xRegion& m, const double& val) 
{
  support = m; 
  clear();
  for(xIter it = support.begin(0); it != support.end(0); ++it) 
  {
    mVertex *v = (mVertex *) *it;
    ls[v] = val;
  }
}

void xLevelSet::reduceSupport(const xRegion& m) 
{
  support = m; 
}


void xLevelSet::clear() {
  ls.clear();
  grad.clear();
  grad_up_to_date = false;
  curv.clear();
  curv_up_to_date = false;
  true_curv.clear();
  true_curv_up_to_date = false;
}




void xLevelSet::accept(xLevelSetModifier& visitor) { 
    visitor.visit(*this, getSupport());
  }
void xLevelSet::accept(xLevelSetModifier& visitor, const xRegion& target) { 
  visitor.visit(*this, target);
}


void xLevelSet::accept(xLevelSetInspector& visitor) const {
  visitor.visit(*this, getSupport());
}
void xLevelSet::accept(xLevelSetInspector& visitor, const xRegion& target) const {
  visitor.visit(*this, target);
}


void xLevelSet::accept(xLevelSetCreator& visitor, xLevelSet& f) const {
  visitor.visit(*this, f);
}


std::vector<double> xLevelSet::getVals(mEntity* e) const {
  const bool debug = xdebug_flag;
  std::vector<double> vals; 
  vals.reserve(e->size(0));
  if (debug) 
    {
      cout << " e in getvals " << endl;
      e->print();
    }
  if (e->getLevel() > 0)
    {
      for (int i = 0; i < e->size(0); i++) {
	mVertex *v = (mVertex*) e->get(0,i);
	lstype::const_iterator itend = ls.end();
	lstype::const_iterator it = ls.find(v);
	if (it!=itend)
	  vals.push_back(it->second);
	else{
	  std::cerr << "can't find vertex "<< v << "in level set " << this;
	  throw 1;
	}
	  
      }
    }
  else {
    lstype::const_iterator itend = ls.end();
    lstype::const_iterator it = ls.find(e);
    if (it!=itend)
      vals.push_back(it->second);
    else{
      std::cerr << "can't find vertex "<< e << "in level set " << this;
      throw 1 ;
    }
  }
  if (debug) 
    {
      cout << " ls values at the nodes " << endl;
      cout << " number of nodes " << vals.size() << endl;
      std::copy(vals.begin(), vals.end(), std::ostream_iterator<double>(cout, " ")); 
      cout << endl; 
    }
  return vals;
}


// get values of the level set for the nodes given in argument
// v is a vector of pointer to mEntity entity
std::vector<double> xLevelSet::getVals(std::vector<mEntity*> &v) const {
  std::vector<double> vals; 
  std::vector<mEntity*>::iterator itv = v.begin();
  std::vector<mEntity*>::const_iterator itve = v.end(); 
  lstype::const_iterator itl ;
  lstype::const_iterator itle = ls.end();

  for (; itv !=  itve; ++itv) 
  {
     itl = ls.find(*itv);
     if (itl!=itle)
         vals.push_back(itl->second);
     else
     {
         mVertex * itvv = (mVertex *) *itv;
         std::cerr << "can't find vertex "<< itvv << "in level set " << this;
         throw 1;
     }
  }
  return vals;
}




xtensor::xVector xLevelSet::getGrad(mEntity* e)  const {
//  ... cette fonction calcule le gradient du xLevelSet sur tous les noeuds du maillage
  const bool debug = xdebug_flag;
  if (debug) cout << " e in getgrad is " << endl;
  
  if (debug) e->print();
  if (debug) cout << " bool grad_up_to_date is " << grad_up_to_date << endl;
  if (!grad_up_to_date) compute_grad();
  auto it = grad.find(e);
  if (it == grad.end()) throw 1;
  xtensor::xVector out(grad.find(e)->second);
  return out;


}


double xLevelSet::getVal(mEntity* e, const Trellis_Util::mPoint& uvw)  const {
//  ... cette fonction calcule EN UN POINT D UN ELEMENT la valeur du champ du xLevelSet
//      a partir des fonctions de forme et de la valeur aux noeuds de l element 
//      de ce meme xLevelSet
  //const bool debug = xdebug_flag;  const int nnodes = e->size(0);
  std::vector<double> vals=getVals(e);
  
  switch  (e->getType()){
  case mEntity::VERTEX :{
    return vals[0];
    break;
  }
  case mEntity::EDGE:{
    const double u = uvw(0);
    return vals[0]*(0.5*(1.-u)) +  vals[1]*(0.5*(1.+u));
    break;
  }
  case mEntity::TRI:{
    const double u = uvw(0);
    const double v = uvw(1);
    return vals[0]*(1.-u-v) +  vals[1]*u +  vals[2]*v;
    break;
  }
  case mEntity::TET:{
    const double u = uvw(0);
    const double v = uvw(1);
    const double w = uvw(2);
    return vals[0]*(1.-u-v-w) +  vals[1]*u +  vals[2]*v + vals[3]*w ;
    break;
  } 
  case mEntity::QUAD:{
    xElement elem(e);
    elem.setUvw(uvw); 
    double out=elem.getInterpoSca(vals);
    return out;
    break;
  }
  case mEntity::HEX:{
    xElement elem(e);
    elem.setUvw(uvw);
    double out =elem.getInterpoSca(vals);
    return out;
    break;
  }  
  default :throw;
  }
}

bool xLevelSet::getVal (mEntity* e, const Trellis_Util::mPoint& uvw, double &val) const{
  try{
    val = getVal( e, uvw);
    return true;
  }
  catch(...){
    val =0.;
    return false;
  }
}

xtensor::xVector xLevelSet::getGrad(mEntity* e, const Trellis_Util::mPoint& uvw)  const {
//  ... cette fonction calcule EN UN POINT D UN ELEMENT la valeur du gradient du xLevelSet
//      a partir des fonctions de forme et de la valeur aux noeuds de l element du gradient
//      de ce meme xLevelSet
  if (!grad_up_to_date) compute_grad();
  // xElement elem(e);
//  ... determination des coordonnees locales du noeud considere
  //elem.xyz2uvw(p);
  //elem.setUvw(uvw);
//  get des gradients aux noeuds
  const int nnodes = e->size(0);
  std::vector<xtensor::xVector> grad_v(nnodes);
  for(int j=0;j< nnodes;j++)
    {
      grad_v[j] =  getGrad( e->get(0,j));
    }
  switch  (e->getType()){
  case mEntity::VERTEX :{
    return grad_v[0];
    break;
  }
  case mEntity::EDGE:{
    const double u = uvw(0);
    return grad_v[0]*(0.5*(1.-u)) +  grad_v[1]*(0.5*(1.+u));
    break;
  }
  case mEntity::TRI:{
    const double u = uvw(0);
    const double v = uvw(1);
    return grad_v[0]*(1.-u-v) +  grad_v[1]*u +  grad_v[2]*v;
    break;
  }
  case mEntity::TET:{
    const double u = uvw(0);
    const double v = uvw(1);
    const double w = uvw(2);
    return grad_v[0]*(1.-u-v-w) +  grad_v[1]*u +  grad_v[2]*v + grad_v[3]*w ;
    break;
  } 
  case mEntity::QUAD:{
    xElement elem(e);
    //  ... determination des coordonnees locales du noeud considere
    elem.setUvw(uvw);
    xtensor::xVector out(elem.getInterpoVec(grad_v));
    return out;
    break;
  }
  case mEntity::HEX:{
    xElement elem(e);
    //  ... determination des coordonnees locales du noeud considere
    elem.setUvw(uvw);
    xtensor::xVector out(elem.getInterpoVec(grad_v));
    return out;
    break;
  }  
  default :throw;
  }
  // The function should have return or throw before this line
  return xfem::xVector();
}

bool xLevelSet::getGrad(mEntity* e, const Trellis_Util::mPoint& uvw, xtensor::xVector & vals)  const {
  try{
    vals = getGrad(e, uvw);
    return true;
  }
  catch(...){
    return false;
  } 
}


xtensor::xTensor2 xLevelSet::getCurv(mEntity* e) const {
  if (!curv_up_to_date) compute_curv();
  return(curv.find(e)->second);
}

xtensor::xTensor2 xLevelSet::getTrueCurv(mEntity* e) const {
  if (!true_curv_up_to_date) compute_true_curv();
  return(true_curv.find(e)->second);
}


//  rentre des valeurs dans grad_region, grad_vertex
void xLevelSet::compute_grad(void) const {
  const bool debug = xdebug_flag;
  if (debug) cout << " support in compute_grad is of size dim=3 " << support.size(3) << endl; 
  for(xIter it = support.begin(); it != support.end(); ++it) {
    mEntity *e = *it;
    if (debug) cout << " compute_grad for element " << endl;
    if (debug) e->print();
    std::vector<double> ls_v = getVals(e);
    if (debug) cout << " nodal level set values " << endl;
    if (debug) std::copy(ls_v.begin(), ls_v.end(), std::ostream_iterator<double>(cout, " ")); 
    xElement elem(e);
    grad[e] = elem.getGradInterpoSca(ls_v);
  }
  for(xIter it = support.begin(0); it != support.end(0); ++it) {
    mEntity *v =  *it;
    grad[v]    = xtensor::xVector(0.,0.,0.);
  }

  eXlibris_types::hash_map<AOMD::mEntity*,int,AOMD::EntityHashKey,AOMD::EntityEqualKey> nbconnect;
 

  for(xIter it = support.begin(); it != support.end(); ++it) {
  //for(mIteratorPtr it(support.getTopIter()); !it->done(); it->next() ) {
    mEntity *e = *it;
    xtensor::xVector g = grad[e];
    for (int j = 0; j < e->size(0); ++j) {
      mEntity* v = e->get(0,j);
      grad[v] += grad[e];
      nbconnect[v] += 1;
    }
  }

  // for vertex on boundary comunication must be done to have the full contribution of all elements to the
  // grad of the vertex : sum grad[v] and nbconnect[v] (reduce) => grad[v]/nbconect[v] is correct
  assert(xtool::workInProgress());

  for(xIter it = support.begin(0); it != support.end(0); ++it) {
  //for(mIteratorPtr it(support.getIter(0)); !it->done(); it->next() ) {
    mEntity *v =  *it;
    grad[v]    /= (double) nbconnect[v];
  }
  grad_up_to_date = true;
  return;
}

//  rentre des valeurs dans grad
void xLevelSet::compute_curv(void) const {
//let grad be the field gradient
//curv is the matrix [grad,1  grad,2 grad,3]
  if (!grad_up_to_date) compute_grad();
  xtensor::xTensor2 ten;
  for(xIter it = support.begin(); it != support.end(); ++it) {
  //for(mIteratorPtr it(support.getTopIter()); !it->done(); it->next() ) {
//  ... determination de l element courant
    mEntity *e = *it;
    std::vector<xtensor::xVector> grad_v;
    for (int j = 0; j < e->size(0); j++) {
      mVertex* v = (mVertex*) e->get(0,j);
      grad_v.push_back(getGrad(v));
    }
    xElement elem(e);
    elem.getGradInterpoVec(grad_v, ten);
//curv must be symmetric it would be good to check
    curv[e] = ten;
  }
  curv_up_to_date = true;
  return;
}

void xLevelSet::compute_true_curv(void) const {
  if (!grad_up_to_date) compute_grad();
  xtensor::xTensor2 ten;
  for(xIter it = support.begin(); it != support.end(); ++it) {
    mEntity *e = *it;
    std::vector<xtensor::xVector> grad_v;
    for (int j = 0; j < e->size(0); j++) {
      mVertex* v = (mVertex*) e->get(0,j);
      grad_v.push_back(getGrad(v).norm());
    }
    xElement elem(e);
    elem.getGradInterpoVec(grad_v, ten);
    true_curv[e] = ten;
  }
  true_curv_up_to_date = true;
  return;
}







void xLevelSet::restrictTo(xRegion new_support, double init_val) 
{ 
   for(xIter it = new_support.begin(0); it != new_support.end(0); ++it) {
   //for(mIteratorPtr it(new_support.getIter(0)); !it->done(); it->next() ) {
    mEntity* v = *it;
    if (!isDefinedAt(v)) ls[v] = init_val;
  }
   for(xIter it = support.begin(0); it != support.end(0); ++it) {
  //for(mIteratorPtr it(support.getIter(0)); !it->done(); it->next() ) {
    mEntity* v = *it;
    if (!new_support.find(v)) { ls.erase(v); }
  }
  support = new_support;
}
       


void xLevelSet::printDebug(void) const
{  
  for(xIter it = support.begin(0); it != support.end(0); ++it) {
  //for(mIteratorPtr it(support.getIter(0)); !it->done(); it->next() ) {
    mVertex *v = (mVertex *) *it;
    Trellis_Util::mPoint p = v->point();
    cout <<"At Point :"<< p << " Value ls is :" << ls.find(v)->second << endl;
  }
}


void xLevelSet::exportMatlab(ostream& fout, const std::string& field_name, 
			  int level) const
{
  fout <<  field_name << "= [ ";
  for(xIter it = support.begin(level); it != support.end(level); ++it) 
  //for(mIteratorPtr it(support.getIter(level)); !it->done() ; it->next())
    {
      mEntity *e = *it;
      fout << ls.find(e)->second << " ";
    }
  fout << "];\n";
  return;
}

void xLevelSet::takeTraceOn(xMesh* mesh, xLevelSet& trace) const
{ 
  const bool debug = xdebug_flag;
  if (debug) cout<<"takeTraceOn : entree\n\n";
  trace.setSupport(mesh);
  std::vector<double> lsv(2);
  //xRegion m_interface = trace.getSupport();
  for(xIter it = mesh->begin(0); it != mesh->end(0); ++it) 
    {
      mVertex *v = (mVertex*) *it;
      mEntity *e = v->getAttachedEntity(xMesh::was_created_by_tag);
      assert(e != 0);
      if (e->getLevel() == 0)  { trace(v) = (*this)((mVertex*) e);}
      else if (e->getLevel() == 1) {
        double r = v->getAttachedDouble(xMesh::r_on_edge_tag);
	//we need to interpolate
	for(int j=0;j<2;j++) {
	  mVertex *ends = (mVertex*) e->get(0,j);	
          lsv[j] = (*this)(ends);
	}
        trace(v) = lsv[0] + r * (lsv[1] - lsv[0]);
      }
      else assert(0);
    }
  if (debug) cout<<"takeTraceOn : sortie\n\n";
  return;
}

/// enables to evaluate x_field on the nodes of mesh_new
/// where mesh_new is necessarily a submesh
void xLevelSet::interpolateTo(xMesh* mesh_new, xLevelSet& lsnew) const
{
  const bool debug = xdebug_flag;
  if (debug) cout<<"interpolateTo : ENTREE\n\n";
  lsnew.setSupport(mesh_new);
  //int dim = mesh_new->dim();
  std::vector<double> lsv(2);
  for(xIter it = mesh_new->begin(0); it!= mesh_new->end(0); ++it)
    {
      mVertex *v = (mVertex*) *it;
      if (debug) cout<<"noeud courant du submesh "<<v->point()<<"****\n";
      if (debug) cout<<" id is " << v->getId() << endl;
      mVertex *vd = (mVertex*) v->getAttachedEntity(xMesh::get_was_duplicated_from_tag());//n
      if(vd == 0) { assert(0); abort(); } 
      if (debug) cout<<"noeud duplicateur "<<vd->point()<<"****\n";
      if (!vd->getAttachedEntity(xMesh::was_created_by_tag))
	{ lsnew(v) = (*this)((mVertex*) vd);
	  if (debug) cout<<"valeur de lsnew en "<<v->point()<<"est "<<lsnew(v)<<"****\n";
	}
      else
	{	      
	  mEntity *e = (mEntity*)vd->getAttachedEntity(xMesh::was_created_by_tag);
	  if (debug){
			cout <<"level de e = "<<e->getLevel()<<"*****\n";
			e->print();
		}
	  if (e->getLevel()==0) { lsnew(v) = (*this)((mVertex*) e);}
	  else //e is then an edge
	    {
	      double r = vd->getAttachedDouble(xMesh::r_on_edge_tag);
	      //we need to interpolate
	      for(int j=0;j<2;j++) {
		mVertex *ends = (mVertex*) e->get(0,j);	
		lsv[j] = (*this)(ends);
	      }
	      lsnew(v) = lsv[0] + r * (lsv[1] - lsv[0]);
	      if (debug) cout<<"valeur de r en  "<<v->point()<<" est "<<r<<"****\n";
	      if (debug) cout<<"valeur de lsv[0] en  "<<v->point()<<" est "<<lsv[0]<<"****\n";
	      if (debug) cout<<"valeur de lsv[1] en  "<<v->point()<<" est "<<lsv[1]<<"****\n";
	      if (debug) cout<<"valeur de lsnew en "<<v->point()<<" est "<<lsnew(v)<<"****\n";
	    }
	}
    }
  if(debug) cout<<"interpolateTo : SORTIE\n\n";
  return;
}
} // end of namespace






// void xLevelSet::exportGmsh(const std::string& name) const {
//   exportGmsh(name, support);
// }

// void xLevelSet::exportGmsh(const std::string& name, xRegion sub_support) const
// {
//   xexport::xExportGmshAscii pexport;
//   exportGeneral(pexport, name, sub_support, xtool::xIdentity<double>());
// }
// void xLevelSet::exportEnsight(const std::string& name) const {
//   exportEnsight(name, support);
// }

// void xLevelSet::exportEnsight(const std::string& name, xRegion sub_support) const
// {
//   xexport::xExportEnsightAscii pexport;
//   exportGeneral(pexport, name, sub_support, xtool::xIdentity<double>());
// }


// //greg
// void xLevelSet::exportGradGmsh(const std::string& name) const {
// 	exportGradGmsh(name, support);
// }


// void xLevelSet::exportGradGmsh(const std::string& name, xRegion sub_support) const
// {

// #ifdef PARALLEL
//   std::stringstream oufilename; 
//   oufilename <<"_"<< ParUtil::Instance()->size() << "_" <<  ParUtil::Instance()->rank()+1;
// 	std::string filename = name + oufilename.str() + ".pos" ;
// #endif


// #ifndef PARALLEL
//   std::string filename = name + ".pos";
// #endif

//   FILE* fout = fopen(filename.c_str(), "w");
//   fprintf(fout, "View \"%s\" { \n", name.c_str());
// 	for(xIter it = sub_support.begin(); it != sub_support.end(); ++it)
// 	{
// 		mEntity *e = *it;
// 		std::vector<Trellis_Util::mPoint> p(e->size(0));
// 		xtensor::xVector lev = getGrad(e);
    
    
// 		for(int j=0;j<e->size(0);j++)
// 		{
// 			mVertex *v = (mVertex*)(e->get(0,j));
// 			p[j] = v->point();
// 		}
// 		switch (e->getLevel()) {
// 			case 0: break;
// 			case 1:
// 				fprintf(fout, "VL(%12.5e, %12.5e, %12.5e,", p[0](0), p[0](1), p[0](2));
// 				fprintf(fout, "%12.5e, %12.5e, %12.5e)",    p[1](0), p[1](1), p[1](2));
// 				fprintf(fout, "{%12.5e, %12.5e, %12.5e,%12.5e, %12.5e, %12.5e};\n", lev(0), lev(1), lev(2),lev(0), lev(1), lev(2));
// 				break;
// 			case 2:
// 				fprintf(fout, "VT(%12.5e, %12.5e, %12.5e,", p[0](0), p[0](1), p[0](2));
// 				fprintf(fout, "%12.5e, %12.5e, %12.5e,",    p[1](0), p[1](1), p[1](2));
// 				fprintf(fout, "%12.5e, %12.5e, %12.5e)",    p[2](0), p[2](1), p[2](2));
// 				fprintf(fout, "{%12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e};\n", lev(0), lev(1), lev(2),lev(0), lev(1), lev(2), lev(0), lev(1), lev(2));
// 				break;
// 			case 3:
// 				fprintf(fout, "VS(%12.5e, %12.5e, %12.5e,", p[0](0), p[0](1), p[0](2));
// 				fprintf(fout, "%12.5e, %12.5e, %12.5e,",    p[1](0), p[1](1), p[1](2));
// 				fprintf(fout, "%12.5e, %12.5e, %12.5e,",    p[2](0), p[2](1), p[2](2));
// 				fprintf(fout, "%12.5e, %12.5e, %12.5e)",    p[3](0), p[3](1), p[3](2));
// 				fprintf(fout, "{%12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e};\n", lev(0), lev(1), lev(2),lev(0), lev(1), lev(2), lev(0), lev(1), lev(2), lev(0), lev(1), lev(2));
// 				break;
// 				default: assert(0); break;
// 		}
// 	}
// 	fprintf(fout, "};\n");
// 	fclose(fout);
// 	return;
// }
