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
#include <limits>
#include "xLevelSetOperators.h"
#include "xLevelSet.h"
#include "xElement.h"
#include "mIterator.h"
#include "mVertex.h"
#include "xDebug.h"
#include "xSubMesh.h"
#include "xLapackInterface.h"

#include "workInProgress.h"

using namespace std;
using namespace lalg;
using namespace AOMD;
using namespace Trellis_Util;

namespace xfem
{

void xGetMinMaxLevelSetInspector::visit(const xLevelSet& f, xRegion target) {
  min_ls = f(*target.begin(0));
  max_ls = min_ls;
  for(xIter it = target.begin(0); it != target.end(0); ++it) {
    mVertex *v = (mVertex*) *it;
    min_ls = std::min(f(v), min_ls);
    max_ls = std::max(f(v), max_ls);
  }
}

double xGetMinMaxLevelSetInspector::getMin() const {return min_ls;} 
double xGetMinMaxLevelSetInspector::getMax() const {return max_ls;} 

xScaleLevelSetModifier::xScaleLevelSetModifier(const double& s ) : scale(s) {}
void xScaleLevelSetModifier::visit(xLevelSet& f, xRegion target) {
  for(xIter it = target.begin(0); it != target.end(0); ++it) {
    mVertex *v = (mVertex*) *it;
    f(v) *= scale;
  }
}

xShiftLevelSetModifier::xShiftLevelSetModifier(const double& s ) : shift(s) {}
 
void xShiftLevelSetModifier::visit(xLevelSet& f, xRegion target) {
  for(xIter it = target.begin(0); it != target.end(0); ++it) {
    mVertex *v = (mVertex*) *it;
    f(v) += shift;
  }
}
  

xShiftLevelSetCreator:: xShiftLevelSetCreator(double _shift):shift(_shift){}
  
void  xShiftLevelSetCreator::visit(const xLevelSet& in, xLevelSet& out){
  const xRegion ls_in_region(in.getSupport().getMesh());
  out.setSupport(ls_in_region);
  for(xIter it = ls_in_region.begin(0); it != ls_in_region.end(0); ++it) {
    mVertex *v = (mVertex*) *it;
    out(v) = in(v) + shift;
  }
}
  
xComplementLevelSetCreator::xComplementLevelSetCreator(){}
  
void  xComplementLevelSetCreator::visit(const xLevelSet& in, xLevelSet& out){
    const xRegion ls_in_region(in.getSupport().getMesh());
    out.setSupport(ls_in_region);
    for(xIter it = ls_in_region.begin(0); it != ls_in_region.end(0); ++it) {
      mVertex *v = (mVertex*) *it;
      out(v) = -in(v);
    }
  }

void xUniformValue::visit(xLevelSet& f, xRegion target) {
//  ... boucle sur les noeuds du maillage global
  for(xIter it = target.begin(0); it != target.end(0); ++it) {
    mVertex *v = (mVertex*) *it;
    f(v) = val;
  }
}

void xNonUniformValue::visit(xLevelSet& f, xRegion target) {
//: autre facon d initialiser le champ des vitesses
   for(xIter it = target.begin(0); it != target.end(0); ++it) {
  //for(mIteratorPtr it(target.getIter(0)); !it->done(); it->next() ) {
    mVertex *v = (mVertex*) *it;
    Trellis_Util::mPoint p = v->point();
//  ... remarque:on ne se sert pas de la variable val mais ca pourrait etre le cas
    f(v) = p(0)/2.;
  }
}

//: procedure modifiee
void xVelocityInitialization::visit(xLevelSet& velo, xRegion target) {

   const bool debug =  xdebug_flag; 
//  ... on va chercher le maillage 1D support du champ de vitesses 1D

   if (debug) cout << " Hello in xVelocityInitialization:" << endl;

  xRegion fs = velofront.getSupport();

//  ... boucle sur les noeuds du maillage global
   for(xIter it = target.begin(0); it != target.end(0); ++it) {
    mVertex *ve = (mVertex*) *it;
    velo(ve) = 0.;
  }
  



  if (debug) cout << "starting init velo " << endl; 
//  ... initialisations diverses
  double   dis1,dis2,mem1;
    //  double   d1,d2,d3,d4,x1,x2,x3;
  mEntity* mem2=0;
//  ... boucle sur les noeuds des elements s appuyant sur le front
  for(xIter itloc = regloc.begin(0); itloc != regloc.end(0); ++itloc ) 
  {
    mVertex* v = (mVertex*) *itloc;
    //  ... on a besoin des coordonnees du noeud
    Trellis_Util::mPoint p = v->point();
    //  ... initialisation (TRES LAIDE) de mem1
    mem1 = -1000.;
    dis1 = -1000.;
    dis2 = -1000.;
    if (target.dim()==3)
    {
      //if(ParUtil::Instance()->rank()==talk) cout<<"a2_1\n";
      //  ON VA CHERCHER LE COUPLE DE POINTS LE PLUS PROCHE DU NOEUD CONSIDERE
      //  ... boucle sur les elements du maillage 1D (ensemble de segments)

      for(xIter itt = fs.begin(); itt != fs.end(); ++itt)
      {

	mEntity *se = (mEntity*) *itt;

        //Greg BEG

	if(se->size(0)==0){
	  //if(ParUtil::Instance()->rank()==0) cout<<"size 1=0\n";
	  //cout<<"Who is Here :"<<ParUtil::Instance()->rank()<<endl;
	  se->print();
	  mVertex *v1 = (mVertex*) se;
	  velo(v) = velofront(v1);
	  cout<<"done !\n";
	  return;
	  
	}
        //Greg END



	//  ... on va chercher les deux noeuds du segment
	mVertex* v1 = (mVertex*) se->get(0,0);
	mVertex* v2 = (mVertex*) se->get(0,1);
	//  ... on a besoin des coordonnees des deux noeuds
	Trellis_Util::mPoint p1 = v1->point();
	Trellis_Util::mPoint p2 = v2->point();

	xtensor::xVector P(p1,p);
	xtensor::xVector AB(p1,p2);
	double mag=AB.mag();
	double alpha=(AB*P)/(mag*mag);
	if (alpha<0) alpha=0;
	if (alpha>1) alpha=1;
	xtensor::xVector PP=P-AB*alpha;
	double dist=PP.mag();
	if ((dist<mem1)||(mem1<0.0)) 
	{
	  mem1=dist;
	  mem2=se;
	  dis2=alpha;
	  dis1=1-alpha;
	}

      }
      
      //  ... determination des deux points du segment
      mVertex* v3 = (mVertex*) mem2->get(0,0);
      mVertex* v4 = (mVertex*) mem2->get(0,1);
      assert ((dis1>=0.) && (dis1<=1.));
      velo(v) = (dis1*velofront(v3)+dis2*velofront(v4))/(dis1+dis2);
    }
    else // cas 2D ... point le + proche 
      {
	for(xIter itt = fs.begin(0); itt != fs.end(0); ++itt)
	  {
	    mVertex *v1 = (mVertex*) *itt;
	    Trellis_Util::mPoint p1 = v1->point();
	    double d = sqrt((p1(0)-p(0))*(p1(0)-p(0))+(p1(1)-p(1))*(p1(1)-p(1))+(p1(2)-p(2))*(p1(2)-p(2)));
	    if (((d) < mem1)||(mem1<0.0))
	      {
		mem2=*itt;
		mem1=d;
	      }
	  }
	velo(v) = velofront(mem2);
      }
  }

  
  
}


bool xPilotError::itersCompleted(const xLevelSet& f, xRegion target) {
  //error;
  double ref;
  if (iter == 0) {
    ref=1.;
    error = max(epsabs, epsrel*ref);
  }
  else {
    f.accept(compute_error, target);
    std::pair<double, double> errorref = compute_error.getError();
    error = errorref.first;
    ref   = errorref.second;
    if (error < (std::max (epsabs, epsrel*ref))) return true;
//: pour l instant on se contente de constater la non convergence
//    if (iter > itmax) {fprintf(stderr, "no convergence"); assert(0);}
    if (iter > itmax) {
      std::cerr <<  "no convergence, error is " << error << std::endl;
      return true;
    }
  }

  compute_error.loadReference(f);

  iter++;
  
  return false;
}


void xEvolveToStationary::visit(xLevelSet& f, xRegion target) {
  //  cout << " Inside lEvolveToStationary target is " << endl;
   pilot.initialize();
   while (!pilot.itersCompleted(f, target) ) {
     //cout << " Inside lEvolveToStationary  loop " << endl;
      double dt = pilot.getVirtualTimeStep();
      std::list<xSpaceOperator*>::iterator it;
      //int i=0;
      for (it = operators.begin(); it !=operators.end(); it++) {
	xSpaceOperator* ope = *it;
        xRungeKutta rung(dt, *ope);
        f.accept(rung,target);
      }
     
   }
   // cout << "iterations " << pilot.getIter() << " error is " << pilot.getLastError() << endl;
}


void xTimeIntegration::visit(xLevelSet& f, xRegion target) {
   pilot.initialize();
   while (!pilot.stepsCompleted() ) {
      double dt = pilot.getTimeStep();
      xRungeKutta rung(dt, soperator);
      f.accept(rung, target);
   }
}

bool  xPilotOneStep::stepsCompleted() {
    if (step != 0) return true;
    step++;
    return false;
}

void xRungeKutta::visit(xLevelSet& ls, xRegion target) {
   xLevelSet temp(ls);
   xLevelSet H;
   ls.accept(space_operator, H);
//  ... boucle sur les noeuds
   for(xIter it = target.begin(0); it != target.end(0); ++it) {
     mVertex *v = (mVertex*) *it;
     temp(v) = ls(v) + dt * H(v);
   }
   temp.accept(space_operator, H);
   
   for(xIter it = target.begin(0); it != target.end(0); ++it) {
     mVertex *v = (mVertex*) *it;
     ls(v) = 0.5 * (ls(v) + temp(v)) + dt/2. * H(v);
   }
}

void xForwardEuler::visit(xLevelSet& ls, xRegion target) {
   xLevelSet temp(ls);
   xLevelSet H;
   ls.accept(space_operator, H);
//  ... boucle sur les noeuds
   for(xIter it = target.begin(0); it != target.end(0); ++it) {
     mVertex *v = (mVertex*) *it;
     ls(v) +=  0.5*dt * H(v);
   }
}
 

void xComputeError::loadReference(const xLevelSet& ref) {
  reference = ref;
}

 std::pair<double, double> xComputeError::getError() const{ return std::make_pair (errnorm, refnorm); }


void xL2norm::visit(const xLevelSet& current, xRegion target) {
 //calcul err
  errnorm = 0.0;
  refnorm = 0.0;
  for(xIter it = target.begin(0); it != target.end(0); ++it) {
    mVertex* v = (mVertex*) *it;
    double r = reference(v);
    double c = current(v);
    errnorm += (r-c)*(r-c);
    refnorm +=  r*r;
  }




  errnorm = sqrt(errnorm);
  refnorm = sqrt(refnorm);
}






void xFitToVertices::visit(xLevelSet& ls, xRegion target)  
{
  const bool debug = xdebug_flag;
  const bool debug_f = false;
 
  if (debug) cout << "Xfit ds visit - 1"<<endl;
  lsf=&ls;
  if (debug) cout << "Xfit ds visit - 2"<<endl;
  targetf=target;
  if (debug_f)
    {
      cout <<"---------------------------------------------------" << endl;
      cout <<"before com" << endl;
      cout <<"----------" << endl;
      cout<<"id " << " " << "ls" << endl;
      cout <<"------" << endl;
      for(xIter itn = targetf.begin(0); itn != targetf.end(0); ++itn)
	{
	  mEntity *n = *itn;
	  cout << n->getId() << " " << (*lsf)(n) << endl;
	}
      cout <<"---------------------------------------------------" << endl;
    }
  if (debug) cout << "Xfit ds visit - 3"<<endl;
  pre_exchange();
  if (debug) cout << "Xfit ds visit - 4"<<endl;
  exchange(targetf.getMesh());
  if (debug) cout << "Xfit ds visit - 5"<<endl;
  post_exchange();
  if (debug) cout << "Xfit ds visit - 6"<<endl;

  if (debug_f)
    {
      cout <<"---------------------------------------------------" << endl;
      cout <<"after com" << endl;
      cout <<"---------" << endl;
      cout<<"id " << " " << "ls" << endl;
      cout <<"------" << endl;
      for(xIter itn = targetf.begin(0); itn != targetf.end(0); ++itn)
	{
	  mEntity *n = *itn;
	  cout << n->getId() << " " << (*lsf)(n) << endl;
	}
      cout <<"---------------------------------------------------" << endl;
    }
  
  return;
}

void xFitToVertices::pre_exchange()
{ 
  const bool debug_pre = false;
  if (debug_pre) cout << "Xfit ds pre_exchange - 3.1"<<endl;
  for(xIter itn = targetf.begin(0); itn != targetf.end(0); ++itn)
    {
      mEntity *n = *itn;
      if (debug_pre) cout << "Xfit ds pre_exchange - 3.1.1"<<endl;
      map_zero_pre[n]=false;
      if (debug_pre) cout << "Xfit ds pre_exchange - 3.1.1"<<endl;
    }

 //loop over the edges in the mesh
  if (debug_pre) cout << "Xfit ds pre_exchange - 3.2"<<endl;
  for(xIter it = targetf.begin(1); it != targetf.end(1); ++it) 
	{
	  //get the value of the level set at the two end points
	  mEntity *e = *it;
	  mVertex *v0 = (mVertex*) e->get(0,0); 
	  mVertex *v1 = (mVertex*) e->get(0,1);
	  if (debug_pre) cout << "Xfit ds pre_exchange - 3.2.1"<<endl;
	  double ls0 = (*lsf)(v0) - fit_value;
	  double ls1 = (*lsf)(v1) - fit_value;  
	  //check if the intersection is inside the segment
	  if (debug_pre) cout << "Xfit ds pre_exchange - 3.2.2"<<endl;
	  if (ls0 * ls1 <= 0.0 && ls0 != 0.0 && ls1 != 0.0) {
	    //get the intersection  r is the solution of l1 + r ( l2-l1) = 0
	    //where r is on the interval ]0,1[1
	    double r = -ls0/(ls1 - ls0);

	    if (debug_pre) cout << "Xfit ds pre_exchange - 3.2.3"<<endl;

	    if (r <= TOL) 
	      { //level-set is attached to the first node
	      //ls(v0) = 0.0;
	      mEntity *nv0 = v0;
	      if (debug_pre) cout << "Xfit ds pre_exchange - 3.2.4"<<endl;
	      map_zero_pre[nv0]=true;
	      if(debug_pre) {cout<<"FitNode="<<endl;
	      v0->print();}

	    }

	    else if (r >= 1.-TOL) 
	      { //level-set is attached to the second node
	      mEntity *nv1 = v1;
	      if (debug_pre) cout << "Xfit ds pre_exchange - 3.2.5"<<endl;
	      map_zero_pre[nv1]=true;
	      if(debug_pre) {cout<<"FitNode="<<endl;
	      v1->print();}
	    }
	  }
	}
  if (debug_pre) cout << "Xfit ds pre_exchange - 3.3" <<endl;
  if (debug_pre)
    { 
        map<mEntity*,bool>::iterator itz = map_zero_pre.begin(), itzEnd = map_zero_pre.end();
	if (itz==itzEnd ) {cout << "MAP EMPTY" << endl;} 
	else 
	  {
	    cout << "map zero: " << endl;
	    for ( ; itz != itzEnd; ++itz) 
	      {
		mEntity *n= (*itz).first;
		cout << n->getId() << " "<< (*itz).second  <<  endl;
	      }
	  }
      }
  if (debug_pre) { cout << "Xfit ds add_bucket - 3.4" <<endl;}
  add_bucket(&map_zero_pre);
  if (debug_pre) { cout << "Xfit ds add_bucket - 3.7" <<endl;}
  return;
}

void xFitToVertices::post_exchange()
{  

  const bool debug_post=false;
  std::vector<bucket_type*>::iterator it = buckets.begin();
  if (debug_post) cout << "Xfit ds post_exchange - 5.1" <<endl;
  for (;it!=buckets.end();++it)
    {
      bucket_type& map_zero_post = **it;
      if (debug_post) cout << "Xfit ds post_exchange - 5.1.1" <<endl;
      bucket_type::iterator itz = map_zero_post.begin(), itzEnd = map_zero_post.end();
      //if (itz==itzEnd ) {cout << "MAP EMPTY" << endl;} 
      for ( ; itz != itzEnd; ++itz) 
	{ 
	  if (debug_post) cout << "Xfit ds post_exchange - 5.1.2" <<endl;
	  mEntity *n=((*itz).first);
	  if ((*itz).second == true) 
	    {
	      if(debug_post) 
		{
		  cout<<"FitNode after com="<<endl;
		  n->print();
		}
	      (*lsf)(n) = fit_value;
	    }
	  if (debug_post) cout << "Xfit ds post_exchange - 5.1.3" <<endl;
	}
    }
  return;
}




void xLsetHamiltonJacobi::byelemcomputation(
					    const std::vector<double > &lsv, 
					    const xtensor::xVector &gradlvlset,
					    const double  &F, const double &f, 
					    const xElementInfo& info ,
					    std::vector<double > & maxalphatild, 
					    double  &sommaxalphatild,
					    double  &somdeltafi) const 
{
  const double vol = info.getVolume();
  const int nnode =  lsv.size();
  vector<double> k(nnode), km(nnode), kp(nnode);
  double sommekmi = 0.;
  double normgls = sqrt(gradlvlset * gradlvlset);
  //  ... CE TEST EST EFFECTUE SI LA NORME DU GRADIENT EST NULLE
  if (normgls!=0.) 
    {
      for (int j = 0; j < nnode; j++) 
	{
	  k[j]  = vol * F/normgls * (gradlvlset * info.getGradFFNode(j));
	  km[j] = min(k[j],0.);
	  kp[j] = max(k[j],0.);
	  somdeltafi += k[j]*lsv[j];
	  sommekmi   += km[j];
	}
    }
  else 
    {
      for (int j = 0; j < nnode; j++)
	{
	  k[j]   = 0.;
	  km[j]  = min(0.,0.);
	  kp[j]  = max(0.,0.);
	}
    }
  vector<double> deltafi(nnode);
  
  //  ... CE TEST EST EFFECTUE SI sommekmi EST NUL
  if (sommekmi!=0.) 
    {
      for(int j=0;j< nnode;j++)
	{
	  double km_int = 0.;
	  for(int k=0; k < nnode; ++k)
	    {
	      km_int +=km[k]*(lsv[j]-lsv[k]);
	    }
	  deltafi[j] = (kp[j]/sommekmi) * km_int;
	}
    }
  else 
    {
      for (int k=0; k< nnode; ++k)
	{
	  deltafi[k] = 0.;
	}
    }
  
  //  ... CE TEST EST EFFECTUE SI LA VITESSE EST NULLE
  if (somdeltafi != 0.)
    {
      for(int j=0; j < nnode;j++)
	{
	  maxalphatild[j] = max(0.,deltafi[j]/somdeltafi);
	}
    }
  else 
    {
      for (int k=0; k< nnode; k++)
	{
	  maxalphatild[k] = 0.;
	}
    } 
  for (int k=0; k< nnode; k++)
    {
      sommaxalphatild += maxalphatild[k];
    }
  return;
}



double xLsetHamiltonJacobi:: evalH(const xLevelSet& ls, mVertex *v){ 
  xRegion s = ls.getSupport();
  const int dim = s.dim();
  const double coef_int = (dim==3)?6./24.:1./3. ;
  xElementInfoManager& info_manager = xElementInfoManagerSingleton::instance();

  starvalue star;
  star.l=star.ll=star.w=star.ww =0.;
  for ( AOMD::mAdjacencyContainer::iter it = v->begin(dim);  it != v->end(dim); ++it){
    mEntity *e = *it;
    const int nnode = e->size(0);
    xElementInfoManager::xElementInfoPtr info = info_manager.getInfo(e);
    double vol = info->getVolume();
    xtensor::xVector gradlvlset = Fop->getGrad(e);
    double F = 0., f=0.;
    vector<double> lsv(nnode);
    for (int j = 0; j < nnode; ++j) 
      {
	mVertex* v = (mVertex*) e->get(0,j); //  determination de l integrale du champ des vitesses
	lsv[j]      = ls(v);
	F   += coef_int * ((*Fop)(v));
	f   += coef_int * ((*fop)(v));      // on integre f (source) sur l'element (membre de droite)
      }
    vector<double>  maxalphatild(nnode);
    double sommaxalphatild = 0., somdeltafi =0.;
    byelemcomputation(lsv, gradlvlset, F, f, *info, maxalphatild, sommaxalphatild, somdeltafi);
    
    //  ... CE TEST EST EFFECTUE POUR LES BORDS A VITESSE RENTRANTE
    int j=0;
    
    while ((mVertex *) e->get(0,j)!= v) ++j;
    if (sommaxalphatild != 0.0) 
      {
	star.l  += (maxalphatild[j]/sommaxalphatild) * (somdeltafi-f*vol);
	star.w  += (maxalphatild[j]/sommaxalphatild) * vol;
	star.ll += somdeltafi;
	star.ww += vol;
      }
    else 
      {
	star.ll += somdeltafi;
	star.ww += vol;
      }
  }
  double H= 0.;
  if (star.w == 0.0) 
    {
      if (rentrant) H = -star.ll/star.ww;
      else H=0.0;
    }
  else
    {
      H = - star.l/star.w;
    }
  //std::cout << H << std::endl;
  return H;
}


//void xLsetHamiltonJacobi::visit(const xLevelSet& ls, xLevelSet &H){ 
// xRegion s = ls.getSupport();
// for (xIter it = s.begin(0); it!=s.end(0); ++it){
//   mVertex * v = (mVertex *)(*it);
//   H (v) = evalH(ls, v);
// }	     
//}


void xLsetHamiltonJacobi::visit(const xLevelSet& ls, xLevelSet& H)
{
  xRegion s = ls.getSupport();
  const int dim = s.dim();
  const double coef_int = (dim==3)?6./24.:1./3. ;
  xElementInfoManager& info_manager = xElementInfoManagerSingleton::instance();
  xVertexValuesStorage< starvalue > starvalues(s);
  
  //loop over the elements of the support loop1 
  for(xIter it = s.begin(); it != s.end(); ++it)
    {
      mEntity *e = *it;
      const int nnode = e->size(0);
      xElementInfoManager::xElementInfoPtr info = info_manager.getInfo(e);
      double vol = info->getVolume();
      xtensor::xVector gradlvlset = Fop->getGrad(e);
      double F = 0., f=0.;
      vector<double> lsv(nnode);
      for (int j = 0; j < nnode; ++j) 
	{
	  mVertex* v = (mVertex*) e->get(0,j); //  determination de l integrale du champ des vitesses
	  lsv[j]      = ls(v);
	  F   += coef_int * ((*Fop)(v));
	  f   += coef_int * ((*fop)(v));      // on integre f (source) sur l'element (membre de droite)
	}
      vector<double>  maxalphatild(nnode);
      double sommaxalphatild = 0., somdeltafi =0.;
      byelemcomputation(lsv, gradlvlset, F, f, *info, maxalphatild, sommaxalphatild, somdeltafi);
      if (sommaxalphatild != 0.0) 
	{
	  for(int j=0;j< nnode; ++j)
	    {
	      starvalue &current = starvalues((mVertex *) e->get(0,j) );
	      current.l  += (maxalphatild[j]/sommaxalphatild) * (somdeltafi-f*vol);
	      current.w  += (maxalphatild[j]/sommaxalphatild) * vol;
	      current.ll += somdeltafi;
	      current.ww += vol;
	    }
	}
      else // ajout eric
	{
	  for(int j=0;j<e->size(0);j++)
	    {  
	      mVertex *v = (mVertex *) e->get(0,j);
	      starvalue &current = starvalues(v);
	      current.ll += somdeltafi;
	      current.ww += vol;
	    }
	}	
    } //end loop over the elements of the support loop1 
  

  
  for(xIter it = s.begin(0); it != s.end(0); ++it)
    {
      mVertex *v = (mVertex *) *it;
      const starvalue &current = starvalues(v);
//  ... on determine l'Hamiltonien au noeud courant
//  ... CE TEST EST EFFECTUE POUR LES BORDS A VITESSE RENTRANTE
      if (current.w == 0.0) 
	{
	  if (rentrant) H(v) = -current.ll/current.ww;
	  else H(v)=0.0;
	}
      else
	{
	  H(v) = - current.l/current.w;
	}
    }
}



double xCFL::getDt(const xLevelSet& velo) {
//  ... on va chercher le maillage
  xRegion s = velo.getSupport();
  //  ... boucle sur les elements pour determiner le pas de temps critique
//  assert(m != 0);
  double vmax;
  double dtc = std::numeric_limits<double>::max();
  
  xElementInfoManager&  info_manager = xElementInfoManagerSingleton::instance();
  double coef;
  switch(s.getMesh()->dim())
    {
    case 3 :
      coef = 1./3.;
      break;
    case 2 :
      coef = 1./2.;
      break;
    default :
      throw;
    }

   for(xIter it = s.begin(); it != s.end(); ++it) {
//  for(it->first(); !it->done(); it->next() ) {  
    //  ... determination de l element courant
    mEntity *e = *it;
    xElementInfoManager::xElementInfoPtr info = info_manager.getInfo(e);
    double h = pow(info->getVolume(), coef);
    vmax = 0.;
    for(int j=0;j<e->size(0);j++)
      {
	mVertex *v = (mVertex*) e->get(0,j);
// eric : fabs car vitesse a prirori signee
	vmax = std::max(fabs(velo(v)), vmax);
      }
//  ... TEST EN CAS DE VITESSE NULLE
    double dt; 
    if (vmax == 0.) dt = 100000000.;
    else            dt = h / vmax / 2.;
    dtc = std::min(dtc, dt);
    
  }
  // dtc have to be reduced on all procs !
  assert(xtool::workInProgress());

  return(dtc);
}

double xCFL::getDt(xRegion s) {
  return ((xCFL::getLength(s))/1./2.);
}

double xCFL::getLength(xRegion s) {
  //  ... boucle sur les elements pour determiner le pas de temps critique
  xElementInfoManager& info_manager = xElementInfoManagerSingleton::instance();
  //mIteratorPtr it(s.getTopIter());
  double hmin = std::numeric_limits<double>::max();
  double coef;
  switch(s.getMesh()->dim())
    {
    case 3 :
      coef = 1./3.;
      break;
    case 2 :
      coef = 1./2.;
      break;
    default:
      abort();
      break;
    }

  for(xIter it = s.begin(); it != s.end(); ++it) {
  //for(it->first(); !it->done(); it->next() ) {  
    //  ... determination de l element courant
    mEntity *e = *it;
    xElementInfoManager::xElementInfoPtr info = info_manager.getInfo(e);
    double h = pow(info->getVolume(), coef);
    hmin = std::min(h, hmin); 
  }
  // hmin have to be reduced on all procs !
  assert(xtool::workInProgress());

  return hmin;
}


void xParallelLevelSetUnifier::visit(xLevelSet& _ls, xRegion target){
  ls =&_ls;
  const xSubMesh &sub = *target.getSubMesh();
  const xPartitionBoundaryInfo &pbi = sub.getPartitionBoundaryInfo();
  DataExchange(pbi, *this );
}

void xParallelLevelSetUnifier::send(mEntity *e, data_type &val) const{
  if (e->getType() == mEntity::VERTEX)
    val = (*ls)(e);
}

void  xParallelLevelSetUnifier::receive(mEntity *e, const std::vector<data_type > &vals)  const{
  if (e->getType() == mEntity::VERTEX){
    (*ls)(e) = std::min(*(std::min_element(vals.begin(), vals.end())), (*ls)(e) );
    // std::cout << vals.size() << std::endl;
  }
}

void xLevelSetSmoother::visit(const xLevelSet & in, xLevelSet &out){
  out.setSupport(in.getSupport());
  for(xIter it = in.getSupport().begin(0); it != in.getSupport().end(0); ++it ) {
    mVertex *v = (mVertex *) *it;
    Trellis_Util::mPoint P0 = v->point();
    int nedge = v->size(1);
    xDenseMatrix M(3,3);
    xCSRVector V(3), X(3);
    for (int i= 0; i < 3 ; ++i){
      V(i) = 0; X(i)=0.;
      for (int j= 0; j < 3 ; ++j) M(i, j) = 0.;
    }
    double ls0 = in(v);
    V(2) += ls0;
    M(2,2) +=1;
    for (int i=0; i < nedge; ++i){
      mEntity * edge =  v->get(1, i);
      mVertex * vi =(mVertex * ) ((edge->get(0, 0)==v) ? edge->get(0,1): edge->get(0,0));
      Trellis_Util::mPoint Pi = vi->point();
      double xi = Pi(0) - P0(0);
      double yi = Pi(1) - P0(1);
      double lsi = in(vi);
      M(0,0) += xi*xi;
      M(0,1) += xi*yi;
      M(0,2) += xi;
      V(0)   += lsi*xi;
      M(1,0) += yi*xi;
      M(1,1) += yi*yi;
      M(1,2) += yi;
      V(1)   += lsi*yi;
      M(2,0) += xi;
      M(2,1) += yi;
      M(2,2) += 1.;
      V(2) += lsi;
    }
    solve(M, V, X);
    out(v) = X(2);
  }
}

xEvalShiftLevelSetModifier::xEvalShiftLevelSetModifier(const xEval<double>& eval, const xEntityFilter filter)
    : eval_shift(eval), filter_upper(filter) 
{
}
 
double xEvalShiftLevelSetModifier::getMaxShift(void)
{
    return max_shift;
}
void xEvalShiftLevelSetModifier::visit(xLevelSet& ls, xRegion target) 
{
  const int dim = target.dim();
  double shift;
  max_shift=0.;
  for(xIter it = target.begin(0); it != target.end(0); ++it) 
  {
        mVertex *v = (mVertex*) *it;
        int nb=v->size(dim);
		mEntity *upper=NULL;
        for (int i=0;i<nb;++i)
        {
		    upper = v->get(dim,i);
            if (filter_upper(upper)) break;
            else upper = NULL;
        }
        if (upper)
        {
		    xGeomElem geomelem(upper);
		    geomelem.setUVWForXYZ(v->point());
		    eval_shift(&geomelem,&geomelem,shift);
            max_shift = std::max(shift,max_shift);
            ls(v) += shift;
        }
#ifndef NDEBUG
        else
        {
            cout<<"Warning xEvalShiftLevelSetModifier : No upper found for vertex "<<v->getId()<<"\n level set on this vertex not shifted\n";
        }
#endif
  }
}
  

} // end of namespace


