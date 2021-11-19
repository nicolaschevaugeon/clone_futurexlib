/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include "xMesh.h"
#include "xRegion.h"
#include "mAOMD.h"
#include "mEdge.h"
#include "mFace.h"
#include "mBuildAdj.h"
#include "mTet.h"
#include "mVertex.h"
#include "mMirrorEntity.h"
#include "mEntity.h"
#include "mEntityContainer.h"
#include "xEntityFilter.h"
#include "xElement.h"
#include "mAOMD.h"
#include "AOMD.h"
#include "AOMD_OwnerManager.h"
#include "mPoint.h"
#include "Mapping.h"
#include "xRegularGrid.h"
#include "xDebug.h"
#include "OctreeCreate2.h"
#include "mEntity.h"
#include "mPrism.h"
#include "mHex.h"
#include "ParUtil.h"
#include "xSubMesh.h"
#ifdef PARALLEL
#include "autopack.h"
#endif


namespace xfem {

  
  using std::cout;
  using std::endl;
  using std::ostream;
  using std::cerr;
  using std::make_pair;
  using namespace AOMD;
  using namespace Trellis_Util;

  unsigned int xMesh::was_created_by_tag          = AOMD_Util::Instance()->lookupMeshDataId("was_created_by");
  unsigned int xMesh::is_the_creator_of_tag       = AOMD_Util::Instance()->lookupMeshDataId("is_the_creator_of");
  unsigned int xMesh::r_on_edge_tag               = AOMD_Util::Instance()->lookupMeshDataId("r_on_edge");
  unsigned int xMesh::is_duplicated_in_tag        = AOMD_Util::Instance()->lookupMeshDataId("is_duplicated_in");
  unsigned int xMesh::was_duplicated_from_tag     = AOMD_Util::Instance()->lookupMeshDataId("was_duplicated_from");
  unsigned int xMesh::is_in_partition_of_tag      = AOMD_Util::Instance()->lookupMeshDataId("is_in_partition_of");
  unsigned int xMesh::partition_tag               = AOMD_Util::Instance()->lookupMeshDataId("partition");
  unsigned int xMesh::is_hanging_on_tag           = AOMD_Util::Instance()->lookupMeshDataId("is_hanging_on");
  unsigned int xMesh::is_hanging_by_tag           = AOMD_Util::Instance()->lookupMeshDataId("is_hanging_by");
  unsigned int xMesh::octree_level_tag            = AOMD_Util::Instance()->lookupMeshDataId("octree_level");
  unsigned int xMesh::down_groupe_tag             = AOMD_Util::Instance()->lookupMeshDataId("down_groupe_tag");
  unsigned int xMesh::bnd_groupe_tag              = AOMD_Util::Instance()->lookupMeshDataId("bnd_groupe_tag");

  unsigned int xMesh::was_created_by_tag2         = AOMD_Util::Instance()->lookupMeshDataId("was_created_by2");
  unsigned int xMesh::is_the_creator_of_tag2      = AOMD_Util::Instance()->lookupMeshDataId("is_the_creator_of2");
  unsigned int xMesh::is_the_copy_of_tag          = AOMD_Util::Instance()->lookupMeshDataId("is_the_copy_of");
  unsigned int xMesh::has_a_copy_tag              = AOMD_Util::Instance()->lookupMeshDataId("has_a_copy_");
  unsigned int xMesh::is_the_element0_of_tag      = AOMD_Util::Instance()->lookupMeshDataId("is_the_element0_of");
  //  unsigned int has_a_copy                         =   AOMD::AOMD_Util::Instance()->lookupMeshDataId ("has_a_copy");
  unsigned int xMesh::refined_elements_tag        = AOMD_Util::Instance()->lookupMeshDataId("refined_elements_tag");  
  unsigned int xMesh::levelset_value_tag          = AOMD_Util::Instance()->lookupMeshDataId("levelset_value_tag"); 
  unsigned int xMesh::cut_edge_tag                = AOMD_Util::Instance()->lookupMeshDataId("cut_edge_tag"); 
  unsigned int xMesh::cut_element_tag             = AOMD_Util::Instance()->lookupMeshDataId("cut_element_tag"); 
  unsigned int xMesh::RH_enrichment_nb_func_tag   = AOMD_Util::Instance()->lookupMeshDataId("RH_enrichment_nb_func_tag"); 
  unsigned int xMesh::RH_enrichment_side_tag      = AOMD_Util::Instance()->lookupMeshDataId("RH_enrichment_side_tag"); 
  unsigned int xMesh::RH_enrichment_status_tag    = AOMD_Util::Instance()->lookupMeshDataId("RH_enrichment_status_tag"); 

  using AOMD::mEdge;
  using AOMD::mFace;
  using AOMD::mTet;
  using AOMD::mAttachableMirrorVertex;

  
  // extension de Aomd pour prendre en compte des
  // sous-domaines dans le maillage

  
  std::pair<Trellis_Util::mPoint, Trellis_Util::mPoint> compute_bounding_box(mEntity * e){

     switch (e->getType()) {
       case mEntity::POINT:
         return make_pair(dynamic_cast<mVertex*>(e)->point(), dynamic_cast<mVertex*>(e)->point());
         break;

       default:
         Trellis_Util::mPoint pmin(dynamic_cast<mVertex*>(e->get(0,0))->point());
         Trellis_Util::mPoint pmax(dynamic_cast<mVertex*>(e->get(0,0))->point());

         for(int i = 1; i < e->size(0); ++i){
             Trellis_Util::mPoint p(dynamic_cast<mVertex* >(e->get(0,i))->point());

             for (int j=0;j<3;j++)
               {
                 if (pmin(j)>p(j)) pmin(j)=p(j);
                 if (pmax(j)<p(j)) pmax(j)=p(j);
               }
           }

         return make_pair(pmin, pmax);
         break;
       }




   }
  
  
  xMesh::~xMesh() {
    for (iter_subsets it = subsetEntities.begin(); it != subsetEntities.end(); ++it) {
      delete it->second;
    }
    if (RegularGrid) delete RegularGrid;
    if (octree)  Octree_Delete(octree);
  }

  mEntity *getSource(mEntity * e){
    mEntity *re = e->getAttachedEntity(xMesh::was_duplicated_from_tag);
    if(re) return( getSource(re));
    else return e;
  }

  
  xMesh::xMesh(int id) : mMesh(id), RegularGrid(0), octree(0) {
	}

  xMesh::xMesh(const string& filename) : mMesh(), RegularGrid(0), octree(0) {
    std::ifstream meshfile(filename.c_str());
    if (!meshfile.is_open()) {
      cout << "can't open meshfile "<< filename << " in xMesh Constructor " << __FILE__ << " " << __LINE__ << endl;
      throw;
    }
    meshfile.close();
    AOMD_Util::Instance()->import(filename.c_str(),this);
    modifyAllState();
    AOMD::ParUtil::Instance()->Barrier( __LINE__, __FILE__ );
    bdryLinkSetup();
    AOMD::ParUtil::Instance()->Barrier( __LINE__, __FILE__ );
  
    //cout << " done!" << endl;
    //Re Classify faces on partition boundary.
    // Those face, for historical reason, are classified at this point on GEntity of level 2
    // They should be classified on the volum. That is  what we do in the following.
    AOMD::AOMD_OwnerManager*  pom = this->theOwnerManager;
    for( AOMD::AOMD_OwnerManager::iter it = pom->begin(); it != pom->end(); ++it ){
      mEntity * e = (*it).first;
      if (e->getLevel() == (dim()-1)){
	int glevel  = GEN_type(e->getClassification());
	if (glevel == (dim()-1)){
	  e->classify (e->get(dim(), 0)->getClassification());
	}
      }  
    }

  
    //findHangingNodes();

  
    std::cout << "Proc " << AOMD::ParUtil::Instance()->rank()+1
	      << " nb nodes "<< size(0) 
	      << " nb edges " << size(1)
	      << " nb faces " << size(2)
	      << " nb regions " << size(3) << endl;
  
  }

  void xMesh::modifyAllState() {
    modifyState(3,2,true);
    modifyState(3,1,true);
    modifyState(3,0,true);
    modifyState(2,1,true);
    modifyState(2,0,true);
    modifyState(1,0,true);
    modifyState(0,1,true);
    modifyState(0,2,true);
    modifyState(0,3,true);
    modifyState(1,2,true);
    modifyState(1,3,true);
    modifyState(2,3,true);
  }


  void xMesh::modifyAllStateFalse() {
    modifyState(3,2,false);
    modifyState(3,1,false); 
    modifyState(2,1,false);
    modifyState(2,0,false);
    modifyState(0,1,false);
    modifyState(0,2,false);
    modifyState(1,2,false);
    modifyState(1,3,false);
    modifyState(2,3,false);
    modifyState(2,1,false);
    modifyState(2,0,false); 
    modifyState(1,0,false);
    modifyState(0,1,false);
    modifyState(0,2,false);
    modifyState(1,2,false);
  }

  

  void xMesh::compute_bounding_box(Trellis_Util::mPoint &min,Trellis_Util::mPoint &max ) const
  {
    iter it=begin(0);
    iter ite=end(0);
    if (it!=ite)
      {
	mVertex* v =(mVertex*) *it;
	min=v->point();
	max=min;
      }
    else
      {
	for (int j=0;j<3;j++)
	  {
	    min(j)=0;
	    max(j)=0;
	  }
      }
    ++it;
    for (;it!=ite;++it)
      {
	mVertex*v =(mVertex*)*it;
	Trellis_Util::mPoint p=v->point();
	for (int j=0;j<3;j++)
	  {
	    if (min(j)>p(j)) min(j)=p(j);
	    if (max(j)<p(j)) max(j)=p(j);
	  }
      }
  }



  std::pair<double,double> xMesh::checkMesh(const int fact) const
  {
    int d=dim();
    int n=1;
    iter it=begin(d);
    iter ite=end(d);
    double vol=xElement(*it).getVolume();
    double max=vol;
    double min=vol;
    double moy=vol;
    for (++it;it!=ite;++it,++n)
    {
       vol=xElement(*it).getVolume();
       if(vol>max) max=vol;
       else if (vol<min) min=vol;
       moy+=vol;
    }
    moy/=n;
    cout<<" "<<endl;
    cout<<"==CheckMesh output ================================================================================================"<<endl;
    cout<<"Volumes statistique : moy="<<moy<<" min="<<min<<" max="<<max<<endl;

    for (it=begin(d);it!=ite;++it)
    {
       vol=xElement(*it).getVolume();
       if (vol*fact<moy)
       {
         cout<<"Element id "<<(*it)->getId()<<" has a volume value of "<<vol<<" wich is "<<fact<<" time lesser then average value "<<moy<<endl;
         cout<<"Its connectivity is :";
         for (int i=0;i<(*it)->size(0);++i)
         {
              cout<<"  "<<((*it)->get(0,i))->getId();
         }
         cout<<endl;
       }
    }
    cout<<"==End CheckMesh output ============================================================================================"<<endl;

    return std::make_pair(min,max);
  }

  //new version
  void xMesh::del(mEntity *e) {
    //mesh->del(e);
    mMesh::del(e);
    //if we are removing a pointer from allEntities, we need
    //to remove it from all the subsets
    for (iter_subsets it = subsetEntities.begin(); 
	 it != subsetEntities.end(); ++it) {
      it->second->del(e);
    }
  }
 

  xSubMesh & xMesh::getSubMesh(const string& sub) const {
    iter_subsets it = subsetEntities.find(sub);
    if (it == subsetEntities.end()) {
      cerr << "The entities set with name " << sub  << " does not exist\n";
      throw;
    }
    return *(it->second);
  }

  /*
  void xMesh::modifyAllState_sub(const string& name ) const{
    getSubMesh(name).modifyAllState();
  }

  void xMesh::updatePartitionBoundary_sub (const string& name) const {
    getSubMesh(name).updatePartitionBoundary();
  }

  void xMesh::updateBoundary_sub (const string& name) const{
    getSubMesh(name).updateBoundary();
  }

  mEntity* xMesh::find_sub(mEntity *e, const string& sub) const
  {
    return getSubMesh(sub).find(e);
  }

  int xMesh::dim_sub(const string& sub) const
  {
    return getSubMesh(sub).dim();
  }
  */
  
  //protected function
 
  
  //NICO2 for debug
  /*
  void xMesh::printSubsetEntities() const
  {
    printf("The subset entities are:\n");
    for (const_iter_subsets it = subsetEntities.begin(); it != subsetEntities.end(); ++it) {
      printf("%s : %d %d %d %d\n", it->first.c_str(), it->second->size(0), 
	     it->second->size(1),  it->second->size(2), it->second->size(3));
      printf("ids of the vertices\n");
      for  (xIter itv = begin_sub(0, it->first); 
	    itv != end_sub(0, it->first); ++itv) (*itv)->print();
      printf("ids of the edges\n");
      for  (xIter itv = begin_sub(1, it->first); 
	    itv != end_sub(1, it->first); ++itv) (*itv)->print();
      printf("ids of the elements\n");
      for  (xIter itv = begin_sub(dim(), it->first); 
	    itv != end_sub(dim(), it->first); ++itv) (*itv)->print();
    }
    return;
    }*/
  
  xSubMesh & xMesh::createSubMesh(const string& sub) const
  {
    deleteSubMesh(sub);
    xSubMesh *subm = new xSubMesh( sub, *this);
    subsetEntities[sub] = subm;
    return *subm;
  }

  xSubMesh & xMesh::createSubMesh(const string& sub, xSubMeshCreator& creator) const
  {
    xSubMesh *subm = &createSubMesh(sub);
    creator.create(*this, sub);
    return *subm;
  }
  
  void xMesh::deleteSubMesh(const string& sub) const
  {
    iter_subsets it = subsetEntities.find(sub);
    if (it != subsetEntities.end()) {
      delete it->second;
      subsetEntities.erase(it);
    }
  }

  void xMesh::deleteAllSubMesh() const{
    for (iter_subsets it = subsetEntities.begin(); it != subsetEntities.end();++it)
      {
	delete it->second;
      }
    subsetEntities.clear();
  }


  int xMesh::dim() const
  {
    //  return mesh->size(3)?3:(mesh->size(2)?2:(mesh->size(1)?1:0));
    return size(3)?3:(size(2)?2:(size(1)?1:0));
  }

  int xMesh::dim_global() const
  {
    int dimloc = dim();
    int dimglob = dimloc;
#ifdef PARALLEL
    MPI_Allreduce(&dimloc,&dimglob,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
    return dimglob;
  }

  

  //  bool xMesh::isHanging(mEntity* e)
  //  {
  //    mEntity*  ent = e->getAttachedEntity("hanging_on");
  //    if (!ent) return false;
  //    return true; 
  //  }

  bool xMesh::isMirror(mEntity* e) 
  {
    //mVertex* v = (mVertex*) e;
    mAttachableMirrorVertex* att = (mAttachableMirrorVertex*) e->getData(mirrorTag);         
    if (!att) return false;
    return true;
  }


  mMirrorVertex* xMesh::getMirrorVertex(mEntity* e)
  {
    //mVertex* v = (mVertex*) e;
    mAttachableMirrorVertex* att = (mAttachableMirrorVertex*) e->getData(mirrorTag);
    if (!att) return 0;
    return (mMirrorVertex*)att->e;
  }

  /*

  //xIter 
  xIter xMesh::begin_sub(int what, const string& sub) const {
    return getSubMesh(sub).begin(what);
  }


  //xIter 
  xIter xMesh::end_sub(int what, const string& sub) const {
    return getSubMesh(sub).end(what);
  }


  //xIter 
  xIterall xMesh::beginall_sub(int what, const string& sub) const {
    return getSubMesh(sub).beginall(what);
  }




  //xIter 
  xIterall xMesh::endall_sub(int what, const string& sub) const {
    return getSubMesh(sub).endall(what);
  }

  int xMesh::size_sub(int what, const string& sub) const {
    return getSubMesh(sub).size(what);
  }



  void xMesh::modifyState_sub (int i , int j , bool state, const string& subset, int with) const
  {
    getSubMesh(subset).modifyState(i,j,state, with);
  }
  */

  // //xIter 
  // xIter xMesh::begin_sub(const string& sub) const {
  //   int dim_mesh = dim();
  //   const mMeshEntityContainer* m = getEntities(sub);
  //   return mLeavesIterator(m->begin(dim_mesh),m->end(dim_mesh));
  // }
  // //xIter 
  // xIter xMesh::end_sub(const string& sub) const {
  //   int dim_mesh = dim();
  //   const mMeshEntityContainer* m = getEntities(sub);
  //   return mLeavesIterator(m->end(dim_mesh),m->end(dim_mesh));
  // }


  //  void xMesh::findHangingNodes() {
  //    int dim_mesh = dim();
  //  cout << "dim_mes " << dim_mesh << endl;
  //    for (iter it = begin(dim_mesh); it != end(dim_mesh); ++it) {
  //    mEntity* e = *it;
  //    e->print();
  //    //printf("elem level raffin %d\n", mesh->getRefinementLevel(e));
  //    //printf("elem depth raffin %d\n", mesh->getRefinementDepth(e));
  //    //printf("nb edges %d\n",  e->size(1));
  //    for (int i = 0; i < e->size(1); ++i) {
  //      mEntity* ed = e->get(1,i);
  //      ed->print();
  //  //    printf("edge depth raffin %d\n", mesh->getRefinementDepth(ed));
  //  //    printf("edge level raffin %d\n", mesh->getRefinementLevel(ed));
  //      if (getRefinementDepth(ed) == 1) {
  //        //printf("edge raffin�\n");
  //        std::list<mEntity*> childs;
  //        ed->getLeaves(childs);      
  //        //printf("nb childs is %d\n", childs.size());
  //        std::list<mEntity*>::iterator ite = childs.begin();
  //        mEdge* ed1 = (mEdge*) *ite;
  //        mEdge* ed2 = (mEdge*) *(++ite);
  //        mVertex* common = ed1->commonVertex(ed2);
  //        if (common->getAttachedEntity("hanging_on") == 0) {
  //  	//printf("hanging node %d added for an edge\n", common->getId());
  //  	//printf("tying nodes %d and %d \n", ed->get(0,0)->getId(), ed->get(0,1)->getId());
  //  	mAttachableEntity * me = new mAttachableEntity;
  //  	me->e = ed;
  //  	common->attachData("hanging_on", me);
  //  	//assert(common->getAttachedEntity("hanging_on") != 0);
  //        }
  //      }
  //    }
  //    if (dim_mesh == 3) {  
  //      for (int i = 0; i < e->size(2); ++i) {
  //        mEntity* fac = e->get(2,i);
  //        //printf("face depth raffin %d\n", mesh->getRefinementDepth(fac));
  //        //printf("face level raffin %d\n", mesh->getRefinementLevel(fac));
  //        if (getRefinementDepth(fac) == 1 && fac->getType() == mEntity::QUAD) {
  //  	std::list<mEntity*> childs;
  //  	fac->getLeaves(childs);      
  //  	std::list<mEntity*>::iterator ite = childs.begin();
  //  	mFace* fac1 = (mFace*) *ite;
  //  	mFace* fac2 = (mFace*) *(++ite);
  //  	mFace* fac3 = (mFace*) *(++ite);
  //  	mVertex* common = fac1->commonVertex(fac2,fac3);
  //  	if (common->getAttachedEntity("hanging_on") == 0) {
  //  	  //printf("hanging node %d added for a face\n", common->getId());
  //  	  //printf("tying nodes %d and %d and %d and %d\n", 
  //  		// fac->get(0,0)->getId(), 
  //  		// fac->get(0,1)->getId(),
  //  		// fac->get(0,2)->getId(),
  //  		// fac->get(0,3)->getId());
  //  	  mAttachableEntity * me = new mAttachableEntity;
  //  	  me->e = fac;
  //  	  common->attachData("hanging_on", me);
  //  	  //assert(common->getAttachedEntity("hanging_on") != 0);
  //  	}
  //        }
  //      }
  //    }
  //  }


  //  //
  //  //if a mirror has hanging in its mirrors, it must be hanging itself
  //  //
  //    for (iter it = begin(0); it != end(0); ++it) {
  //      mEntity* e = *it;
  //      bool mirror  = xMesh::isMirror(e);
  //      bool hanging = xMesh::isHanging(e);
  //      if (mirror && hanging) // we must propagate the hanging info to the mirrors 
  //        {
  //          mEntity* master = e->getAttachedEntity("hanging_on");
  //  	mAttachableMirrorVertex* att = (mAttachableMirrorVertex*) e->getData(mirrorTag);   
  //  	mMirrorVertex* mv = (mMirrorVertex*)att->e;
  //  	for (int i = 0; i < mv->nbCopies(); ++i) 
  //  	  {
  //  	    mVertex* v = mv->getCopy(i);
  //  	    if (!v->getAttachedEntity("hanging_on"))
  //  	      {
  //  		mAttachableEntity * me = new mAttachableEntity;
  //  		me->e = master;
  //  		v->attachData("hanging_on", me);
  //  	      }
  //  	  }
  //        }
  //    }
  
  //    //Printf("nb of hanging nodes %d \n", hanging.size()); 
  //    for (iter it = begin(0); it != end(0); ++it) {
  //      cout << (*it)->getId() << " hanging " << ((*it)->getAttachedEntity("hanging_on") != 0) << endl;
  //    }  
  //  }


  //mEdge *xMesh::getEdge(mVertex *v1, mVertex *v2) {return mesh->getEdge(v1, v2);}


  struct lessUpToEpsilon {

    //given x1 and x2 it returns true
    //if x1 and x2 differ by more than epsilon
    //and x1 < x2
    //for vertices v1 and v2 it returns true
    //if the distance between v1 and v2 differ by more than epsilon
    //and x1 < x2 (alphabetique sur X, y, z) 
    bool operator()(mVertex* v1, mVertex* v2) const
    {
      double eps = 1.e-06;
      return v1->point().lexicographicLessThan(v2->point(), eps);
      //      double eps = 1.e-06;
      //      double x1 = v1->point()(0);
      //      double y1 = v1->point()(1);
      //      double z1 = v1->point()(2);
      //      double x2 = v2->point()(0);
      //      double y2 = v2->point()(1);
      //      double z2 = v2->point()(2);
      //      if (std::sqrt((x2-x1)*(y2-y1)*(z2-z1)) < eps) return false;
      //      if (x1 < x2) return true;
      //      if (y1 < y2) return true;
      //      if (y1 < y2) return true;
      //      return false;
    }
  };

  bool eps_equal(const double& x1, const double& x2)
  {
    double eps = 1.e-6;
    if (fabs(x1-x2) < eps) return true;
    return false;
  }



void xMesh::periodicAssociation() {
  const bool debug = xdebug_flag; 
  cout << "begin association" << endl;
  std::multimap<int,int> vertexAssociations;
  /// select the nodes on the boundary
  std::set<mVertex*, lessUpToEpsilon> boundary;
  typedef std::set<mVertex*, lessUpToEpsilon>::iterator bnd_iterator;
  Trellis_Util::mPoint pmin,pmax ;
  compute_bounding_box(pmin,pmax);
  std::cout << pmax << " " << pmin << std::endl;
  for(iter it = begin(0); it != end(0); ++it)
    {
      mVertex *v = (mVertex*)(*it);
      double x = v->point()(0);
      double y = v->point()(1);
      double z = v->point()(2);
      bool onboundary;
      onboundary =(eps_equal(x, pmax(0)) || eps_equal(x, pmin(0)) ||
		   eps_equal(y, pmax(1)) || eps_equal(y, pmin(1))) ;
      if (getDim()==3)
	onboundary = (onboundary||eps_equal(z, pmax(2)) || eps_equal(z, pmin(2)) );
      if (onboundary) boundary.insert(v);
    } 
  std::cout << "size bound" << boundary.size() << std::endl;
  
  //group 
  for(bnd_iterator it = boundary.begin(); it != boundary.end(); ++it)
    {
      mVertex *v = (mVertex*)(*it);
      double x = v->point()(0); double y = v->point()(1); double z = v->point()(2);

	if (eps_equal(x, pmax(0)))  
	{
	  mVertex vf(-1, Trellis_Util::mPoint(pmin(0), y, z), 0);
          bnd_iterator itf = boundary.find(&vf);
	  if (itf != boundary.end())  {
	    if (debug) cout << "the nodes " << v->getId() << " and " 
			    << (*itf)->getId() << " are opposite on x" << endl;
            if (v->getId() < (*itf)->getId()) 
	      vertexAssociations.insert(make_pair(v->getId(), (*itf)->getId()));  
            else
	      vertexAssociations.insert(make_pair((*itf)->getId(), v->getId()));  
	  }
	}
      
      //mat ->taille boite en dure ...
      if (eps_equal(y, pmax(1))) 
	//      if (eps_equal(y, 1.0)) 
	//mat
	{
	  mVertex vf(-1, Trellis_Util::mPoint(x, pmin(1), z), 0);
          bnd_iterator itf = boundary.find(&vf);
	  if (itf != boundary.end()) 
	    {

	      if (debug) cout << "the nodes " << v->getId() << " and " 
			      << (*itf)->getId() << " are opposite on x" << endl;
	      if (v->getId() < (*itf)->getId()) 
		vertexAssociations.insert(make_pair(v->getId(), (*itf)->getId()));  
	      else
		vertexAssociations.insert(make_pair((*itf)->getId(), v->getId()));  
	    }
	  }
      
      if (getDim()==3)
      if (eps_equal(z, pmax(2))) 
	{
	  mVertex vf(-1, Trellis_Util::mPoint(x, y, pmin(2)), 0);
          bnd_iterator itf = boundary.find(&vf);
	  if (itf != boundary.end()) 
	    {
	      if (debug) cout << "the nodes " << v->getId() << " and " 
			      << (*itf)->getId() << " are opposite on z" << endl;
	      if (v->getId() < (*itf)->getId()) 
		vertexAssociations.insert(make_pair(v->getId(), (*itf)->getId()));  
	      else
		vertexAssociations.insert(make_pair((*itf)->getId(), v->getId()));  
	    }
	}
    }
  
  AOMD_Util::resolveAssociations(this, vertexAssociations);
  try{
    bdryLinkSetup(); 
  }
  catch(...){
    cout << "Some problem in AOMD::mMesh::bdryLinkSetup(), catch in" << __FILE__<<":"<<__LINE__<< endl;
  }
  if (debug) {
    cout << "we check if evey node on the boundary has a mirror" << endl;
    cout << "if not it means we dot have a periodic mesh" << endl;
    for(bnd_iterator it = boundary.begin(); it != boundary.end(); ++it)
      { 
	std::list<mEntity*> mirrors;
	lookupForMirrorEntities(*it, mirrors);
	cout << "mirrors.size() is " << mirrors.size() << endl;
	assert(!mirrors.empty()); 
 
      }
  }
}
  

  //if the entity e has a mirror the support is bigger!!
  void xMesh::lookupSupportBasic(mEntity* e, std::set<mEntity*> &l)
  {
    int n = e->getLevel();
    const int d = dim();
    if (n == d)     { 
      l.insert(e);
      return;
    }
    for (int i = 0; i < e->size(d); ++i) {
      mEntity* edim = e->get(d, i);
      l.insert(edim);
    }
    return;
  }


//KO
  void xMesh::lookupSupport(mEntity* e, std::set<mEntity*> &l){
    lookupSupportBasic(e, l);
    std::list<mEntity*> mirrors;
    bool domirror =false;
    if (size(3)){
      if (e->getLevel()<3) domirror =true;
      else if (size(2)){
        if (e->getLevel()<2) domirror =true;
      }
    }
  else if (size(1))
    if (e->getLevel()<1) domirror =true;

    if (domirror){
      try    {
	lookupForMirrorEntities(e, mirrors);
      }
      catch(...)     {
	cout << "Some problem in lookupForMirrorEntities(), catch in" << __FILE__<<":"<<__LINE__<< endl;
	cout << "e "<< e << "type  "<< e->getType() << endl;
      }
      
      std::list<mEntity*>::const_iterator it;
      for (it = mirrors.begin() ; it != mirrors.end(); ++it){
	std::set<mEntity*> lmir;
	lookupSupportBasic(*it, lmir);
	l.insert(lmir.begin(), lmir.end());
      }
    }
    return;
  }   
  
  
  
  
  void xMesh::getPartition(mEntity* e, xPartition& partition, xEntityFilter filter) 
{
	const bool debug = xdebug_flag;
	if (debug) cout << " In getPartition  e is ";
	if (debug)  e->print();
	xAttachableMesh *am = (xAttachableMesh*)e->getData(partition_tag);
	if (!am)     {
		if (debug) cout << " pas de am " << endl;
		if (debug) cout << " filter is " << filter(e) << endl;
		if (filter(e)) 
		{
		    partition.insert(e);
		} 
		return;
	}
	xMesh* m =  am->m;
	if (debug) cout << " am existe " << endl;
	for (xIter it = m->begin(m->dim()); it != m->end(m->dim()); ++it)
	{
		mEntity *e = *it;
		if (debug) cout << " e in partition " << endl;
		if (debug) e->print();
		xMesh::getPartition(*it, partition, filter);
	}
}


void xMesh::createOctree(void) const
{
	if (octree)  Octree_Delete(octree);
	octree=OctreeCreate2(begin(dim()),
			 end(dim()),5);
  }

  void xMesh::locateElementOctree(const Trellis_Util::mPoint& p, std::set<mEntity*>& elts) const
  { 
    if(!octree) createOctree();
    void * ptr;
    double px[3];
    px[0] = p(0); 
    px[1] = p(1);
    px[2] = p(2); 

    ptr= Octree_Search (&px[0], octree);
    if (ptr) elts.insert((mEntity*)ptr);
  } 

 
  void xMesh::locateElement(const Trellis_Util::mPoint& p, std::set<mEntity*>& elts)
  {
    if (!RegularGrid)
      {
	RegularGrid=new xRegularGrid(this);
      }
    RegularGrid->GetElementsForPoint(p, elts);
  }

  //  void xMesh::lookupSupport(mEntity* e, std::set<mEntity*> &l)
  //  {
  //    int n = e->getLevel();
  //    if (n == dim()) 
  //      {
  //        std::list<mEntity*> lis;
  //        e->getLeaves(lis);
  //        l.insert(lis.begin(), lis.end());
  //      }
  //    else if (n < dim() && n != 0)
  //      {
  //        for (int i = 0; i < e->size(dim()); ++i)
  //  	{
  //  	  mEntity* edim = e->get(dim(), i);
  //  	  if (getRefinementDepth(edim) == 0) l.insert(edim);
  //  	  else   
  //  	    for (int j = 0; j < e->size(n); ++j)
  //  	      {
  //  		lookupSupport(e->get(n, j), l);
  //  	      }
  //  	}
  //      }
  //    else //n = 0
  //      for (int i = 0; i < e->size(1); ++i) lookupSupport(e->get(1, i)  ,l);
  //    return;
  //  }   


  /*

  void xMesh::export_all(const string& filename, int formatoption, int merging) const
  {
   
    export_sub(filename, "all", formatoption, merging);
  }

  void xMesh::export_sub(const string& filename, const string& subset, int formatoption, int merging) const{
    getSubMesh(subset).exportGmsh(filename, formatoption, merging);
  }
*/



  //////////////////////////////////////////// SIMPLEXIFY /////////////////////////////////////
  void xMesh:: Simplexify(xEntityFilter filter){
    list<mEntity*> toTetrize;
    enum AOMD::mEntity::mType _t;
    // build the list of element to divide into triangles or tetrahedra
    for(xIter it = this->begin(2);it != this->end(2); ++it) {
      mEntity* _pmE =(*it);
      switch( _t = _pmE->getType() )  {
      case mEntity::QUAD : { if(filter(_pmE)){toTetrize.push_back(_pmE) ;} break;}
      default: {}
      }  
    }
    for(xIter it = this->begin(3);it != this->end(3); ++it) {
      mEntity* _pmE =(*it);
      switch( _t = _pmE->getType() )  {
      case mEntity::HEX :     { if(filter(_pmE)) {toTetrize.push_back(_pmE) ;} break;}
      case mEntity::PRISM :   { if(filter(_pmE)) {toTetrize.push_back(_pmE) ;} break;}
      case mEntity::PYRAMID : { if(filter(_pmE)) {toTetrize.push_back(_pmE) ;} break;}
      default: {}
      }
    }
    // devide the elements of the list
    for(list<mEntity*>::iterator it = toTetrize.begin() ; it!=toTetrize.end() ; ++it) {
      mEntity* _pmE =(*it);
      switch( _t = _pmE->getType() )  {
      case mEntity::QUAD :    { QuadToTri( _pmE )   ; DEL(_pmE); break;}	      
      case mEntity::HEX :     { HexToTet( _pmE )    ; DEL(_pmE); break;}	      
      case mEntity::PRISM :   { PrismToTet( _pmE )  ; DEL(_pmE); break;}	      
      case mEntity::PYRAMID : { PyramidToTet( _pmE ); DEL(_pmE); break;}	      
      default:{}
      }
    }
  }


  void xMesh:: PyramidToTet(mEntity* _pmE){
    vector<int> connect;
    connect.reserve(5) ;
    connect.clear();
    for (int i=0 ; i<5 ;i++){ connect.push_back(  _pmE->get(0,i)->getId() );   }   
    int selection[] = {0,1,2,3};   vector<int> index;   index.assign(selection,selection+4);  
    int imin = IndexOfMinAmong(connect,index) ;
    vector<int> new_connect;
    new_connect.reserve(5) ;
    new_connect.clear();
    for (int i=0 ; i<4 ;i++){  new_connect[(i-imin+4)%4]= connect[i];   }
    int N0=new_connect[0];
    int N1=new_connect[1];
    int N2=new_connect[2];
    int N3=new_connect[3];
    int N4=connect[4];
    createTetWithVertices(N0,N2,N3,N4,_pmE->getClassification() );
    createTetWithVertices(N0,N1,N2,N4,_pmE->getClassification() );
  }

  void xMesh:: PrismToTet(mEntity* _pmE){
    vector<int> connect;
    connect.reserve(6) ;
    connect.clear();
    for (int i=0 ; i<6 ;i++){ connect.push_back(  _pmE->get(0,i)->getId() );  }   
    int imin = min_element(connect.begin(), connect.end()) - connect.begin();
    vector<int> new_connect;
    new_connect.reserve(6);
    new_connect.clear();
    if(imin<3){
      for (int i=0 ; i<3 ;i++){  new_connect[(i-imin+3)%3]= connect[i];    }
      for (int i=3 ; i<6 ;i++){  new_connect[(i-imin+3)%3+3]= connect[i];  }
    }
    if(2<imin){
      for (int i=3 ; i<6 ;i++){ new_connect[(-i+imin+3)%3]= connect[i]; }
      for (int i=0 ; i<3 ;i++){ new_connect[(-i+imin+3)%3+3]= connect[i]; }
    }
    int N0=new_connect[0];
    int N1=new_connect[1];
    int N2=new_connect[2];
    int N3=new_connect[3];
    int N4=new_connect[4];
    int N5=new_connect[5];
    createTetWithVertices(N3,N5,N4,N0,_pmE->getClassification() );
    mEntity* _pmP;
    _pmP = createHexWithVertices(N4,N5,N2,N1,N0,N0,N0,N0,_pmE->getClassification() );
    PyramidToTet( _pmP);
    DEL(_pmP);
  }

  void xMesh:: HexToTet(mEntity* _pmE){
    vector<int> connect;
    connect.reserve(8) ;
    connect.clear();
    for (int i=0 ; i<8 ;i++){  connect.push_back(  _pmE->get(0,i)->getId() );  }   
    int imin = min_element(connect.begin(), connect.end()) - connect.begin();
    vector<int> new_connect;
    new_connect.reserve(8);
    new_connect.clear();
    if(imin<4){
      for (int i=0 ; i<4 ;i++){  new_connect[(i-imin+4)%4]= connect[i];    }
      for (int i=4 ; i<8 ;i++){  new_connect[(i-imin+4)%4+4]= connect[i];  }
    }
    if(3<imin){
      for (int i=4 ; i<8 ;i++){ new_connect[(-i+imin+4)%4]= connect[i];    }
      for (int i=0 ; i<4 ;i++){ new_connect[(-i+imin+4)%4+4]= connect[i];  }
    }
    int N0=new_connect[0];
    int N1=new_connect[1];
    int N2=new_connect[2];
    int N3=new_connect[3];
    int N4=new_connect[4];
    int N5=new_connect[5];
    int N6=new_connect[6];
    int N7=new_connect[7];
    mEntity* _pmP;
    int selection[] = {2,3,6,7};   vector<int> index;   index.assign(selection,selection+4);  
    imin = IndexOfMinAmong(new_connect,index) ;
    if (imin==6 || imin==3){
      _pmP = createPrismWithVertices(N0,N5,N1,N3,N6,N2,_pmE->getClassification() );
      PrismToTet( _pmP);
      DEL(_pmP);
      _pmP =  createPrismWithVertices(N5,N0,N4,N6,N3,N7,_pmE->getClassification() );
      PrismToTet( _pmP);
      DEL(_pmP);
    }
    else {
      int selection[] = {4,5,6,7};   vector<int> index;   index.assign(selection,selection+4);  
      imin = IndexOfMinAmong(new_connect,index) ;
      if (imin==4 || imin==6){
	_pmP = createPrismWithVertices(N4,N6,N5,N0,N2,N1,_pmE->getClassification() );
	PrismToTet( _pmP);
	DEL(_pmP);
	_pmP = createPrismWithVertices(N4,N7,N6,N0,N3,N2,_pmE->getClassification() );
	PrismToTet( _pmP);
	DEL(_pmP);
      }
      else {
	int selection[] = {5,6,2,1};   vector<int> index;   index.assign(selection,selection+4);  
	imin = IndexOfMinAmong(new_connect,index) ;
	if (imin==1 || imin==6){
	  _pmP = createPrismWithVertices(N5,N6,N1,N4,N7,N0,_pmE->getClassification() );
	  PrismToTet( _pmP);
	  DEL(_pmP);
	  _pmP = createPrismWithVertices(N1,N6,N2,N0,N7,N3,_pmE->getClassification() );
	  PrismToTet( _pmP);
	  DEL(_pmP);
	}
	else {
	  createTetWithVertices(N5,N4,N7,N0,_pmE->getClassification() );
	  createTetWithVertices(N0,N7,N5,N2,_pmE->getClassification() );
	  createTetWithVertices(N5,N0,N2,N1,_pmE->getClassification() );
	  createTetWithVertices(N0,N7,N2,N3,_pmE->getClassification() );
	  createTetWithVertices(N5,N2,N7,N6,_pmE->getClassification() );
	}
      }
    }
  }

  int xMesh::IndexOfMinAmong(vector<int>& iValue , vector<int>& index){
    int imin=0;
    int iValMin=iValue[ index[imin] ];
    for (size_t i=1 ; i<index.size() ;i++){
      if(iValue[ index[i] ]<iValMin ){
	imin= i;
	iValMin=iValue[index[imin]];
      }
    }
    return index[imin];
  }

  void xMesh:: QuadToTri(mEntity* _pmE){
    vector<int> connect;
    connect.clear();
    for (int i=0 ; i<4 ;i++){	   
      connect.push_back(  _pmE->get(0,i)->getId() );
    }  
    int imin = min_element(connect.begin(), connect.end()) - connect.begin();
    vector<int> new_connect;
    new_connect.reserve(4);
    new_connect.clear();
    for (int i=0 ; i<4 ;i++){	   
      new_connect[(i-imin+4)%4]= connect[i];
    }  
    createFaceWithVertices( new_connect[0], new_connect[2], new_connect[3],_pmE->getClassification());
    createFaceWithVertices( new_connect[0], new_connect[1], new_connect[2],_pmE->getClassification());
  }
  //////////////////////////// END OF DEFINITIONS OF FUNCTIONS FOR SIMPLEXIFY ///////////////////







  //////////////////////////////// END OF DEFINITIONS OF XMESH FUNCTIONS ///////////////////

  /*
  void xAllCreator::create(const xMesh* m,const string& name)
  {
    for (int d = 0; d < 4; d++) 
      {
	for (xIter it=m->begin(d); it!=m->end(d); ++it ) 
	  {
	    mEntity* e = *it;
	    m->add_sub(e, name);
	  }
      }
    return;
  }
  */

  /*
  void xAddLayerCreator::create(const xMesh* m, const string& name)
  {
    for (int d = 0; d < 4; d++) 
      {
	for(xIter it=m->begin_sub(d,initial); it!=m->end_sub(d,initial); ++it)
	  {
	    m->add_sub(*it, name);
	  }
      }
    //now we add a couple of layers
    //based on the layer parameter
    const int d_sub = m->dim();
 
    for (int l = 1; l <= layers; l++) 
      {
	std::set<mEntity*> extra_layers;
	for(xIter it=m->begin_sub(dim_growth,name); it!=m->end_sub(dim_growth,name); ++it) 
	  {
	    mEntity* e = *it;
	    for (int j = 0; j < e->size(d_sub); j++) 
	      {
		mEntity* ee = e->get(d_sub, j);
		if (filter(ee)) extra_layers.insert(ee);
	      }
	  }
	for (std::set<mEntity*>::const_iterator it = extra_layers.begin(); 
	     it != extra_layers.end(); it++) 
	  {
	    mEntity* e = *it; 
	    m->add_sub(e, name);
	
	  }
	if (d_sub > 2)
	  {
	    m->modifyState_sub(3,2,true,name);
	    m->modifyState_sub(3,1,true,name);
	    m->modifyState_sub(3,0,true,name);
	  }
	if (d_sub > 1)
	  {
	    m->modifyState_sub(2,1,true,name);
	    m->modifyState_sub(2,0,true,name);
	  }
      }

    return;
  }

  int  xAddLayerModifier::tag() const{
    return 789;
  };

  void *  xAddLayerModifier::AP_alloc_and_fill_buffer (mEntity *e, AOMD_SharedInfo &si, int tag){
#ifdef PARALLEL
    const unsigned int  visited_tag  =  AOMD::AOMD_Util::Instance()->lookupMeshDataId ("visted_entity");
    int buffer_size = sizeof( mEntity* ) +  sizeof(int);
    void* buf = AP_alloc( si.pid(), this->tag(), buffer_size );
    mEntity** ebuf = reinterpret_cast< mEntity ** >( buf );
    *(ebuf++) = si.getRemotePointer();
    int* ibuf = reinterpret_cast< int* >( ebuf );
    *ibuf = e->getAttachedInt(visited_tag);
    return buf;
#endif
  };

  void  xAddLayerModifier::receiveData (int pid, void *buf){
#ifdef PARALLEL
    const unsigned int  visited_tag  =  AOMD::AOMD_Util::Instance()->lookupMeshDataId ("visted_entity");
    mEntity** ebuf = reinterpret_cast< mEntity** >( buf );
    mEntity*  e    = *(ebuf++);
    int * tag = reinterpret_cast< int* >(ebuf );
    if (*tag){
      e->attachInt(visited_tag, 1);
      received.insert(e);
    }    
#endif
  };

  void xAddLayerModifier::modify(const xMesh* m, const string& name) const
  {
    //  const bool debug = xdebug_flag;
    const bool debug = true;
    if (layers == 0) return; 
    received.clear();
    //now we add a couple of layers
    //based on the layer parameter
    m->modifyAllState_sub(name);
    const int d_sub = m->dim();
 
    Container  tocheckold;
    Container tochecknew;
    Container *tocheck1 = &tocheckold;
    Container *tocheck2 = &tochecknew;
    Container *tmp = 0;

  
    const unsigned int  visited_tag  = AOMD::AOMD_Util::Instance()->newMeshDataId("visted_entity");
  
    for(xIter it=m->begin_sub(dim_growth,name); it!=m->end_sub(dim_growth,name); ++it){
      (*it)->attachInt(visited_tag, 1); 
      tocheckold.insert(*it);
    }

    xMesh* m_ = const_cast <xMesh *> (m);
    m_->exchangeDataOnPartBdrys(const_cast< xAddLayerModifier& >( *this), true);
    // std::inserter( tocheckold, tocheckold.begin() );
    //std::inserter( tocheckold);
    //  for some reason inserter does not work:
    //std::copy(received.begin(), received.end(), std::inserter( tocheckold, tocheckold.begin() ) );
    //considerering using tr1::unordered_set instead of hashset
    // std::tr1::unordered_set<AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey, std::allocator
    //   <AOMD::mEntity*> > ga;
    //  POST HACK BEG
    // eXlibris_types::hash_set<AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey, std::allocator
    //   <AOMD::mEntity*> > ga;
    //  POST HACK END
    for (Container::const_iterator it= received.begin(); it!= received.end(); ++it)
      tocheckold.insert(*it);
    if (debug) std::cout << "received " << received.size() << std::endl;
    received.clear();
  
    for (int l = 1; l <= layers; l++) 
      {
	if (debug) std::cout << "proc " << ParUtil::Instance()->rank()+1 << " start layer " << l << std::endl;
	if (debug) std::cout << "proc " << ParUtil::Instance()->rank()+1 <<
		     " left to treat " << tocheck1->size();
	std::set<mEntity*> extra_layers;
	for(Container::iterator it= tocheck1->begin(); it!= tocheck1->end() ; ++it) 
	  {
	    mEntity* e = *it;
	  
	    for (int j = 0; j < e->size(d_sub); j++) 
	      {
		mEntity* ee = e->get(d_sub, j);
		if (filter(ee)) {
		  m->add_sub(ee, name);
		  for (int k=0; k < ee->size(dim_growth); ++k){
		    mEntity* eee = ee->get(dim_growth, k);
		    if (!(eee->getAttachedInt(visited_tag))){
		      tocheck2->insert(eee);
		      eee->attachInt(visited_tag, 1);
		    }
		  }
		}
	      }
	  }
	ParUtil::Instance()->Barrier(33, "b1");
	m_->exchangeDataOnPartBdrys(const_cast<xAddLayerModifier &>(*this), true);
	if (debug) std::cout << "proc " << ParUtil::Instance()->rank()+1
			     << " end layer " << l << " received " << received.size() << std::endl;
	ParUtil::Instance()->Barrier(33, "b2");
	for (Container::const_iterator it = received.begin(); it!= received.end(); ++it)
	  tocheck2->insert(*it);
	received.clear();

	tmp = tocheck1;
	tocheck1 = tocheck2;
	tocheck2 = tmp;
	tocheck2->clear();
      }
    if (debug) std::cout << "proc " << ParUtil::Instance()->rank()+1 << " finallyzing." << std::endl;
    m->modifyAllState_sub(name);
  
    for_each(m->begin_sub(dim_growth,name),m->end_sub(dim_growth,name), bind2nd( mem_fun(&mEntity::deleteData), visited_tag));
    AOMD::AOMD_Util::Instance()->deleteMeshDataId (visited_tag);


    if (debug) std::cout << "proc " << ParUtil::Instance()->rank()+1 << " end." << std::endl;
    return;
  }
  */




void xMesh::getPartitionEntity(mEntity* e, std::set<mEntity*>& partition) 
{
  const bool debug = xdebug_flag;
  if (debug) cout << " In getPartition2  e is ";
  if (debug)  e->print();
  xAttachableMesh *am = (xAttachableMesh*)e->getData(partition_tag);
  if (!am) 
    {
      if (debug) cout << " pas de am " << endl;
      partition.insert(e);
      return;
    }
  if (am)
    {
      xMesh* m =  am->m;
      if (debug) cout << " am existe " << endl;
      for (xIter it = m->begin(m->dim()); it != m->end(m->dim()); ++it)
	{
	  mEntity *e = *it;
	  if (debug) cout << " e in partition " << endl;
	  if (debug) e->print();
	  xMesh::getPartitionEntity(*it, partition);
	}
    }
}


void xMesh::getPartitionMesh(xMesh & partitionMesh) {
  xIter it =    begin(dim());
  xIter itend =  end(dim());
    
  while (it!=itend){
  void * pointerToData =   (*it)->getData(xMesh::partition_tag);
  if (!pointerToData) {
    std::cout << "No data attached to entity "<< (*it) << " with tag xMesh::partition_tag" << std::endl;
    throw;
  }
  xMesh *subpat =  (( xAttachableMesh*)
		   (pointerToData ))->m;
  xCopyMesh(subpat, &partitionMesh);
  ++it;
  }
  return;
}

void xMesh::getOctreeMesh(xMesh & octreeMesh) 
{
  for(xIter it=begin(dim()); it!=end(dim()); ++it)
    {
      mEntity* e = (mEntity*) *it;
      xAttachableMesh* octree_att = ( xAttachableMesh*) e->getData(xMesh::refined_elements_tag);
      std::cout << "(*it) is : "<< (*it)<< std::endl;
      if (!octree_att) {
	std::cout << "No data attached to entity "<< (*it) << " with tag xMesh::refined_elements_tag" << std::endl;
      }
      else{
	std::cout << " octree mesh tied to (*it) "<< std::endl;
	xMesh *octree = octree_att->m;
	xCopyMesh(octree, &octreeMesh);
      }
    }
  return;
}

/*mVertex * xMesh::createVertexUnique (const double x, const double y, const double z , pGEntity classif){
  Trellis_Util::mPoint p(x, y, z);
  std::map<Trellis_Util::mPoint, mVertex *>::iterator it =  pointVertexMap.find(p);
  if (it != pointVertexMap.end()) return it->second;
  mVertex *newV = createVertex(x, y, z, classif);
  pointVertexMap.insert(std::make_pair(p, newV ) );
  return newV;
  }*/
  

  void xCopyMesh(xMesh * pin ,  xMesh * pout, bool clean_tag_on_source_mesh) 
  {
    const bool debug = xdebug_flag;
    if (pin == pout) {
      std::cout << "warnning : copying a mesh onto itself" << std::endl;
      return;
    }
    //  should check first if in and out are not  0 ....
    xMesh & in =*pin;
    xMesh & out =*pout;
    
    xIter it = in.begin(0);
    xIter ite = in.end(0);
    while (it!=ite){
      mVertex *vin = (mVertex * ) (*it);
      if (debug)  cout << " vin is : " << vin << endl;
      if (debug)  vin->print();
      mVertex *vout = out.createVertex(vin->point()(0),vin->point()(1),vin->point()(2), vin->getClassification());
      if (debug) cout << " vout is : " << vout << endl; 
      if (debug)  vout->print();
      vin->attachEntity(xMesh::has_a_copy_tag, vout);
      vout->attachEntity(xMesh::is_the_copy_of_tag,vin);
      ++it;
    }
    
    it = in.begin(1);
    ite = in.end(1);
    while (it!=ite){
      mEdge *ein = (mEdge * ) (*it);
      mVertex *v0in = (mVertex *) ein->get(0,0);
      mVertex *v1in = (mVertex *) ein->get(0,1);
      mVertex *v0out = (mVertex * ) v0in->getAttachedEntity(xMesh::has_a_copy_tag);
      mVertex *v1out = (mVertex * ) v1in->getAttachedEntity(xMesh::has_a_copy_tag);
      mEdge* eout = out.createEdge(v0out, v1out, ein->getClassification());
      ein->attachEntity(xMesh::has_a_copy_tag,eout);
      eout->attachEntity(xMesh::is_the_copy_of_tag,ein);
      ++it;
    }
    
    it = in.begin(2);
    ite = in.end(2);
    while (it!=ite){
      mFace *fin = (mFace * ) (*it);
      int nnode = fin->size(0);
      switch (nnode){
      case 3:{
	mVertex *v0in =(mVertex *) fin->get(0,0);
	mVertex *v1in =(mVertex *) fin->get(0,1);
	mVertex *v2in = (mVertex *)fin->get(0,2);
	mVertex *v0out = (mVertex * ) v0in->getAttachedEntity(xMesh::has_a_copy_tag);
	mVertex *v1out = (mVertex * ) v1in->getAttachedEntity(xMesh::has_a_copy_tag);
	mVertex *v2out = (mVertex * ) v2in->getAttachedEntity(xMesh::has_a_copy_tag);
	mFace* fout = out.createFaceWithVertices(v0out, v1out, v2out, fin->getClassification());
	fin->attachEntity(xMesh::has_a_copy_tag,fout);
	fout->attachEntity(xMesh::is_the_copy_of_tag,fin);
	break;
      }
      case 4:{
	cout << "xCopyMesh , quadrangles : PAS ENCORE PROGRAMME ! " << endl;
	throw 1;
      }
      default:{
	throw 1;
      }
      }
      ++it;
    }
    
    it = in.begin(3);
    ite = in.end(3);
    while (it!=ite){
      mTet *tetin = (mTet * ) (*it);
      int nnode = tetin->size(0);
      switch (nnode){
      case 4:{
	mVertex *v0in =(mVertex *) tetin->get(0,0);
	mVertex *v1in =(mVertex *) tetin->get(0,1);
	mVertex *v2in = (mVertex *) tetin->get(0,2);
	mVertex *v3in = (mVertex *) tetin->get(0,3);
	mVertex *v0out = (mVertex * ) v0in->getAttachedEntity(xMesh::has_a_copy_tag);
	mVertex *v1out = (mVertex * ) v1in->getAttachedEntity(xMesh::has_a_copy_tag);
	mVertex *v2out = (mVertex * ) v2in->getAttachedEntity(xMesh::has_a_copy_tag);
	mVertex *v3out = (mVertex * ) v3in->getAttachedEntity(xMesh::has_a_copy_tag);
	mTet* tetout = out.createTetWithVertices(v0out, v1out, v2out, v3out,tetin->getClassification());
	tetin->attachEntity(xMesh::has_a_copy_tag,tetout);
	tetout->attachEntity(xMesh::is_the_copy_of_tag,tetin);
	break;
      }
      case 5:{
	cout << "xCopyMesh ,3D : PAS ENCORE PROGRAMME ! " << endl;
	throw 1;
      }
      default:{
	throw 1;
      }
      }
      ++it;
    }
    if(clean_tag_on_source_mesh){
      for (int i =0; i< 4 ; ++i){
	it = in.begin(i);
	ite = in.end(i);
	while (it!=ite){
	  mVertex *vin = (mVertex * ) (*it);
	  vin->deleteData(xMesh::has_a_copy_tag);
	  ++it;
	}
      }
    }
    
    return;
  }
  
  void getSupportBoundary(mEntity *in, std::set<mEntity *> &supportBoundary ){
    supportBoundary.clear();
    int nface = in->size(2);
    if (in->getLevel() == 0){
      for (int i=0; i < nface; ++i){
	mEntity * currentface = in->get(2, i);
	for (int k=0; k <3; ++k){
	  mEntity *currentedge = currentface->get(1, k);
	  if (currentedge->size(2) == 1) supportBoundary.insert(currentedge);
	  else {
	    mEntity * v0 = currentedge->get(0, 0);
	    mEntity * v1 = currentedge->get(0, 1);
	    if (v0 != in && v1 != in) supportBoundary.insert(currentedge); 
	  }
	}
      }
    }
    else if(in->getLevel() == 1){
      if (nface == 1) supportBoundary.insert(in);
      for (int i=0; i < nface; ++i){
	mEntity * currentface = in->get(2, i);
	for (int k=0; k <3; ++k){
	  mEntity *currentedge = currentface->get(1, k);
	  if (currentedge != in) supportBoundary.insert(currentedge);
	}
      }
    }
    else{
      cout << "Warning : not  a node  nor an edge in get support boundary " << endl;
    }
  }


  
// methode which construct a set of constitutive entity of entity e with level egale to d
void getSubentities(mEntity* e, const int d, std::vector<mEntity*> &l)
{
    // local variable
    int i;

    // get level of entity e
    int n = e->getLevel();

    // if the constitutive entity searched are of same level as e, or of higher level
    // this call don't have much sens
    // return imediatly leaving the caller with it's vector
    // => caller have the responsability of checking the corectness of l
    if (d >= n) return;

    // number of Constitutive entity of entity e
    n= e->size(d);

    // loop on constitutive entity
    for (i = 0; i < n; ++i) 
	l.push_back(e->get(d,i));

    // end
    return;
}


xMesh* xMesh::refineInNewMesh(int n, bool debut)//ne pas utiliser le 2nd param
{
  n--;
  xMesh* mesh_result=new xMesh;//Nouveau mesh (resultat)
  
  
  //Recopie des noeuds grossiers (avec les meme Ids) du raffinement precedent
  //--------------------------------------------------------------------------
  for(xIter it = begin(0); it!= end(0); it++)//Pour chaque noeud
  {
    mVertex* vOld=(mVertex*)(*it);//On recupere le noeud...
    Trellis_Util::mPoint p=vOld->point();//...puis son point
    mVertex* vNew = mesh_result->createVertex(vOld->getId(),p(0),p(1),p(2),vOld->getClassification());//On le cree dans le nouveau mesh
    
    //On les tags pour qu'ils se rappellent de quels noeuds il viennent
    if(debut)//Si c'est le premier remaillage
      vNew->attachEntity(was_duplicated_from_tag,vOld);//On lui dit qui est son ancetre
    else//Si ce n'est pas le premier remaillge
    {
      if(vOld->getAttachedEntity(was_created_by_tag))//On verifie si c'est un noeud originel
        vNew->attachEntity(was_duplicated_from_tag,vOld->getAttachedEntity(was_duplicated_from_tag));//On transmet le tag
    }
  }
  
  
  map<mEdge*, mVertex*>m_nodes;//Map des nouveaux noeuds
    
  //on cree les points en milieu  d'aretes
  //----------------------------------------
  for(xIter it=begin(2); it!=end(2);it++)//Pour chaque face
  {
    mFace* fOld=(mFace*)(*it);//On recupere la face grossiere
    
    //On recupere les 3 noeuds parents
    mVertex* v0=(mVertex*)fOld->get(0,0);
    mVertex* v1=(mVertex*)fOld->get(0,1);
    mVertex* v2=(mVertex*)fOld->get(0,2);
    Trellis_Util::mPoint p0=v0->point();
    Trellis_Util::mPoint p1=v1->point();
    Trellis_Util::mPoint p2=v2->point();
    
    //On recupere les edges (pour verifier les doublons)
    mEdge* e01=getEdge(v0,v1);
    mEdge* e12=getEdge(v1,v2);
    mEdge* e20=getEdge(v2,v0);
    
    //On cree les noeuds milieux s'ils n'existent pas
    mVertex *v01,*v12,*v20;
    mEdge *e0n, *e1n;
    if(m_nodes.find(e01)==m_nodes.end())//Si le noeud n'existe pas on le cree
                        {v01=mesh_result->createVertex(0.5*(p0(0)+p1(0)),
                                                       0.5*(p0(1)+p1(1)),
                                                       0.5*(p0(2)+p2(2)),
                                                       getEdge(v0,v1)->getClassification());//Nouveau noeud milieu
                        e0n = mesh_result->createEdge(v0->getId() ,v01->getId(), getEdge(v0,v1)->getClassification());//On cree les edge avec la bonne classification
                        e1n = mesh_result->createEdge(v01->getId() ,v1->getId(), getEdge(v0,v1)->getClassification());//On cree les edge avec la bonne classification
                        m_nodes[e01]=v01;

                        //On tag l'edge pour qu'elle connaisse sa maman
                        if(debut){
                            e0n->attachEntity(was_created_by_tag2,e01);
                            e1n->attachEntity(was_created_by_tag2,e01);
                        }else{
                            e0n->attachEntity(was_created_by_tag2,e01->getAttachedEntity(was_created_by_tag2));
                            e1n->attachEntity(was_created_by_tag2,e01->getAttachedEntity(was_created_by_tag2));
                        }
                        }
                        else v01=m_nodes[e01];//si le noeud existe deja, on le recupere

                if(m_nodes.find(e12)==m_nodes.end())//Si le noeud n'existe pas
                        {v12=mesh_result->createVertex(0.5*(p1(0)+p2(0)),
                                                       0.5*(p1(1)+p2(1)),
                                                       0.5*(p1(2)+p2(2)),
                                                       getEdge(v1,v2)->getClassification());
                        e0n = mesh_result->createEdge(v1->getId() ,v12->getId(), getEdge(v1,v2)->getClassification());//On cree les edge avec la bonne classification
                        e1n = mesh_result->createEdge(v12->getId() ,v2->getId(), getEdge(v1,v2)->getClassification());//On cree les edge avec la bonne classification
                        m_nodes[e12]=v12;

                        //On tag l'edge pour qu'elle connaisse sa maman
                        if(debut){
                            e0n->attachEntity(was_created_by_tag2,e12);
                            e1n->attachEntity(was_created_by_tag2,e12);
                        }else{
                            e0n->attachEntity(was_created_by_tag2,e12->getAttachedEntity(was_created_by_tag2));
                            e1n->attachEntity(was_created_by_tag2,e12->getAttachedEntity(was_created_by_tag2));
                        }
                        }
                        else  v12=m_nodes[e12];//si le noeud existe deja, on le recupere

                if(m_nodes.find(e20)==m_nodes.end())//Si le noeud n'existe pas
                        {v20=mesh_result->createVertex(0.5*(p2(0)+p0(0)),
                                                       0.5*(p2(1)+p0(1)),
                                                       0.5*(p2(2)+p0(2)),
                                                       getEdge(v2,v0)->getClassification());
                        e0n = mesh_result->createEdge(v2->getId() ,v20->getId(), getEdge(v2,v0)->getClassification());//On cree les edge avec la bonne classification
                        e1n = mesh_result->createEdge(v20->getId() ,v0->getId(), getEdge(v2,v0)->getClassification());//On cree les edge avec la bonne classification
                        m_nodes[e20]=v20;

                        //On tag l'edge pour qu'elle connaisse sa maman
                        if(debut){
                            e0n->attachEntity(was_created_by_tag2,e20);
                            e1n->attachEntity(was_created_by_tag2,e20);
                        }else{
                            e0n->attachEntity(was_created_by_tag2,e20->getAttachedEntity(was_created_by_tag2));
                            e1n->attachEntity(was_created_by_tag2,e20->getAttachedEntity(was_created_by_tag2));
                        }
                        }
                        else v20=m_nodes[e20];//si le noeud existe deja, on le recupere


                // On cree les nouvelles faces
                // -------------------------------
                
                mFace* f0=mesh_result->createFaceWithVertices (v0->getId() , v01->getId(), v20->getId(),fOld->getClassification());
                mFace* f1=mesh_result->createFaceWithVertices (v1->getId() , v12->getId(), v01->getId(),fOld->getClassification());
                mFace* f2=mesh_result->createFaceWithVertices (v2->getId() , v20->getId(), v12->getId(),fOld->getClassification());
                mFace* f3=mesh_result->createFaceWithVertices (v01->getId() , v12->getId(), v20->getId(),fOld->getClassification());
                
                
                //On attach l'ancienne face
                // -------------------------------
                if(debut)//Si c'est le premier remaillage, on tag la face parente
                        {
                          f0->attachEntity(was_created_by_tag2,fOld);
                          f1->attachEntity(was_created_by_tag2,fOld);
                          f2->attachEntity(was_created_by_tag2,fOld);
                          f3->attachEntity(was_created_by_tag2,fOld);
                        }
                        else//sinon, on transmet le tag initial
                        {
                          f0->attachEntity(was_created_by_tag2,fOld->getAttachedEntity(was_created_by_tag2));
                          f1->attachEntity(was_created_by_tag2,fOld->getAttachedEntity(was_created_by_tag2));
                          f2->attachEntity(was_created_by_tag2,fOld->getAttachedEntity(was_created_by_tag2));
                          f3->attachEntity(was_created_by_tag2,fOld->getAttachedEntity(was_created_by_tag2));
                        }
                        
                        //On attach les noeuds a l'ancienne face
                        // --------------------------------------
                        (mesh_result->getVertex(v0->getId()))->attachEntity(was_created_by_tag2,f0->getAttachedEntity(was_created_by_tag2));//V0
                        (mesh_result->getVertex(v1->getId()))->attachEntity(was_created_by_tag2,f0->getAttachedEntity(was_created_by_tag2));//V1
                        (mesh_result->getVertex(v2->getId()))->attachEntity(was_created_by_tag2,f0->getAttachedEntity(was_created_by_tag2));//V2
                        v01->attachEntity(was_created_by_tag2,f0->getAttachedEntity(was_created_by_tag2));
                        v12->attachEntity(was_created_by_tag2,f0->getAttachedEntity(was_created_by_tag2));
                        v20->attachEntity(was_created_by_tag2,f0->getAttachedEntity(was_created_by_tag2));
                        
                        
                        
  }//Fin de la boucle de remaillage -----------------
  
  classifyUnclassifiedVerices(mesh_result);
  mesh_result->modifyAllState();
  
  
  
  
  // RE-REMAILLAGE
  // ----------------
  if(n>0){//Si il faut re-remailler
  xMesh* mesh_result2;//On creer un nouveau maillage
  mesh_result2=mesh_result->refineInNewMesh(n, false);//On raffine dans ce nouveau maillage
  delete mesh_result;//On supprime l'ancien
  return mesh_result2;
  }
  else return mesh_result;//Si c'est le remaillage maximum
  
}


xSubMesh* xMesh::setBoundary()
{  
  bool debug=false;
  if(debug) cout<<"\nSet boundary *********\n  Dim = "<<dim()<<endl;
  if(boundary) deleteSubMesh("boundary");
  //boundary = new xSubMesh("boundary",*this);
  boundary = &(createSubMesh("boundary"));
  
  if(dim()==0)
  {
    cout<<"Warning : coded but not tested !!\n";
    boundary->add(*(this->begin(0)));
  }
  
  if(dim()==1)
  {
    cout<<"Warning : coded but not tested !!\n";
    for(xIter it=begin(0);it!=end(0);it++)//pour chaque noeud
    {
      mEntity* v=*it;//on recupere le noeud
      if(v->size(1)<=1) boundary->add(v);
    }
  }
  
  if(dim()==2)
  {
    for(xIter it=begin(1);it!=end(1);it++)//pour chaque edge
    {
      mEntity* e=*it;//on recupere les edge
      if(debug) {cout<<"  - Edge : ";e->print();}
      if(debug) cout<<"     Size upper : "<<e->size(2)<<endl;
        if(e->size(2)<=1){
          boundary->add(e);
          if(debug)cout<<"      -> Added\n";
        }
    }
  }
  
  if(dim()==3)
  {
    cout<<"Warning : coded but not tested !!\n";
    for(xIter it=begin(2);it!=end(2);it++)//pour chaque face
    {
      mEntity* f=*it;//on recupere le noeud
      if(f->size(3)<=1)
        boundary->add(f);
    }
  }
  
  
  if(debug) cout<<"End set boundary *********\n";
  boundary->modifyAllState();
  return boundary;
}


xSubMesh* xMesh::getBoundary()
{
  if(!boundary) setBoundary();
  //cout<<"adresse boundary : "<<boundary<<endl;
  return boundary;
}


bool xMesh::isOnBoundary(mEntity* e)
{
  getBoundary();
  return (boundary->find(e));
}

set<mVertex*> xMesh::getNeighbors(mVertex* v)
{
  set<mVertex*> neighbors;
  for(int i=0;i<v->size(1);i++)//Pour chaque arete
    {
      mEntity* e=v->get(1,i);
      if(e->get(0,0)!=v) neighbors.insert((mVertex*)e->get(0,0));
      else neighbors.insert((mVertex*)e->get(0,1));
    }
    return neighbors;
}

} // end of namespace
