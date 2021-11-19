/**************************************************************************** 

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written 
   and maintained by the Scientific Computation Research Center (SCOREC) at 
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*****************************************************************************/

#include <cassert>
#include <cstdio>
#include <iostream>
#include "mEntity.h"
#include "mException.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mEntityContainer.h"
#include "mAOMD.h"
#include "mCommonBdry.h"


#include "modeler.h"

#include <cstdlib>  
#include <cmath> 
#include <algorithm> 
#include <vector> 
#include <string> 

#ifdef HAVE_TOTALVIEW
#include "tv_data_display.h"
#endif

using std::cout;
using std::sort;
using std::set;
using std::list;
using std::string;
using std::vector;

namespace AOMD {

#if defined(DMUM) || defined(FLEXDB)
   mMesh* mEntity::g_mesh=0;
#endif

  // create an empty entity
  mEntity::mEntity()
  : theFastData(nullptr),theClassification((pGEntity)nullptr),
    theCommonBdry((mCommonBdry*)nullptr),iD(0)
  {

    for(int i=0;i<4;i++)theAdjacencies[i] = (mAdjacencyContainer*)nullptr;
  }
  // destructor
  mEntity::~mEntity()
  {
    for(int i=0;i<4;i++)
      {
	if(theAdjacencies[i])
	  delete theAdjacencies[i];
      }
    if(theFastData)delete theFastData;
  }


  Trellis_Util::mPoint mEntity::getCentroid()
  {
    //    const bool debug = xdebug_flag;
    //    if (debug) std::cout<<"dim de e : "<<this->getLevel()<<"\n";
    Trellis_Util::mPoint g;
    if (this->getLevel()==0){
      mVertex *v = (mVertex*) this->get(0,0);
      g =  v->point();      
    }
    else
      {
	for (int i=0; i<this->size(0); ++i) {
	  mVertex *v = (mVertex*) this->get(0,i);
	  g += v->point();
	}
	g *= 1./(double) this->size(0);
      }
    return g;
  }


  double mEntity::gSize ()
  {
    double lmin;
    //  std::cout << e->size(1) << " edges\n";

    if(getType() == TRI)
      {
	Trellis_Util::mPoint p0 = ((mVertex*)get(0,0))->point();
	Trellis_Util::mPoint p1 = ((mVertex*)get(0,1))->point();
	Trellis_Util::mPoint p2 = ((mVertex*)get(0,2))->point();
	Trellis_Util::mVector v1(p1,p0);
	Trellis_Util::mVector v2(p2,p0);
	Trellis_Util::mVector v3(p2,p1);
	
	Trellis_Util::mVector v1xv2 = v1 % v2;

	double semi_perimeter = 0.5 * (v1.normValue() + v2.normValue() + v3.normValue());
	double area = v1xv2.normValue();

	return area / semi_perimeter;
	
      }

    if(getType() == TET)
      {
	Trellis_Util::mPoint p0 = ((mVertex*)get(0,0))->point();
	Trellis_Util::mPoint p1 = ((mVertex*)get(0,1))->point();
	Trellis_Util::mPoint p2 = ((mVertex*)get(0,2))->point();
	Trellis_Util::mPoint p3 = ((mVertex*)get(0,3))->point();
	Trellis_Util::mVector v1(p1,p0);
	Trellis_Util::mVector v2(p2,p0);
	Trellis_Util::mVector v3(p3,p0);
	Trellis_Util::mVector v4(p1,p2);
	Trellis_Util::mVector v5(p1,p3);
	Trellis_Util::mVector v1v2 = v1 % v2;
	Trellis_Util::mVector v1v3 = v1 % v3;
	Trellis_Util::mVector v2v3 = v2 % v3;
	Trellis_Util::mVector v4v5 = v4 % v5;
	double V = v1v2 * v3;
	double A1 = sqrt (v1v2 * v1v2);
	double A2 = sqrt (v1v3 * v1v3);
	double A3 = sqrt (v2v3 * v2v3);
	double A4 = sqrt (v4v5 * v4v5);
	return fabs(V) / (A1+A2+A3+A4);
      }

    for(int i=0;i<size(1);i++)
      {
	Trellis_Util::mPoint p1 = ((mEdge*)get(1,i))->vertex(0)->point();
	Trellis_Util::mPoint p2 = ((mEdge*)get(1,i))->vertex(1)->point();
	Trellis_Util::mVector v(p2,p1);
	double l = sqrt (v(0)*v(0) + v(1)*v(1) + v(2)*v(2));
	//      std :: cout << "edge " << i << " length " << l << "\n";
	if(!i)lmin = l;
	else lmin = (l<lmin)?l:lmin;
      }
    return lmin;
  }

// end of gEntity functions
//*****************************************

  int mEntity::getNbTemplates(int what) const
  {
    throw new mException (__LINE__,__FILE__,"no template in mEntity base class");
  }

  mEntity* mEntity::getTemplate(int,int,int)const 
  {
    throw new mException (__LINE__,__FILE__,"no template in mEntity base class");
  }

  mEntity * mEntity::find(mEntity*me)const
  {
    mEntity *ret;
    int dim = me->getLevel();
    if(theAdjacencies[dim])
      if(ret = theAdjacencies[dim]->find(me))return ret;
    return (mEntity *)nullptr;
  }

  void mEntity::add(mEntity*me)
  {
    int dim = me->getLevel();
    if(!theAdjacencies[dim]) 
      theAdjacencies[dim] = new mAdjacencyContainer;
    theAdjacencies[dim]->add(me);
  }

  void mEntity::appendUnique(mEntity*me)
  {
    int dim = me->getLevel();
    if(!theAdjacencies[dim]) 
      theAdjacencies[dim] = new mAdjacencyContainer;
    theAdjacencies[dim]->appendUnique(me);
  }

  void mEntity::del(mEntity*me)
  {
    int dim = me->getLevel();
    if(theAdjacencies[dim])theAdjacencies[dim]->del(me);
  }

  void mEntity::classify(pGEntity g)
  {
    theClassification = g;
    if(parent())parent()->classify(g);
  }

  void mEntity::setCommonBdry(mCommonBdry* g)
  {
    theCommonBdry = g;
    if(parent())parent()->setCommonBdry(g);
  }

  pGEntity mEntity::getClassification() const
  {
    return theClassification;
  }

  mCommonBdry* mEntity::getCommonBdry() const
  {
    return theCommonBdry;
  }

  void mEntity::deleteAdjacencies(int what)
  {
    if(theAdjacencies[what])
      {
	delete theAdjacencies[what];
	theAdjacencies[what] = (mAdjacencyContainer*)nullptr;
      }
  }

  mAdjacencyContainer::iter mEntity::begin(int what)
  {
    if(!theAdjacencies[what])throw new mException (__LINE__,__FILE__, "Adjacencies Not Created");
    return theAdjacencies[what]->begin();
  }

  mAdjacencyContainer::iter mEntity::end(int what)
  {
    if(!theAdjacencies[what])throw new mException (__LINE__,__FILE__, "Adjacencies Not Created");
    return theAdjacencies[what]->end();
  }

  // are two mesh entities equal ?
  bool mEntity::equal (mEntity *other) const
  {
    if(this == other)return true;
    if(getId() != other->getId())return false;
    int dim = getLevel();
    if(dim != other->getLevel())return false;
    if(dim == 0)return getId() == other->getId();
  
    if(dim == 1)
      {
	int v11 = get(0,0)->getId();
	int v12 = get(0,1)->getId();
	int v21 = other->get(0,0)->getId();
	int v22 = other->get(0,1)->getId();
	if(v11 == v21 && v12 == v22)return true;
	if(v11 == v22 && v12 == v21)return true;
	return false;
      }
  
    int v1[100],v2[100];
    int s1 = size(0);
    int s2 = other->size(0);
    if(s1 != s2)return false;
    int i;
    for(i=0;i<s1;i++)
      {
	v1[i]=get(0,i)->getId();
	v2[i]=other->get(0,i)->getId();
      }
    sort(v1,v1+s1);
    sort(v2,v2+s2);
    for(i=0 ; i<s1 ;i++)
      {
	if(v1[i] != v2[i])return false;
      }
    return true;
  }

  // compare two entities
  bool mEntity::lessthan (mEntity *other) const
  {
    /* first do some fast checks**/
    int dim = getLevel();
    int other_dim = other->getLevel();
    if(dim != other_dim)return dim < other_dim;
    int id = getId();
    int other_id =  other->getId();
    if(dim == 0)return id < other_id;
    if(id < other_id)return true;
    if(id > other_id)return false;
  
    /*then do the whole comparison on vertices*/
    std::vector<mVertex*> thisVertices;
    std::vector<mVertex*> otherVertices;

    {  for(int i=0;i<getNbTemplates(0);i++)thisVertices.push_back((mVertex*)get(0,i));}
    {  for(int i=0;i<other->getNbTemplates(0);i++)otherVertices.push_back((mVertex*)other->get(0,i));}
    if(thisVertices.size() < otherVertices.size())return true;
    if(thisVertices.size() > otherVertices.size())return false;
    std::sort(thisVertices.begin(),thisVertices.end());
    std::sort(otherVertices.begin(),otherVertices.end());
    for( int i=0;i<thisVertices.size();i++)
      {      
	int id1 = thisVertices[i]->getId();
	int id2 = otherVertices[i]->getId();
	if(id1 < id2)return true;
	if(id1 > id2)return false;
      }  
    return false;
  }


  void mEntity::print() const
  {
  }

  void mEntity::setParent (mEntity* e)
  {
    unsigned int x = AOMD_Util::Instance()->getParent();
    mAttachableEntity *ai = (mAttachableEntity *)getData(x);
    if(!ai)
      {
	ai = new mAttachableEntity;
	attachData(x,ai);
      }
    ai->e = e;
  }

  mEntity* mEntity::parent ()
  {
    unsigned int x = AOMD_Util::Instance()->getParent();
    mAttachableEntity *ai = (mAttachableEntity *)getData(x);
    if(!ai)return nullptr;
    return ai->e;
  }

  mEntity* mEntity::root () 
  {
    mEntity *p = parent();
    if(!p)return this;
    else return p->root();
  }

  void mEntity::deleteParent ()
  {
    unsigned int x = AOMD_Util::Instance()->getParent();
    mAttachableEntity *ai = (mAttachableEntity *)getData(x);
    if(!ai)return;
    deleteData(x);
  }

  void recur_get_leaves(mEntity * e, list<mEntity*> &leaves, int n)
  {
    for(int i=0;i<e->size(n);i++)
      {
	mEntity *sub = e->get(n,i);
	if(sub->isAdjacencyCreated(n))recur_get_leaves(sub,leaves,n);
	else leaves.push_back(sub);
      }
  }

  void mEntity::getLeaves(list<mEntity*> &leaves)
  {
    int n = getLevel();
    if (isAdjacencyCreated(n))recur_get_leaves(this,leaves,n);
    else leaves.push_back(this);
  }

  void recur_get_all(mEntity * e, list<mEntity*> &leaves, int n)
  {
    for(int i=0;i<e->size(n);i++)
      {
	mEntity *sub = e->get(n,i);
	leaves.push_back(sub);
	if(sub->isAdjacencyCreated(n))recur_get_all(sub,leaves,n);
      }
  }

  void mEntity::getAllSubTree(list<mEntity*> &leaves)
  {
    int n = getLevel();
    leaves.push_back(this);
    if (isAdjacencyCreated(n))recur_get_all(this,leaves,n);
  }


  /* 
     This function gets all sub-entities of all dimensions for the element. 
     I'm trying to make that efficient
  */

  void mEntity:: getAllEntitiesOfAllDimensions(std::set<mEntity*> &the_whole_family)
  {
    int n = getLevel();
    the_whole_family.insert(this);
    list<mEntity *> family;
    getAllSubTree(family);
  
    for(std::list<mEntity*>::iterator it = family.begin() ; it != family.end() ; ++it)
      {
	mEntity *e = *it;
	the_whole_family.insert(e);
	for(int dim = 1; dim<n ;dim++)
	  {
	    if(e->isAdjacencyCreated(dim))
	      {

		for(int j=0; j<e->size(dim);j++)
		  {
		    mEntity *sub = e->get(dim,j);
		    list<mEntity *> subfamily;
		    if(the_whole_family.find(sub) == the_whole_family.end())
		      {
			sub->getAllSubTree(subfamily);
			sub->getAllEntitiesOfAllDimensions (the_whole_family);
			for(std::list<mEntity*>::iterator it2 = subfamily.begin() ; it2 != subfamily.end() ; ++it2)
			  {
			    mEntity *sub_f = *it2;
			  
			    the_whole_family.insert(sub_f);
			    for(int k=0; k<sub_f->size(0);k++)
			      the_whole_family.insert(sub_f->get(0,k));
			  }		    
		      }
		  }
	      }
	  }
      }
  }

  void mEntity:: getAllEntitiesOfAllDimensions(std::set<mEntity*,EntityLessThanKey> &the_whole_family)
  {
    int n = getLevel();
    the_whole_family.insert(this);
    list<mEntity *> family;
    getAllSubTree(family);
  
    for(std::list<mEntity*>::iterator it = family.begin() ; it != family.end() ; ++it)
      {
	mEntity *e = *it;
	the_whole_family.insert(e);
	for(int dim = 1; dim<n ;dim++)
	  {
	    if(e->isAdjacencyCreated(dim))
	      {

		for(int j=0; j<e->size(dim);j++)
		  {
		    mEntity *sub = e->get(dim,j);
		    list<mEntity *> subfamily;
		    if(the_whole_family.find(sub) == the_whole_family.end())
		      {
			sub->getAllSubTree(subfamily);
			sub->getAllEntitiesOfAllDimensions (the_whole_family);
			for(std::list<mEntity*>::iterator it2 = subfamily.begin() ; it2 != subfamily.end() ; ++it2)
			  {
			    mEntity *sub_f = *it2;
			  
			    the_whole_family.insert(sub_f);
			    for(int k=0; k<sub_f->size(0);k++)
			      the_whole_family.insert(sub_f->get(0,k));
			  }		    
		      }
		  }
	      }
	  }
      }
  }


  void mEntity::getHigherOrderUpward (int dim, mAdjacencyContainer &ents) const
  {
    if(isAdjacencyCreated(dim))
      {
	for(int i=0;i<size(dim);i++)ents.appendUnique(get(dim,i));
      }
    else 
      {
	int myDim = getLevel();
	for(int DIM = myDim+1; DIM!=dim ,DIM <= 3; DIM++)
	  {
	    if(isAdjacencyCreated(DIM))
	      {
		for(int i=0;i<size(DIM);i++)get(DIM,i)->getHigherOrderUpward(dim,ents);
		break;
	      }
	  }
      }
  }

  /**
     compute use of an entity
  */

  int mEntity :: getUse ( mEntity *e ) const
  {
    if (getLevel() == 2 && e->getLevel() == 1)
      {
	mEdge *ed = (mEdge*)e;
	mVertex *v1 = ed->vertex(0);
	mVertex *v2 = ed->vertex(1);
	int nbVert = size(0);
	for(int i=0;i<nbVert;i++)
	  {
	    mVertex *vv1 = (mVertex*)get(0,i);
	    mVertex *vv2 = (mVertex*)get(0,(i+1)%nbVert);
	    if(vv1 == v1 && vv2 == v2)return 1;
	    if(vv1 == v2 && vv2 == v1)return 0;
	  }
	throw 1;
      }

    int tetfacedir [] = {1,0,0,1};
    int dim = e->getLevel();
    assert (dim < getLevel());

    mEntity *v1 = e->get (0,0);
    mEntity *v2 = e->get (0,1);

    for(int i=0 ; i< getNbTemplates (dim) ; i++)
      {
	if ( get(dim,i) == e )
	  {
	    mEntity *other = getTemplate(i, dim, 0);
	    for(int j=0; j<other->size(0);j++)
	      {
		if(other->get(0,j) == v1)
		  {
		    mEntity *vtest = other->get(0,(j+1)%other->size(0));
		    if(vtest == v2)
		      {
			if(other != get(dim,i))delete other;
			if(dim == 2)return ( 1 - tetfacedir[i] );
			return 1;
		      }
		    else
		      {
			if(other != get(dim,i))delete other;
			if(dim == 2)return ( tetfacedir[i] );
			return 0;
		      }
		  }
	      }
	    delete other;
	  }
      }
    throw 1;
  }

#ifdef TSTT_
  mEntity::mTopology mEntity::getTopo()
  {
    int mDim;
    mDim=getLevel();
    switch(mDim)
      {
      case 0:
	return POINT;
	break;
      case 1:
	return LINE_SEGMENT;
	break;
      case 2:
	if     (size(0) == 3) return TRIANGLE;
	else if (size(0) == 4) return QUADRILATERAL;
	else if (size(0) >= 5) return POLYGON;
	break;
      case 3:
	if      (size(0) == 4 && size(2) == 4) return TETRAHEDRON;
	else if (size(0) == 8 && size(2) == 6) return HEXAHEDRON;
	else if (size(0) == 6 && size(2) == 5) return PRISM_;
	else if (size(0) == 5 && size(2) == 5) return PYRAMID_;
	else if (size(0) > 8)                  return POLYHEDRON;
	break;
      default:
	return UNDEFINED;
      }
  }
#endif

string mEntity::getUid() const
  {

    char cuid[256];
    //int uid;
    std::vector<int> uids;
    switch(getType())
    {
      case VERTEX : 
	uids.push_back(getId());
	break;
      default: 
        for (int i=0; i<size(0);++i)
          uids.push_back(get(0,i)->getId());
	std::sort(uids.begin(), uids.end());
	break;
    }
    
    switch(getType())
      {
      case VERTEX : 
	sprintf (cuid,"V%d",uids[0]);
	break;
      case EDGE : 
	sprintf (cuid,"E%d_%d",uids[0],uids[1]);
	break;
      case TRI : 
	sprintf (cuid,"F%d_%d_%d",uids[0],uids[1],uids[2]);
	break;  
      case QUAD : 
	sprintf (cuid,"F%d_%d_%d_%d",uids[0],uids[1],uids[2],uids[3]);
	break;
      case TET : 
	sprintf (cuid,"R%d_%d_%d_%d",uids[0],uids[1],uids[2],uids[3]);
	break;
      case HEX : 
	sprintf (cuid,"R%d_%d_%d_%d_%d_%d_%d_%d",uids[0],uids[1],uids[2],uids[3],
  		 				 uids[4],uids[5],uids[6],uids[7]);
	break;
    }

    return string(cuid);
  }


}

#ifdef HAVE_TOTALVIEW
int TV_ttf_display_type (  const AOMD::mEntity *e )
{
  return AOMD::mEntity::TV_ttf_display_type_ (e);
}
// Really durty here but realy kool : resetting definition of certain xfem attachable
// derived class so that they get recognized by dynamic_cast test. Clearly if they
// change in xfem output from tv display will be impredictible.
// This avoid introducing dependancy of xfem in AOMD
namespace xfem
{
class xAttachableChar : public AOMD::mAttachableData { 
    public: 
       virtual ~xAttachableChar(){}
       char c;
};
class xAttachableString : public AOMD::mAttachableData 
{ 
    public: 
       virtual ~xAttachableString(){}
       std::string s;
};
class xAttachableEntitiesSet : public AOMD::mAttachableData { 
    public: 
       virtual ~xAttachableEntitiesSet(){}
       std::set<AOMD::mEntity *> vect;
};
class xAttachableEntitiesVector : public AOMD::mAttachableData { 
    public: 
       virtual ~xAttachableEntitiesVector(){}
       std::vector<AOMD::mEntity *> vect;
};
}
// Same as above with some xcut class
namespace xcut
{
class xElemCutData : public AOMD::mAttachableData 
{
  public : 
    int key_code;
    int idx_code;
    std::vector<AOMD::mVertex*> cut_nodes;
    std::vector<AOMD::mEntity*> iso_elem ;
    std::vector<AOMD::mVertex*> warp_node ;
    std::vector< std::vector<AOMD::mVertex*> > warp_def ;
  private:
};
}

namespace AOMD {
int mEntity::TV_ttf_display_type_ (  const AOMD::mEntity *e )
{
    int   ret_val = TV_ttf_format_ok;
    int   err=TV_ttf_ec_ok;
    char name[512];

    // id
    if (e->theAdjacencies[0])
    {
        int i,n=e->theAdjacencies[0]->size();

        for(i=0;i<n&&err==TV_ttf_ec_ok ;++i)
        {
            sprintf(name, "node[%d]", i);
            err = TV_ttf_add_row ( name, TV_ttf_type_int, &((*(e->theAdjacencies[0]))[i]->iD) );
        }
    }
    else
    {
            err = TV_ttf_add_row ( "iD", TV_ttf_type_int, &(e->iD) );
    }
    if(err != TV_ttf_ec_ok) return TV_ttf_format_raw;

    // give acces to attachable
    if(std::distance(e->begin_attachdata(),e->end_attachdata()) )
    {
        mAttachableDataContainer::citer_attachdata it=e->begin_attachdata(); 
        mAttachableDataContainer::citer_attachdata ite=e->end_attachdata(); 
        for(;it!=ite&&err==TV_ttf_ec_ok ;++it)
        {
            const mAttachableDataContainer::info &a = *it;
            sprintf(name, "tag %ld", a.first);
            if (mAttachableInt * at = dynamic_cast<mAttachableInt *>(a.second))
                err = TV_ttf_add_row ( name, "in", &(at->i) );
            else if (mAttachableIntVector * at = dynamic_cast< mAttachableIntVector *>(a.second))
                err = TV_ttf_add_row ( name, "std::vector<int>", &(at->v) );
            else if (mAttachableDouble * at = dynamic_cast< mAttachableDouble *>(a.second))
                err = TV_ttf_add_row ( name, "double", &(at->d) );
            else if (mAttachableEntity * at = dynamic_cast< mAttachableEntity *>(a.second))
                err = TV_ttf_add_row ( name, "AOMD::mEntity *", &(at->e) );
            else if (mAttachablePoint * at = dynamic_cast< mAttachablePoint *>(a.second))
                err = TV_ttf_add_row ( name, "Trellis_Util::mpoint", &(at->p) );
            else if (mAttachablePointVector * at = dynamic_cast< mAttachablePointVector *>(a.second))
                err = TV_ttf_add_row ( name, "std::vector<Trellis_Util::mpoint>", &(at->v) );
            else if (mAttachable_mVector * at = dynamic_cast< mAttachable_mVector *>(a.second))
                err = TV_ttf_add_row ( name, "Trellis_Util::mVector", &(at->v) );
            else if (mAttachable_mTensor2 * at = dynamic_cast< mAttachable_mTensor2 *>(a.second))
                err = TV_ttf_add_row ( name, "Trellis_Util::mTensor2", &(at->t) );
            else if (xfem::xAttachableChar * at = dynamic_cast< xfem::xAttachableChar *>(a.second))
                err = TV_ttf_add_row ( name, "char", &(at->c) );
            else if (xfem::xAttachableString * at = dynamic_cast< xfem::xAttachableString *>(a.second))
                err = TV_ttf_add_row ( name, "std::string", &(at->s) );
            else if (xfem::xAttachableEntitiesSet * at = dynamic_cast< xfem::xAttachableEntitiesSet *>(a.second))
                err = TV_ttf_add_row ( name, "std::set<AOMD::mEntity *>", &(at->vect) );
            else if (xfem::xAttachableEntitiesVector * at = dynamic_cast< xfem::xAttachableEntitiesVector *>(a.second))
                err = TV_ttf_add_row ( name, "std::vector<AOMD::mEntity *>", &(at->vect) );
            else if (xcut::xElemCutData * at = dynamic_cast< xcut::xElemCutData *>(a.second))
            {
                sprintf(name, "tag %ld - key-code", a.first);
                err = TV_ttf_add_row ( name, "int", &(at->key_code) );
                sprintf(name, "tag %ld - idx-code", a.first);
                err = TV_ttf_add_row ( name, "int", &(at->idx_code) );
                sprintf(name, "tag %ld - full", a.first);
                err = TV_ttf_add_row ( name, "xcut::xElemCutData ", at );
            }
            else
                err = TV_ttf_add_row ( name, "AOMD::mAttachableData *", &(a.second) );
        }
        if(err != TV_ttf_ec_ok) return TV_ttf_format_raw;
    }

    // give acce to adjancy
    err = TV_ttf_add_row ( "theAdjancies", "AOMD::mAdjacencyContainer *[4]", &(e->theAdjacencies) );
    if(err != TV_ttf_ec_ok) return TV_ttf_format_raw;
        
    // classification
    err = TV_ttf_add_row ( "theClassification", "AOMD::pGEntity", &(e->theClassification) );
    if(err != TV_ttf_ec_ok) return TV_ttf_format_raw;

    // commonBdry
    err = TV_ttf_add_row ( "theCommonBdry", "AOMD::mCommonBdry *", &(e->theCommonBdry) );
    if(err != TV_ttf_ec_ok) return TV_ttf_format_raw;

  return ret_val;


}
}
#endif


